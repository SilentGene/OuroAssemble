#!/usr/bin/env python3
"""Write the final OuroAssemble tables, FASTAs and a single self-contained HTML report."""

from __future__ import annotations

import argparse
import datetime as dt
import html
import json
from pathlib import Path

from _common import n50, read_fasta, read_tsv, setup_logging, to_float, write_fasta, write_tsv
from summarise_classification import COLUMNS as CLASS_COLUMNS

CLUSTER_COLUMNS = ["cluster_id", "cluster_size", "representative", "representative_id"]
COLUMNS = CLASS_COLUMNS + CLUSTER_COLUMNS
FASTA_FOR = {
    "chromosome": "chromosomes.fa",
    "plasmid": "plasmids.fa",
    "virus": "viruses.fa",
    "large_circular": "other_circular.fa",
    "other_circular": "other_circular.fa",
}
MERGE_NUMERIC = {"circular_pass", "pct_reads_mapped", "pct_of_merged", "dropped_redundant",
                 "circular"}  # fmt: skip
STATUS_COLOURS = {"PASS": "var(--pass)", "PENDING": "var(--pending)", "FAIL": "var(--fail)"}


# -----------------------------------------------------------------------------
# data
# -----------------------------------------------------------------------------
def assembly_stats(contig_tables: list[str], rows: list[dict]) -> list[dict]:
    out = []
    for path in contig_tables:
        contigs = read_tsv(path)
        asm = contigs[0]["assembler"] if contigs else Path(path).parent.name
        lengths = [int(c["length"]) for c in contigs]
        cands = [r for r in rows if r["assembler"] == asm]
        out.append(
            {
                "assembler": asm,
                "contigs": len(contigs),
                "total_bp": sum(lengths),
                "n50": n50(lengths),
                "longest": max(lengths, default=0),
                "circular": sum(c["circular"] == "yes" for c in contigs),
                "candidates": len(cands),
                "PASS": sum(r["status"] == "PASS" for r in cands),
                "PENDING": sum(r["status"] == "PENDING" for r in cands),
                "FAIL": sum(r["status"] == "FAIL" for r in cands),
            }
        )
    return out


def merge_comparison(asm_stats, selection, mapping) -> tuple[list[dict], list[dict]]:
    """(assembler vs merged comparison rows, per-assembler contribution rows)."""
    pct = {r["assembly"]: r["pct_reads_mapped"] for r in mapping}
    kept = [r for r in selection if r["kept"] == "yes"]
    comparison = [
        {
            "assembly": s["assembler"],
            "contigs": s["contigs"],
            "total_bp": s["total_bp"],
            "n50": s["n50"],
            "circular_pass": s["PASS"],
            "pct_reads_mapped": pct.get(s["assembler"], ""),
        }
        for s in asm_stats
    ]
    lengths = [int(r["length"]) for r in kept]
    comparison.append(
        {
            "assembly": "merged",
            "contigs": len(kept),
            "total_bp": sum(lengths),
            "n50": n50(lengths),
            "circular_pass": sum(r["tier"] == "0" for r in kept),
            "pct_reads_mapped": pct.get("merged", ""),
        }
    )
    total = sum(lengths) or 1
    contribution = []
    for s in asm_stats:
        mine = [r for r in kept if r["assembler"] == s["assembler"]]
        bp = sum(int(r["length"]) for r in mine)
        contribution.append(
            {
                "assembly": s["assembler"],
                "contigs": len(mine),
                "total_bp": bp,
                "circular": sum(r["topology"] == "circular" for r in mine),
                "pct_of_merged": round(100 * bp / total, 2),
                "dropped_redundant": sum(
                    r["assembler"] == s["assembler"] and r["kept"] != "yes" for r in selection
                ),
            }
        )
    return comparison, contribution


def merge_rows(classification: list[dict], clusters: list[dict]) -> list[dict]:
    by_id = {c["contig_id"]: c for c in clusters}
    for row in classification:
        row.update({k: by_id.get(row["contig_id"], {}).get(k, "") for k in CLUSTER_COLUMNS})
    return classification


def write_fastas(rows, seqs, outdir: Path, sample: str) -> dict[str, int]:
    buckets: dict[str, list] = {name: [] for name in set(FASTA_FOR.values())}
    buckets["pending_review.fa"] = []
    for r in rows:
        if r.get("representative") != "yes":
            continue
        target = FASTA_FOR.get(r["category"])
        if target is None:
            continue
        if r["status"] == "PENDING":
            target = "pending_review.fa"
        tax = r.get("gtdb_classification") or r.get("genomad_taxonomy") or ""
        header = (
            f"{r['contig_id']} sample={sample} cluster={r['cluster_id']} "
            f"category={r['category']} status={r['status']} length={r['length']}"
        )
        if tax:
            header += f" taxonomy={tax.replace(' ', '_')}"
        buckets[target].append((header, seqs[r["contig_id"]]))
    return {name: write_fasta(recs, outdir / name) for name, recs in buckets.items()}


# -----------------------------------------------------------------------------
# HTML helpers
# -----------------------------------------------------------------------------
def esc(value) -> str:
    return html.escape("" if value is None else str(value))


def fmt_int(value) -> str:
    v = to_float(value)
    return f"{int(v):,}" if v is not None else esc(value)


def table(rows: list[dict], columns: list[tuple[str, str]], table_id: str, numeric=()) -> str:
    head = "".join(
        f'<th data-col="{i}" data-num="{int(k in numeric)}">{esc(label)}</th>'
        for i, (k, label) in enumerate(columns)
    )
    body = []
    for r in rows:
        cells = []
        for k, _label in columns:
            v = r.get(k, "")
            if k == "status" and v in STATUS_COLOURS:
                cells.append(
                    f'<td><span class="pill" style="--c:{STATUS_COLOURS[v]}">{v}</span></td>'
                )
            elif k in numeric:
                as_int = k.endswith(("length", "bp", "n50", "longest", "contigs"))
                shown = fmt_int(v) if as_int else esc(v)
                cells.append(f'<td class="num" data-v="{esc(v)}">{shown}</td>')
            else:
                cells.append(f"<td>{esc(v)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    if not body:
        body.append(f'<tr><td colspan="{len(columns)}" class="empty">none</td></tr>')
    return (
        f'<div class="tablewrap"><table id="{table_id}" class="sortable">'
        f"<thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"
    )


def status_chart(stats: list[dict]) -> str:
    """Horizontal stacked bars of PASS / PENDING / FAIL candidates per assembler (inline SVG)."""
    if not stats:
        return ""
    width, bar_h, gap, left = 640, 26, 12, 120
    peak = max([s["candidates"] for s in stats] + [1])
    height = len(stats) * (bar_h + gap) + 30
    parts = [
        f'<svg viewBox="0 0 {width} {height}" class="chart" role="img" '
        f'aria-label="circular candidates per assembler by junction status">'
    ]
    for i, s in enumerate(stats):
        y = i * (bar_h + gap)
        parts.append(
            f'<text x="{left - 8}" y="{y + bar_h * 0.68}" text-anchor="end" class="lbl">'
            f"{esc(s['assembler'])}</text>"
        )
        x = left
        for status in ("PASS", "PENDING", "FAIL"):
            n = s[status]
            if not n:
                continue
            w = (width - left - 60) * n / peak
            parts.append(
                f'<rect x="{x:.1f}" y="{y}" width="{max(w, 1):.1f}" height="{bar_h}" rx="3" '
                f'fill="{STATUS_COLOURS[status]}"><title>{status}: {n}</title></rect>'
            )
            if w > 24:
                parts.append(
                    f'<text x="{x + w / 2:.1f}" y="{y + bar_h * 0.68}" text-anchor="middle" '
                    f'class="inbar">{n}</text>'
                )
            x += w
        parts.append(
            f'<text x="{x + 6:.1f}" y="{y + bar_h * 0.68}" class="lbl">{s["candidates"]}</text>'
        )
    ly = len(stats) * (bar_h + gap) + 12
    for j, status in enumerate(("PASS", "PENDING", "FAIL")):
        lx = left + j * 110
        parts.append(
            f'<rect x="{lx}" y="{ly}" width="12" height="12" rx="2" '
            f'fill="{STATUS_COLOURS[status]}"/>'
            f'<text x="{lx + 18}" y="{ly + 10}" class="lbl">{status}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


CSS = """
:root{--bg:#fbfbfa;--fg:#1d1d1b;--muted:#6b6b66;--card:#ffffff;--line:#e4e3df;--accent:#2f5d8a;
--pass:#2e8b57;--pending:#d49a1a;--fail:#c0453b;--th:#f3f2ee}
@media (prefers-color-scheme:dark){:root{--bg:#161615;--fg:#ecebe6;--muted:#a3a29b;--card:#1f1f1d;
--line:#34332f;--accent:#7fb0e0;--pass:#4fb37a;--pending:#e0b04a;--fail:#e06a5f;--th:#262624}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);
font:15px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:1200px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:1.7rem;margin:0 0 4px}h3{font-size:1rem;margin:18px 0 8px}
h2{font-size:1.2rem;margin:36px 0 10px;
border-bottom:1px solid var(--line);padding-bottom:6px}.sub{color:var(--muted);margin:0}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin-top:20px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 14px}
.card .v{font-size:1.5rem;font-weight:650}.card .k{color:var(--muted);font-size:.85rem}
.tablewrap{overflow-x:auto;border:1px solid var(--line);border-radius:8px;background:var(--card)}
table{border-collapse:collapse;width:100%;font-size:.86rem}
th,td{padding:6px 10px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}
th{background:var(--th);cursor:pointer;position:sticky;top:0;user-select:none}
th:after{content:" \\2195";color:var(--muted);font-size:.75em}
td.num{text-align:right;font-variant-numeric:tabular-nums}td.empty{color:var(--muted);text-align:center}
.pill{display:inline-block;padding:1px 8px;border-radius:999px;color:#fff;background:var(--c);
font-size:.78rem;font-weight:600}
.chart{width:100%;max-width:760px;height:auto}.chart .lbl{fill:var(--fg);font-size:12px}
.chart .inbar{fill:#fff;font-size:11px;font-weight:600}
input.filter{margin:0 0 8px;padding:6px 10px;border:1px solid var(--line);border-radius:6px;
background:var(--card);color:var(--fg);width:100%;max-width:360px}
details{margin-top:8px}summary{cursor:pointer;color:var(--accent)}
pre{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px;
overflow:auto;font-size:.8rem}
dl.legend{display:grid;grid-template-columns:max-content 1fr;gap:4px 14px;font-size:.9rem}
dl.legend dt{font-weight:600}
"""

JS = """
document.querySelectorAll('table.sortable th').forEach(th=>{th.addEventListener('click',()=>{
 const t=th.closest('table'),tb=t.tBodies[0],i=+th.dataset.col,num=th.dataset.num==='1';
 const asc=th.dataset.dir!=='asc';t.querySelectorAll('th').forEach(h=>delete h.dataset.dir);
 th.dataset.dir=asc?'asc':'desc';
 const val=r=>{const c=r.cells[i];if(!c)return'';const v=c.dataset.v??c.textContent;
  return num?(parseFloat(v)||-Infinity):v.toLowerCase()};
 [...tb.rows].sort((a,b)=>{const x=val(a),y=val(b);return (x>y?1:x<y?-1:0)*(asc?1:-1)})
  .forEach(r=>tb.appendChild(r));});});
document.querySelectorAll('input.filter').forEach(inp=>{inp.addEventListener('input',()=>{
 const q=inp.value.toLowerCase();document.querySelectorAll('#'+inp.dataset.table+' tbody tr')
 .forEach(r=>{r.style.display=r.textContent.toLowerCase().includes(q)?'':'none'})})});
"""


def card(value, label) -> str:
    return (
        f'<div class="card"><div class="v">{esc(value)}</div>'
        f'<div class="k">{esc(label)}</div></div>'
    )


def merged_section(comparison, contribution, numeric) -> str:
    if not comparison:
        return ""
    comp_cols = [
        ("assembly", "assembly"),
        ("contigs", "contigs"),
        ("total_bp", "total bp"),
        ("n50", "N50"),
        ("circular_pass", "validated circles"),
        ("pct_reads_mapped", "% reads mapped"),
    ]
    contrib_cols = [
        ("assembly", "assembler"),
        ("contigs", "contigs kept"),
        ("total_bp", "bp kept"),
        ("circular", "circles kept"),
        ("pct_of_merged", "% of merged bp"),
        ("dropped_redundant", "dropped as redundant"),
    ]
    return f"""
<h2>Merged assembly</h2>
<p class="sub">One non-redundant assembly built from all assemblers by contig selection only
(no contigs are cut or joined): validated circles first, then linear contigs ranked by
HiFi read support, length and assembler priority; a contig is dropped when it is already
covered by kept contigs. <code>final/merged_assembly.fa</code>,
provenance in <code>final/merged_assembly.tsv</code>.</p>
{table(comparison, comp_cols, "mergecmp", numeric)}
<h3>Contribution of each assembler</h3>
{table(contribution, contrib_cols, "mergecontrib", numeric)}
"""


def build_html(
    sample, read_stats, asm_stats, rows, fasta_counts, cfg, comparison=(), contribution=()
) -> str:
    reps = [r for r in rows if r.get("representative") == "yes"]
    final = [r for r in reps if r["status"] == "PASS"]
    by_cat = {}
    for r in final:
        by_cat[r["category"]] = by_cat.get(r["category"], 0) + 1
    filtered = read_stats[-1] if read_stats else {}

    cards = "".join(
        [
            card(fmt_int(filtered.get("num_seqs", "")), "reads after filtering"),
            card(fmt_int(filtered.get("sum_len", "")), "bases after filtering"),
            card(len(rows), "circular candidates"),
            card(sum(r["status"] == "PASS" for r in rows), "junction PASS"),
            card(by_cat.get("chromosome", 0), "complete chromosomes"),
            card(by_cat.get("plasmid", 0), "plasmids"),
            card(by_cat.get("virus", 0), "viruses"),
            card(sum(r["status"] == "PENDING" for r in reps), "pending review"),
        ]
        + ([card(fmt_int(comparison[-1]["n50"]), "merged assembly N50")] if comparison else [])
    )
    read_cols = [
        ("file", "file"),
        ("num_seqs", "reads"),
        ("sum_len", "bases"),
        ("min_len", "min"),
        ("avg_len", "mean"),
        ("max_len", "max"),
        ("N50", "N50"),
    ]
    rs_rows = [{**r, "file": Path(r.get("file", "")).name} for r in read_stats]
    asm_cols = [
        ("assembler", "assembler"),
        ("contigs", "contigs"),
        ("total_bp", "total bp"),
        ("n50", "N50"),
        ("longest", "longest"),
        ("circular", "circular"),
        ("candidates", "candidates"),
        ("PASS", "PASS"),
        ("PENDING", "PENDING"),
        ("FAIL", "FAIL"),
    ]
    final_cols = [
        ("cluster_id", "cluster"),
        ("contig_id", "representative"),
        ("assembler", "assembler"),
        ("length", "length"),
        ("category", "category"),
        ("status", "status"),
        ("spanning_reads", "spanning"),
        ("median_depth", "depth"),
        ("checkm2_completeness", "completeness"),
        ("checkm2_contamination", "contamination"),
        ("gtdb_classification", "GTDB taxonomy"),
        ("length_ratio", "length / ref"),
        ("genomad_class", "geNomad"),
        ("genomad_score", "score"),
        ("cluster_size", "n assemblers"),
    ]
    cand_cols = [
        ("contig_id", "contig"),
        ("assembler", "assembler"),
        ("orig_length", "raw length"),
        ("overlap_len", "trimmed overlap"),
        ("length", "length"),
        ("depth_at_junction", "depth at J"),
        ("unique_depth_at_junction", "unique depth"),
        ("spanning_reads", "spanning"),
        ("span_ratio", "span ratio"),
        ("clipped_at_junction", "clipped at J"),
        ("low_mapq_frac", "low-MAPQ frac"),
        ("n_suspicious_sites", "suspicious sites"),
        ("status", "status"),
        ("reason", "reason"),
        ("category", "category"),
        ("cluster_id", "cluster"),
    ]
    numeric = {
        "length", "orig_length", "num_seqs", "sum_len", "min_len", "avg_len", "max_len", "N50",
        "contigs", "total_bp", "n50", "longest", "circular", "candidates", "PASS", "PENDING",
        "FAIL", "spanning_reads", "median_depth", "checkm2_completeness",
        "checkm2_contamination", "length_ratio", "genomad_score", "cluster_size", "overlap_len",
        "depth_at_junction", "unique_depth_at_junction", "span_ratio", "clipped_at_junction",
        "low_mapq_frac", "n_suspicious_sites",
    }  # fmt: skip
    fasta_list = "".join(
        f"<li><code>final/{esc(k)}</code> — {v} sequences</li>"
        for k, v in sorted(fasta_counts.items())
    )
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>OuroAssemble report</title><style>{CSS}</style></head><body><main>
<h1>OuroAssemble report — {esc(sample)}</h1>
<p class="sub">Circular genomes from PacBio HiFi metagenome assemblies ·
generated {dt.datetime.now().strftime("%Y-%m-%d %H:%M")}</p>
<div class="cards">{cards}</div>

<h2>Final circular genomes (dereplicated, junction PASS)</h2>
<input class="filter" data-table="final" placeholder="filter rows…">
{table(final, final_cols, "final", numeric)}
<ul>{fasta_list}</ul>

{merged_section(list(comparison), list(contribution), numeric | MERGE_NUMERIC)}

<h2>Pending review</h2>
<p class="sub">Representatives whose junction or contig could not be confirmed (low coverage,
repeat at the junction, or suspicious sites elsewhere). Inspect them in IGV with
<code>05_validation/&lt;assembler&gt;/igv_batch.txt</code>.</p>
{table([r for r in reps if r["status"] == "PENDING"], final_cols, "pending", numeric)}

<h2>Junction validation per assembler</h2>
{status_chart(asm_stats)}
{table(asm_stats, asm_cols, "asm", numeric)}

<h2>All circular candidates</h2>
<input class="filter" data-table="cands" placeholder="filter rows…">
{table(rows, cand_cols, "cands", numeric)}

<h2>Reads</h2>
{table(rs_rows, read_cols, "reads", numeric)}

<h2>Definitions</h2>
<dl class="legend">
<dt>PASS</dt><dd>&ge; min_spanning reads (MAPQ &ge; min_mapq) cross the circular junction with
&ge; flank bp anchored on both sides, spanning/unique depth &ge; min_span_ratio, and no other
suspicious sites in the contig.</dd>
<dt>PENDING</dt><dd>Low coverage at the junction, a repeat at the junction (many low-MAPQ
reads), or a supported junction with suspicious sites elsewhere.</dd>
<dt>FAIL</dt><dd>Enough coverage but the junction is not crossed by reads.</dd>
<dt>chromosome</dt><dd>&ge; min_chrom_len and CheckM2 completeness &ge; threshold and
contamination &le; threshold.</dd>
<dt>length / ref</dt><dd>contig length divided by the genome size of the GTDB-Tk closest
reference.</dd>
</dl>
<details><summary>Parameters</summary><pre>{esc(json.dumps(cfg, indent=2))}</pre></details>
</main><script>{JS}</script></body></html>
"""


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--read-stats", required=True)
    ap.add_argument("--contigs", nargs="+", required=True)
    ap.add_argument("--classification", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--config-json", default="{}")
    ap.add_argument("--merge-selection", help="08_merge/selection.tsv")
    ap.add_argument("--mapping-stats", help="08_merge/mapping_stats.tsv")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    log = setup_logging()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = merge_rows(read_tsv(args.classification), read_tsv(args.clusters))
    rows.sort(key=lambda r: (r.get("cluster_id") or "Z", r["contig_id"]))
    write_tsv(rows, outdir / "circular_genomes.tsv", COLUMNS)

    seqs = {name: seq for name, _d, seq in read_fasta(args.fasta)}
    fasta_counts = write_fastas(rows, seqs, outdir, args.sample)

    asm_stats = assembly_stats(args.contigs, rows)
    comparison, contribution = [], []
    if args.merge_selection and args.mapping_stats:
        comparison, contribution = merge_comparison(
            asm_stats, read_tsv(args.merge_selection), read_tsv(args.mapping_stats)
        )
    page = build_html(
        args.sample,
        read_tsv(args.read_stats),
        asm_stats,
        rows,
        fasta_counts,
        json.loads(args.config_json),
        comparison,
        contribution,
    )
    (outdir / "report.html").write_text(page)
    log.info("report written; final FASTAs: %s", fasta_counts)


if __name__ == "__main__":
    main()
