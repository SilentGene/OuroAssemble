import merge_assemblies as ma
from eval_junctions import Aln, Params, scan

PRIORITY = ["metaflye", "metamdbg", "hifiasm_meta", "hicanu"]


def contig(cid, asm, length, circular="no"):
    return {"contig_id": cid, "assembler": asm, "length": str(length), "circular": circular}


def test_build_pool_tiers_and_redundant_circles():
    tables = {
        "metaflye": [contig("metaflye__c1", "metaflye", 2_000_000, "yes"),
                     contig("metaflye__l1", "metaflye", 50_000)],
        "metamdbg": [contig("metamdbg__c1", "metamdbg", 2_000_100, "yes"),
                     contig("metamdbg__l1", "metamdbg", 80_000),
                     contig("metamdbg__tiny", "metamdbg", 500)],
        "hicanu": [contig("hicanu__f1", "hicanu", 30_000, "yes")],
    }  # fmt: skip
    classification = {
        "metaflye__c1": {"length": "2000000", "status": "PASS", "selected": "yes",
                         "category": "chromosome", "n_suspicious_sites": "0"},
        "metamdbg__c1": {"length": "2000000", "status": "PASS", "selected": "yes",
                         "category": "chromosome", "n_suspicious_sites": "0"},
        "hicanu__f1": {"length": "29000", "status": "FAIL", "selected": "no",
                       "category": "unvalidated", "n_suspicious_sites": "1"},
    }  # fmt: skip
    clusters = {
        "metaflye__c1": {"representative": "yes", "cluster_id": "C0001"},
        "metamdbg__c1": {"representative": "no", "cluster_id": "C0001"},
    }
    scans = {
        "metamdbg__l1": {"n_suspicious_sites": "2"},
        "metaflye__l1": {"n_suspicious_sites": "0"},
    }
    pool = ma.build_pool(tables, scans, classification, clusters, 1000, True, PRIORITY)
    ids = [r["contig_id"] for r in pool]
    # circle representative first; non-representative copy and tiny contig left out
    assert ids[0] == "metaflye__c1"
    assert "metamdbg__c1" not in ids and "metamdbg__tiny" not in ids
    # failed circle enters as linear, using its trimmed length
    failed = next(r for r in pool if r["contig_id"] == "hicanu__f1")
    assert failed["tier"] == 2 and failed["category"] == "circular_unverified"
    assert failed["length"] == 29000
    # among linear contigs fewer suspicious sites beats length
    assert ids.index("metaflye__l1") < ids.index("metamdbg__l1")
    assert [r["rank"] for r in pool] == list(range(1, len(pool) + 1))


def test_greedy_select_drops_contained_keeps_variants(tmp_path):
    pool = [
        {"rank": "1", "contig_id": "A", "length": "100000"},
        {"rank": "2", "contig_id": "B", "length": "40000"},  # inside A at 99.9 %
        {"rank": "3", "contig_id": "C", "length": "40000"},  # inside A but only 97 % (strain)
        {"rank": "4", "contig_id": "D", "length": "40000"},  # half overlaps A
    ]
    paf = tmp_path / "x.paf"
    paf.write_text(
        "B\t40000\t0\t40000\t+\tA\t100000\t10000\t50000\t39960\t40000\t60\n"
        "C\t40000\t0\t40000\t+\tA\t100000\t50000\t90000\t38800\t40000\t60\n"
        "D\t40000\t0\t20000\t+\tA\t100000\t80000\t100000\t20000\t20000\t60\n"
    )
    hits = ma.parse_paf(str(paf), 0.99, 500)
    out = {r["contig_id"]: r for r in ma.greedy_select(pool, hits, 0.9)}
    assert out["A"]["kept"] and out["C"]["kept"] and out["D"]["kept"]
    assert not out["B"]["kept"] and out["B"]["redundant_with"] == "A"
    assert out["D"]["covered_frac"] == 0.5


def test_lower_ranked_cannot_remove_higher(tmp_path):
    # A (rank 1) is contained in B (rank 2): A is kept first, so B stays unless covered
    pool = [{"rank": "1", "contig_id": "A", "length": "10000"},
            {"rank": "2", "contig_id": "B", "length": "100000"}]  # fmt: skip
    paf = tmp_path / "x.paf"
    paf.write_text("A\t10000\t0\t10000\t+\tB\t100000\t0\t10000\t10000\t10000\t60\n")
    out = {
        r["contig_id"]: r for r in ma.greedy_select(pool, ma.parse_paf(str(paf), 0.99, 500), 0.9)
    }
    assert out["A"]["kept"] and out["B"]["kept"] and out["B"]["covered_frac"] == 0.1


def test_scan_flags_clip_pileup():
    length = 50_000
    alns = [Aln(f"t{i}", max(0, s), min(length, s + 12000), 60, [(0, 12000)])
            for i, s in enumerate(range(-10000, length, 1500))]  # fmt: skip
    alns += [Aln(f"c{i}", 14000, 20000, 60, [(0, 6000), (4, 3000)]) for i in range(6)]
    res = scan(alns, length, Params())
    assert res["n_clip_sites"] == 1 and res["n_low_cov_regions"] == 0
