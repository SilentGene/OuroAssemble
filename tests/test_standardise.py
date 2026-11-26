import gzip

import standardise_assembly as sa
from _common import read_fasta, read_tsv


def test_metamdbg(tmp_path, make_seq):
    f = tmp_path / "contigs.fasta.gz"
    with gzip.open(f, "wt") as fh:
        fh.write(f">ctg1 length=3000 coverage=12.5 circular=yes\n{make_seq(3000)}\n")
        fh.write(f">ctg2 length=2000 coverage=3 circular=no\n{make_seq(2000)}\n")
        fh.write(f">ctg3 length=500 coverage=3 circular=yes\n{make_seq(500)}\n")
    contigs = sa.parse_metamdbg(str(f))
    assert [c.circular for c in contigs] == [True, False, True]
    assert contigs[0].depth == 12.5
    rows = sa.standardise("metamdbg", contigs, 1000)
    assert [r["candidate"] for r in rows] == [True, False, False]
    assert rows[0]["contig_id"] == "metamdbg__ctg1"


def test_hifiasm_gfa_and_rescue(tmp_path, make_seq):
    gfa = tmp_path / "asm.p_ctg.gfa"
    gfa.write_text(
        f"H\tVN:Z:1.0\nS\ts1.ctg000001c\t{make_seq(2000)}\tLN:i:2000\trd:i:30\n"
        f"S\ts2.ctg000002l\t{make_seq(1500)}\tLN:i:1500\trd:i:8\nA\tfoo\n"
    )
    rescue = tmp_path / "asm.rescue.fa"
    rescue.write_text(f">s9.ctg000009c\n{make_seq(1200)}\n")
    contigs = sa.parse_hifiasm_meta(str(gfa), str(rescue))
    assert [(c.orig_id, c.circular, c.source) for c in contigs] == [
        ("s1.ctg000001c", True, "main"),
        ("s2.ctg000002l", False, "main"),
        ("s9.ctg000009c", True, "rescue"),
    ]
    assert contigs[0].depth == 30


def test_metaflye(tmp_path, make_seq):
    fa = tmp_path / "assembly.fasta"
    fa.write_text(f">contig_1\n{make_seq(5000)}\n>contig_2\n{make_seq(4000)}\n")
    info = tmp_path / "assembly_info.txt"
    info.write_text(
        "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
        "contig_1\t5000\t20\tY\tN\t1\t*\t1\n"
        "contig_2\t4000\t5\tN\tY\t2\t*\t2\n"
    )
    contigs = sa.parse_metaflye(str(fa), str(info))
    assert [c.circular for c in contigs] == [True, False]
    assert contigs[1].repeat == "Y"


def test_hicanu_and_main(tmp_path, make_seq):
    fa = tmp_path / "asm.contigs.fasta"
    fa.write_text(
        f">tig00000001 len=3000 reads=50 class=contig suggestRepeat=no suggestBubble=no "
        f"suggestCircular=yes trim=0-3000\n{make_seq(3000)}\n"
        f">tig00000002 len=2500 reads=10 class=contig suggestRepeat=no suggestBubble=no "
        f"suggestCircular=no\n{make_seq(2500)}\n"
    )
    out = tmp_path / "out"
    out.mkdir()
    sa.main(
        [
            "--assembler", "hicanu", "--fasta", str(fa), "--min-circular-len", "1000",
            "--out-fasta", str(out / "c.fa"), "--out-tsv", str(out / "c.tsv"),
            "--out-candidates", str(out / "cand.fa"),
        ]
    )  # fmt: skip
    rows = read_tsv(out / "c.tsv")
    assert [r["circular"] for r in rows] == ["yes", "no"]
    assert "trim=0-3000" in rows[0]["notes"]
    assert [n for n, _d, _s in read_fasta(out / "cand.fa")] == ["hicanu__tig00000001"]
