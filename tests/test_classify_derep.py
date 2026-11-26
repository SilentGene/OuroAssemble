import gzip

import dereplicate as dr
import summarise_classification as sc
from _common import read_tsv, write_tsv


def test_genome_sizes_stream(tmp_path):
    meta = tmp_path / "bac120_metadata.tsv.gz"
    with gzip.open(meta, "wt") as fh:
        fh.write("accession\tambiguous_bases\tgenome_size\n")
        fh.write("RS_GCF_000001.1\t0\t4000000\nGB_GCA_000002.1\t0\t2500000\n")
    sizes = sc.genome_sizes([str(meta)], {"GCF_000001.1", "GCA_000002.1", "GCF_9"})
    assert sizes == {"GCF_000001.1": 4000000, "GCA_000002.1": 2500000}


def test_genomad_calls_prefers_higher_score(tmp_path):
    p = tmp_path / "p.tsv"
    v = tmp_path / "v.tsv"
    p.write_text("seq_name\tplasmid_score\nx\t0.9\ny\t0.6\n")
    v.write_text("seq_name\tvirus_score\ttaxonomy\ny|provirus_1_100\t0.8\tViruses;Caudo\n")
    calls = sc.genomad_calls(str(p), str(v))
    assert calls["x"]["genomad_class"] == "plasmid"
    assert calls["y"]["genomad_class"] == "virus"


def test_categorise():
    base = {"selected": "yes", "length": "3000000"}
    assert (
        sc.categorise(
            {**base, "checkm2_completeness": "95", "checkm2_contamination": "1"}, 500000, 90, 5
        )
        == "chromosome"
    )
    assert (
        sc.categorise(
            {**base, "checkm2_completeness": "60", "checkm2_contamination": "1"}, 500000, 90, 5
        )
        == "large_circular"
    )
    assert (
        sc.categorise(
            {"selected": "yes", "length": "5000", "genomad_class": "plasmid"}, 500000, 90, 5
        )
        == "plasmid"
    )
    assert sc.categorise({"selected": "yes", "length": "5000"}, 500000, 90, 5) == "other_circular"
    assert sc.categorise({"selected": "no", "length": "5000"}, 500000, 90, 5) == "unvalidated"


def test_small_edges_rotation_invariant(make_seq):
    a = make_seq(8000)
    b = a[3000:] + a[:3000]  # same circle, different start
    c = make_seq(8000)
    edges = dr.small_edges({"a": a, "b": b, "c": c}, 20000, 99, 90)
    assert edges == [("a", "b")] or edges == [("b", "a")]


def test_cluster_representative_ranking():
    rows = [
        {
            "contig_id": "metaflye__1",
            "assembler": "metaflye",
            "status": "PENDING",
            "length": "2000000",
            "checkm2_completeness": "99",
            "checkm2_contamination": "0",
            "spanning_reads": "50",
            "n_suspicious_sites": "1",
        },
        {
            "contig_id": "metamdbg__1",
            "assembler": "metamdbg",
            "status": "PASS",
            "length": "2000000",
            "checkm2_completeness": "98",
            "checkm2_contamination": "0",
            "spanning_reads": "20",
            "n_suspicious_sites": "0",
        },
        {
            "contig_id": "hicanu__1",
            "assembler": "hicanu",
            "status": "PASS",
            "length": "2000000",
            "checkm2_completeness": "98",
            "checkm2_contamination": "0",
            "spanning_reads": "20",
            "n_suspicious_sites": "0",
        },
        {"contig_id": "hicanu__2", "assembler": "hicanu", "status": "PASS", "length": "5000"},
    ]
    edges = [("metaflye__1", "metamdbg__1"), ("hicanu__1", "metamdbg__1")]
    out = dr.cluster(rows, edges, ["metaflye", "metamdbg", "hifiasm_meta", "hicanu"])
    reps = {r["cluster_id"]: r["representative_id"] for r in out}
    assert reps == {"C0001": "metamdbg__1", "C0002": "hicanu__2"}
    assert {r["cluster_size"] for r in out if r["cluster_id"] == "C0001"} == {3}


def test_skani_edges(tmp_path):
    f = tmp_path / "s.tsv"
    write_tsv(
        [
            {
                "ANI": 99.9,
                "Align_fraction_ref": 98,
                "Align_fraction_query": 97,
                "Ref_name": "a",
                "Query_name": "b",
            },
            {
                "ANI": 99.9,
                "Align_fraction_ref": 50,
                "Align_fraction_query": 97,
                "Ref_name": "a",
                "Query_name": "c",
            },
            {
                "ANI": 95.0,
                "Align_fraction_ref": 98,
                "Align_fraction_query": 97,
                "Ref_name": "b",
                "Query_name": "c",
            },
        ],
        f,
        [
            "Ref_file",
            "Query_file",
            "ANI",
            "Align_fraction_ref",
            "Align_fraction_query",
            "Ref_name",
            "Query_name",
        ],
    )
    assert dr.skani_edges(str(f), 99, 90) == [("a", "b")]
    assert read_tsv(f)[0]["Ref_name"] == "a"
