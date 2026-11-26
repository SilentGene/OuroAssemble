import eval_junctions as ej

L, J = 100_000, 50_000
P = ej.Params(min_spanning=3, min_depth=5)


def full(name, start, end, mapq=60):
    return ej.Aln(name, start, end, mapq, [(0, end - start)])


def tiling(skip_junction=False, step=2000, readlen=15000):
    """reads tiling the whole contig (overhanging reads are clipped at the ends)"""
    alns = []
    for i, s in enumerate(range(-readlen + step, L, step)):
        start, end = max(0, s), min(L, s + readlen)
        if skip_junction and start < J < end:
            continue
        alns.append(full(f"t{i}", start, end))
    return alns


def test_pass_when_junction_spanned():
    alns = [full(f"r{i}", J - 8000 + i * 200, J + 7000 + i * 200) for i in range(10)]
    alns += tiling()
    res = ej.evaluate(alns, L, J, P)
    assert res.metrics["status"] == "PASS", res.metrics
    assert res.metrics["spanning_reads"] >= 10


def test_fail_when_reads_clipped_at_junction():
    alns = tiling(skip_junction=True)
    for i in range(8):  # reads end at J with a long clip
        alns.append(ej.Aln(f"l{i}", J - 9000 + i * 100, J, 60, [(0, 9000 - i * 100), (4, 5000)]))
        alns.append(ej.Aln(f"r{i}", J, J + 9000 - i * 100, 60, [(4, 5000), (0, 9000 - i * 100)]))
    res = ej.evaluate(alns, L, J, P)
    assert res.metrics["status"] == "FAIL", res.metrics
    assert res.metrics["clipped_at_junction"] >= 8


def test_pending_low_coverage():
    alns = [full("a", J - 5000, J + 5000), full("b", J - 4000, J + 6000)]
    res = ej.evaluate(alns, L, J, P)
    assert res.metrics["status"] == "PENDING"
    assert res.metrics["reason"] == "low_coverage_at_junction"


def test_pending_repeat_low_mapq():
    alns = tiling(skip_junction=True)
    alns += [full(f"m{i}", J - 6000 + i * 100, J + 6000 + i * 100, mapq=0) for i in range(10)]
    res = ej.evaluate(alns, L, J, P)
    assert res.metrics["status"] == "PENDING"
    assert res.metrics["reason"] == "repeat_at_junction"


def test_internal_clip_pileup_makes_pending():
    alns = [full(f"r{i}", J - 8000 + i * 200, J + 7000 + i * 200) for i in range(10)]
    alns += tiling()
    site = 20_000
    alns += [ej.Aln(f"c{i}", site - 6000, site, 60, [(0, 6000), (4, 3000)]) for i in range(6)]
    res = ej.evaluate(alns, L, J, P)
    assert res.metrics["n_clip_sites"] == 1
    assert res.metrics["status"] == "PENDING"
    assert res.sites[0][0] == site


def test_clips_at_contig_ends_ignored():
    aln = ej.Aln("x", 0, 5000, 60, [(4, 3000), (0, 5000)])
    assert ej.clip_events(aln, L, 200) == []
    aln = ej.Aln("y", L - 5000, L, 60, [(0, 5000), (4, 3000)])
    assert ej.clip_events(aln, L, 200) == []
