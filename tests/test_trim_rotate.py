import trim_rotate as tr


def test_exact_overlap(make_seq):
    core = make_seq(5000)
    seq = core + core[:300]
    assert tr.exact_overlap(seq, 1000) == 300
    assert tr.exact_overlap(core, 1000) < 20


def test_approx_overlap_with_mismatch(make_seq):
    core = make_seq(20000)
    dup = list(core[:2000])
    dup[1000] = "A" if dup[1000] != "A" else "C"  # 1 mismatch in 2 kb = 99.95 %
    seq = core + "".join(dup)
    ov = tr.find_overlap(seq, 50000, 100, 0.999)
    assert ov.method == "minimap2"
    assert abs(ov.length - 2000) <= 5


def test_no_overlap(make_seq):
    ov = tr.find_overlap(make_seq(20000), 50000, 100, 0.99)
    assert ov.length == 0 and ov.method == "none"


def test_rotation_moves_seam_inside(make_seq):
    seq = make_seq(1001)
    rotated, j = tr.rotate(seq)
    assert len(rotated) == len(seq)
    # the base before J is the old last base, the base at J the old first base
    assert rotated[j - 1] == seq[-1] and rotated[j] == seq[0]
    assert rotated[j:] + rotated[:j] == seq


def test_process_trims_then_rotates(make_seq):
    core = make_seq(4000)
    ov, trimmed, rotated, j = tr.process(core + core[:500], 50000, 100, 0.99)
    assert ov.length == 500 and trimmed == core
    assert rotated[j:] + rotated[:j] == core


def test_diverged_tandem_repeat_not_collapsed(make_seq):
    unit = list(make_seq(3600))
    copy2 = unit[:]
    for i in range(0, 3600, 300):  # ~99.7 % identity, like Buchnera pTrp
        copy2[i] = "A" if copy2[i] != "A" else "C"
    ov = tr.find_overlap("".join(unit + copy2), 50000, 100, 0.999)
    assert ov.length == 0


def test_hicanu_double_copy_circle_trimmed(make_seq):
    unit = make_seq(7786)
    seq = unit + unit[:7451]  # HiCanu-style ~2 copies of a small plasmid
    ov, trimmed, _rot, _j = tr.process(seq, 50000, 100, 0.999)
    assert ov.length == 7451 and trimmed == unit


def test_low_identity_overlap_not_trimmed(make_seq):
    core = list(make_seq(20000))
    dup = core[:2000]
    for i in range(0, 2000, 100):  # 1 % divergence
        dup[i] = "A" if dup[i] != "A" else "C"
    ov = tr.find_overlap("".join(core + dup), 50000, 100, 0.999)
    assert ov.length == 0
