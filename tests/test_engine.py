"""Tests for the NEO-CRISPR engine, using small synthetic genomes.

Run with:  pytest
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import engine  # noqa: E402

SPACER = "ACGTACGTACGTCCGGTTAA"  # 20 nt, GC = 50%


# --- sequence helpers --------------------------------------------------------

def test_clean_strips_header_whitespace_and_marks_unknown():
    raw = ">seq1 desc\nacgt ACGT\n\tNNxy"
    assert engine.clean_sequence(raw) == "ACGTACGTNNNN"


def test_clean_empty():
    assert engine.clean_sequence("") == ""
    assert engine.clean_sequence(">only header\n") == ""


def test_reverse_complement():
    assert engine.reverse_complement("ACGT") == "ACGT"
    assert engine.reverse_complement("AAAC") == "GTTT"
    assert engine.reverse_complement("ACGTN") == "NACGT"


def test_gc_and_heuristic_bounds():
    assert engine.gc_percent("GGGGCCCCGGGGCCCCGGGG") == 100
    assert engine.gc_percent("AAAAAAAAAAAAAAAAAAAA") == 0
    assert engine.heuristic_score("") == 0
    assert 0 <= engine.heuristic_score(SPACER) <= 100
    # 3' G and mid-range GC both rewarded.
    assert engine.heuristic_score("ACGTACGTACGTACGTACGG") == 90


# --- guide discovery on both strands ----------------------------------------

def test_find_sites_forward():
    seq = SPACER + "TGG"  # N20 + NGG
    sites = engine.find_guide_sites(seq)
    plus = [s for s in sites if s.strand == "+"]
    assert len(plus) == 1
    assert plus[0].spacer == SPACER
    assert plus[0].pam == "TGG"
    assert (plus[0].start, plus[0].end) == (0, 23)


def test_find_sites_reverse_strand():
    # A minus-strand guide appears as CCN + 20 nt on the forward strand.
    target = SPACER + "AGG"
    forward = engine.reverse_complement(target)  # starts with CCT...
    sites = engine.find_guide_sites(forward)
    minus = [s for s in sites if s.strand == "-"]
    assert len(minus) == 1
    assert minus[0].spacer == SPACER
    assert minus[0].pam == "AGG"
    # Highlight coordinates are always forward-strand.
    assert (minus[0].start, minus[0].end) == (0, 23)


# --- off-target scanning -----------------------------------------------------

@pytest.fixture
def genome(tmp_path):
    def _make(fasta_text, gtf_text=None):
        d = tmp_path / "sp"
        d.mkdir(exist_ok=True)
        fa = d / "genome.fa"
        fa.write_text(fasta_text)
        gtf = None
        if gtf_text is not None:
            gtf = d / "annotation.gtf"
            gtf.write_text(gtf_text)
        engine.clear_cache()
        return str(fa), (str(gtf) if gtf else None)
    return _make


def test_on_target_excluded_single_perfect(genome):
    # One perfect site in the genome == the on-target => excluded => 0 offs.
    fa, _ = genome(">c\n" + "GGGG" + SPACER + "TGG" + "GGGG\n")
    offs = engine.scan_off_targets([SPACER], fa)[0]
    assert offs == []


def test_duplicate_perfect_site_is_offtarget(genome):
    # Two perfect sites => one is the on-target, the other is a real off-target.
    block = SPACER + "TGG"
    fa, _ = genome(">c\n" + block + "AAAA" + block + "\n")
    offs = engine.scan_off_targets([SPACER], fa)[0]
    assert len(offs) == 1
    assert offs[0].mm == 0


def test_mismatch_offtarget_detected_with_seed_count(genome):
    on = SPACER + "TGG"
    # 2 mismatches in the PAM-proximal seed, NGG PAM.
    off = "ACGTACGTACGTCCGGTTTT" + "AGG"
    fa, _ = genome(">c\n" + on + "CCCC" + off + "\n")
    offs = engine.scan_off_targets([SPACER], fa)[0]
    assert len(offs) == 1
    assert offs[0].mm == 2
    assert offs[0].seed_mm == 2
    assert offs[0].pam == "AGG"
    assert offs[0].strand == "+"


def test_nag_offtarget_detected(genome):
    on = SPACER + "TGG"
    off = SPACER + "TAG"  # perfect protospacer but weak NAG PAM elsewhere
    fa, _ = genome(">c\n" + on + "TTTT" + off + "\n")
    offs = engine.scan_off_targets([SPACER], fa)[0]
    pams = sorted(o.pam for o in offs)
    assert "TAG" in pams  # the NAG site is reported


def test_both_strands_scanned(genome):
    on = SPACER + "TGG"
    # Same protospacer present on the minus strand (store its revcomp on fwd).
    minus_copy = engine.reverse_complement(SPACER + "CGG")
    fa, _ = genome(">c\n" + on + "TTTT" + minus_copy + "TTTT\n")
    offs = engine.scan_off_targets([SPACER], fa)[0]
    assert any(o.strand == "-" and o.mm == 0 for o in offs)


def test_mismatch_cutoff_respected(genome):
    on = SPACER + "TGG"
    far = "TTTTTTTTTTTTTTTTTTTT" + "TGG"  # 16 mismatches
    fa, _ = genome(">c\n" + on + "CCCC" + far + "\n")
    offs = engine.scan_off_targets([SPACER], fa, max_mismatch=3)[0]
    assert offs == []


def test_multiple_scaffolds(genome):
    on = SPACER + "TGG"
    off = "ACGTACGTACGTCCGGTTTT" + "TGG"  # 2 mm, on scaffold 2
    fa, _ = genome(">c1\n" + on + "\n>c2\nGGGG" + off + "GGGG\n")
    offs = engine.scan_off_targets([SPACER], fa)[0]
    assert len(offs) == 1
    assert offs[0].chrom == "c2"


# --- GTF annotation ----------------------------------------------------------

def test_gtf_annotation(genome):
    on = SPACER + "TGG"
    off = "ACGTACGTACGTCCGGTTTT" + "TGG"
    fa, gtf = genome(
        ">c\n" + on + "CCCC" + off + "\n",
        'c\tsrc\tgene\t1\t200\t.\t+\t.\tgene_name "MYGENE";\n',
    )
    offs = engine.scan_off_targets([SPACER], fa)[0]
    index = engine.load_gtf(gtf)
    assert len(offs) == 1
    assert engine.annotate(index, offs[0].chrom, offs[0].pos) == "MYGENE"
    # A position outside any gene returns None.
    assert engine.annotate(index, "c", 10_000) is None


def test_gtf_falls_back_to_exon(genome):
    _, gtf = genome(">c\nACGT\n", 'c\tsrc\texon\t1\t50\t.\t+\t.\tgene_id "EX1";\n')
    index = engine.load_gtf(gtf)
    assert engine.annotate(index, "c", 10) == "EX1"
