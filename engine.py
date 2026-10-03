"""
NEO-CRISPR bioinformatics engine.

This module holds the scientific core of the tool, separated from the web
server so it can be tested and reused independently.

What it does:

* Cleans pasted input (strips FASTA headers, whitespace, lowercases, keeps
  only A/C/G/T, marks anything else as N).
* Finds SpCas9 target sites (N20 + NGG) on BOTH strands of the input.
* Scans the WHOLE genome (every scaffold) on BOTH strands for off-targets,
  tolerating NGG and the weaker NAG off-target PAM.
* Encodes each genome exactly once and caches it (keyed by path + mtime), then
  scans every candidate guide against it in a single batched tensor pass on the
  GPU if one is available, otherwise the CPU.
* Excludes each guide's own on-target site from its off-target list.
* Reports, for every off-target, its locus, strand, aligned sequence, total
  mismatches, seed-region (PAM-proximal) mismatches, PAM, and -- when a GTF is
  present -- the gene/feature it falls in.
* Computes a transparent HEURISTIC guide score (not a validated efficiency
  predictor; see ``heuristic_score``).
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import torch

# --- constants ---------------------------------------------------------------

SPACER_LEN = 20
PAM_LEN = 3
TARGET_LEN = SPACER_LEN + PAM_LEN  # 23
SEED_LEN = 10  # PAM-proximal nucleotides that matter most for specificity
DEFAULT_MAX_MISMATCH = 3

# A=0 C=1 G=2 T=3; anything else (incl. N) -> 4 so it never matches a real base.
BASE_MAP = {"A": 0, "C": 1, "G": 2, "T": 3}
UNKNOWN_BASE = 4
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

# Scan the genome in blocks of this many candidate PAM windows at a time, so
# peak memory stays bounded no matter how large the genome or how many guides.
_WINDOW_CHUNK = 200_000

_DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def device_name() -> str:
    return "GPU (CUDA)" if _DEVICE.type == "cuda" else "CPU"


# --- sequence helpers --------------------------------------------------------

def clean_sequence(raw: str) -> str:
    """Normalise pasted text into a bare nucleotide string.

    Removes FASTA header lines (``>...``), all whitespace, and uppercases.
    Any character that is not A/C/G/T becomes ``N`` so it is preserved in the
    coordinate map but never counts as a match.
    """
    if not raw:
        return ""
    lines = [ln for ln in raw.splitlines() if not ln.lstrip().startswith(">")]
    seq = "".join(lines)
    seq = re.sub(r"\s+", "", seq).upper()
    return "".join(c if c in BASE_MAP else "N" for c in seq)


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def encode(seq: str) -> torch.Tensor:
    """Encode a nucleotide string to an int8 tensor (A0 C1 G2 T3, else 4)."""
    return torch.tensor(
        [BASE_MAP.get(b, UNKNOWN_BASE) for b in seq],
        dtype=torch.int8,
        device=_DEVICE,
    )


# --- candidate guide discovery -----------------------------------------------

@dataclass
class GuideSite:
    """A candidate SpCas9 target on the input sequence.

    ``start``/``end`` are always expressed in forward-strand coordinates of the
    cleaned input, so the UI can highlight the binding region regardless of
    which strand the guide targets.
    """

    spacer: str          # 20 nt protospacer, 5'->3'
    target: str          # 23 nt protospacer + PAM, 5'->3'
    pam: str             # 3 nt PAM
    strand: str          # '+' or '-'
    start: int           # forward-coord start of the 23-mer (0-based)
    end: int             # forward-coord end (exclusive)


def find_guide_sites(seq: str) -> List[GuideSite]:
    """Find all N20-NGG target sites on both strands of ``seq``.

    Forward strand: match ``[ACGT]{20}[ACGT]GG`` directly.
    Reverse strand: a guide on the minus strand appears on the forward strand
    as ``CCN`` followed by 20 nt; its spacer is the reverse complement of that
    downstream 20 nt. We report forward coordinates for highlighting.
    """
    sites: List[GuideSite] = []

    fwd = re.compile(r"(?=([ACGT]{20}[ACGT]GG))")
    for m in fwd.finditer(seq):
        target = m.group(1)
        start = m.start()
        sites.append(
            GuideSite(
                spacer=target[:SPACER_LEN],
                target=target,
                pam=target[SPACER_LEN:],
                strand="+",
                start=start,
                end=start + TARGET_LEN,
            )
        )

    # Reverse strand: CC[ACGT] then 20 nt on the forward sequence.
    rev = re.compile(r"(?=(CC[ACGT][ACGT]{20}))")
    for m in rev.finditer(seq):
        chunk = m.group(1)              # 23 nt on forward strand
        start = m.start()
        target = reverse_complement(chunk)  # 5'->3' on the minus strand
        sites.append(
            GuideSite(
                spacer=target[:SPACER_LEN],
                target=target,
                pam=target[SPACER_LEN:],
                strand="-",
                start=start,
                end=start + TARGET_LEN,
            )
        )

    return sites


def heuristic_score(spacer: str) -> int:
    """A transparent rule-of-thumb score in [0, 100].

    This is NOT a validated on-target efficiency predictor (e.g. Rule Set 2 /
    DeepSpCas9). It only rewards a mid-range GC content and a 3' G, which are
    weak positive signals. Treat it as a tie-breaker, not a ground truth.
    """
    if not spacer:
        return 0
    gc = (spacer.count("G") + spacer.count("C")) / len(spacer) * 100
    score = 50
    if 40 <= gc <= 60:
        score += 25
    if spacer[-1] == "G":
        score += 15
    return max(0, min(100, score))


def gc_percent(spacer: str) -> int:
    if not spacer:
        return 0
    return round((spacer.count("G") + spacer.count("C")) / len(spacer) * 100)


# --- genome caching ----------------------------------------------------------

@dataclass
class _EncodedGenome:
    """Both strands of every scaffold, encoded once and kept in memory."""

    mtime: float
    # chrom name -> (forward tensor, reverse-complement tensor)
    chroms: Dict[str, Tuple[torch.Tensor, torch.Tensor]] = field(default_factory=dict)


_GENOME_CACHE: Dict[str, _EncodedGenome] = {}


def _encode_rc_tensor(fwd: torch.Tensor) -> torch.Tensor:
    """Reverse-complement an encoded tensor. A/C/G/T (0..3) map to 3-x; 4 stays."""
    rc = torch.flip(fwd, dims=[0]).clone()
    known = rc != UNKNOWN_BASE
    rc[known] = 3 - rc[known]
    return rc


def load_genome(fasta_path: str) -> _EncodedGenome:
    """Load + encode every scaffold of a FASTA, caching by path and mtime."""
    from pyfaidx import Fasta

    mtime = os.path.getmtime(fasta_path)
    cached = _GENOME_CACHE.get(fasta_path)
    if cached and cached.mtime == mtime:
        return cached

    genome = Fasta(fasta_path)
    enc = _EncodedGenome(mtime=mtime)
    for name in genome.keys():
        seq = str(genome[name][:]).upper()
        fwd = encode(seq)
        enc.chroms[name] = (fwd, _encode_rc_tensor(fwd))
    _GENOME_CACHE[fasta_path] = enc
    return enc


def clear_cache() -> None:
    _GENOME_CACHE.clear()


# --- off-target scanning -----------------------------------------------------

@dataclass
class OffTarget:
    loc: str
    chrom: str
    pos: int          # 1-based forward-strand coordinate of the 23-mer start
    strand: str
    seq: str          # 23-mer as found (protospacer + PAM, 5'->3' on its strand)
    mm: int           # total mismatches in the 20 nt protospacer
    seed_mm: int      # mismatches within the PAM-proximal seed region
    pam: str
    gene: Optional[str] = None


def _pam_mask(windows: torch.Tensor) -> torch.Tensor:
    """Boolean mask of windows whose PAM is NGG or NAG (position 21==G, 22 in {A,G})."""
    p1 = windows[:, SPACER_LEN + 1]   # the 'G' of NGG / 'A' of NAG
    p2 = windows[:, SPACER_LEN + 2]   # final G
    return (p2 == BASE_MAP["G"]) & ((p1 == BASE_MAP["G"]) | (p1 == BASE_MAP["A"]))


def _decode(window: torch.Tensor) -> str:
    inv = ("A", "C", "G", "T", "N")
    return "".join(inv[int(b)] for b in window)


def _scan_strand(
    chrom: str,
    strand: str,
    chrom_tensor: torch.Tensor,
    spacers: torch.Tensor,     # [S, 20] int8
    max_mismatch: int,
    chrom_len: int,
    per_guide: List[List[OffTarget]],
) -> None:
    """Scan one strand of one chromosome for all spacers at once, chunked."""
    if chrom_tensor.numel() < TARGET_LEN:
        return

    seed_cols = torch.arange(SPACER_LEN - SEED_LEN, SPACER_LEN, device=_DEVICE)
    windows_all = chrom_tensor.unfold(0, TARGET_LEN, 1)  # [W, 23], a view

    n_windows = windows_all.shape[0]
    for begin in range(0, n_windows, _WINDOW_CHUNK):
        windows = windows_all[begin : begin + _WINDOW_CHUNK]  # [w, 23]
        pam_ok = _pam_mask(windows)
        if not bool(pam_ok.any()):
            continue
        idx = torch.nonzero(pam_ok, as_tuple=False).squeeze(1)  # window offsets
        protospacers = windows[idx, :SPACER_LEN].to(torch.int16)  # [p, 20]

        # [S, p, 20] mismatch comparison, summed over the 20 nt.
        diff = protospacers.unsqueeze(0) != spacers.unsqueeze(1).to(torch.int16)
        mm = diff.sum(dim=2)                       # [S, p]
        seed_mm = diff[:, :, seed_cols].sum(dim=2)  # [S, p]

        hit_s, hit_p = torch.where(mm <= max_mismatch)
        for s, p in zip(hit_s.tolist(), hit_p.tolist()):
            win_offset = int(idx[p])
            fwd_pos0 = _to_forward_pos(strand, begin + win_offset, chrom_len)
            window = windows[win_offset]
            per_guide[s].append(
                OffTarget(
                    loc=f"{chrom}:{fwd_pos0 + 1}({strand})",
                    chrom=chrom,
                    pos=fwd_pos0 + 1,
                    strand=strand,
                    seq=_decode(window),
                    mm=int(mm[s, p]),
                    seed_mm=int(seed_mm[s, p]),
                    pam=_decode(window[SPACER_LEN:]),
                )
            )


def _to_forward_pos(strand: str, window_start: int, chrom_len: int) -> int:
    """Map a 0-based window start on the scanned strand to forward coordinates."""
    if strand == "+":
        return window_start
    # On the reverse-complement tensor, index w corresponds to a 23-mer whose
    # forward-strand start is chrom_len - (w + TARGET_LEN).
    return chrom_len - (window_start + TARGET_LEN)


def scan_off_targets(
    spacers: List[str],
    fasta_path: str,
    max_mismatch: int = DEFAULT_MAX_MISMATCH,
) -> List[List[OffTarget]]:
    """Return, for each input spacer, its off-target hits across the whole genome.

    The intended on-target site (exactly one perfect, zero-mismatch hit) is
    removed from each list; any *additional* perfect matches are retained,
    because duplicated perfect sites are genuine specificity problems.
    """
    results: List[List[OffTarget]] = [[] for _ in spacers]
    if not spacers:
        return results

    genome = load_genome(fasta_path)
    spacer_tensor = torch.stack([encode(s[:SPACER_LEN]) for s in spacers])  # [S,20]

    for chrom, (fwd, rc) in genome.chroms.items():
        chrom_len = fwd.shape[0]
        _scan_strand(chrom, "+", fwd, spacer_tensor, max_mismatch, chrom_len, results)
        _scan_strand(chrom, "-", rc, spacer_tensor, max_mismatch, chrom_len, results)

    cleaned: List[List[OffTarget]] = []
    for hits in results:
        hits.sort(key=lambda h: (h.mm, h.seed_mm, h.chrom, h.pos))
        # Drop one perfect match: the guide's own intended cut site.
        for i, h in enumerate(hits):
            if h.mm == 0:
                hits.pop(i)
                break
        cleaned.append(hits)
    return cleaned


# --- GTF annotation ----------------------------------------------------------

_GTF_CACHE: Dict[str, Tuple[float, Dict[str, List[Tuple[int, int, str]]]]] = {}


def load_gtf(gtf_path: str) -> Dict[str, List[Tuple[int, int, str]]]:
    """Parse a GTF into per-chromosome sorted (start, end, gene) intervals.

    Cached by path + mtime. Uses pandas for the heavy parse. Only ``gene``
    features are indexed (falling back to ``exon`` if no gene features exist),
    which is enough to say which gene an off-target lands in.
    """
    import pandas as pd

    mtime = os.path.getmtime(gtf_path)
    cached = _GTF_CACHE.get(gtf_path)
    if cached and cached[0] == mtime:
        return cached[1]

    cols = ["seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute"]
    df = pd.read_csv(
        gtf_path, sep="\t", comment="#", header=None, names=cols,
        dtype={"seqname": str}, engine="c",
    )

    feature = "gene" if (df["feature"] == "gene").any() else "exon"
    df = df[df["feature"] == feature]

    name_re = re.compile(r'gene_name "([^"]+)"')
    id_re = re.compile(r'gene_id "([^"]+)"')

    def gene_of(attr: str) -> str:
        m = name_re.search(attr) or id_re.search(attr)
        return m.group(1) if m else "?"

    index: Dict[str, List[Tuple[int, int, str]]] = {}
    for seqname, start, end, attr in zip(df["seqname"], df["start"], df["end"], df["attribute"]):
        index.setdefault(str(seqname), []).append((int(start), int(end), gene_of(attr)))
    for seqname in index:
        index[seqname].sort()

    _GTF_CACHE[gtf_path] = (mtime, index)
    return index


def annotate(index: Dict[str, List[Tuple[int, int, str]]], chrom: str, pos: int) -> Optional[str]:
    """Return the gene overlapping ``pos`` (1-based) on ``chrom``, if any."""
    intervals = index.get(chrom)
    if not intervals:
        return None
    for start, end, gene in intervals:
        if start <= pos <= end:
            return gene
        if start > pos:
            break
    return None
