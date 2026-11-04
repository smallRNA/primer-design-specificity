"""
qPCR primer design with in silico specificity screening.

Generates candidate primer pairs across a target region, filters them on the
thermodynamic and structural criteria that matter at the bench, then screens
every survivor against the full background genome for off-target priming.

Melting temperatures use SantaLucia (1998) unified nearest-neighbour
thermodynamics with salt correction, not the 2(A+T)+4(G+C) rule of thumb,
which is unusable above ~20 nt.
"""

from dataclasses import dataclass
from itertools import product

import numpy as np

# SantaLucia 1998 unified parameters: dH (kcal/mol), dS (cal/mol/K)
NN = {
    "AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "AC": (-8.4, -22.4), "AG": (-7.8, -21.0),
    "TA": (-7.2, -21.3), "TT": (-7.9, -22.2), "TC": (-8.2, -22.2), "TG": (-8.5, -22.7),
    "CA": (-8.5, -22.7), "CT": (-7.8, -21.0), "CC": (-8.0, -19.9), "CG": (-10.6, -27.2),
    "GA": (-8.2, -22.2), "GT": (-8.4, -22.4), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9),
}
INIT = {"A": (2.3, 4.1), "T": (2.3, 4.1), "C": (0.1, -2.8), "G": (0.1, -2.8)}
COMP = str.maketrans("ACGT", "TGCA")

R = 1.987  # cal/mol/K


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)


def tm(seq, primer_nM=250.0, na_M=0.05):
    """Nearest-neighbour Tm with Owczarzy-style salt correction, in C."""
    dh, ds = INIT[seq[0]]
    dh2, ds2 = INIT[seq[-1]]
    dh += dh2
    ds += ds2
    for i in range(len(seq) - 1):
        d_h, d_s = NN[seq[i:i + 2]]
        dh += d_h
        ds += d_s

    dh *= 1000.0                      # kcal -> cal
    ct = primer_nM * 1e-9 / 4.0       # non-self-complementary duplex
    tm_k = dh / (ds + R * np.log(ct))
    # salt correction on 1/Tm
    tm_k = 1 / (1 / tm_k + (4.29 * gc(seq) - 3.95) * 1e-5 * np.log(na_M)
                + 9.40e-6 * np.log(na_M) ** 2)
    return tm_k - 273.15


def max_self_complement(seq, min_run=4):
    """Longest self-complementary run: a proxy for hairpin / self-dimer risk."""
    rc = revcomp(seq)
    best = 0
    for size in range(min_run, len(seq) + 1):
        hit = any(seq[i:i + size] in rc for i in range(len(seq) - size + 1))
        if hit:
            best = size
        else:
            break
    return best


def three_prime_clamp(seq):
    """G or C count in the last 5 bases. 1-3 is the usual target.

    Too few and the 3' end breathes; too many and it anchors to mismatched
    sites, which is a specificity problem rather than a stability one.
    """
    return sum(b in "GC" for b in seq[-5:])


def ends_in_gc_run(seq):
    return seq[-3:] in ("GGG", "CCC") or seq[-1] not in "ACGT"


@dataclass
class Primer:
    seq: str
    start: int
    strand: str

    @property
    def tm(self):
        return tm(self.seq)

    @property
    def gc(self):
        return gc(self.seq)


def candidates(genome, region, lengths=range(18, 26)):
    """Every subsequence of allowed length starting inside the target region."""
    lo, hi = region
    out = []
    for start in range(lo, hi):
        for n in lengths:
            sub = genome[start:start + n]
            if len(sub) == n and set(sub) <= set("ACGT"):
                out.append(Primer(sub, start, "+"))
                out.append(Primer(revcomp(sub), start, "-"))
    return out


def passes_filters(p, tm_range=(58, 62), gc_range=(0.40, 0.60), max_hairpin=4):
    return (
        tm_range[0] <= p.tm <= tm_range[1]
        and gc_range[0] <= p.gc <= gc_range[1]
        and 1 <= three_prime_clamp(p.seq) <= 3
        and not ends_in_gc_run(p.seq)
        and max_self_complement(p.seq) <= max_hairpin
        and "GGGG" not in p.seq
        and "CCCC" not in p.seq
        and "AAAA" not in p.seq
        and "TTTT" not in p.seq
    )


def count_sites(genome, primer_seq, max_mismatch=2, seed_len=8):
    """Count genome positions where this primer could prime.

    A primer only extends if its 3' end is matched, so the last `seed_len`
    bases are required to match exactly and mismatches are allowed only
    further upstream. Counting full-length approximate matches instead would
    badly overstate the off-target load.
    """
    rc = revcomp(primer_seq)
    n = len(primer_seq)
    hits = 0

    for probe, _ in ((primer_seq, "+"), (rc, "-")):
        seed = probe[-seed_len:] if probe is primer_seq else probe[:seed_len]
        pos = genome.find(seed)
        while pos != -1:
            if probe is primer_seq:
                s = pos + seed_len - n
            else:
                s = pos
            if 0 <= s <= len(genome) - n:
                window = genome[s:s + n]
                mism = sum(a != b for a, b in zip(window, probe))
                if mism <= max_mismatch:
                    hits += 1
            pos = genome.find(seed, pos + 1)
    return hits


def pair_primers(fwd, rev, amplicon=(70, 150), max_tm_diff=1.0):
    pairs = []
    for f, r in product(fwd, rev):
        size = r.start + len(r.seq) - f.start
        if amplicon[0] <= size <= amplicon[1] and abs(f.tm - r.tm) <= max_tm_diff:
            cross = max_self_complement(f.seq + r.seq)
            pairs.append({"fwd": f, "rev": r, "amplicon": size, "cross_dimer": cross})
    return pairs
