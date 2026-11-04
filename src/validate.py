"""
Is the specificity screen calibrated?

Every forward primer comes back unique against the 16.5 kb mitochondrial
genome. That is not evidence the screen works - a 20-mer has 4^20 possible
sequences, so uniqueness in a genome this small is close to guaranteed either
way. A test that cannot fail reports nothing.

So the screen is checked against a background where the right answer is known
in advance. Random sequence is generated with the same base composition as the
real genome, at sizes from 10 kb to 10 Mb, and the observed off-target counts
are compared with the number expected by chance.

For a seed of length k requiring an exact 3' match, scanned on both strands of
a background of length N, the expected number of seed hits is

    E = 2 * N * prod(p_base) over the seed

Observed counts tracking E means the screen is finding what is there and no
more. Observed well above E would mean it is over-calling.
"""

import numpy as np
import pandas as pd

from primers import revcomp
from design import read_fasta

SEED = 0
SIZES = [10_000, 100_000, 1_000_000, 10_000_000]
SEED_LENS = [8, 10, 12]
N_PRIMERS = 40


def base_freqs(seq):
    return np.array([seq.count(b) for b in "ACGT"], dtype=float) / len(seq)


def random_background(n, freqs, rng):
    return "".join(rng.choice(list("ACGT"), size=n, p=freqs))


def expected_hits(seed, n, freqs):
    """Expected exact seed occurrences in n bases of random sequence, both strands."""
    p = {b: f for b, f in zip("ACGT", freqs)}
    prob = np.prod([p[b] for b in seed])
    return 2 * n * prob


def count_seed_hits(background, seed):
    """Exact occurrences of the seed and its reverse complement."""
    total = 0
    for probe in (seed, revcomp(seed)):
        pos = background.find(probe)
        while pos != -1:
            total += 1
            pos = background.find(probe, pos + 1)
    return total


def main():
    rng = np.random.default_rng(SEED)
    genome = read_fasta("data/MT-human.fa")
    freqs = base_freqs(genome)
    print("genome base composition A/C/G/T: " +
          ", ".join(f"{f:.3f}" for f in freqs))

    pairs = pd.read_csv("results/pairs.csv")
    primers = pairs["fwd"].drop_duplicates().head(N_PRIMERS).tolist()
    print(f"primers tested: {len(primers)}\n")

    rows = []
    for n in SIZES:
        background = random_background(n, freqs, rng)
        for k in SEED_LENS:
            obs, exp = [], []
            for p in primers:
                seed = p[-k:]
                obs.append(count_seed_hits(background, seed))
                exp.append(expected_hits(seed, n, freqs))
            rows.append(
                {
                    "background_bp": n,
                    "seed_len": k,
                    "mean_observed": float(np.mean(obs)),
                    "mean_expected": float(np.mean(exp)),
                    "frac_unique": float(np.mean([o <= 1 for o in obs])),
                }
            )
            print(f"  {n:>10,} bp  seed {k:>2}  observed {np.mean(obs):8.2f}  "
                  f"expected {np.mean(exp):8.2f}")

    df = pd.DataFrame(rows)
    df["ratio"] = df.mean_observed / df.mean_expected
    df.to_csv("results/calibration.csv", index=False)

    print("\nObserved / expected:")
    print(df.pivot(index="background_bp", columns="seed_len",
                   values="ratio").round(3).to_string())
    print("\nFraction of primers still unique:")
    print(df.pivot(index="background_bp", columns="seed_len",
                   values="frac_unique").round(3).to_string())


if __name__ == "__main__":
    main()
