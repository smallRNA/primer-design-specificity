"""Design and screen qPCR primer pairs against the human mitochondrial genome."""

import pandas as pd

from primers import (candidates, count_sites, pair_primers, passes_filters,
                     three_prime_clamp, max_self_complement)

GENOME = "data/MT-human.fa"
# MT-CO1 coding region in the revised Cambridge reference sequence.
TARGET = (5904, 7445)
REPORT_REGION = (6000, 6400)   # window searched for forward primers


def read_fasta(path):
    with open(path) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">")).upper()


def main():
    genome = read_fasta(GENOME)
    print(f"genome: {len(genome):,} bp   target region: {TARGET[0]}-{TARGET[1]}")

    cands = candidates(genome, REPORT_REGION)
    print(f"candidate oligos generated: {len(cands):,}")

    kept = [p for p in cands if passes_filters(p)]
    print(f"passing thermodynamic and structural filters: {len(kept):,} "
          f"({100*len(kept)/len(cands):.1f}%)")

    fwd = [p for p in kept if p.strand == "+"]
    rev = [p for p in kept if p.strand == "-"]
    pairs = pair_primers(fwd, rev)
    print(f"valid pairs (70-150 bp amplicon, dTm <= 1 C): {len(pairs):,}\n")

    rows = []
    for p in pairs:
        f, r = p["fwd"], p["rev"]
        f_sites = count_sites(genome, f.seq)
        r_sites = count_sites(genome, r.seq)
        rows.append(
            {
                "fwd": f.seq,
                "rev": r.seq,
                "fwd_tm": round(f.tm, 2),
                "rev_tm": round(r.tm, 2),
                "fwd_gc": round(f.gc, 3),
                "rev_gc": round(r.gc, 3),
                "amplicon_bp": p["amplicon"],
                "cross_dimer": p["cross_dimer"],
                "fwd_clamp": three_prime_clamp(f.seq),
                "rev_clamp": three_prime_clamp(r.seq),
                "fwd_hairpin": max_self_complement(f.seq),
                "fwd_sites": f_sites,
                "rev_sites": r_sites,
                "specific": f_sites == 1 and r_sites == 1,
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv("results/pairs.csv", index=False)

    n_spec = int(df.specific.sum())
    print(f"pairs with a single genomic binding site for both primers: "
          f"{n_spec:,} / {len(df):,} ({100*n_spec/len(df):.1f}%)")
    print(f"off-target load: fwd max {df.fwd_sites.max()} sites, "
          f"rev max {df.rev_sites.max()} sites\n")

    best = (
        df[df.specific & (df.cross_dimer <= 4)]
        .assign(tm_diff=lambda t: (t.fwd_tm - t.rev_tm).abs())
        .sort_values(["tm_diff", "amplicon_bp"])
        .head(5)
    )
    print("Top candidate pairs:")
    cols = ["fwd", "rev", "fwd_tm", "rev_tm", "amplicon_bp", "cross_dimer"]
    print(best[cols].to_string(index=False))
    best.to_csv("results/top_pairs.csv", index=False)

    # How much work each filter stage does.
    stages = pd.DataFrame(
        [
            {"stage": "candidate oligos", "n": len(cands)},
            {"stage": "pass Tm/GC/structure", "n": len(kept)},
            {"stage": "form valid pairs", "n": len(pairs)},
            {"stage": "single-site both primers", "n": n_spec},
            {"stage": "and no cross-dimer", "n": int(((df.specific) & (df.cross_dimer <= 4)).sum())},
        ]
    )
    stages.to_csv("results/funnel.csv", index=False)
    print("\nFilter funnel:")
    print(stages.to_string(index=False))

    # Does the specificity screen discriminate at all on this background?
    # At full stringency every primer is unique, which may say more about a
    # 16.5 kb genome than about the screen. Degrade stringency until it bites.
    print("\nStringency sweep (fraction of forward primers that are unique):")
    sample = [p.seq for p in fwd[:150]]
    sweep = []
    for seed in (12, 10, 8, 6, 5, 4):
        for mism in (0, 2, 4, 6):
            counts = [count_sites(genome, s, max_mismatch=mism, seed_len=seed)
                      for s in sample]
            sweep.append(
                {
                    "seed_len": seed,
                    "max_mismatch": mism,
                    "mean_sites": round(float(sum(counts) / len(counts)), 3),
                    "frac_unique": round(sum(c == 1 for c in counts) / len(counts), 3),
                }
            )
    sweep = pd.DataFrame(sweep)
    sweep.to_csv("results/stringency_sweep.csv", index=False)
    print(sweep.pivot(index="seed_len", columns="max_mismatch",
                      values="frac_unique").to_string())


if __name__ == "__main__":
    main()
