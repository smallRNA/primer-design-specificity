"""Filter funnel and specificity calibration."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

funnel = pd.read_csv("results/funnel.csv")
cal = pd.read_csv("results/calibration.csv")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))

# --- funnel --------------------------------------------------------------
ax1.barh(funnel.stage, funnel.n, color="#1f3a5f")
ax1.set_xscale("log")
ax1.invert_yaxis()
ax1.set_xlabel("count (log scale)")
ax1.set_title("Design funnel\n6,400 oligos to 1,104 usable pairs", fontsize=10)
for i, v in enumerate(funnel.n):
    ax1.text(v * 1.15, i, f"{v:,}", va="center", fontsize=8.5)
ax1.set_xlim(right=funnel.n.max() * 4)
ax1.spines[["top", "right"]].set_visible(False)

# --- calibration ---------------------------------------------------------
markers = {8: "o", 10: "s", 12: "^"}
for k, grp in cal.groupby("seed_len"):
    ax2.plot(grp.background_bp, grp.frac_unique, marker=markers[k], lw=1.4,
             label=f"{k}-base 3' seed")
ax2.axvline(16569, ls="--", lw=1, color="#c0392b")
ax2.text(17500, 0.45, "mtDNA\n16.5 kb", color="#c0392b", fontsize=8)
ax2.set_xscale("log")
ax2.set_xlabel("background genome size (bp)")
ax2.set_ylabel("fraction of primers still unique")
ax2.set_title("Specificity vs background size\nthe screen only bites on large genomes",
              fontsize=10)
ax2.legend(fontsize=8, frameon=False)
ax2.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
fig.savefig("results/primers.png", dpi=160)
print("wrote results/primers.png")
