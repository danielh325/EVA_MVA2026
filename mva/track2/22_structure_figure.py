"""Figure for the structural context of p.Asn1002Lys, from script 21's JSON.

Two panels, because the result has two halves and one of them is a failure:

  left   burial of the measured alleles and the query. The honest message is
         that the control FAILS - Y155C is tolerated and almost completely
         buried - so burial alone cannot call this variant.
  right  how far the query sits from the positions whose behaviour is known.
         It is 22-35 A from the structural pocket, so the K795A precedent
         does not transfer, and nearest to L1012, which is measured
         destabilising.

Reads 21_structure_stability_results.json. Run: python3 -I 22_structure_figure.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "../../reports/figures"
SRC = "21_structure_stability_results.json"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#a3a29c"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
GRID = "#e6e5e1"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 10,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": GRID,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
})

d = json.load(open(SRC))
rows = d["rows"]
by = {r["pos"]: r for r in rows}
COLOUR = {"destabilised": ORANGE, "tolerated": AQUA, "UNKNOWN": BLUE}

fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(12.6, 5.5), gridspec_kw={"width_ratios": [1.15, 1]})

# ---------------------------------------------------------------- left
order = [795, 909, 1012, 155, 550, 1002]
labels, vals, cols = [], [], []
for p in order:
    r = by[p]
    tag = {"destabilised": "measured\ndestabilising",
           "tolerated": "measured\ntolerated",
           "UNKNOWN": "never assayed"}[r["class"]]
    labels.append(f"{r['variant']}\n{tag}")
    vals.append(r["rsa_percent"])
    cols.append(COLOUR[r["class"]])

bars = ax1.bar(range(len(order)), vals, color=cols, width=.62,
               edgecolor=SURFACE, linewidth=2)
ax1.axhline(d["domain_median_rsa"], color=MUTED, lw=1.2, ls=(0, (4, 3)),
            zorder=0)
ax1.text(-.42, d["domain_median_rsa"] + 2.0,
         f"domain median {d['domain_median_rsa']:.1f}%", ha="left",
         fontsize=8.6, color=INK2)

for i, (b, v) in enumerate(zip(bars, vals)):
    ax1.text(b.get_x() + b.get_width() / 2, v + 1.6, f"{v:.1f}%",
             ha="center", fontsize=9.6, fontweight="bold", color=INK)

ax1.set_xticks(range(len(order)))
ax1.set_xticklabels(labels, fontsize=8.4, color=INK2, linespacing=1.45)
ax1.set_ylabel("relative solvent accessibility  (buried → exposed)")
ax1.set_ylim(0, 96)
ax1.grid(axis="y", color=GRID, lw=.8)
ax1.set_axisbelow(True)
ax1.set_title("Burial does not separate the measured alleles",
              loc="left", fontsize=12.5, fontweight="bold", color=INK, pad=30)
ax1.text(0, 1.015,
         "The control fails: Tyr155Cys is tolerated and almost fully buried.\n"
         "So burial alone cannot call this variant — reported, not dropped.",
         transform=ax1.transAxes, fontsize=9.2, color=INK2, va="bottom",
         linespacing=1.4)

# call out the control failure
ax1.annotate("tolerated, yet more buried than\ntwo of the destabilising alleles",
             xy=(3, by[155]["rsa_percent"] + 2.5), xytext=(2.25, 30),
             fontsize=8.6, color=INK2, ha="center", va="bottom",
             arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))

# ---------------------------------------------------------------- right
dist = d["query_cb_distances"]
lm = [(1012, "Leu1012 — measured destabilising"),
      (882, "Asp882 — catalytic proton acceptor"),
      (909, "Ile909 — measured destabilising"),
      (911, "Asp911 — catalytic triad"),
      (795, "Lys795 — pocket, measured destabilising")]
lm = [(p, t) for p, t in lm if str(p) in dist or p in dist]
vals2 = [dist.get(str(p), dist.get(p)) for p, _ in lm]
names2 = [t for _, t in lm]
cols2 = [ORANGE if p in (795, 909, 1012) else MUTED for p, _ in lm]

y = range(len(lm))
b2 = ax2.barh(list(y), vals2, color=cols2, height=.58,
              edgecolor=SURFACE, linewidth=2)
for bb, v in zip(b2, vals2):
    ax2.text(v + .7, bb.get_y() + bb.get_height() / 2, f"{v:.1f} Å",
             va="center", fontsize=9.6, fontweight="bold", color=INK)

ax2.axvline(8, color=AQUA, lw=1.3, ls=(0, (4, 3)), zorder=0)
ax2.text(8.6, -.46, "8 Å — contact range", fontsize=8.6,
         color=AQUA, va="bottom", ha="left")

ax2.set_yticks(list(y))
ax2.set_yticklabels(names2, fontsize=8.8, color=INK2)
ax2.invert_yaxis()
ax2.set_xlabel("Cβ–Cβ distance from Asn1002  (Å)")
ax2.set_xlim(0, 42)
ax2.set_ylim(len(lm) - .45, -.72)
ax2.grid(axis="x", color=GRID, lw=.8)
ax2.set_axisbelow(True)
ax2.set_title("Asn1002 is not in the structural pocket",
              loc="left", fontsize=12.5, fontweight="bold", color=INK, pad=30)
ax2.text(0, 1.015,
         "22–35 Å from the pocket residues, so the Lys795Ala precedent does\n"
         "not transfer. Nearest landmark is Leu1012, 14 Å, in the same lobe.",
         transform=ax2.transAxes, fontsize=9.2, color=INK2, va="bottom",
         linespacing=1.4)

def _ord(v):
    n = int(round(v))
    suf = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suf}"


nb = d["query_neighbours_8A"]
fig.text(0.5, 0.012,
         f"Asn1002 itself: {by[1002]['rsa_percent']}% RSA "
         f"({_ord(by[1002]['burial_percentile_in_domain'])} burial percentile "
         f"in the domain), pLDDT {by[1002]['plddt']}, "
         f"{nb['hydrophobic']} of {nb['n']} neighbours within 8 Å hydrophobic. "
         f"Asn→Lys places a longer, positively charged side chain there.",
         ha="center", fontsize=9, color=INK2)

fig.tight_layout(rect=(0, 0.045, 1, 1))
os.makedirs(OUT, exist_ok=True)
path = f"{OUT}/t2-structure-context.png"
fig.savefig(path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  wrote {path}")
