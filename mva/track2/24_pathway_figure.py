"""The MVA biochemical network, with every candidate placed on it.

One picture answering the question a reader asks after §5: where does each
candidate actually act, and why are the exclusions excluded? Three bands,
read top to bottom, because that is the causal order:

  1  WHY BubR1 IS SCARCE   the lesion. Both alleles, and the degradation
                           route that sets the level of the one that
                           survives. Candidates 1, 2a, 3, 7, 8, 9 act here.
  2  WHAT BubR1 DOES       the checkpoint itself, and the four proteins in
                           it that ARE druggable - which is exactly §6's
                           exclusion list, verified in PrimeKG by `23`.
  3  WHAT GOES WRONG       the consequence. Mosaic aneuploidy, its
                           proteostatic cost, and tumour risk. Candidates
                           4, 5, 10 act here.

Nothing in this figure is a new claim: every arrow restates mechanism cited
in §2, §5 or §6. It exists because eleven candidates spread over three
causal layers is hard to hold in the head from prose alone.

Run: python3 -I 24_pathway_figure.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = "../../reports/figures"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#a3a29c"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
GRID = "#e6e5e1"
BAND = "#f4f3ef"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "text.color": INK,
})

fig, ax = plt.subplots(figsize=(13.4, 11.8))
ax.set_xlim(0, 100); ax.set_ylim(-13, 118); ax.axis("off")

PAD_HEAD = 3.0      # band top -> header baseline
PAD_SUB = 6.6       # band top -> subtitle baseline
PAD_BODY = 11.0     # band top -> first row of boxes may start


def band(top, h, label, sub):
    ax.add_patch(FancyBboxPatch((1.2, top - h), 97.6, h,
                                boxstyle="round,pad=0.4", fc=BAND, ec="none",
                                zorder=0))
    ax.text(3.0, top - PAD_HEAD, label, fontsize=11.5, fontweight="bold",
            color=INK, va="center")
    ax.text(3.0, top - PAD_SUB, sub, fontsize=9, color=INK2, va="center")


def node(x, y, w, h, text, fc=SURFACE, ec=MUTED, fs=9.2, bold=False,
         tc=INK, lw=1.3):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3",
                                fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold" if bold else "normal",
            linespacing=1.4, zorder=3)


def arrow(x1, y1, x2, y2, col=MUTED, lw=1.5, style="-|>", rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=13, color=col, lw=lw,
                                 connectionstyle=f"arc3,rad={rad}", zorder=1))


def cand(x, y, nums, note, col=BLUE):
    ax.text(x, y, nums, fontsize=9.4, fontweight="bold", color=col,
            ha="left", va="center", zorder=4)
    ax.text(x + 5.2, y, note, fontsize=8.5, color=INK2, ha="left",
            va="center", zorder=4)


ax.text(1.2, 116, "Where each candidate acts on this child's biology",
        fontsize=15.5, fontweight="bold", color=INK, va="center")

# ============================================================ band 1
B1 = 110.0
band(B1, 41, "1 \u00b7 WHY BubR1 IS SCARCE \u2014 the lesion",
     "Both alleles, and the degradation route that sets the level of the one that survives")
r1a, r1b = B1 - PAD_BODY - 7.5, B1 - PAD_BODY - 18.5
node(4, r1a, 20, 7.5, "nonsense allele\np.Leu737Ter", ec=ORANGE)
node(4, r1b, 20, 7.5, "missense allele\np.Asn1002Lys", ec=ORANGE)
node(30, r1a, 15, 7.5, "NMD\n(predicted)", fs=8.8)
node(30, r1b, 15, 7.5, "unstable fold", fs=8.8)
node(51, r1a, 15, 7.5, "no protein", fs=8.8, tc=INK2)
node(51, r1b, 15, 7.5, "proteasome", fs=8.8)
node(74, (r1a + r1b) / 2, 22, 9.0, "too little\nBubR1", fc="#eaf1fb",
     ec=BLUE, bold=True, fs=11.5)
arrow(24, r1a + 3.75, 30, r1a + 3.75); arrow(45, r1a + 3.75, 51, r1a + 3.75)
arrow(24, r1b + 3.75, 30, r1b + 3.75); arrow(45, r1b + 3.75, 51, r1b + 3.75)
arrow(66, r1a + 3.75, 74, (r1a + r1b) / 2 + 6.5, rad=-0.1)
arrow(66, r1b + 3.75, 74, (r1a + r1b) / 2 + 2.5, rad=0.1)

k668 = r1b - 7.0
node(30, k668, 36, 5.0, "K668 acetylation:  CBP writes it \u2192 degradation"
     "   \u00b7   SIRT2 removes it \u2192 stability", fs=8.3, ec=GRID)
arrow(52, k668 + 5.0, 56, r1b, col=MUTED, lw=1.1)
cand(4, k668 + 2.5, "3", "readthrough \u2014 rescues the nonsense allele",
     col=AQUA)
cand(70, k668 + 5.5, "7 \u00b7 8 \u00b7 9", "stabilise the fold / raise\nchaperone capacity")
cand(70, k668 + 0.5, "1 \u00b7 2a", "act on the K668 switch")

# ============================================================ band 2
B2 = 66.0
band(B2, 36, "2 \u00b7 WHAT BubR1 DOES \u2014 the checkpoint",
     "The only proteins here carrying drug edges in PrimeKG are the four this report excludes (\u00a76, script 23)")
r2 = B2 - PAD_BODY - 7.5
node(4, r2, 17, 7.5, "unattached\nkinetochore", fs=8.8)
node(25, r2, 13, 7.5, "MPS1 / TTK", fs=8.8, ec=ORANGE)
node(42, r2, 14, 7.5, "BUB1 \u00b7 BUB3", fs=8.8)
node(60, r2 - 1.0, 15, 9.5, "BUBR1", fc="#eaf1fb", ec=BLUE, bold=True, fs=12.5)
node(79, r2, 17, 7.5, "MCC\ninhibits APC/C", fs=8.8)
for a, b in ((21, 25), (38, 42), (56, 60), (75, 79)):
    arrow(a, r2 + 3.75, b, r2 + 3.75)

r2b = r2 - 13.5
node(42, r2b, 20, 7.5, "PP2A-B56\nvia the KARD", fs=8.8)
node(68, r2b, 18, 7.5, "Aurora B\nerror correction", fs=8.8, ec=ORANGE)
arrow(64, r2 - 1.0, 56, r2b + 7.5, rad=0.15)
arrow(62, r2b + 3.75, 68, r2b + 3.75, style="-[", col=MUTED)
ax.text(65, r2b - 1.4, "opposes", fontsize=8, color=INK2, ha="center")

ax.text(4, r2 - 5.5, "\u2717 EXCLUDED on mechanism (\u00a76)", fontsize=9,
        fontweight="bold", color=ORANGE, va="center")
ax.text(4, r2 - 10.5, "MPS1/TTK \u00b7 Aurora B \u00b7 PLK1\n"
        "KIF11/Eg5 \u00b7 CENPE \u00b7 HSP90i", fontsize=8.6, color=INK2,
        va="center", linespacing=1.6)
ax.text(4, r2 - 15.5, "these are the druggable ones", fontsize=8,
        color=MUTED, style="italic", va="center")

# ============================================================ band 3
B3 = 27.0
band(B3, 35, "3 \u00b7 WHAT GOES WRONG \u2014 the consequence",
     "Where the only approved candidate with direct human-cell evidence acts")
r3 = B3 - PAD_BODY - 7.5
node(4, r3, 18, 7.5, "checkpoint\nleaks", fs=8.8)
node(26, r3, 19, 7.5, "missegregation", fs=8.8)
node(49, r3 - 1.0, 20, 9.5, "mosaic\naneuploidy", fc="#fdeee7", ec=ORANGE,
     bold=True, fs=11)
node(75, r3 + 4.0, 21, 7.0, "tumour\npredisposition", fs=8.6)
node(75, r3 - 5.5, 21, 7.5, "proteostatic +\nmitochondrial load", fs=8.6)
arrow(22, r3 + 3.75, 26, r3 + 3.75); arrow(45, r3 + 3.75, 49, r3 + 3.75)
arrow(69, r3 + 5.0, 75, r3 + 7.0, rad=-0.1)
arrow(69, r3 + 2.5, 75, r3 - 1.5, rad=0.1)
arrow(85, (r1a + r1b) / 2, 85, B2 + 0.5, col=BLUE, lw=1.6)

cand(4, r3 - 7.5, "10", "metformin \u2014 restores autophagy\u2013mitophagy flux;\nthe only route that does not wait on E0", col=AQUA)
cand(4, r3 - 13.5, "4 \u00b7 5", "hydroxychloroquine, senolytics \u2014 treat the\nconsequence; 4 now looks pointed the wrong way", col=ORANGE)

ax.text(50, -11.0,
        "Every arrow restates mechanism cited in \u00a72, \u00a75 or \u00a76. "
        "Numbers are candidate IDs from \u00a75.1. Nothing here is a new claim.",
        ha="center", fontsize=8.6, color=MUTED, style="italic")

fig.tight_layout()
os.makedirs(OUT, exist_ok=True)
path = f"{OUT}/t2-pathway-map.png"
fig.savefig(path, dpi=185, bbox_inches="tight")
plt.close(fig)
print(f"  wrote {path}")
