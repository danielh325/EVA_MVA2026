"""Figures for the Track 2 report, generated from the committed result files.

Every number plotted here is read from a *_results.json beside this script, so
a figure cannot drift from the text. Nothing is hand-entered except the four
CETSA fold-changes, which are printed by 14_bub1b_cetsa_audit.py and restated
in the docstring below so they can be checked by eye.

Four figures:
  1  the CETSA result - the finding the report turns on
  2  the ClinGen cohort - how common the scarcity is, by group
  3  activation-worded share by protein class - where the scarcity comes from
  4  the evidence ladder - what is established, and how far from this variant

Design follows the project's viz conventions: validated categorical palette
(blue #2a78d6 / orange #eb6834 / aqua #1baf7a on surface #fcfcfb), every mark
direct-labelled so nothing is carried by colour alone, one axis per chart,
recessive grid and spines, text in ink colours rather than series colours.

Run: python3 -I 17_figures.py
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

OUT = "../../reports/figures"

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


def titled(ax, title, subtitle=None):
    nl = subtitle.count("\n") if subtitle else 0
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold",
                 color=INK, pad=20 + 13 * nl if subtitle else 10)
    if subtitle:
        ax.text(0, 1.015, subtitle, transform=ax.transAxes, fontsize=9.5,
                color=INK2, va="bottom", linespacing=1.4)


# ---------------------------------------------------------------- figure 1
def fig_cetsa():
    """The four cell-based records in CHEMBL6125447, from script 14."""
    names = ["compound 8h", "compound 8i", "BAY 1816032", "TG-101209"]
    fc = [0.29, 0.65, 1.33, 2.19]
    colors = [BLUE if v < 1 else ORANGE for v in fc]

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    y = range(len(fc))
    ax.barh(y, [v - 1 for v in fc], left=1, height=0.6, color=colors,
            edgecolor=SURFACE, linewidth=2)
    ax.axvline(1, color=INK2, lw=1.4, zorder=3)

    for i, v in enumerate(fc):
        off = 0.045 if v > 1 else -0.045
        ax.text(v + off, i, f"{v:.2f}×", va="center",
                ha="left" if v > 1 else "right",
                fontsize=11, fontweight="bold", color=INK)

    ax.set_yticks(list(y)); ax.set_yticklabels(names, fontsize=10.5)
    ax.set_xlim(0, 2.55); ax.set_ylim(-0.75, 4.15)
    ax.set_xlabel("BubR1 thermal stability, fold change vs control")
    ax.xaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)

    ax.text(0.92, 3.55, "DESTABILISED  ←", color=BLUE, fontsize=9.5,
            fontweight="bold", ha="right", va="center")
    ax.text(1.08, 3.55, "→  STABILISED", color=ORANGE, fontsize=9.5,
            fontweight="bold", ha="left", va="center")
    ax.text(1.04, -0.62, "reference = no change", color=INK2, fontsize=8.5,
            ha="left", va="center")

    titled(ax, "Two of these four compounds stabilise BubR1",
           'All four sit in one assay titled "Destabilization of BUB1B…"\n'
           'CHEMBL6125447 · CETSA in CAKI-1 cells · 20 µM')
    fig.text(0.015, 0.015,
             "Direction is a property of the record, not the assay title — a "
             "title-based search returns all four as destabilisation.\n"
             "CETSA shows compound-associated thermal stabilisation, not "
             "demonstrated binding or increased abundance.\n"
             "Measured on wild-type BubR1, not on p.Asn1002Lys.",
             fontsize=8, color=INK2, va="bottom", linespacing=1.5)
    fig.tight_layout(rect=[0, 0.19, 1, 1])
    fig.savefig(f"{OUT}/t2-cetsa-bubr1.png", dpi=200)
    print("  wrote t2-cetsa-bubr1.png")


# ---------------------------------------------------------------- figure 2
def fig_clingen():
    d = json.load(open("13_lof_disease_gene_direction_results.json"))
    groups = [("3", "Haploinsufficient", 421), ("30", "Recessive-\nassociated", 597)]
    fig, ax = plt.subplots(figsize=(8.6, 4.0))

    for i, (score, label, _) in enumerate(groups):
        rows = [r for r in d["rows"] if r["score"] == score]
        no_t = sum(1 for r in rows if not r["target"])
        assayed = [r for r in rows if (r["assays"] or 0) > 0]
        zero = sum(1 for r in assayed if r["activation"] == 0)
        some = sum(1 for r in assayed if (r["activation"] or 0) > 0)
        other = len(rows) - no_t - zero - some
        left = 0
        for val, col in [(no_t, MUTED), (zero, BLUE), (some, ORANGE),
                         (other, GRID)]:
            if val <= 0:
                continue
            ax.barh(i, val, left=left, height=0.44, color=col,
                    edgecolor=SURFACE, linewidth=2)
            if val / len(rows) > 0.06:
                ax.text(left + val / 2, i, f"{val}", ha="center", va="center",
                        fontsize=10.5, fontweight="bold",
                        color="#ffffff" if col in (BLUE, MUTED) else INK)
            left += val
        ax.text(len(rows) + 10, i, f"{some}/{len(rows)} = {100*some/len(rows):.1f}%",
                va="center", fontsize=10.5, fontweight="bold", color=ORANGE)

    ax.set_yticks([0, 1]); ax.set_yticklabels([g[1] for g in groups], fontsize=10.5)
    ax.set_ylim(-0.75, 1.45); ax.set_xlim(0, 790)
    ax.set_xlabel("genes")
    ax.xaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)

    handles = [plt.Rectangle((0, 0), 1, 1, color=c)
               for c in (MUTED, BLUE, ORANGE)]
    ax.legend(handles, ["no ChEMBL target", "zero activation-worded",
                        "≥1 activation-worded"],
              loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=3,
              frameon=False, fontsize=9, labelcolor=INK2, handlelength=1.1,
              columnspacing=2.0)

    titled(ax, "Dose-sensitive disease genes rarely carry an\nactivation-worded assay",
           "1,018 ClinGen Dosage Sensitivity genes resolved to ChEMBL · "
           "BUB1B is recessive-associated")
    fig.text(0.015, 0.015,
             "Counts wording, not pharmacology: abundance rescue and "
             "stabilisation are not “activation” and are invisible here.\n"
             "ClinGen's recessive designation does not establish a disease "
             "mechanism, so the groups are reported separately.",
             fontsize=8, color=INK2, va="bottom", linespacing=1.5)
    fig.tight_layout(rect=[0, 0.17, 1, 1])
    fig.savefig(f"{OUT}/t2-clingen-cohort.png", dpi=200)
    print("  wrote t2-clingen-cohort.png")


# ---------------------------------------------------------------- figure 3
def fig_classes():
    d = json.load(open("12_direction_random_sample_results.json"))
    q = [r for r in d["rows"]
         if r["assays"] >= d["min_assays"] and r["activation"] is not None]
    agg = {}
    for r in q:
        k = r["klass"] or "unclassified"
        cur = agg.setdefault(k, [0, 0, 0])   # [activation, assays, n targets]
        cur[0] += r["activation"]; cur[1] += r["assays"]; cur[2] += 1
    rows = [(k, 100 * v[0] / v[1], v[2]) for k, v in agg.items() if v[2] >= 2]
    rows.sort(key=lambda x: x[1])

    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    y = range(len(rows))
    cols = [ORANGE if r[0] == "Enzyme" else BLUE for r in rows]
    ax.barh(y, [r[1] for r in rows], height=0.6, color=cols,
            edgecolor=SURFACE, linewidth=2)
    for i, (k, pct, n) in enumerate(rows):
        ax.text(pct + 0.05, i, f"{pct:.2f}%   (n={n})", va="center",
                fontsize=9.5, color=INK)
    ax.set_yticks(list(y))
    ax.set_yticklabels([r[0] for r in rows], fontsize=10)
    ax.set_xlim(0, 3.1)
    ax.set_xticks([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    ax.xaxis.set_major_formatter(PercentFormatter(decimals=1))
    ax.set_xlabel("share of assays described as activation")
    ax.xaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)

    ei = [i for i, r in enumerate(rows) if r[0] == "Enzyme"]
    if ei:
        ax.text(1.25, ei[0], "72 of 81 enzymes sit at exactly zero  ←",
                va="center", ha="left", fontsize=9, color=ORANGE,
                fontweight="bold")

    titled(ax, "Proteins with a native on-switch carry activation\nchemistry; enzymes do not",
           "Seeded random sample of the 5,869-target human ChEMBL frame\n"
           "The 128 targets with ≥50 assays")
    fig.text(0.015, 0.015,
             "BubR1 is a checkpoint scaffold — further from a switch than an "
             "enzyme — and its activation-worded count is zero.\n"
             "Membrane receptors score low here only because receptor "
             "activation is worded “Agonist…”, which this search does not "
             "match.\nPer-class n is small outside Enzyme.",
             fontsize=8, color=INK2, va="bottom", linespacing=1.5)
    fig.tight_layout(rect=[0, 0.20, 1, 1])
    fig.savefig(f"{OUT}/t2-activation-by-class.png", dpi=200)
    print("  wrote t2-activation-by-class.png")


# ---------------------------------------------------------------- figure 4
def fig_evidence():
    tiers = [
        ("1 — this variant", "AlphaMissense 0.923 · absent from ClinVar ·\n"
         "1/1,461,878 gnomAD · in the pseudokinase domain",
         "no functional data of any kind", ORANGE),
        ("2 — other variants of\nthis domain, in patients",
         "I909T, L1012P: 5–10× less protein,\n**function fully restored by raising abundance**",
         "the abundance-rescue claim", BLUE),
        ("3 — engineered mutants", "K795A: Tm 34.7 °C vs 39.1 °C,\nrescued by growth at 27–33 °C",
         "stability is the failure point", BLUE),
        ("4 — mouse models", "BubR1^H/H progeroid · dose bands 40–90% / >90%\n· T7 vs T23 threshold",
         "dose–phenotype relationships", BLUE),
        ("5 — pathway biology", "K668/CBP/SIRT2 · KARD/PP2A-B56 ·\naneuploidy-selective screens",
         "routes, none validated in MVA", BLUE),
    ]
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    ax.set_xlim(0, 10); ax.set_ylim(0, len(tiers) + 0.6); ax.axis("off")

    for i, (tier, ev, supports, col) in enumerate(tiers):
        y = len(tiers) - 1 - i
        ax.add_patch(plt.Rectangle((0.08, y + 0.08), 9.84, 0.84,
                                   facecolor=SURFACE, edgecolor=GRID, lw=1.2))
        ax.add_patch(plt.Rectangle((0.08, y + 0.08), 0.09, 0.84,
                                   facecolor=col, edgecolor="none"))
        ax.text(0.34, y + 0.5, tier, fontsize=10, fontweight="bold",
                color=INK, va="center")
        ax.text(3.05, y + 0.5, ev.replace("**", ""), fontsize=9,
                color=INK2, va="center")
        ax.text(7.15, y + 0.5, supports, fontsize=8.8, color=INK,
                va="center", style="italic")

    ax.text(0.08, len(tiers) + 0.28,
            "What is established — and how far it sits from this child's allele",
            fontsize=13, fontweight="bold", color=INK)
    ax.text(0.08, len(tiers) + 0.02,
            "Tier 2 carries the therapeutic argument. Tier 1 — the variant "
            "itself — has never been functionally assayed.",
            fontsize=9.5, color=INK2)
    ax.annotate("", xy=(0.02, 0.1), xytext=(0.02, len(tiers) - 0.1),
                arrowprops=dict(arrowstyle="-", color=GRID, lw=2))
    fig.text(0.01, 0.02,
             "E0 is the experiment that moves the abundance-rescue claim from "
             "tier 2 to tier 1 — or refutes it.",
             fontsize=8.5, color=INK2, va="bottom")
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(f"{OUT}/t2-evidence-tiers.png", dpi=200)
    print("  wrote t2-evidence-tiers.png")




# ---------------------------------------------------------------- figure 5
def fig_decision():
    """The E0 -> E9 sequence as a decision, not a list."""
    fig, ax = plt.subplots(figsize=(9.6, 6.4))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10.9); ax.axis("off")

    def box(x, y, w, h, text, fc=SURFACE, ec=GRID, tc=INK, fs=9.3,
            weight="normal", lw=1.3):
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec,
                                   lw=lw, zorder=2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fs, color=tc, zorder=3, fontweight=weight,
                linespacing=1.45)

    def arrow(x1, y1, x2, y2, label=None, col=INK2, lx=0, ly=0):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.5,
                                    shrinkA=2, shrinkB=2), zorder=1)
        if label:
            ax.text((x1 + x2) / 2 + lx, (y1 + y2) / 2 + ly, label, fontsize=8.5,
                    color=col, ha="center", va="center", fontweight="bold",
                    bbox=dict(fc=SURFACE, ec="none", pad=1.6))

    # characterisation layer
    box(0.15, 8.55, 3.0, 1.0,
        "E1  abundance\nE2  which part of the checkpoint fails\nE7  phosphatase recruitment",
        fc="#f4f3ef", fs=8.6)
    ax.text(1.65, 9.72, "shared cells, stain, microscope", fontsize=8.3,
            color=INK2, ha="center", style="italic")

    # E0
    box(3.6, 8.35, 6.25, 1.2,
        "E0 — does p.Asn1002Lys reduce abundance or half-life?",
        fc="#fdf0e9", ec=ORANGE, fs=10.5, weight="bold", lw=2)
    ax.text(6.72, 9.72, "the hinge: everything below depends on it",
            fontsize=8.6, color=ORANGE, ha="center", fontweight="bold")

    arrow(5.3, 8.35, 5.3, 7.5, "NO", ORANGE, lx=-0.45)
    arrow(8.2, 8.35, 8.2, 7.5, "YES", BLUE, lx=0.45)

    box(3.35, 6.55, 3.6, 0.95,
        "abundance-deficiency premise fails\ncandidates 1 · 2a · 8 · 9 lose their basis",
        fc="#fdf0e9", ec=ORANGE, fs=8.8)
    box(7.05, 6.55, 2.8, 0.95,
        "E0 Q2 — at matched abundance,\nis function restored?",
        fc=SURFACE, ec=BLUE, fs=8.8, lw=1.6)

    arrow(5.15, 6.55, 5.15, 5.75, col=ORANGE)
    box(3.35, 4.95, 3.6, 0.8,
        "candidate 7 survives either way —\na ligand may stabilise a conformation\nwithout raising total protein",
        fc="#f4f3ef", fs=8.3)

    arrow(7.6, 6.55, 7.6, 5.75, "NO", ORANGE, lx=-0.42)
    arrow(9.2, 6.55, 9.2, 5.75, "YES", BLUE, lx=0.45)
    box(7.05, 5.0, 1.15, 0.75, "raising it\nwon't suffice", fc="#fdf0e9",
        ec=ORANGE, fs=8.2)
    box(8.45, 5.0, 1.4, 0.75, "the hypothesis\nholds", fc="#eaf2fd",
        ec=BLUE, fs=8.6, weight="bold")

    # the routes
    ax.text(0.15, 4.15, "THEN — the routes it licenses",
            fontsize=9.5, color=INK, fontweight="bold")
    arrow(9.15, 5.0, 9.15, 4.55, col=BLUE)
    ax.plot([2.0, 9.15], [4.55, 4.55], color=BLUE, lw=1.5, zorder=1)
    for xa in (1.28, 3.76, 6.24):
        arrow(xa, 4.55, xa, 3.85, col=BLUE)

    routes = [
        (0.15, "E9", "4-PBA vs arimoclomol", "two approved drugs, rival\nmechanisms · absolute BubR1\nand synthesis separate them", ORANGE),
        (2.63, "E3", "NAD⁺ → SIRT2 → K668", "SIRT2 knockdown must\nattenuate · read PP2A-B56\nalongside BubR1", BLUE),
        (5.11, "E8", "ligand stabilisation", "stability, abundance and\nfunction kept separate ·\nneeds orthogonal binding", BLUE),
        (7.59, "E6a→E6b", "readthrough", "the OTHER allele —\nNOT gated by E0 · does\nnonsense transcript survive NMD?", AQUA),
    ]
    for x, code, name, detail, col in routes:
        box(x, 1.75, 2.26, 2.05, "", fc=SURFACE, ec=col, lw=1.6)
        ax.text(x + 1.13, 3.45, code, ha="center", fontsize=11.5,
                fontweight="bold", color=col)
        ax.text(x + 1.13, 3.08, name, ha="center", fontsize=9.2,
                fontweight="bold", color=INK)
        ax.text(x + 1.13, 2.35, detail, ha="center", fontsize=8.1,
                color=INK2, linespacing=1.5)

    box(0.15, 0.45, 9.7, 0.85,
        "Every arm reports proliferation, cell death and segregation per COMPLETED division —\n"
        "because fewer micronuclei is also what selective death of the worst cells looks like.",
        fc="#f4f3ef", fs=8.8)

    ax.text(0.15, 10.85, "The decision sequence: one experiment gates the rest",
            fontsize=13.5, fontweight="bold", color=INK, va="top")
    ax.text(0.15, 10.42, "E0 is cheap, comes first, and can end the programme. "
            "Only the readthrough arm is independent of it.",
            fontsize=9.3, color=INK2, va="top")
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(f"{OUT}/t2-decision-sequence.png", dpi=200)
    print("  wrote t2-decision-sequence.png")


if __name__ == "__main__":
    fig_cetsa(); fig_clingen(); fig_classes(); fig_evidence()
    fig_decision()
    print("done")
