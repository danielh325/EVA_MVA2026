"""BUB1B has no ligand. Which proteins next to it in the mechanism do, and which way do they act?

The report's mechanism sets BubR1 abundance through K668 acetylation (written
by CBP/p300, erased by SIRT2), NAD+ supply (NAMPT, NMNAT1), degradation
(APC/C-Cdh1 via FZR1), folding (HSP90) and, for the nonsense allele, NMD
(SMG1, UPF1). For each, the useful drug has a required direction: an inhibitor
of a writer, an activator of an eraser or of NAD+ synthesis.

A graph edge that only says "drug touches protein" cannot answer that.
OptimusKG labels many drug-gene edges with a mechanism (inhibitor, agonist,
activator ...), so this lists every ligand with its label and counts how many
point the needed way. Where the label is only "target", the direction is
unknown and is counted as such, not guessed.

The protein list follows texdata's 07_alternative_targets.py (MVA Hackathon
2026); the direction check and the implementation are ours.
"""
from __future__ import annotations

import polars as pl

import mvakg

# OptimusKG mechanism labels that say which way a drug moves its target. Every
# other label (target, binding_agent, modulator, substrate, enzyme ...) does not.
RAISE = ("activator", "agonist", "partial_agonist", "positive_modulator",
         "positive_allosteric_modulator", "stabiliser", "opener")
LOWER = ("inhibitor", "antagonist", "allosteric_antagonist", "blocker",
         "negative_modulator", "negative_allosteric_modulator", "inverse_agonist",
         "degrader", "rnai_inhibitor", "antisense_inhibitor", "disrupting_agent")
# gene -> (role in the mechanism, direction a useful drug must push it)
TARGETS = {
    "CREBBP":   ("writes the K668 acetyl mark", LOWER),
    "EP300":    ("writes the K668 acetyl mark", LOWER),
    "SIRT2":    ("erases the K668 acetyl mark", RAISE),
    "NAMPT":    ("rate-limiting step of NAD+ salvage", RAISE),
    "NMNAT1":   ("NAD+ synthesis", RAISE),
    "FZR1":     ("APC/C-Cdh1 adaptor that degrades BubR1", LOWER),
    "HSP90AA1": ("chaperone for a fold-defective client", RAISE),
    "SMG1":     ("NMD kinase; blocking it spares the nonsense transcript", LOWER),
    "UPF1":     ("NMD core helicase", LOWER),
}
SHOWN = 8


def classify(kinds: pl.Expr, useful: tuple[str, ...]) -> pl.Expr:
    wrong = LOWER if useful is RAISE else RAISE
    return (pl.when(kinds.is_in(list(useful))).then(pl.lit("useful direction"))
              .when(kinds.is_in(list(wrong))).then(pl.lit("wrong direction"))
              .otherwise(pl.lit("direction unknown")))


def main() -> None:
    nodes = mvakg.optimuskg_nodes()
    _, links = mvakg.optimuskg_tables(nodes)

    summary = []
    for gene, (role, useful) in TARGETS.items():
        found = mvakg.ligands(links, mvakg.ENSEMBL[gene])
        found = found.with_columns(direction=classify(pl.col("kind"), useful))
        counts = dict(found.group_by("direction").agg(pl.col("drug").n_unique()).iter_rows())
        n_drugs = found["drug"].n_unique()
        print(f"\n{gene} - {role}; needs {'a raiser' if useful is RAISE else 'a lowerer'}")
        print(f"  {n_drugs} ligand(s)")
        for direction in ("useful direction", "direction unknown", "wrong direction"):
            rows = found.filter(pl.col("direction") == direction)
            if rows.is_empty():
                continue
            shown = [f"{d} ({k})" for d, k in rows.select("drug", "kind").head(SHOWN).iter_rows()]
            more = f" ... +{rows.height - SHOWN}" if rows.height > SHOWN else ""
            print(f"    {direction:17s} {counts[direction]:3d}: {', '.join(shown)}{more}")
        summary.append((gene, n_drugs, counts.get("useful direction", 0)))

    print("\nsummary: ligands / with a label pointing the useful way")
    for gene, n_drugs, n_useful in summary:
        print(f"  {gene:8s} {n_drugs:3d} / {n_useful}")


if __name__ == "__main__":
    main()
