"""Which drugs sit one protein-protein interaction away from the MVA genes, in PrimeKG?

None of BUB1B, BUB1, BUB3, CEP57 or TRIP13 has a drug edge of its own, so the
nearest mechanism-anchored question PrimeKG can answer is one step out: which
drugs act on a protein that interacts with a protein the disease disrupts?
The answer is a set, not a ranking, and every row names the interactor it came
through so the claim can be checked by hand.

Two things are reported separately because they change what a hit means:

  * the PrimeKG edge kind. A drug linked to an interactor only as its enzyme,
    carrier or transporter is reached through how it is cleared, not through
    what it does, so those hits are counted apart from target hits;
  * interactors whose inhibition weakens the spindle checkpoint further
    (mvakg.CHECKPOINT_WEAKENING). Drugs reached through them are flagged.

The question follows texdata's 05_mechanism_anchored.py (MVA Hackathon 2026);
the implementation and the edge-kind split are ours.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

import mvakg

TOP = 25


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--kg", type=Path, default=mvakg.PRIMEKG_CSV)
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).with_name("05_mechanism_anchored.tsv"))
    args = ap.parse_args()

    interactions, links = mvakg.primekg_tables(args.kg)
    print(f"PrimeKG: {interactions.height // 2:,} interaction pairs, "
          f"{links.height:,} drug-protein edges")
    print(links.group_by("kind").len().sort("len", descending=True)
               .select(pl.format("  {}: {}", "kind", "len")).to_series().str.join("\n").item())

    seeds = {g: g for g in mvakg.MVA_GENES}
    reach = mvakg.reach_through_interactors(interactions, links, seeds)
    n_partners = mvakg.interactor_counts(interactions, seeds)

    print("\nper MVA gene:")
    per_gene = reach.group_by("mva_gene").agg(
        drugs=pl.col("drug_id").n_unique(),
        target_drugs=pl.col("drug_id").filter(~pl.col("adme")).n_unique())
    for gene, drugs, target_drugs in per_gene.sort("mva_gene").iter_rows():
        print(f"  {gene:7s} {n_partners[gene]:4d} interactors  {drugs:4d} drugs "
              f"({target_drugs} through a non-ADME edge)")

    mvakg.summarise_reach(reach)

    by_drug = (reach.group_by("drug")
                    .agg(paths=pl.struct("mva_gene", "interactor").n_unique(),
                         via=pl.format("{}-{}", "mva_gene", "interactor").unique().sort(),
                         kinds=pl.col("kind").unique().sort(),
                         flag=pl.col("flag").any())
                    .sort(["paths", "drug"], descending=[True, False]))
    print(f"\nmost-connected reachable drugs (top {TOP}):")
    for drug, paths, via, kinds, flag in by_drug.head(TOP).iter_rows():
        mark = "  [!] via a checkpoint-weakening protein" if flag else ""
        print(f"  {drug:30.30s} {paths:2d} paths  [{'/'.join(kinds)}]  "
              f"via {', '.join(via[:4])}{mark}")

    flagged = by_drug.filter(pl.col("flag"))["drug"].to_list()
    print(f"\nflagged ({len(flagged)}): {', '.join(flagged)}")

    reach.select("mva_gene", "interactor", "drug", "kind", "flag").write_csv(
        args.out, separator="\t")
    print(f"\nwrote {args.out.name} ({reach.height} rows)")


if __name__ == "__main__":
    main()
