"""Repeat 05's question on OptimusKG, and check where the candidate drugs land.

OptimusKG is about five times the size of PrimeKG and labels drug-protein edges
by mechanism (inhibitor, agonist, activator ...) as well as target / enzyme /
transporter. So the same one-interaction reach is asked again, and the drugs the
report discusses are looked up one by one:

  * reachable   the drug is in the set, with the interactor(s) it came through;
  * present     the drug is in the graph but no MVA-gene interactor reaches it;
  * absent      no drug node of that name or synonym.

"Present but unreachable" is the informative case: the graph knows the drug and
still has no mechanistic route from it to an MVA gene.

The question follows texdata's 06_optimuskg_reach.py (MVA Hackathon 2026); the
implementation, the mechanism-kind split and the three-way status are ours.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

import mvakg

WATCH = (
    "nicotinamide riboside", "nicotinamide", "niacin", "acadesine",
    "chloroquine", "hydroxychloroquine", "dasatinib", "quercetin",
    "tanespimycin", "metformin", "cambinol", "reversine",
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).with_name("06_optimuskg_reach.tsv"))
    args = ap.parse_args()

    nodes = mvakg.optimuskg_nodes()
    interactions, links = mvakg.optimuskg_tables(nodes)
    print(f"OptimusKG: {interactions.height // 2:,} interaction pairs, "
          f"{links.height:,} drug-gene edges")

    seeds = {g: mvakg.ENSEMBL[g] for g in mvakg.MVA_GENES}
    reach = mvakg.reach_through_interactors(interactions, links, seeds)
    n_partners = mvakg.interactor_counts(interactions, seeds)

    print("\nper MVA gene:")
    per_gene = dict(reach.group_by("mva_gene")
                         .agg(pl.col("drug_id").n_unique()).iter_rows())
    for gene in mvakg.MVA_GENES:
        print(f"  {gene:7s} {n_partners[gene]:4d} interactors  "
              f"{per_gene.get(gene, 0):4d} drugs")

    mvakg.summarise_reach(reach)
    print("\nedge kinds among reachable drug rows:")
    for kind, n in reach.group_by("kind").len().sort("len", descending=True).iter_rows():
        print(f"  {kind:32s} {n}")

    print("\nwatch list:")
    ids = mvakg.find_by_name(nodes, "DRG", WATCH)
    for name in WATCH:
        drug_id = ids[name]
        if drug_id is None:
            print(f"  {name:22s} absent")
            continue
        hits = reach.filter(pl.col("drug_id") == drug_id)
        if hits.is_empty():
            print(f"  {name:22s} present ({drug_id}), not reachable")
            continue
        via = sorted({f"{g}-{i} ({k})" for g, i, k in
                      hits.select("mva_gene", "interactor", "kind").iter_rows()})
        mark = "  [!] via a checkpoint-weakening protein" if hits["flag"].any() else ""
        print(f"  {name:22s} reachable via {', '.join(via)}{mark}")

    reach.select("mva_gene", "interactor", "drug", "drug_id", "kind", "flag").write_csv(
        args.out, separator="\t")
    print(f"\nwrote {args.out.name} ({reach.height} rows)")


if __name__ == "__main__":
    main()
