"""Does OptimusKG fill the two gaps that limited PrimeKG and TxGNN?

PrimeKG has no drug edge on any MVA gene, and the NAD+ precursors are missing
from TxGNN's annotated drug set (see 02). OptimusKG is a newer and larger graph,
so before anything is run on it, three lookups:

  1. edges and drug edges on each MVA gene, and on the checkpoint and
     acetylation proteins the report discusses;
  2. whether each NAD+ precursor and each probe drug has a drug node;
  3. whether MVA itself has a disease node.

Lookups go through fixed Ensembl ids and through drug names *and* synonyms;
OptimusKG keeps names inside a JSON column, so a lookup on a plain name column
returns zero for everything.

The question follows texdata's 04_optimuskg_coverage.py (MVA Hackathon 2026);
the implementation is ours.
"""
from __future__ import annotations

import polars as pl

import mvakg

CONTEXT_GENES = ("AURKB", "PLK1", "TTK", "CENPE", "SIRT2", "CREBBP")
NAD_PRECURSORS = ("nicotinamide riboside", "nicotinamide mononucleotide",
                  "nicotinamide", "niacin", "nicotinic acid", "nadh")
PROBES = ("dasatinib", "quercetin", "chloroquine", "hydroxychloroquine",
          "metformin", "tanespimycin", "navitoclax", "fisetin", "acadesine",
          "vorinostat", "romidepsin", "reversine")
SHOWN = 6


def gene_table(title: str, genes, degree, links) -> None:
    print(f"-- {title}")
    for gene in genes:
        gene_id = mvakg.ENSEMBL[gene]
        drugs = mvakg.ligands(links, gene_id)["drug"].unique().sort().to_list()
        example = ", ".join(drugs[:SHOWN]) if drugs else "none"
        print(f"  {gene:8s} edges={degree[gene_id]:5d}  drug edges={len(drugs):3d}  {example}")


def lookup_table(title: str, nodes, label: str, names) -> None:
    print(f"\n=== {title}")
    for name, node_id in mvakg.find_by_name(nodes, label, names).items():
        print(f"  {name:42s} {node_id or 'NOT FOUND'}")


def main() -> None:
    nodes = mvakg.optimuskg_nodes()
    _, links = mvakg.optimuskg_tables(nodes)
    n_edges = mvakg.optimuskg_edges().select(pl.len()).collect().item()
    print(f"OptimusKG: {nodes.height:,} nodes, {n_edges:,} edges")
    print(nodes.group_by("label").len().sort("len", descending=True)
               .select(pl.format("  {} {}", "label", "len")).to_series().str.join("\n").item())

    genes = mvakg.MVA_GENES + CONTEXT_GENES
    degree = mvakg.optimuskg_degree(mvakg.ENSEMBL[g] for g in genes)
    print("\n=== edges and drug edges per gene")
    gene_table("MVA genes", mvakg.MVA_GENES, degree, links)
    gene_table("checkpoint and acetylation context", CONTEXT_GENES, degree, links)

    lookup_table("NAD+ precursors", nodes, "DRG", NAD_PRECURSORS)
    lookup_table("probe drugs", nodes, "DRG", PROBES)
    lookup_table("MVA disease node", nodes, "DIS",
                 (mvakg.MVA_DISEASE,))


if __name__ == "__main__":
    main()
