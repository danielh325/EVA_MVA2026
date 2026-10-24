"""Shared graph access for the Track 2 knowledge-graph scripts 01-07.

Every script in that range asks a question of one of two graphs, PrimeKG (the
graph TxGNN was trained on) or OptimusKG (a larger, newer one). This module
reduces both to the same two tables, so one function answers each question on
either graph:

    interactions  a_id, b_id                              protein-protein, both directions
    drug_links    gene_id, gene, drug_id, drug, kind      every drug-protein edge

`kind` is the graph's own label for the drug-protein edge, lower-cased:
PrimeKG has target / enzyme / carrier / transporter, OptimusKG has those plus
mechanism labels such as inhibitor, agonist and activator. Enzyme, carrier and
transporter edges record how a drug is metabolised or moved, not what it does,
so they are kept but marked (`ADME_KINDS`).

Data location: set MVA_DATA to a folder holding

    primekg/kg.csv                       Harvard Dataverse, doi:10.7910/DVN/IXA7BM
    txgnn-data/                          TxData folder (kg.csv, node.csv, edges.csv)
    txgnn-ckpt/TxGNNExplorer/            TxGNN Explorer bundle (model + GraphMask gates)
    work/                                outputs that are too large to commit

OptimusKG is read from the `optimuskg` package's download cache.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl

DATA = Path(os.environ.get("MVA_DATA", "data"))
PRIMEKG_CSV = DATA / "primekg" / "kg.csv"
TXGNN_DATA = DATA / "txgnn-data"
TXGNN_CKPT = DATA / "txgnn-ckpt" / "TxGNNExplorer"
WORK = DATA / "work"

MVA_DISEASE = "mosaic variegated aneuploidy syndrome"
# PrimeKG groups MONDO 13582, 9759, 54736 and 141 under this one disease node.
MVA_PRIMEKG_ID = "13582_9759_54736_141"

MVA_GENES = ("BUB1B", "BUB1", "BUB3", "CEP57", "TRIP13")

# Inhibiting any of these weakens the spindle assembly checkpoint further, which
# is the lesion in MVA. A drug reached through one is flagged, never dropped.
CHECKPOINT_WEAKENING = frozenset({
    "AURKA", "AURKB", "PLK1", "TTK", "CENPE", "KIF11",
    "BUB1", "BUB1B", "MAD2L1", "CDC20",
})

ADME_KINDS = frozenset({"enzyme", "carrier", "transporter"})

# OptimusKG keys genes by Ensembl id; symbols are only inside the JSON properties
# and are not unique, so lookups go through these fixed ids.
ENSEMBL = {
    "BUB1B": "ENSG00000156970", "BUB1": "ENSG00000169679",
    "BUB3": "ENSG00000154473", "CEP57": "ENSG00000166037",
    "TRIP13": "ENSG00000071539", "AURKA": "ENSG00000087586",
    "AURKB": "ENSG00000178999", "PLK1": "ENSG00000166851",
    "TTK": "ENSG00000112742", "CENPE": "ENSG00000138778",
    "KIF11": "ENSG00000138160", "MAD2L1": "ENSG00000164109",
    "CDC20": "ENSG00000117399", "SIRT2": "ENSG00000068903",
    "CREBBP": "ENSG00000005339", "EP300": "ENSG00000100393",
    "NAMPT": "ENSG00000105835", "NMNAT1": "ENSG00000173614",
    "FZR1": "ENSG00000105325", "HSP90AA1": "ENSG00000080824",
    "SMG1": "ENSG00000157106", "UPF1": "ENSG00000005007",
}


# --------------------------------------------------------------------- PrimeKG

def _primekg(path: Path = PRIMEKG_CSV) -> pl.LazyFrame:
    if not path.exists():
        raise SystemExit(f"PrimeKG not found at {path}; set MVA_DATA (see mvakg.py)")
    return pl.scan_csv(path, infer_schema=False)


def primekg_tables(path: Path = PRIMEKG_CSV) -> tuple[pl.DataFrame, pl.DataFrame]:
    """(interactions, drug_links) from PrimeKG, keyed by gene symbol.

    PrimeKG's gene names are HGNC symbols and unique, so the symbol is the id.
    """
    kg = _primekg(path)
    ppi = kg.filter(pl.col("relation") == "protein_protein").select(
        a_id=pl.col("x_name"), b_id=pl.col("y_name"))
    gene_is_x = pl.col("x_type") == "gene/protein"
    links = kg.filter(pl.col("relation") == "drug_protein").select(
        gene_id=pl.when(gene_is_x).then(pl.col("x_name")).otherwise(pl.col("y_name")),
        drug_id=pl.when(gene_is_x).then(pl.col("y_id")).otherwise(pl.col("x_id")),
        drug=pl.when(gene_is_x).then(pl.col("y_name")).otherwise(pl.col("x_name")),
        kind=pl.col("display_relation").str.to_lowercase(),
    ).with_columns(gene=pl.col("gene_id"))
    return _symmetric(ppi).collect(), links.unique().collect()


def primekg_drug_degree(path: Path = PRIMEKG_CSV) -> pl.DataFrame:
    """drug_id, degree: distinct neighbours of each drug node in PrimeKG."""
    kg = _primekg(path)
    sides = [
        kg.filter(pl.col(f"{s}_type") == "drug")
          .select(drug_id=pl.col(f"{s}_id"), nb=pl.col(f"{o}_id") + "|" + pl.col(f"{o}_type"))
        for s, o in (("x", "y"), ("y", "x"))
    ]
    return (pl.concat(sides).unique()
              .group_by("drug_id").agg(degree=pl.len()).collect())


# ------------------------------------------------------------------- OptimusKG

def _optimuskg_file(name: str) -> Path:
    """Path to an OptimusKG table, downloading it on first use.

    Harvard Dataverse answers 403 to the python-requests User-Agent and 200 to
    the same URL from curl, so the download sends a curl-style agent.
    """
    import optimuskg
    import requests.utils

    requests.utils.default_user_agent = lambda name="curl": "curl/8.5.0"
    return optimuskg.get_file(name)


def optimuskg_nodes() -> pl.DataFrame:
    """id, label, name, synonyms for every OptimusKG node.

    `name` is the HGNC symbol for genes and the preferred name otherwise.

    `synonyms` is the lower-cased JSON text of the synonym list. Drug synonyms
    are strings but gene synonyms are objects, so the list is not decoded;
    find_by_name matches a quoted name inside it instead.
    """
    props = pl.col("properties")
    return (pl.scan_parquet(_optimuskg_file("nodes.parquet"))
              .select(
                  "id", "label",
                  name=pl.coalesce(props.str.json_path_match("$.symbol"),
                                   props.str.json_path_match("$.name")),
                  synonyms=props.str.json_path_match("$.synonyms").str.to_lowercase(),
              )
              .collect())


def optimuskg_edges() -> pl.LazyFrame:
    return pl.scan_parquet(_optimuskg_file("edges.parquet")).select(
        "from", "to", "label", "relation")


def optimuskg_tables(nodes: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """(interactions, drug_links) from OptimusKG, keyed by Ensembl / ChEMBL id."""
    edges = optimuskg_edges()
    ppi = edges.filter(pl.col("label") == "GEN-GEN").select(
        a_id=pl.col("from"), b_id=pl.col("to"))
    names = nodes.lazy().select("id", "name")
    links = (edges.filter(pl.col("label") == "DRG-GEN")
                  .select(drug_id=pl.col("from"), gene_id=pl.col("to"),
                          kind=pl.col("relation").str.to_lowercase())
                  .join(names.rename({"id": "gene_id", "name": "gene"}), on="gene_id", how="left")
                  .join(names.rename({"id": "drug_id", "name": "drug"}), on="drug_id", how="left")
                  .unique())
    return _symmetric(ppi).collect(), links.collect()


def optimuskg_degree(ids) -> dict[str, int]:
    """Number of OptimusKG edges touching each id, in either direction."""
    ids = list(ids)
    edges = optimuskg_edges()
    ends = pl.concat([edges.select(id=pl.col("from")), edges.select(id=pl.col("to"))])
    counts = dict(ends.filter(pl.col("id").is_in(ids)).group_by("id").len().collect().iter_rows())
    return {i: counts.get(i, 0) for i in ids}


def ligands(drug_links: pl.DataFrame, gene_id: str) -> pl.DataFrame:
    """drug, kind rows for every drug edge on one gene, sorted by drug name."""
    return (drug_links.filter(pl.col("gene_id") == gene_id)
                      .select(drug=pl.coalesce("drug", "drug_id"), kind="kind")
                      .unique().sort("drug", "kind"))


def find_by_name(nodes: pl.DataFrame, label: str, names) -> dict[str, str | None]:
    """Map each name to the id of an OptimusKG node with that name or synonym."""
    pool = nodes.filter(pl.col("label") == label)
    found = {}
    for name in names:
        key = name.lower()
        hit = pool.filter((pl.col("name").str.to_lowercase() == key)
                          | pl.col("synonyms").str.contains(f'"{key}"', literal=True))
        found[name] = hit["id"][0] if hit.height else None
    return found


# ------------------------------------------------------------------- questions

def _symmetric(pairs: pl.LazyFrame) -> pl.LazyFrame:
    flipped = pairs.select(a_id=pl.col("b_id"), b_id=pl.col("a_id"))
    return pl.concat([pairs, flipped]).unique()


def reach_through_interactors(interactions: pl.DataFrame, drug_links: pl.DataFrame,
                              seeds: dict[str, str]) -> pl.DataFrame:
    """Drugs one protein-protein interaction away from each seed gene.

    seeds maps a display symbol to the id used in `interactions`. Returns one
    row per (seed, interactor, drug, kind) with `adme` and `flag` columns.
    """
    seed_df = pl.DataFrame({"mva_gene": list(seeds), "a_id": list(seeds.values())})
    partners = seed_df.join(interactions, on="a_id").rename({"b_id": "gene_id"})
    return (partners.join(drug_links, on="gene_id")
                    .select("mva_gene", interactor="gene", drug="drug",
                            drug_id="drug_id", kind="kind")
                    .with_columns(adme=pl.col("kind").is_in(list(ADME_KINDS)),
                                  flag=pl.col("interactor").is_in(list(CHECKPOINT_WEAKENING)))
                    .unique()
                    .sort("mva_gene", "interactor", "drug", "kind"))


def interactor_counts(interactions: pl.DataFrame, seeds: dict[str, str]) -> dict[str, int]:
    ids = list(seeds.values())
    counts = (interactions.filter(pl.col("a_id").is_in(ids))
                          .group_by("a_id").agg(n=pl.col("b_id").n_unique()))
    by_id = dict(counts.iter_rows())
    return {sym: by_id.get(i, 0) for sym, i in seeds.items()}


def summarise_reach(reach: pl.DataFrame) -> None:
    """Print the counts every reachability script reports."""
    drugs = reach["drug_id"].n_unique()
    pharm = reach.filter(~pl.col("adme"))["drug_id"].n_unique()
    flagged = reach.filter(pl.col("flag"))["drug_id"].n_unique()
    print(f"\n{drugs} distinct drugs reachable in two hops "
          f"({reach.height} seed-interactor-drug-kind rows)")
    print(f"  {pharm} through at least one non-ADME edge "
          f"(kind not in {', '.join(sorted(ADME_KINDS))}); "
          f"{drugs - pharm} only through metabolism or transport")
    print(f"  {flagged} reached through a checkpoint-weakening interactor")
