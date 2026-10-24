"""Which graph paths does TxGNN's GraphMask say connect MVA to a drug?

GraphMask, released with the TxGNN checkpoint, gives every PrimeKG edge two
gate values (one per message-passing layer) saying how much the trained model
relied on it. An explanation is a short path MVA -> ... -> drug whose edges
the model leaned on.

Method (exact, no beam):
  * edge weight = mean of the two layer gates; the graph is treated as
    undirected and a pair keeps its strongest direction;
  * anatomy and exposure nodes and drug-drug edges are left out (1.4M
    tissue-expression edges and 2.5M drug-drug edges say nothing about a
    mechanism and would otherwise fill every path);
  * every path of 2 to 4 edges is scored by the geometric mean of its
    weights, so lengths compare fairly. Paths are built from two halves
    (MVA side, drug side) that meet in the middle; for each meeting node only
    the best TOP_HALVES halves are kept, separately for halves that do and do
    not pass a gene. The top paths are therefore exact up to TOP_HALVES;
  * only simple paths that pass through at least one gene are reported.

Raw gates (--raw-gates) put every drug through CYP3A4, a drug-metabolising
enzyme that is near every drug in PrimeKG. The default mode removes that by
using PrimeKG's own edge typing: a drug-protein edge is kept only if PrimeKG
calls the protein the drug's *target*, not its enzyme, carrier or
transporter. Drug-side-effect edges are also left out.

The question follows texdata's 03_txgnn_explain.py (MVA Hackathon 2026),
which used a beam search and a gene-name ADME list; this method and code are ours.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

import mvakg

GATES = mvakg.TXGNN_CKPT / "graphmask_output_indication.pkl"
LEFT_OUT_TYPES = {"anatomy", "exposure"}
LEFT_OUT_RELATIONS = {"drug_drug"}
CORRECTED_LEFT_OUT = {"drug_effect", "rev_drug_effect"}
GENE = "gene/protein"
TOP_HALVES = 25
DEFAULT_DRUGS = ("Dasatinib", "Hydroxychloroquine", "Chloroquine", "Paclitaxel")


def target_pairs(kg: Path) -> set[tuple[str, str]]:
    """(drug id, gene id) pairs that PrimeKG labels drug -> *target*."""
    rows = (pl.scan_csv(kg, infer_schema=False)
              .filter((pl.col("relation") == "drug_protein")
                      & (pl.col("display_relation") == "target")
                      & (pl.col("x_type") == "drug"))
              .select("x_id", "y_id").unique().collect())
    return set(rows.iter_rows())


def load_edges(gates: Path, kg: Path, raw: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Undirected weighted edge list (u, v, w, rel) and node table (type, name)."""
    g = pd.read_pickle(gates)
    keep = (~g["x_type"].isin(LEFT_OUT_TYPES) & ~g["y_type"].isin(LEFT_OUT_TYPES)
            & ~g["relation"].isin(LEFT_OUT_RELATIONS))
    if not raw:
        keep &= ~g["relation"].isin(CORRECTED_LEFT_OUT)
    g = g[keep]
    if not raw:
        g = drop_non_target_drug_edges(g, target_pairs(kg))

    x = g["x_type"] + "|" + g["x_idx"].astype(str)
    y = g["y_type"] + "|" + g["y_idx"].astype(str)
    codes, keys = pd.factorize(pd.concat([x, y], ignore_index=True))
    n = len(g)
    edges = pd.DataFrame({"u": codes[:n], "v": codes[n:],
                          "w": (g["indication_layer1_att"].to_numpy(float)
                                + g["indication_layer2_att"].to_numpy(float)) / 2,
                          "rel": g["relation"].str.removeprefix("rev_").to_numpy()})
    both = pd.concat([edges, edges.rename(columns={"u": "v", "v": "u"})], ignore_index=True)
    both = (both[both["u"] != both["v"]]
                .sort_values("w", ascending=False)
                .drop_duplicates(["u", "v"]))

    names = pd.concat([g["x_name"], g["y_name"]], ignore_index=True)
    types = pd.concat([g["x_type"], g["y_type"]], ignore_index=True)
    nodes = (pd.DataFrame({"code": codes, "type": types, "name": names})
               .drop_duplicates("code").set_index("code").sort_index())
    return both.reset_index(drop=True), nodes


def drop_non_target_drug_edges(g: pd.DataFrame, targets: set) -> pd.DataFrame:
    is_dp = g["relation"].isin({"drug_protein", "rev_drug_protein"})
    forward = g["relation"] == "drug_protein"
    drug = np.where(forward, g["x_id"], g["y_id"])
    gene = np.where(forward, g["y_id"], g["x_id"])
    gene = pd.Series(gene, index=g.index).astype(str).str.removesuffix(".0")
    is_target = pd.Series(list(zip(drug, gene)), index=g.index).isin(targets)
    return g[~is_dp | is_target]


def halves(edges: pd.DataFrame, source: int, is_gene: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One- and two-edge paths out of `source`.

    one: end, w1, rel1          two: mid, end, w1, w2, rel1, rel2, gene
    `gene` says whether mid or end is a gene; ends keep their best TOP_HALVES
    halves per gene/no-gene group.
    """
    one = edges[edges["u"] == source].rename(columns={"v": "end", "w": "w1", "rel": "rel1"})
    one = one.drop(columns="u")
    two = one.rename(columns={"end": "mid"}).merge(
        edges.rename(columns={"u": "mid", "v": "end", "w": "w2", "rel": "rel2"}), on="mid")
    two = two[two["end"] != source]
    two["s"] = two["w1"] * two["w2"]
    two["gene"] = is_gene[two["mid"].to_numpy()] | is_gene[two["end"].to_numpy()]
    two = (two.sort_values("s", ascending=False)
              .groupby(["end", "gene"], sort=False).head(TOP_HALVES))
    return one, two


def join_paths(d1, d2, r1, r2, disease: int, drug: int, is_gene: np.ndarray) -> pd.DataFrame:
    """All simple 2-4 edge paths disease -> drug, as node and relation tuples."""
    rows = []
    # 2 edges: disease - m - drug
    p = d1.merge(r1, on="end", suffixes=("_d", "_r"))
    rows.append(pd.DataFrame({
        "nodes": list(zip([disease] * len(p), p["end"], [drug] * len(p))),
        "rels": list(zip(p["rel1_d"], p["rel1_r"])),
        "prod": p["w1_d"] * p["w1_r"], "hops": 2,
        "gene": is_gene[p["end"].to_numpy()]}))
    # 3 edges: disease - a - m - drug
    p = d2.merge(r1, on="end", suffixes=("", "_r"))
    p = p[p["mid"] != drug]
    rows.append(pd.DataFrame({
        "nodes": list(zip([disease] * len(p), p["mid"], p["end"], [drug] * len(p))),
        "rels": list(zip(p["rel1"], p["rel2"], p["rel1_r"])),
        "prod": p["s"] * p["w1_r"], "hops": 3, "gene": p["gene"].to_numpy()}))
    # 4 edges: disease - a - m - b - drug
    p = d2.merge(r2, on="end", suffixes=("_d", "_r"))
    p = p[(p["mid_d"] != p["mid_r"]) & (p["mid_d"] != drug) & (p["mid_r"] != disease)]
    rows.append(pd.DataFrame({
        "nodes": list(zip([disease] * len(p), p["mid_d"], p["end"], p["mid_r"], [drug] * len(p))),
        "rels": list(zip(p["rel1_d"], p["rel2_d"], p["rel2_r"], p["rel1_r"])),
        "prod": p["s_d"] * p["s_r"], "hops": 4,
        "gene": (p["gene_d"] | p["gene_r"]).to_numpy()}))
    paths = pd.concat(rows, ignore_index=True)
    paths["score"] = paths["prod"] ** (1 / paths["hops"])
    return paths.sort_values("score", ascending=False)


def render(nodes_tuple, rels, nodes: pd.DataFrame) -> str:
    names = [nodes.at[n, "name"] for n in nodes_tuple]
    parts = [names[0]]
    for rel, name in zip(rels, names[1:]):
        parts.append(f" -[{rel}]- {name}")
    return "".join(parts)


def node_by_name(nodes: pd.DataFrame, name: str, node_type: str) -> int | None:
    hit = nodes[(nodes["type"] == node_type) & (nodes["name"].str.lower() == name.lower())]
    return int(hit.index[0]) if len(hit) else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drugs", nargs="+", default=list(DEFAULT_DRUGS))
    ap.add_argument("--gates", type=Path, default=GATES)
    ap.add_argument("--kg", type=Path, default=mvakg.PRIMEKG_CSV)
    ap.add_argument("--raw-gates", action="store_true",
                    help="keep every drug-protein and drug-side-effect edge")
    ap.add_argument("--top", type=int, default=6)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    mode = "raw gates" if args.raw_gates else "target edges only, no side-effect edges"
    edges, nodes = load_edges(args.gates, args.kg, args.raw_gates)
    is_gene = (nodes["type"] == GENE).to_numpy()
    print(f"GraphMask gates, {mode}: {len(nodes):,} nodes, {len(edges) // 2:,} edges")

    disease = node_by_name(nodes, mvakg.MVA_DISEASE, "disease")
    d1, d2 = halves(edges, disease, is_gene)
    print(f"MVA: {len(d1)} neighbours")

    found = []
    for name in args.drugs:
        drug = node_by_name(nodes, name, "drug")
        if drug is None:
            print(f"\n{name}: not a drug node in the gate table")
            continue
        r1, r2 = halves(edges, drug, is_gene)
        paths = join_paths(d1, d2, r1, r2, disease, drug, is_gene)
        with_gene = paths[paths["gene"]]
        print(f"\n{name}: {len(with_gene):,} simple paths through a gene "
              f"({len(paths):,} in all); top {args.top}:")
        if with_gene.empty:
            print("    none")
        for _, p in with_gene.head(args.top).iterrows():
            text = render(p["nodes"], p["rels"], nodes)
            print(f"    {p['score']:.4f}  {p['hops']} edges  {text}")
            found.append((name, round(p["score"], 4), p["hops"], text))

    if args.out:
        pd.DataFrame(found, columns=["drug", "score", "edges", "path"]) \
          .to_csv(args.out, sep="\t", index=False)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
