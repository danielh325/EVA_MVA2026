"""Rank TxGNN's zero-shot scores for MVA, and test them for a node-degree artefact.

MVA has no drug edge in PrimeKG, so TxGNN scores it purely from its neighbours.
On such a disease, a link predictor can rank drugs by how connected they are
rather than by anything about the disease. Two checks:

  * over all 7,957 scored drugs, Spearman rho between score and log degree;
  * inside the drugs TxGNN itself treats as clinically annotated for this
    relation (its "Ranked List"), the same rho, plus a residual rank after a
    linear fit of score on log degree.

If the annotated-list rho is near 0, the ranking inside that list is not a
degree ranking and its order can be read. The probe drugs are the report's
hypotheses and its contraindicated checkpoint-weakeners.

Reads the raw dumps written by 01. The analysis follows texdata's
02_txgnn_rank.py (MVA Hackathon 2026); the implementation is ours.
"""
from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import mvakg

PROBES = {
    "NAD+ / SIRT2 (stabilise BubR1)": ("Nicotinamide", "Nicotinamide riboside",
                                       "Nicotinamide mononucleotide", "NADH",
                                       "Niacin", "Nicotinic acid"),
    "senolytic": ("Dasatinib", "Quercetin", "Navitoclax", "Fisetin"),
    "aneuploidy-selective stress": ("Metformin", "Chloroquine", "Hydroxychloroquine",
                                    "Tanespimycin", "Geldanamycin", "Acadesine"),
    "contraindicated (weakens the checkpoint)": ("Paclitaxel", "Docetaxel", "Vincristine",
                                                 "Vinblastine", "Eribulin", "Volasertib",
                                                 "Alisertib", "Barasertib"),
}
TOP = 15


def load_scores(raw_path: Path, names: dict[str, str]) -> tuple[pd.DataFrame, list[str]]:
    raw = pd.read_pickle(raw_path)
    key = mvakg.MVA_PRIMEKG_ID
    scores = pd.DataFrame({"score": pd.Series(raw["prediction"][key], dtype=float),
                           "label": pd.Series(raw["label"][key])})
    if scores["label"].any():
        raise SystemExit(f"{raw_path.name}: MVA has a labelled drug, so this is not zero-shot")
    scores["drug"] = scores.index.map(names)
    return scores.rename_axis("drug_id").reset_index(), raw["result"]["Ranked List"][key]


def with_ranks(scores: pd.DataFrame) -> tuple[pd.DataFrame, float, float]:
    x = np.log10(scores["degree"] + 1)
    slope, intercept = np.polyfit(x, scores["score"], 1)
    rho = spearmanr(x, scores["score"]).statistic
    ranked = scores.assign(residual=scores["score"] - (slope * x + intercept),
                           rank=scores["score"].rank(ascending=False, method="first").astype(int))
    ranked["residual_rank"] = ranked["residual"].rank(ascending=False, method="first").astype(int)
    return ranked.sort_values("rank"), rho, slope


def print_probes(ranked: pd.DataFrame) -> None:
    by_name = ranked.set_index("drug")
    n = len(ranked)
    for group, drugs in PROBES.items():
        print(f"  {group}")
        missing = [d for d in drugs if d not in by_name.index]
        for d in drugs:
            if d in by_name.index:
                r = by_name.loc[d]
                print(f"      {d:28s} rank {r['rank']:5d} of {n} ({r['rank'] / n:6.1%})"
                      f"  residual rank {r['residual_rank']:5d}  degree {r['degree']:5d}")
        if missing:
            print(f"      not in the annotated list: {', '.join(missing)}")


def analyse(relation: str, raw_path: Path, names, degree, out_dir: Path) -> None:
    scores, annotated = load_scores(raw_path, names)
    scores = scores.merge(degree, on="drug_id", how="left").fillna({"degree": 0})
    scores["degree"] = scores["degree"].astype(int)

    rho_all = spearmanr(np.log10(scores["degree"] + 1), scores["score"]).statistic
    ranked, rho, slope = with_ranks(scores[scores["drug"].isin(annotated)])
    print(f"\n===== {relation}")
    print(f"all {len(scores):,} scored drugs: rho(score, log10 degree) = {rho_all:+.3f}")
    print(f"  top {TOP} over all scored drugs:")
    for _, r in scores.nlargest(TOP, "score").iterrows():
        print(f"      {r['drug']:32s} score {r['score']:+7.3f}  degree {r['degree']:5d}")
    print(f"{len(ranked):,} annotated drugs: rho = {rho:+.3f}, slope = {slope:+.3f}")
    print(f"  top {TOP} annotated by score:")
    for _, r in ranked.head(TOP).iterrows():
        print(f"    {r['rank']:4d}  {r['drug']:32s} score {r['score']:+7.3f}"
              f"  residual rank {r['residual_rank']:4d}  degree {r['degree']:5d}")
    print_probes(ranked)

    out = out_dir / f"ranked_{relation}.tsv"
    ranked[["rank", "drug_id", "drug", "score", "degree", "residual", "residual_rank"]] \
        .to_csv(out, sep="\t", index=False, float_format="%.5f")
    print(f"  wrote {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--work", type=Path, default=mvakg.WORK / "txgnn")
    ap.add_argument("--ckpt", type=Path, default=mvakg.TXGNN_CKPT)
    ap.add_argument("--kg", type=Path, default=mvakg.PRIMEKG_CSV)
    args = ap.parse_args()

    with open(args.ckpt / "name_mapping.pkl", "rb") as fh:
        names = pickle.load(fh)["id2name_drug"]
    degree = mvakg.primekg_drug_degree(args.kg).to_pandas()
    for relation in ("indication", "contraindication"):
        analyse(relation, args.work / f"raw_{relation}.pkl", names, degree, args.work)


if __name__ == "__main__":
    main()
