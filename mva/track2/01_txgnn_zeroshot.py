"""Score every drug against MVA with the released TxGNN checkpoint, zero-shot.

MVA (PrimeKG node 13582_9759_54736_141) has no indication or contraindication
edge, so whatever TxGNN predicts for it comes from the disease's neighbours
(its genes, phenotypes and related diseases), not from a known treatment.
This script only runs the model and saves its raw output; 02 ranks it and 03
explains it.

The check we wanted the output to pass: drugs that weaken the spindle
checkpoint (taxanes, vinca alkaloids, eribulin, PLK1/Aurora/TTK inhibitors)
belong on the contraindication side for a checkpoint-deficiency disease. If
TxGNN ranks them high for *indication*, its MVA prediction is following
"cancer-like neighbourhood", not mechanism.

TxGNN (Huang et al., Nat Med 2024) is used through its own API. The setup
follows texdata's 01_txgnn_zeroshot.py (MVA Hackathon 2026); the code is ours.
"""
from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import pandas as pd

import mvakg

RELATIONS = ("indication", "contraindication")
SPLIT_SEED = 1


def restore_dataframe_append() -> None:
    """TxGNN calls DataFrame.append, which pandas 2 removed; give it pd.concat."""
    if hasattr(pd.DataFrame, "append"):
        return

    def append(self, other, ignore_index=False, verify_integrity=False, sort=False):
        if isinstance(other, (dict, pd.Series)):
            other = pd.DataFrame([other])
        return pd.concat([self, other], ignore_index=ignore_index,
                         verify_integrity=verify_integrity, sort=sort)

    pd.DataFrame.append = append


def disease_index(edges: pd.DataFrame, node_id: str) -> float:
    """TxGNN's internal index for a PrimeKG disease id, read from its edge table."""
    for side in ("x", "y"):
        hit = edges.loc[(edges[f"{side}_type"] == "disease") & (edges[f"{side}_id"] == node_id),
                        f"{side}_idx"]
        if len(hit):
            return hit.iloc[0]
    raise SystemExit(f"{node_id} is not a disease node in TxGNN's graph")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", type=Path, default=mvakg.TXGNN_DATA)
    ap.add_argument("--ckpt", type=Path, default=mvakg.TXGNN_CKPT)
    ap.add_argument("--out", type=Path, default=mvakg.WORK / "txgnn")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    restore_dataframe_append()
    from txgnn import TxData, TxEval, TxGNN  # imported after the pandas patch

    print("reading the PrimeKG split", flush=True)
    data = TxData(data_folder_path=str(args.data))
    data.prepare_split(split="full_graph", seed=SPLIT_SEED)

    print("reading the TxGNN checkpoint", flush=True)
    model = TxGNN(data=data, weight_bias_track=False, proj_name="mva", exp_name="mva",
                  device="cpu")
    model.load_pretrained(str(args.ckpt))
    model.best_model = model.model  # TxEval scores with best_model

    idx = disease_index(data.df, mvakg.MVA_PRIMEKG_ID)
    print(f"MVA ({mvakg.MVA_PRIMEKG_ID}) is disease index {idx}", flush=True)

    evaluator = TxEval(model=model)
    for relation in RELATIONS:
        raw = evaluator.eval_disease_centric(disease_idxs=[idx], relation=relation,
                                             save_result=False, verbose=False, return_raw=True)
        scores = raw["prediction"][mvakg.MVA_PRIMEKG_ID]
        labelled = sum(raw["label"][mvakg.MVA_PRIMEKG_ID].values())
        out = args.out / f"raw_{relation}.pkl"
        with open(out, "wb") as fh:
            pickle.dump(raw, fh)
        print(f"{relation}: {len(scores):,} drugs scored, {labelled} with a known "
              f"{relation} edge to MVA; wrote {out.name}", flush=True)


if __name__ == "__main__":
    main()
