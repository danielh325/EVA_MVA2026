"""Does reading the record instead of the title CHANGE THE REPURPOSING DECISION?

Scripts 15 and 16 answered narrower questions. 15 counted assays whose title
direction does not describe every record filed under them (a flag rate). 16
asked whether the measured endpoint is the target itself or a downstream
proxy (an interpretability audit). Neither answers the question that actually
matters to a repurposing pipeline:

    If I am searching for compounds that RAISE OR STABILISE a protein, how
    often does reading the record change whether a compound is kept or
    discarded -- and in which direction is the mistake?

That is the question this script answers, because it is the one our own
BUB1B finding is an instance of. ChEMBL assay CHEMBL6125447 is titled
"Destabilization of BUB1B ... assessed as reduction in BUB1B protein level".
A stabilisation-seeking search reading titles DISCARDS it. Two of the four
records filed under it show the protein becoming MORE thermally stable. The
compounds we built candidate 7 on were discarded by our own first pass.

THE DECISION MODEL, prespecified.

We simulate one concrete, realistic query -- the one this dossier ran:
"find compounds that move target X UP". Two pipelines answer it.

  TITLE PIPELINE   reads the assay description only.
      leading word is an UP word   -> KEEP
      leading word is a DOWN word  -> DISCARD
      no direction                 -> ABSTAIN

  RECORD PIPELINE  reads the activity record: the sign of the fold-change
      value against the reference, outside a noise band.
      value > 1 + NOISE            -> KEEP
      value < 1 - NOISE            -> DISCARD
      inside the band              -> ABSTAIN

A DECISION CHANGE is a record where the two pipelines disagree:

  FALSE NEGATIVE (the costly one) title DISCARD, record KEEP.
      A compound that moves the protein the way the patient needs, thrown
      away unread. This is the BUB1B/TG-101209 case exactly.
  FALSE POSITIVE                  title KEEP, record DISCARD.
      A compound promoted into a shortlist that the data does not support.

RESTRICTED TO CLASS A, and this restriction is what makes the comparison
fair. Using script 16's prespecified rule, Class A is "the measured endpoint
IS the target itself"; Class B is "the endpoint is a downstream proxy". In
Class B a title and a value may CORRECTLY point opposite ways -- inhibiting
MAGL raises 2-AG -- so a disagreement there is not an error and is reported
separately, never as a decision change.

HELD-OUT EVALUATION. The rules above were written against the sample drawn
with SEED_EXPLORE in scripts 15 and 16. They are frozen here and re-run
unchanged on a DIFFERENT sample drawn with SEED_HELDOUT. Both are reported.
If the held-out numbers match the exploratory ones, the rules were not tuned
to a sample; if they diverge, that is reported too, not hidden.

WHAT THIS IS NOT. This is RULE-BASED INTERNAL ADJUDICATION, not independent
review: the rules were written by the same system they are testing, and the
per-record sheet is written to disk precisely so a third party can disagree.
It is not a pharmacological error rate for ChEMBL. It covers fold-change
records only; FC conventions are not uniform (a minority report "fold
reduction", where a larger number means a larger DECREASE), so the false-
negative figure is an UPPER bound on genuine decision change. Abstentions are
coverage gaps, not errors, and are counted separately throughout.

Read-only. No API key, no challenge data. Run: python3 -I this_script.py
"""
import collections
import csv
import json
import random
import re
import sys
import time
import urllib.request

UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}
ACT = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

SEED_EXPLORE = 20261022      # the sample scripts 15/16 developed the rules on
SEED_HELDOUT = 20261024      # rules frozen, re-run unchanged on this one
PAGES = 24
NOISE = 0.20                 # values within +-20% of reference are no direction

DOWN = ("inhibition", "destabilization", "destabilisation", "decrease",
        "reduction", "downregulation", "suppression", "displacement",
        "antagonis", "degradation", "depletion", "blockade")
UP = ("activation", "stabilization", "stabilisation", "increase",
      "upregulation", "induction", "potentiation", "enhancement",
      "agonis", "stimulation")

STOP = {"the", "a", "an", "human", "rat", "mouse", "recombinant", "and",
        "in", "at", "by", "of", "cells", "cell", "expressed", "assessed"}

SYSTEM_CUT = re.compile(r"\b(in|from|expressed|cultured|transfected|"
                        r"incubated|treated|induced|by|using|after|with)\b")


def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.load(r)
        except Exception:
            if i == tries - 1:
                return None
            time.sleep(2 * (i + 1))
    return None


# ---- the frozen rules, lifted unchanged from scripts 15 and 16 ----------
def title_direction(desc):
    w = re.split(r"[^a-z]+", desc.lower().strip())[0]
    if w.startswith(DOWN):
        return "down"
    if w.startswith(UP):
        return "up"
    return None


def target_tokens(desc):
    head = desc.lower().split("assessed as")[0]
    m = re.search(r"\bof\s+(.+)", head) or re.search(r"\bat\s+(.+)", head)
    if not m:
        return set()
    frag = m.group(1)
    cut = SYSTEM_CUT.search(frag)
    if cut:
        frag = frag[:cut.start()]
    toks = re.split(r"[^a-z0-9]+", frag)
    return {t for t in toks if len(t) > 2 and t not in STOP}


def endpoint_class(desc):
    low = desc.lower()
    i = low.find("assessed as")
    if i < 0:
        return "no-endpoint-clause"
    ep = low[i + len("assessed as"):]
    tgt = target_tokens(desc)
    if not tgt:
        return "ambiguous-no-target-token"
    return "A" if (tgt & set(re.split(r"[^a-z0-9]+", ep))) else "B"


# ---- the decision model ------------------------------------------------
def title_decision(desc):
    """What a title-reading pipeline does with this record, for an
    'find compounds that move the target UP' query."""
    d = title_direction(desc)
    return {"up": "KEEP", "down": "DISCARD"}.get(d, "ABSTAIN")


def record_decision(value):
    if value > 1 + NOISE:
        return "KEEP"
    if value < 1 - NOISE:
        return "DISCARD"
    return "ABSTAIN"


def fetch(seed, label):
    total = (get(f"{ACT}?standard_type=FC&limit=1") or {}).get(
        "page_meta", {}).get("total_count", 0)
    rng = random.Random(seed)
    offsets = sorted(rng.sample(range(0, max(total - 1000, 1000), 1000), PAGES))
    recs = []
    for n, off in enumerate(offsets, 1):
        d = get(f"{ACT}?standard_type=FC&limit=1000&offset={off}")
        if not d:
            print(f"  [{label}] page {n}/{len(offsets)} FAILED - not counted "
                  f"as zero", flush=True)
            continue
        for a in d.get("activities", []):
            v, desc = a.get("standard_value"), a.get("assay_description") or ""
            if v is None or not desc:
                continue
            try:
                recs.append({"assay": a.get("assay_chembl_id"),
                             "molecule": a.get("molecule_chembl_id"),
                             "value": float(v), "desc": desc,
                             "comment": a.get("activity_comment") or ""})
            except (TypeError, ValueError):
                continue
        print(f"  [{label}] page {n}/{len(offsets)} · {len(recs)} records",
              flush=True)
        time.sleep(0.1)
    return total, recs


def evaluate(recs, label):
    st = collections.Counter()
    rows = []
    for r in recs:
        td = title_decision(r["desc"])
        st["records"] += 1
        if td == "ABSTAIN":
            st["title abstains (no direction stated)"] += 1
            continue
        st["title states a decision"] += 1

        cls = endpoint_class(r["desc"])
        st[f"class {cls}"] += 1
        if cls != "A":
            # endpoint is a proxy or unidentifiable - a disagreement here is
            # not necessarily an error, so it is NOT counted as one
            st["not adjudicable (proxy or unidentifiable endpoint)"] += 1
            continue
        st["adjudicable (endpoint IS the target)"] += 1

        rd = record_decision(r["value"])
        if rd == "ABSTAIN":
            st["record abstains (inside noise band)"] += 1
            continue
        st["both pipelines decide"] += 1

        if td == rd:
            verdict = "AGREE"
        elif td == "DISCARD" and rd == "KEEP":
            verdict = "FALSE NEGATIVE - title discards what the record keeps"
        else:
            verdict = "FALSE POSITIVE - title keeps what the record discards"
        st[verdict] += 1
        rows.append({"set": label, "assay": r["assay"],
                     "molecule": r["molecule"], "value": r["value"],
                     "title_decision": td, "record_decision": rd,
                     "endpoint_class": cls, "verdict": verdict,
                     "activity_comment": r["comment"],
                     "description": r["desc"]})
    return st, rows


def report(st, label):
    n = st["records"]
    dec = st["both pipelines decide"]
    fn = st["FALSE NEGATIVE - title discards what the record keeps"]
    fp = st["FALSE POSITIVE - title keeps what the record discards"]
    ag = st["AGREE"]
    print(f"\n=== {label}: {n:,} fold-change records ===")
    for k in ("title abstains (no direction stated)",
              "title states a decision",
              "not adjudicable (proxy or unidentifiable endpoint)",
              "adjudicable (endpoint IS the target)",
              "record abstains (inside noise band)",
              "both pipelines decide"):
        if st[k]:
            print(f"  {st[k]:>6}  ({100*st[k]/n:>4.1f}%)  {k}")
    if not dec:
        print("  no adjudicable records - nothing to report")
        return None
    print(f"\n  Of the {dec} records where BOTH pipelines reach a decision:")
    print(f"    {ag:>6}  ({100*ag/dec:>4.1f}%)  agree")
    print(f"    {fn:>6}  ({100*fn/dec:>4.1f}%)  FALSE NEGATIVE - discarded "
          f"unread, record moves the wanted way")
    print(f"    {fp:>6}  ({100*fp/dec:>4.1f}%)  FALSE POSITIVE - kept on a "
          f"title the record does not support")
    print(f"\n  >> Record-level reading CHANGES THE DECISION for "
          f"{fn+fp} of {dec} ({100*(fn+fp)/dec:.1f}%).")
    # the number that matters to a repurposing search
    disc = st["AGREE"] + fn  # agreements where title said DISCARD are unknown
    print(f"  >> Of records a title-reading search DISCARDS, the record shows "
          f"{fn} moving the wanted way.")
    return {"decided": dec, "agree": ag, "false_negative": fn,
            "false_positive": fp,
            "decision_change_pct": round(100*(fn+fp)/dec, 2)}


print(__doc__.split("Read-only")[0][:0] or "", end="")
print("Decision-impact benchmark: does reading the record change the call?\n")

out = {"noise_band": NOISE, "pages": PAGES, "sets": {}}
all_rows = []
for seed, label in ((SEED_EXPLORE, "exploratory (rules written here)"),
                    (SEED_HELDOUT, "HELD-OUT (rules frozen)")):
    total, recs = fetch(seed, label.split()[0])
    if not recs:
        print(f"  [{label}] no records retrieved - reporting as FAILURE, "
              f"not as zero", flush=True)
        out["sets"][label] = {"status": "retrieval failed"}
        continue
    st, rows = evaluate(recs, label)
    summary = report(st, label)
    out["sets"][label] = {"seed": seed, "fc_records_total": total,
                          "stats": dict(st), "summary": summary}
    all_rows.extend(rows)

# ---- do the two samples agree? the point of holding one out ------------
a = out["sets"].get("exploratory (rules written here)", {}).get("summary")
b = out["sets"].get("HELD-OUT (rules frozen)", {}).get("summary")
if a and b:
    d = abs(a["decision_change_pct"] - b["decision_change_pct"])
    print(f"\n=== held-out check ===")
    print(f"  decision-change rate  exploratory {a['decision_change_pct']}% "
          f"vs held-out {b['decision_change_pct']}%  (difference {d:.2f} pp)")
    print("  The rules were written on the exploratory sample and run "
          "unchanged on the\n  held-out one. A small difference means they "
          "were not tuned to a sample;\n  it does NOT make them correct - "
          "only stable.")
    out["held_out_difference_pp"] = round(d, 2)

if all_rows:
    with open("18_decision_impact_sample.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)
    print(f"\n  Per-record sheet: 18_decision_impact_sample.csv "
          f"({len(all_rows)} adjudicated records)")

with open("18_decision_impact_results.json", "w") as fh:
    json.dump(out, fh, indent=1)

# ---- the case this was built from --------------------------------------
print("\n=== the record this benchmark exists because of ===")
print("  CHEMBL6125447  'Destabilization of BUB1B ... assessed as reduction")
print("  in BUB1B protein level'  -> a title-reading search DISCARDS it.")
print("  Two of its four records are above the reference, both annotated")
print("  'Increase In Heat Resistance Observed'. Candidate 7 is built on")
print("  compounds our own first pass threw away. See 14_bub1b_cetsa_audit.py")
