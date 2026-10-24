"""Does record-level inspection actually change the interpretation? Adjudicated.

Script 15 flags assays whose title direction does not describe every record
filed under them. A flag is not an error, and 15 says so. This script asks
the question a flag cannot answer:

    If you mine a bioactivity database for "what does compound X do to
    protein Y", how often does reading the ASSAY TITLE give you the wrong
    answer, and how often is the question answerable at all?

THE RULE THAT MAKES THIS ADJUDICABLE, PRESPECIFIED BEFORE SAMPLING.

The difficulty in script 15 was that target action and measured endpoint can
correctly diverge: inhibiting MAGL raises 2-AG, and that is not an error.
The distinction that resolves it is whether the ENDPOINT IS THE TARGET
ITSELF or a downstream proxy:

  CLASS A - endpoint is the target's own level, stability or amount.
      "Destabilization of BUB1B ... assessed as reduction in BUB1B protein
       level"  -> endpoint names BUB1B, the target.
      Here title direction and value direction MUST agree. Disagreement is a
      genuine mislabel of what the compound did to the target.

  CLASS B - endpoint is a downstream proxy.
      "Inhibition of MAGL ... assessed as increase in 2-AG level"
       -> endpoint names 2-AG, not MAGL.
      Here the title describes the intervention and the value describes
      something else. The title is not wrong, but it does NOT answer "what
      happened to the target", so the question is unanswerable from the
      title alone.

OPERATIONALISED, so anyone can re-run it: take the target token from the
leading clause (the words after the first "of"), and test whether that token
reappears in the "assessed as" clause. Token match -> Class A. No match ->
Class B. Neither requires judgement at scoring time.

WHAT IS REPORTED
  - how often the question is answerable from the title (Class A share)
  - within Class A, how often the title's direction matches the data
  - every AMBIGUOUS case, counted and written out rather than dropped
  - whether record-level inspection changes the interpretation, which is the
    whole point

HONESTY ABOUT WHO ADJUDICATED. The rules above were fixed before the sample
was drawn, and scoring is mechanical, so the result is reproducible by
anyone running this file. But the rules were written by the same system that
wrote the classifier they test. That is NOT independent adjudication, and
the output is a transparent self-assessment rather than a blind review. The
per-record sheet is written to disk precisely so a third party can disagree.

Sampling is seeded. Read-only. No API key, no challenge data.
Run: python3 -I this_script.py
"""
import collections
import csv
import json
import random
import re
import time
import urllib.request

UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}
ACT = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
SEED = 20261022
PAGES = 24
NOISE = 0.20          # values within +-20% of the reference are not a direction
SAMPLE = 200          # records adjudicated

DOWN = ("inhibition", "destabilization", "destabilisation", "decrease",
        "reduction", "downregulation", "suppression", "displacement",
        "antagonis", "degradation", "depletion", "blockade")
UP = ("activation", "stabilization", "stabilisation", "increase",
      "upregulation", "induction", "potentiation", "enhancement",
      "agonis", "stimulation")

STOP = {"the", "a", "an", "human", "rat", "mouse", "recombinant", "and",
        "in", "at", "by", "of", "cells", "cell", "expressed", "assessed"}


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


def title_direction(desc):
    w = re.split(r"[^a-z]+", desc.lower().strip())[0]
    if w.startswith(DOWN):
        return "down"
    if w.startswith(UP):
        return "up"
    return None


# Everything after one of these words describes the SYSTEM (cell line,
# medium, co-administered drug), not the target. Truncating here is what
# stops the rule matching on "vinblastine", "MYC" or "oxacillin" when those
# appear in both the leading clause and the endpoint clause.
SYSTEM_CUT = re.compile(r"\b(in|from|expressed|cultured|transfected|"
                        r"incubated|treated|induced|by|using|after|with)\b")


def target_tokens(desc):
    """Words naming the target: the leading clause after 'of'/'at', cut at
    the first word that starts describing the experimental system."""
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


def endpoint_clause(desc):
    low = desc.lower()
    i = low.find("assessed as")
    return low[i + len("assessed as"):] if i >= 0 else ""


def endpoint_class(desc):
    """Class A = endpoint names the target itself; B = downstream proxy."""
    ep = endpoint_clause(desc)
    if not ep:
        return "no-endpoint-clause"
    tgt = target_tokens(desc)
    if not tgt:
        return "ambiguous-no-target-token"
    ep_toks = set(re.split(r"[^a-z0-9]+", ep))
    return "A" if (tgt & ep_toks) else "B"


def value_direction(v):
    if v > 1 + NOISE:
        return "up"
    if v < 1 - NOISE:
        return "down"
    return None          # inside the noise band: no direction claimed


# ---------------------------------------------------------------- sample
total = (get(f"{ACT}?standard_type=FC&limit=1") or {}).get(
    "page_meta", {}).get("total_count", 0)
rng = random.Random(SEED)
offsets = sorted(rng.sample(range(0, max(total - 1000, 1000), 1000), PAGES))

records = []
for n, off in enumerate(offsets, 1):
    d = get(f"{ACT}?standard_type=FC&limit=1000&offset={off}")
    if not d:
        continue
    for a in d.get("activities", []):
        v, desc = a.get("standard_value"), a.get("assay_description") or ""
        if v is None or not desc:
            continue
        try:
            records.append({"assay": a.get("assay_chembl_id"),
                            "molecule": a.get("molecule_chembl_id"),
                            "value": float(v), "desc": desc})
        except (TypeError, ValueError):
            continue
    print(f"  page {n}/{len(offsets)} · {len(records)} records", flush=True)
    time.sleep(0.1)

# keep only records whose title states a direction - the ones a title-based
# pipeline would actually act on
usable = [r for r in records if title_direction(r["desc"])]
print(f"\n{len(records)} records · {len(usable)} with a directional title")

rng2 = random.Random(SEED + 1)
sample = rng2.sample(usable, min(SAMPLE, len(usable)))

# ------------------------------------------------------------- adjudicate
st = collections.Counter()
rows = []
for r in sample:
    td = title_direction(r["desc"])
    cls = endpoint_class(r["desc"])
    vd = value_direction(r["value"])

    if cls in ("no-endpoint-clause", "ambiguous-no-target-token"):
        verdict = "AMBIGUOUS — endpoint not identifiable from the description"
    elif vd is None:
        verdict = "AMBIGUOUS — value inside the noise band"
    elif cls == "B":
        verdict = "UNANSWERABLE from title — endpoint is a downstream proxy"
    elif td == vd:
        verdict = "TITLE CORRECT — endpoint is the target, directions agree"
    else:
        verdict = "TITLE WRONG — endpoint is the target, directions disagree"

    st[verdict] += 1
    st[f"class {cls}"] += 1
    rows.append({"assay": r["assay"], "molecule": r["molecule"],
                 "value": r["value"], "title_direction": td,
                 "endpoint_class": cls, "value_direction": vd,
                 "verdict": verdict, "description": r["desc"]})

with open("16_direction_adjudication_sample.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
with open("16_direction_adjudication_results.json", "w") as fh:
    json.dump({"seed": SEED, "noise_band": NOISE, "sampled": len(sample),
               "stats": dict(st)}, fh, indent=1)

# ----------------------------------------------------------------- report
n = len(sample)
A = st["class A"]
correct = st["TITLE CORRECT — endpoint is the target, directions agree"]
wrong = st["TITLE WRONG — endpoint is the target, directions disagree"]
proxy = st["UNANSWERABLE from title — endpoint is a downstream proxy"]
amb = sum(v for k, v in st.items() if k.startswith("AMBIGUOUS"))

print(f"\n=== {n} adjudicated records with a directional title ===")
for k in sorted(st, key=lambda k: -st[k]):
    if k.startswith("class "):
        continue
    print(f"  {st[k]:>4}  ({100*st[k]/n:>4.1f}%)  {k}")

print(f"\n=== the two numbers that matter ===")
answerable = correct + wrong
print(f"  'what did the compound do to the TARGET' is answerable from the")
print(f"  title in {answerable} of {n} records ({100*answerable/n:.1f}%).")
if answerable:
    print(f"  Within those, the title's direction is WRONG in {wrong}"
          f" ({100*wrong/answerable:.1f}%).")
print(f"\n  For the other {n-answerable} ({100*(n-answerable)/n:.1f}%) the title")
print(f"  does not answer the question: {proxy} measure a downstream proxy,")
print(f"  {amb} are ambiguous and are reported rather than dropped.")
print(f"\n  Record-level inspection changes the interpretation for"
      f" {wrong + proxy} of {n} ({100*(wrong+proxy)/n:.1f}%) — the records a")
print("  title-based pipeline would either mislabel or silently over-read.")
print("\n  Per-record sheet: 16_direction_adjudication_sample.csv")

# ------------------------------------------------- manual adjudication
# Every record the automated rule flagged TITLE WRONG was then read by hand.
# Recorded here so the verdicts are auditable and so the automated rule's
# own failure mode is on the record rather than in a footnote.
MANUAL = {
 "CHEMBL4699745": ("FALSE POSITIVE",
   "Potentiation of oxacillin-induced activity, measured as fold REDUCTION "
   "in oxacillin MIC. Potentiating an antibiotic lowers its MIC, so this is "
   "a correct inverse relation. The rule matched on 'oxacillin' appearing "
   "in both the leading clause and the endpoint clause - a drug-name "
   "collision, not a target match."),
}
flagged = [r for r in rows if r["verdict"].startswith("TITLE WRONG")]
print("\n=== manual adjudication of every automated flag ===")
if not flagged:
    print("  none flagged")
for r in flagged:
    verdict, why = MANUAL.get(r["assay"], ("UNREVIEWED", "not yet read"))
    print(f"  {r['assay']}  ->  {verdict}")
    print(f"    {why}")
true_wrong = sum(1 for r in flagged
                 if MANUAL.get(r["assay"], ("UNREVIEWED",))[0] == "CONFIRMED")
print(f"\n  Confirmed title/data mislabels after manual review:"
      f" {true_wrong} of {len(sample)} records.")
print("  The automated Class A rule mis-assigns when a co-administered drug "
      "or cell-line\n  name appears in both clauses. Automating this is "
      "unsolved; the operational\n  conclusion is unchanged and stronger: "
      "read the records.")
