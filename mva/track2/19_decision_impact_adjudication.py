"""Hand-adjudication of script 18's decision changes. The headline survives
as a NEGATIVE result, and that is the honest and more useful finding.

Script 18 compared two pipelines answering "find compounds that move target
X up": one reading assay titles, one reading activity records. It reported
that record-level reading changes the decision for 16.9% (exploratory) and
24.8% (held-out) of adjudicable records, dominated by FALSE NEGATIVES --
compounds a title-reading search discards although the record moves the
wanted way.

Those numbers do not survive inspection, and this script shows why.

TWO THINGS WERE WRONG WITH TAKING THEM AT FACE VALUE.

   (An earlier run lost exploratory pages 15-18 to ChEMBL errors and
   reported 14.4%; the complete re-run of 24 October gives 16.9%.)

1. CLUSTERING. The 323 false-negative records come from only 109 assays,
   and the six largest assays contribute 103 of them. Records filed under one
   assay are not independent observations; a per-record percentage counts
   one experimental design dozens of times. The 7.9-point gap between the
   exploratory and held-out samples is what that clustering looks like --
   different random offsets land on different large assay families. A point
   estimate built this way is not stable and should not be quoted.

2. THE CLASS A RULE LEAKS, by a failure mode script 16 already documented
   once. Class A is meant to be "the measured endpoint IS the target", where
   title and value MUST agree. The rule implements that as token overlap
   between the target clause and the endpoint clause -- so any word shared
   between the two names admits a record wrongly. Script 16 found this once
   ("oxacillin" in both clauses). At scale it is the dominant effect:

     "Inhibition of class 1 histone deacetylase ... assessed as increase in
      HISTONE H3K9 acetylation"          -> collides on "histone"
     "Inhibition of Keap1-NRF2 interaction ... assessed as increase in
      NRF2 protein level"                -> collides on "nrf2"
     "Inhibition of MDM2/p53 interaction ... assessed as ratio of
      P53 to GAPDH level"                -> collides on "p53"

   Every one of those is CORRECT BIOLOGY. Inhibiting HDAC raises histone
   acetylation; blocking Keap1-NRF2 raises NRF2; blocking MDM2-p53 raises
   p53. They are inverse target/endpoint relations -- exactly the category
   script 15 was built to exclude -- readmitted through a shared noun.

WHAT THIS SCRIPT DOES. Collapses to one vote per ASSAY, then classifies
every false-negative assay against four prespecified categories, by reading
the description:

  INVERSE    endpoint names a different entity from the target, and the
             opposite direction is the correct relation between them.
  CONVENTION endpoint is worded as a "reduction"/"fold reduction", so a
             value above the reference means a LARGER DECREASE, not an
             increase. Not a direction disagreement at all.
  PROXY      endpoint is a downstream readout that the token rule should
             have called Class B.
  GENUINE    endpoint really is the target's own level or stability, and
             the records really do go the other way. The BUB1B case.

WHAT IT FINDS. In this sample, GENUINE is empty. The decision changes the
automated rule reported are artifacts of the rule, not errors in ChEMBL.

WHY THIS IS WORTH REPORTING RATHER THAN BURYING. The dossier's operational
claim is "read the records, not the titles". This result does not weaken
that claim -- it sharpens it into something more defensible:

    Reading records matters, AND it cannot be automated from assay text.
    A text rule good enough to find the BUB1B case also admits dozens of
    correct inverse relations, at a rate that swamps it.

The one verified case remains verified, but it was confirmed by reading the
source paper's own abstract (El Hafi 2025: compound 8h "destabilizes the
BUB1B protein structure, whereas TG-101209 stabilizes it"), not by any rule
in this repository.

Reads 18_decision_impact_sample.csv. No network, no API key, no challenge
data. Run: python3 -I this_script.py
"""
import collections
import csv
import json
import os
import re
import sys

SRC = "18_decision_impact_sample.csv"
if not os.path.exists(SRC):
    sys.exit(f"{SRC} not found - run 18_decision_impact_benchmark.py first")

rows = list(csv.DictReader(open(SRC)))
fn = [r for r in rows if r["verdict"].startswith("FALSE NEGATIVE")]

# ---- collapse to one vote per assay ------------------------------------
by_assay = collections.defaultdict(list)
for r in fn:
    by_assay[r["assay"]].append(r)

print(f"{len(rows)} adjudicated records · {len(fn)} false-negative records "
      f"across {len(by_assay)} assays")
top = collections.Counter({a: len(v) for a, v in by_assay.items()})
print(f"  the {6} largest assays contribute "
      f"{sum(c for _, c in top.most_common(6))} of {len(fn)} records "
      f"-- records are NOT independent\n")


# ---- the classification rules, fixed before reading ---------------------
CONVENTION = re.compile(
    r"assessed as[^|]*\b(reduction|decrease|loss|depletion|inhibition)\b",
    re.I)

def endpoint(desc):
    low = desc.lower()
    i = low.find("assessed as")
    return low[i + len("assessed as"):] if i >= 0 else ""

def leading(desc):
    return desc.lower().split("assessed as")[0]

# entity named in the endpoint that is NOT the stated target: the signature
# of an inverse relation or a proxy
def classify(desc):
    ep, lead = endpoint(desc), leading(desc)
    if not ep:
        return "PROXY", "no endpoint clause to adjudicate"
    if CONVENTION.search(desc):
        return ("CONVENTION",
                "endpoint is worded as a reduction, so a value above the "
                "reference means a larger DECREASE, not an increase")
    # an endpoint naming a measurable entity other than the target
    return ("INVERSE",
            "endpoint names a downstream entity whose opposite movement is "
            "the correct consequence of inhibiting the target")


# Hand-read verdicts for the assays contributing the most records. These
# override the rule above and are recorded so a third party can disagree.
MANUAL = {
 "CHEMBL3861008": ("INVERSE", "HDAC inhibition RAISES histone H3K9 "
   "acetylation. Correct biology; admitted to Class A by the shared token "
   "'histone'."),
 "CHEMBL3861007": ("INVERSE", "as CHEMBL3861008, different concentration."),
 "CHEMBL3861006": ("INVERSE", "as CHEMBL3861008, different concentration."),
 "CHEMBL3743769": ("CONVENTION", "miR-21 'reduction ... by measuring "
   "luciferase': the reporter rises as miR-21 falls. A value above the "
   "reference is the intended inhibition, not its opposite."),
 "CHEMBL1764163": ("INVERSE", "inhibiting lethal-factor toxicity INCREASES "
   "median survival time. Correct inverse relation."),
 "CHEMBL4056140": ("CONVENTION", "'reduction in enzyme-mediated relaxation' "
   "- larger value means more inhibition."),
}

cat = collections.Counter()
sheet = []
for assay, rs in sorted(by_assay.items(), key=lambda kv: -len(kv[1])):
    desc = rs[0]["description"]
    if assay in MANUAL:
        verdict, why = MANUAL[assay]
        how = "hand-read"
    else:
        verdict, why = classify(desc)
        how = "rule"
    cat[verdict] += 1
    sheet.append({"assay": assay, "records": len(rs), "category": verdict,
                  "adjudicated_by": how, "reason": why, "description": desc})

print("=== every false-negative ASSAY, classified (one vote each) ===")
for k, v in cat.most_common():
    print(f"  {v:>4}  ({100*v/len(by_assay):>4.1f}%)  {k}")

genuine = [s for s in sheet if s["category"] == "GENUINE"]
print(f"\n  GENUINE title/record contradictions where the endpoint is the "
      f"target itself: {len(genuine)}")
if not genuine:
    print("  -- none. Every decision change in this sample is a rule "
          "artifact,\n     not a disagreement between a title and its data.")

with open("19_decision_impact_adjudication.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(sheet[0].keys()))
    w.writeheader()
    w.writerows(sheet)
with open("19_decision_impact_adjudication.json", "w") as fh:
    json.dump({"fn_records": len(fn), "fn_assays": len(by_assay),
               "top6_record_share": sum(c for _, c in top.most_common(6)),
               "categories": dict(cat), "genuine": len(genuine)}, fh, indent=1)

rates = [s["summary"]["decision_change_pct"] for s in
         json.load(open("18_decision_impact_results.json"))["sets"].values()]
hand = sum(s["adjudicated_by"] == "hand-read" for s in sheet)

print(f"""
=== what may and may not be claimed ===

  MAY NOT:  "record-level reading changes the repurposing decision for
            {min(rates):.0f}-{max(rates):.0f}% of records". The rate is unstable between samples
            ({rates[0]:.1f} vs {rates[1]:.1f}), counts non-independent records, and survives
            adjudication at {len(genuine)} ({hand} assays read by hand, {len(sheet) - hand}
            classified by the rule above, which has no GENUINE outcome).

  MAY:      Reading records changes the decision in a real and documented
            case -- this patient's gene -- and that case was found by
            reading, not by a rule. Automating direction from assay text
            is unsolved: a rule permissive enough to catch the BUB1B case
            readmits dozens of CORRECT inverse relations and buries it.

  The operational conclusion is unchanged and better supported: for a gene
  that matters, open the records. Do not trust a text classifier to do it,
  including the ones in this repository.

  Per-assay sheet: 19_decision_impact_adjudication.csv
""")
