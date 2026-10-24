"""What title-based direction classification can and cannot be blamed for.

Scripts 10-12 classify a bioassay's DIRECTION from the leading words of its
description. Script 14 found that method failing on this patient's gene:
ChEMBL assay CHEMBL6125447 is titled "Destabilization of BUB1B ... assessed
as reduction in BUB1B protein level ... (Rvb = 1)" and contains two
compounds whose values are ABOVE the reference - i.e. the protein became
more thermally stable, not less.

A FIRST VERSION OF THIS SCRIPT OVERCLAIMED, and the error is worth stating
because it is the same class of error the report is about.

It counted any assay whose leading word disagreed with a later "assessed as"
clause as a misclassification, and reported ~11% as an error rate. That is
wrong. Target action and assay endpoint CAN correctly move in opposite
directions:

    "Inhibition of MAGL ... assessed as increase in 2-AG level"
        -> MAGL degrades 2-AG. Inhibiting it RAISES 2-AG. Correct.
    "Inhibition of Pgp ... assessed as increase in doxorubicin accumulation"
        -> Blocking efflux RAISES intracellular drug. Correct.
    "Antagonist activity ... assessed as fold increase in EC50 of agonist"
        -> That is what antagonism does to an agonist's EC50. Correct.

None of those is a classification failure. They show that TARGET ACTIVITY
and MEASURED ENDPOINT are different variables, which is a real difficulty
for any wording-based classifier but is not an error rate.

WHAT THIS SCRIPT NOW MEASURES, which is narrower and defensible.

The BUB1B case is different in kind: the leading word ("Destabilization")
and the endpoint clause ("reduction in BUB1B protein level") name the SAME
property in the SAME direction, so the title is internally unambiguous - and
the recorded values still go both ways. That is a genuine contradiction
between a title and the data filed under it.

So the benchmark restricts to assays where:

  1. the leading word states a direction, AND
  2. the "assessed as" clause either agrees with it or names no direction
     (so the title is not describing an inverse target/endpoint relation), AND
  3. the assay has at least two fold-change records, AND
  4. those records fall on both sides of the reference by more than a noise
     band (default 20%), so near-1 jitter is not counted.

Everything else is reported separately and NOT as error:

  - no direction in the leading word          -> ABSTENTION, not error
  - leading word vs endpoint clause disagree  -> candidate inverse relation,
                                                 not adjudicated here

STILL NOT CLAIMED. This is not a pharmacological error rate for ChEMBL. It
covers fold-change records only; FC conventions are not uniform (a minority
report "fold reduction", where a larger number means a larger decrease), so
the residual figure remains an UPPER bound on genuine title/data conflict.
Establishing a true error rate needs hand-adjudicated labels with target
action, endpoint, and the sign of their relationship recorded separately -
which is the right next piece of work and is not done here.

Read-only. No API key, no challenge data. Run: python3 -I this_script.py
"""
import collections
import json
import random
import re
import time
import urllib.request

UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}
ACT = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
SEED = 20261022
PAGES = 24
NOISE = 0.20          # ignore values within +-20% of the reference

DOWN = ("inhibition", "destabilization", "destabilisation", "decrease",
        "reduction", "downregulation", "suppression", "displacement",
        "antagonis", "degradation", "depletion", "blockade")
UP = ("activation", "stabilization", "stabilisation", "increase",
      "upregulation", "induction", "potentiation", "enhancement",
      "agonis", "stimulation")

READOUT_UP = re.compile(r"\b(increase|induction|upregulation|accumulation|"
                        r"stabilization|stabilisation|enhancement|rise|potentiation|sensitization|sensitisation|restoration|recovery)\\b")
READOUT_DOWN = re.compile(r"\b(decrease|reduction|loss|depletion|"
                          r"downregulation|destabilization|destabilisation)\b")


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


def lead_direction(desc):
    w = re.split(r"[^a-z]+", desc.lower().strip())[0]
    if w.startswith(DOWN):
        return "down"
    if w.startswith(UP):
        return "up"
    return None


def endpoint_direction(desc):
    """Direction named in the 'assessed as' clause, if any."""
    low = desc.lower()
    cut = low.find("assessed as")
    if cut < 0:
        return None
    tail = low[cut:]
    up, dn = READOUT_UP.search(tail), READOUT_DOWN.search(tail)
    if up and not dn:
        return "up"
    if dn and not up:
        return "down"
    return None


total = (get(f"{ACT}?standard_type=FC&limit=1") or {}).get(
    "page_meta", {}).get("total_count", 0)
print(f"FC activity records in ChEMBL: {total:,}", flush=True)

rng = random.Random(SEED)
offsets = sorted(rng.sample(range(0, max(total - 1000, 1000), 1000),
                            min(PAGES, max(total // 1000, 1))))

by_assay = collections.defaultdict(list)
desc_of = {}
for n, off in enumerate(offsets, 1):
    d = get(f"{ACT}?standard_type=FC&limit=1000&offset={off}")
    if not d:
        continue
    for a in d.get("activities", []):
        aid, v = a.get("assay_chembl_id"), a.get("standard_value")
        if not aid or v is None:
            continue
        try:
            by_assay[aid].append(float(v))
        except (TypeError, ValueError):
            continue
        desc_of[aid] = a.get("assay_description") or ""
    print(f"  page {n}/{len(offsets)} · {len(by_assay)} assays", flush=True)
    time.sleep(0.1)

st = collections.Counter()
inverse, conflict = [], []

for aid, vals in by_assay.items():
    desc = desc_of[aid]
    lead = lead_direction(desc)
    st["assays"] += 1

    if lead is None:
        st["abstain: no direction in leading word"] += 1
        continue
    st["leading word states a direction"] += 1

    ep = endpoint_direction(desc)
    if ep and ep != lead:
        # target action and endpoint point opposite ways - frequently CORRECT
        st["title names an inverse target/endpoint relation"] += 1
        inverse.append((aid, lead, ep, desc[:110]))
        continue

    st["title internally unambiguous"] += 1
    if len(vals) < 2:
        st["  (single record - cannot conflict)"] += 1
        continue

    hi = [v for v in vals if v > 1 + NOISE]
    lo = [v for v in vals if v < 1 - NOISE]
    if hi and lo:
        st["title's direction does not describe every record"] += 1
        conflict.append((aid, lead, len(hi), len(lo), desc[:110]))
        # a lone outlier is likely an inactive compound in a screen; two or
        # more on the minority side is a substantive split
        if min(len(hi), len(lo)) >= 2:
            st["  ...with >=2 records on each side"] += 1

with open("15_direction_classifier_benchmark_results.json", "w") as fh:
    json.dump({"seed": SEED, "noise_band": NOISE, "fc_records_total": total,
               "pages": len(offsets), "stats": dict(st),
               "inverse_examples": inverse[:100],
               "conflicts": conflict[:200]}, fh, indent=1)

n = st["assays"]
lead_n = st["leading word states a direction"]
unamb = st["title internally unambiguous"]
testable = unamb - st["  (single record - cannot conflict)"]
conf = st["title's direction does not describe every record"]
conf2 = st["  ...with >=2 records on each side"]

print(f"\n=== {n} assays with fold-change records ===")
print(f"  abstain (no direction in leading word) : {st['abstain: no direction in leading word']:>5}"
      f" ({100*st['abstain: no direction in leading word']/n:.1f}%)  <- coverage gap, NOT error")
print(f"  leading word states a direction        : {lead_n:>5}")
print(f"    of which title names an inverse")
print(f"    target/endpoint relation             : {st['title names an inverse target/endpoint relation']:>5}"
      f" ({100*st['title names an inverse target/endpoint relation']/lead_n:.1f}%)  <- often CORRECT biology")
print(f"    of which title is unambiguous        : {unamb:>5}")
print(f"      testable (>=2 records)             : {testable:>5}")
if testable:
    print(f"\n=== the defensible number ===")
    print(f"  title's stated direction does not describe every compound"
          f" filed under it: {conf} of {testable} ({100*conf/testable:.1f}%)")
    print(f"  ...with >=2 records on each side (not a lone outlier):"
          f" {conf2} ({100*conf2/testable:.1f}%)")
    print(f"  — upper bound: FC conventions are not uniform, and a minority")
    print(f"    side can be inactive compounds rather than opposite effects.")

print("\n  examples NOT counted as error (inverse target/endpoint, often correct):")
for aid, lead, ep, d in inverse[:5]:
    print(f"    [{lead}->{ep}] {d}")
print("\n  examples counted as conflict (title unambiguous, data both ways):")
for aid, lead, hi, lo, d in conflict[:6]:
    print(f"    [{lead}] {hi} above / {lo} below — {d}")
