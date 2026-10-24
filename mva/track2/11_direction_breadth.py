"""Is the 42:1 asymmetry systematic, or did we pick a convenient panel?

10_direction_audit.py measures the direction of medicinal chemistry globally
and across a 16-target panel chosen for this disease. A panel chosen by the
authors of a report is exactly the thing a reviewer should distrust, so this
repeats the measurement on targets selected without reference to MVA: the
first 600 human SINGLE PROTEIN targets the ChEMBL API returns, keeping the
200 that carry at least 50 assays - 117,870 assays in total.

Result (2026-10-22): the MEDIAN human protein target has ZERO assays
describing activation. 117 of 200 (58.5%) have none at all, 164 (82.0%) are
under 1%, and the pooled ratio is 63.4:1 - steeper than the global figure.

The activation-rich tail is the thesis rather than an exception to it:
Nav1.1 (18.7%), leukotriene A-4 hydrolase (16.4%), caspase-3 (15.1%),
RXR-alpha (8.3%), vitamin D3 receptor (6.5%), PPAR-gamma (5.0%), glucagon
receptor (3.5%). Ion channels, nuclear receptors and a zymogen - proteins
with a native on-switch chemistry can reach. A mitotic checkpoint scaffold
has no such switch.

LIMITATION, stated rather than worked around: this is the first page of an
unsorted API listing, not a random draw. It is a convenience sample, and the
report calls it one. Like 10_direction_audit.py it counts STATED direction
from curator-written free text, and uses prefix matching so that activation
counts are conservative.

Read-only. No API key, no challenge data. Run: python3 -I this_script.py
"""
import json
import time
import urllib.request

UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}
A = "https://www.ebi.ac.uk/chembl/api/data/assay.json"
T = "https://www.ebi.ac.uk/chembl/api/data/target.json"


def g(u, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(u, headers=UA)
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception:
            if i == tries - 1:
                return None
            time.sleep(2 * (i + 1))
    return None


def cnt(u):
    d = g(u)
    return None if d is None else d.get("page_meta", {}).get("total_count")


# collect human SINGLE PROTEIN targets
targets = []
offset = 0
while len(targets) < 600 and offset < 2000:
    d = g(f"{T}?target_type=SINGLE%20PROTEIN&organism=Homo%20sapiens"
          f"&limit=200&offset={offset}")
    if not d:
        break
    batch = d.get("targets", [])
    if not batch:
        break
    for t in batch:
        targets.append((t["target_chembl_id"], t.get("pref_name")))
    offset += 200
    time.sleep(0.25)

print(f"collected {len(targets)} human SINGLE PROTEIN targets", flush=True)

rows = []
checked = 0
for tid, name in targets:
    tot = cnt(f"{A}?target_chembl_id={tid}&limit=1")
    time.sleep(0.12)
    checked += 1
    if not tot or tot < 50:          # only well-studied targets
        continue
    act = cnt(f"{A}?target_chembl_id={tid}&description__istartswith=Activation&limit=1")
    time.sleep(0.12)
    inh = cnt(f"{A}?target_chembl_id={tid}&description__istartswith=Inhibition&limit=1")
    time.sleep(0.12)
    rows.append({"target": tid, "name": name, "assays": tot,
                 "activation": act, "inhibition": inh})
    if len(rows) % 20 == 0:
        print(f"  {len(rows)} qualifying targets ({checked} checked)", flush=True)
    if len(rows) >= 200:
        break

with open("broad_sample.json", "w") as fh:
    json.dump(rows, fh, indent=1)

n = len(rows)
zero = sum(1 for r in rows if not r["activation"])
under1 = sum(1 for r in rows
             if r["assays"] and (r["activation"] or 0) / r["assays"] < 0.01)
tot_a = sum(r["activation"] or 0 for r in rows)
tot_i = sum(r["inhibition"] or 0 for r in rows)

print(f"\n=== {n} human targets with >=50 assays ===")
print(f"  targets with ZERO activation assays : {zero} ({100*zero/n:.1f}%)")
print(f"  targets with <1% activation assays  : {under1} ({100*under1/n:.1f}%)")
print(f"  pooled activation assays            : {tot_a}")
print(f"  pooled inhibition assays            : {tot_i}")
if tot_a:
    print(f"  pooled inhibition:activation ratio  : {tot_i/tot_a:.1f} : 1")

rows.sort(key=lambda r: -((r["activation"] or 0) / r["assays"]))
print("\n  most activation-rich targets:")
for r in rows[:12]:
    pct = 100 * (r["activation"] or 0) / r["assays"]
    print(f"    {pct:5.1f}%  {r['activation']:>5}/{r['assays']:<6} {r['name']}")
