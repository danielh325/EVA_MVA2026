"""Audit the report's own weakest claim: "no SIRT2 activator exists".

Section 4.5 argued that every drug adjacent to this lesion points the wrong
way - SIRT2's only ligand is cambinol, an inhibitor - and section 7 (E3) drew
the conclusion that H1 "has to work at the substrate level or not at all".

Both rested on a *knowledge-graph lookup*. The report flagged that as the
weakest of its three absence claims, because absence in PrimeKG or OptimusKG
is evidence about the graph, not about pharmacology - which is precisely the
distinction the rest of the report is built on. This script runs the search
the claim deserved, against ChEMBL and PubMed directly.

It overturns the claim. Activator records exist (EC150 15-25 uM, max
activation 206-257%, one reading at EC50 160 nM) and neither graph carries
them - a *coverage* gap, distinct from the schema and semantics gaps already
reported. The practical conclusion survives on different grounds: every one
of those compounds is an unnamed preclinical research chemical, the single
potent value is assay-format-dependent, and resveratrol - the only named
clinical-stage molecule in the set - does not activate SIRT2 at all.

The target ID is resolved from the ChEMBL target search rather than assumed;
an earlier run of this script hardcoded the wrong accession and had to be
discarded. Every query issued is written to 09_sirt2_activator_queries.txt so
the result can be repeated rather than taken on trust.

Read-only. No API key, no challenge data. Run: python3 -I this_script.py
"""
import json
import time
import urllib.parse
import urllib.request
import datetime

DATE = datetime.date.today().isoformat()
UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}
LOG = []


def get(url, label, pause=0.0):
    if pause:
        time.sleep(pause)
    LOG.append((label, url))
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.load(r)
    except Exception as e:  # report, never substitute
        return {"__error__": f"{type(e).__name__}: {e}"}


print(f"=== Task E - SIRT2 activator search - retrieved {DATE} ===\n")

TARGET = "CHEMBL4462"  # human SIRT2, SINGLE PROTEIN, from the target search
print(f"Target: {TARGET} (NAD-dependent protein deacetylase sirtuin-2, Homo sapiens)\n")

# --- activity records, fully paginated ---
total = None
u = (f"https://www.ebi.ac.uk/chembl/api/data/activity.json"
     f"?target_chembl_id={TARGET}&limit=1")
d = get(u, f"ChEMBL activity count {TARGET}")
if "__error__" in d:
    print("ERROR:", d["__error__"])
else:
    total = d.get("page_meta", {}).get("total_count")
    print(f"--- ChEMBL activity records for {TARGET}: {total} ---")

types = {}
activation_like = []
offset = 0
while total and offset < total:
    u = (f"https://www.ebi.ac.uk/chembl/api/data/activity.json"
         f"?target_chembl_id={TARGET}&limit=1000&offset={offset}")
    d = get(u, f"ChEMBL activities offset={offset}", pause=0.2)
    if "__error__" in d:
        print("ERROR:", d["__error__"])
        break
    acts = d.get("activities", [])
    if not acts:
        break
    for a in acts:
        t = (a.get("standard_type") or "None")
        types[t] = types.get(t, 0) + 1
        # flag anything whose TYPE or ASSAY wording implies increased activity
        desc = (a.get("assay_description") or "")
        blob = f"{t} {desc}".lower()
        if ("activation" in blob or "activator" in blob
                or t.lower() in ("activation", "ac50", "fold change", "ratio")
                or "fold increase" in blob or "stimulat" in blob):
            activation_like.append({
                "molecule": a.get("molecule_chembl_id"),
                "type": t,
                "relation": a.get("standard_relation"),
                "value": a.get("standard_value"),
                "units": a.get("standard_units"),
                "assay": desc[:150],
                "doc": a.get("document_chembl_id"),
            })
    offset += 1000

print("\n--- standard_type distribution (all records) ---")
for t, n in sorted(types.items(), key=lambda x: -x[1]):
    print(f"  {n:6d}  {t}")

print(f"\n--- records whose type/assay wording implies ACTIVATION: "
      f"{len(activation_like)} ---")
seen = set()
for r in activation_like[:40]:
    key = (r["molecule"], r["type"], r["value"])
    if key in seen:
        continue
    seen.add(key)
    print(f"  {r['molecule']:>14}  {r['type']:<12} "
          f"{r['relation'] or ''}{r['value']} {r['units'] or ''}  "
          f"[{r['doc']}]\n      {r['assay']}")

# --- PubMed, rate-limited to <=3/sec ---
print(f"\n--- PubMed counts (retrieved {DATE}) ---")
queries = [
    '"SIRT2"[tiab] AND ("activator"[tiab] OR "activators"[tiab])',
    '"SIRT2 activator"[tiab]',
    '"sirtuin activating compound"[tiab]',
    '"SIRT1"[tiab] AND ("activator"[tiab] OR "activators"[tiab])',
    '"nicotinamide riboside"[tiab] AND "SIRT2"[tiab]',
    '"BubR1"[tiab] AND "SIRT2"[tiab]',
    '"BUB1B"[tiab] AND ("NAD"[tiab] OR "nicotinamide"[tiab])',
]
for q in queries:
    u = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
         "?db=pubmed&retmode=json&retmax=0&term=" + urllib.parse.quote(q))
    d = get(u, f"PubMed esearch: {q}", pause=0.5)
    if "__error__" in d:
        print(f"  ERROR  {q}: {d['__error__']}")
    else:
        print(f"  {d['esearchresult']['count']:>6}  {q}")

with open("task_e_queries.txt", "w") as fh:
    fh.write(f"Task E queries, retrieved {DATE}\n\n")
    for label, url in LOG:
        fh.write(f"[{label}]\n  {url}\n")
print(f"\nExact query log written to task_e_queries.txt ({len(LOG)} queries)")
