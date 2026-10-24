"""The assay title said "Destabilization". Two of its records say the opposite.

This script exists because it refutes a claim this report made twice.

10_direction_audit.py counts the direction of medicinal chemistry by reading
what assay DESCRIPTIONS say, and on that basis the report asserted:

    "Not one molecule ever measured against BubR1 raises it. The only
     cell-based compound is an explicit destabiliser."

Both sentences are wrong, and the method that produced them is the reason.

ChEMBL assay CHEMBL6125447 is a cellular thermal shift assay (CETSA) in
CAKI-1 cells, titled:

    "Destabilization of BUB1B in human CAKI-1 cells assessed as reduction in
     BUB1B protein level at 20 uM preincubated for 1 hrs followed by
     incubation under 45 degreeC for 3 mins by CETSA analysis (Rvb = 1)"

The title states one direction. The four records inside it do not agree with
it: two report fold-change BELOW 1 (destabilised) and two report fold-change
ABOVE 1, each carrying the curator annotation "Increase In Heat Resistance
Observed". The source paper (El Hafi et al., Eur J Med Chem 2025;285:117247,
PMID 39818011) says so in its own abstract: compound 8h "destabilizes the
BUB1B protein structure, whereas TG-101209 stabilizes it".

THE METHODOLOGICAL POINT, which is the part that generalises:

    Direction is a property of the RECORD, not of the assay title. An assay
    titled for one direction can contain both. Any audit that counts
    description text - including this report's own 42:1 analysis - will
    misclassify every such record, and will do so silently.

So this script reads direction at the record level: standard_type,
standard_relation, standard_value, activity_comment. That is the corrected
procedure, and it is what found the result the title had hidden.

WHAT THIS DOES NOT SHOW, stated here so the report cannot drift:
  - CETSA measures target engagement and thermal stabilisation of the FOLD.
    It is not a measurement of steady-state protein abundance. A ligand can
    stabilise a fold and inhibit the protein at the same time.
  - The measurement is on WILD-TYPE BubR1 in a renal carcinoma line. Nothing
    here was measured on p.Asn1002Lys.
  - TG-101209 is a cytotoxic JAK2 inhibitor. It is a tool compound and a
    chemical starting point, not a proposed therapy for anyone.

Read-only. No API key, no challenge data. Run: python3 -I this_script.py
"""
import json
import time
import urllib.request

UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}
BASE = "https://www.ebi.ac.uk/chembl/api/data"
BUB1B = "CHEMBL4295998"


def get(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.load(r)
        except Exception as exc:
            if i == tries - 1:
                print(f"  ! failed: {url} ({exc})")
                return None
            time.sleep(2 * (i + 1))
    return None


d = get(f"{BASE}/activity.json?target_chembl_id={BUB1B}&limit=200")
if d is None:
    raise SystemExit("could not retrieve BUB1B activities")

acts = d["activities"]
print(f"=== BUB1B ({BUB1B}) ===")
print(f"activity records: {d['page_meta']['total_count']}\n")

by_type = {}
for a in acts:
    by_type.setdefault(a.get("standard_type"), []).append(a)
print("records by standard_type:")
for k, v in sorted(by_type.items(), key=lambda kv: -len(kv[1])):
    print(f"  {str(k):<14} {len(v)}")

# ---------------------------------------------------------------- the IC50s
print("\n=== IC50 records (report stated the range as 9,670-67,300 nM) ===")
ic = sorted((float(a["standard_value"]), a) for a in by_type.get("IC50", [])
            if a.get("standard_value"))
for v, a in ic:
    print(f"  {a['molecule_chembl_id']:<16} {v:>10,.0f} nM   "
          f"doc={a.get('document_chembl_id')}")
if ic:
    print(f"  -> n={len(ic)}, true range {ic[0][0]:,.0f} - "
          f"{ic[-1][0]:,.0f} nM")
    print("  NOTE: the upper bound is 100,000 nM, not 67,300. The report's")
    print("        stated range omitted the weakest compound.")

# --------------------------------------------------- the direction reversal
print("\n=== cell-based CETSA records: direction read PER RECORD ===")
print("    (the assay TITLE says 'Destabilization' for all four)\n")
fc = by_type.get("FC", [])
for a in sorted(fc, key=lambda a: float(a["standard_value"])):
    v = float(a["standard_value"])
    verdict = "STABILISED  (heat resistance UP)" if v > 1 else "destabilised"
    print(f"  {a['molecule_chembl_id']:<16} FC={v:<6} {verdict}")
    if a.get("activity_comment"):
        print(f"      activity_comment: {a['activity_comment']!r}")
print(f"\n  assay: {fc[0]['assay_chembl_id'] if fc else 'n/a'}")
if fc:
    print(f"  title: {fc[0].get('assay_description')}")

up = [a for a in fc if float(a["standard_value"]) > 1]
print(f"\n  -> {len(up)} of {len(fc)} cell-based compounds RAISED BubR1's")
print("     thermal stability. The assay title says the opposite, and an")
print("     audit that reads titles counts all four as destabilisation.")

# ------------------------------------------------------------- the molecules
print("\n=== identity of the stabilising compounds ===")
for a in up:
    m = get(f"{BASE}/molecule/{a['molecule_chembl_id']}.json")
    time.sleep(0.2)
    if not m:
        continue
    syns = [s.get("molecule_synonym") for s in (m.get("molecule_synonyms") or [])]
    print(f"  {m['molecule_chembl_id']:<16} pref_name={m.get('pref_name')} "
          f"max_phase={m.get('max_phase')}")
    print(f"      synonyms: {syns[:8] if syns else 'none in ChEMBL'}")
    st = m.get("molecule_structures") or {}
    print(f"      smiles: {(st.get('canonical_smiles') or '')[:80]}")

# ----------------------------------------------------------------- the paper
doc = get(f"{BASE}/document/CHEMBL6121040.json")
if doc:
    print("\n=== source publication ===")
    print(f"  {doc.get('title')}")
    print(f"  {doc.get('journal')} {doc.get('year')};"
          f"{doc.get('volume')}:{doc.get('first_page')}  "
          f"PMID {doc.get('pubmed_id')}  doi:{doc.get('doi')}")
    ab = doc.get("abstract") or ""
    key = "whereas TG-101209 stabilizes it"
    if key in ab:
        i = ab.index(key)
        print(f"\n  the authors say it themselves:")
        print(f"    \"...{ab[max(0, i - 190):i + len(key)]}.\"")
