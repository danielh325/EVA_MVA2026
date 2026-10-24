"""Is the directional asymmetry real, and does it track what the protein IS?

11_direction_breadth.py measured the inhibition:activation asymmetry across
200 human targets, and reported an honest limitation: those targets were the
first page of an unsorted API listing. A convenience sample cannot support a
claim about "the median human protein target", so this script replaces it.

Two changes, and the second is the interesting one.

1. A SEEDED SIMPLE RANDOM SAMPLE over the full sampling frame. ChEMBL holds
   5,869 human SINGLE PROTEIN targets; this enumerates all of them by paging
   the API, then draws from the complete list with a fixed seed. Sampling
   from an enumerated frame is both reproducible and unbiased by listing
   order, and it costs ~6 requests rather than one per draw.

2. STRATIFICATION BY ChEMBL PROTEIN CLASS. The report's explanation for the
   asymmetry is mechanistic: you can mimic a ligand at a receptor and you
   cannot easily make an enzyme go faster, so proteins with a native
   on-switch should carry activation chemistry and proteins without one
   should not. That is a testable prediction about target CLASS, not a
   property of any one protein. If the explanation is right, activation
   share should rank: ion channels and nuclear receptors high, enzymes and
   structural/scaffold proteins at zero.

   This converts the 42:1 number from an observation into a hypothesis with
   a held-out prediction attached.

Method notes, because they decide whether the number means anything:
  - Counting is by STATED direction in curator-written assay descriptions,
    using prefix matching (description__istartswith). An earlier version used
    icontains, which scored "Inhibition of Aurora B activation" as an
    activation assay. Prefix matching makes every activation count here
    CONSERVATIVE: the true asymmetry is if anything larger.
  - API failures are recorded as failures, NOT as zeros. The previous script
    tested `if not tot`, which silently merged a failed request with a target
    that genuinely has no assays. That is the difference between a missing
    measurement and a measurement of zero, and a script that conflates them
    will report a cleaner result than it earned.
  - Targets are kept at >=50 assays so that a share is meaningful; the
    unfiltered distribution is reported alongside it rather than discarded.

Read-only. No API key, no challenge data. Run: python3 -I this_script.py
"""
import json
import random
import time
import urllib.request

UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}
BASE = "https://www.ebi.ac.uk/chembl/api/data"
A = f"{BASE}/assay.json"
T = f"{BASE}/target.json"

SEED = 20261022          # the retrieval date, so the draw is reproducible
DRAW = 600               # random targets to draw from the frame
MIN_ASSAYS = 50          # threshold for "a share is meaningful"


def get(url, tries=4):
    """Fetch JSON. Returns None only after genuine repeated failure."""
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.load(r)
        except Exception:
            if i == tries - 1:
                return None
            time.sleep(2 * (i + 1))
    return None


def count(url):
    """Total_count for a query, or None if the request failed.

    None and 0 are different answers and are kept different.
    """
    d = get(url)
    if d is None:
        return None
    return d.get("page_meta", {}).get("total_count")


# ---------------------------------------------------------------- class tree
print("building ChEMBL protein class hierarchy...", flush=True)
classes = {}
off = 0
while True:
    d = get(f"{BASE}/protein_classification.json?limit=1000&offset={off}")
    if not d:
        break
    batch = d.get("protein_classifications", [])
    if not batch:
        break
    for c in batch:
        classes[c["protein_class_id"]] = c
    off += 1000
    if off >= d.get("page_meta", {}).get("total_count", 0):
        break
print(f"  {len(classes)} protein classes", flush=True)


def top_level(class_id):
    """Walk a class id up to its level-1 ancestor (Enzyme, Receptor, ...)."""
    seen = set()
    cur = class_id
    while cur is not None and cur not in seen:
        seen.add(cur)
        c = classes.get(cur)
        if c is None:
            return None
        if c.get("class_level") == 1:
            return c.get("pref_name")
        cur = c.get("parent_id")
    return None


# ------------------------------------------------------------ sampling frame
# Enumerate the WHOLE frame first, then sample from it. Drawing by random
# offset means one API round trip per draw; paging the full list costs ~6
# requests and gives a true simple random sample from an enumerated frame,
# which is also the stronger claim.
frame = []
off = 0
while True:
    d = get(f"{T}?target_type=SINGLE%20PROTEIN&organism=Homo%20sapiens"
            f"&limit=1000&offset={off}")
    if not d:
        raise SystemExit("could not enumerate the sampling frame")
    batch = d.get("targets", [])
    if not batch:
        break
    frame.extend(batch)
    total = d.get("page_meta", {}).get("total_count", 0)
    print(f"  enumerated {len(frame)}/{total}", flush=True)
    off += 1000
    if off >= total:
        break

FRAME = len(frame)
print(f"sampling frame: {FRAME} human SINGLE PROTEIN targets", flush=True)

rng = random.Random(SEED)
sample = rng.sample(frame, min(DRAW, FRAME))
print(f"drew {len(sample)} at random (seed {SEED})", flush=True)

rows = []
failures = []
for n, t in enumerate(sample, 1):
    tid = t["target_chembl_id"]

    gene, comp_id = None, None
    for comp in t.get("target_components", []):
        comp_id = comp.get("component_id")
        for syn in comp.get("target_component_synonyms", []):
            if syn.get("syn_type") == "GENE_SYMBOL":
                gene = syn.get("component_synonym")
                break
        break

    tot = count(f"{A}?target_chembl_id={tid}&limit=1")
    time.sleep(0.1)
    if tot is None:
        failures.append({"target": tid, "stage": "assay_total"})
        continue

    row = {"target": tid, "name": t.get("pref_name"), "gene": gene,
           "assays": tot, "activation": None, "inhibition": None,
           "klass": None}

    if tot >= MIN_ASSAYS:
        row["activation"] = count(
            f"{A}?target_chembl_id={tid}&description__istartswith=Activation&limit=1")
        time.sleep(0.1)
        row["inhibition"] = count(
            f"{A}?target_chembl_id={tid}&description__istartswith=Inhibition&limit=1")
        time.sleep(0.1)
        if comp_id is not None:
            cd = get(f"{BASE}/target_component/{comp_id}.json")
            if cd:
                pcs = cd.get("protein_classifications") or []
                if pcs:
                    row["klass"] = top_level(pcs[0].get("protein_classification_id"))
            time.sleep(0.1)

    rows.append(row)
    if n % 50 == 0:
        q = sum(1 for r in rows if r["assays"] >= MIN_ASSAYS)
        print(f"  {n}/{len(sample)} profiled, {q} qualifying, "
              f"{len(failures)} failed", flush=True)

with open("12_direction_random_sample_results.json", "w") as fh:
    json.dump({"seed": SEED, "frame": FRAME, "drawn": len(sample),
               "min_assays": MIN_ASSAYS, "failures": failures,
               "rows": rows}, fh, indent=1)

# ----------------------------------------------------------------- reporting
print(f"\n=== random draw: {len(rows)} targets retrieved, "
      f"{len(failures)} failed ===")

allc = sorted(r["assays"] for r in rows)
med_all = allc[len(allc) // 2] if allc else 0
print(f"  assays per target, all drawn : median {med_all}, "
      f"max {allc[-1] if allc else 0}")
print(f"  targets with 0 assays        : {sum(1 for a in allc if a == 0)}")

q = [r for r in rows if r["assays"] >= MIN_ASSAYS
     and r["activation"] is not None and r["inhibition"] is not None]
n = len(q)
print(f"\n=== {n} targets with >={MIN_ASSAYS} assays ===")
if n:
    shares = sorted(100 * r["activation"] / r["assays"] for r in q)
    zero = sum(1 for r in q if r["activation"] == 0)
    under1 = sum(1 for r in q if r["activation"] / r["assays"] < 0.01)
    ta = sum(r["activation"] for r in q)
    ti = sum(r["inhibition"] for r in q)
    print(f"  median activation share   : {shares[n // 2]:.2f}%")
    print(f"  25th / 75th percentile    : {shares[n // 4]:.2f}% / "
          f"{shares[3 * n // 4]:.2f}%")
    print(f"  ZERO activation assays    : {zero} ({100 * zero / n:.1f}%)")
    print(f"  <1% activation assays     : {under1} ({100 * under1 / n:.1f}%)")
    print(f"  pooled activation         : {ta}")
    print(f"  pooled inhibition         : {ti}")
    if ta:
        print(f"  pooled inhibition:activation : {ti / ta:.1f} : 1")

    # ---------------------------------------------- the held-out prediction
    print(f"\n=== by ChEMBL top-level protein class ===")
    print("  (prediction: native on-switch high, scaffolds/enzymes zero)")
    by = {}
    for r in q:
        by.setdefault(r["klass"] or "unclassified", []).append(r)
    stats = []
    for k, rs in by.items():
        ta = sum(r["activation"] for r in rs)
        tt = sum(r["assays"] for r in rs)
        z = sum(1 for r in rs if r["activation"] == 0)
        stats.append((100 * ta / tt if tt else 0, k, len(rs), ta, tt, z))
    stats.sort(reverse=True)
    print(f"  {'class':<26} {'n':>4} {'act%':>7} {'zero':>6}  act/assays")
    for pct, k, cnt_, ta, tt, z in stats:
        print(f"  {k:<26} {cnt_:>4} {pct:>6.2f}% {z:>5}/{cnt_}  {ta}/{tt}")
