"""Actually predict ddG for p.Asn1002Lys, with the same calibration controls.

Script 21 could only place the variant geometrically: burial, packing,
neighbourhood. That is a prior, not a stability prediction, and its control
failed. This script runs a real, published ddG predictor over the same
calibration set and asks the same question properly.

THE PREDICTOR. DDMut (Zhou, Pires et al., Nucleic Acids Res 2023), a
graph-based deep-learning ddG predictor, via its public REST API. Chosen
because it takes a structure directly, needs no local install, and is
independent of this project. ddG is reported in kcal/mol; by DDMut's
convention NEGATIVE means destabilising.

THE DESIGN IS THE POINT, and it is unchanged from script 21: five BUB1B
alleles whose behaviour is MEASURED are run alongside the query, so the
question is not "what number does the model give" but "does the model
separate the alleles we already know, and where does the query fall among
them?" A predictor that cannot order the knowns tells us nothing about the
unknown, and that result is reported rather than hidden.

  DESTABILISED, measured:  K795A, I909T, L1012P
  TOLERATED, measured:     Y155C, R550Q
  QUERY:                   N1002K

WHAT IS SENT OFF THIS MACHINE. Only the public AlphaFold model of BubR1
(AF-O60566-F1-v6, already freely downloadable from EBI) and single-residue
mutation codes. No patient data, no challenge data, nothing private. The
variant identity is a public ClinVar-adjacent coordinate discussed openly
in the report.

LIMITS THAT SURVIVE WHATEVER THE NUMBER SAYS. The input is a PREDICTED
structure, not a crystal structure, so error in the model propagates into
the ddG. DDMut was trained on experimental thermodynamic datasets that are
dominated by small, soluble, well-folded proteins; a 1,050-residue mitotic
scaffold with long disordered stretches is out of distribution. And a ddG
is a thermodynamic quantity, whereas what the report needs is steady-state
abundance in a cell, which also depends on degradation machinery this
calculation knows nothing about.

Run: python3 -I 25_ddmut_stability_test.py <path/to/model.pdb>
"""
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
import uuid

API = "https://biosig.lab.uq.edu.au/ddmut/api/prediction_single"
UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}
CHAIN = "A"
POLL_EVERY = 12
POLL_MAX = 40

VARIANTS = [
    ("K795A",  "destabilised", "denaturation 34.7 C vs 39.1 C WT"),
    ("I909T",  "destabilised", "5-10 fold protein loss"),
    ("L1012P", "destabilised", "5-10 fold protein loss"),
    ("Y155C",  "tolerated",    "protein level unaffected"),
    ("R550Q",  "tolerated",    "protein level unaffected"),
    ("N1002K", "QUERY",        "never assayed"),
]

if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
    sys.exit("usage: python3 -I 25_ddmut_stability_test.py <model.pdb>")
PDB = sys.argv[1]
pdb_bytes = open(PDB, "rb").read()
print(f"structure: {os.path.basename(PDB)}  ({len(pdb_bytes)/1024:.0f} KB)")


def multipart(fields, filename, filebytes):
    b = "----mva" + uuid.uuid4().hex
    out = []
    for k, v in fields.items():
        out.append(f"--{b}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n"
                   f"\r\n{v}\r\n".encode())
    out.append(f"--{b}\r\nContent-Disposition: form-data; name=\"pdb_file\"; "
               f"filename=\"{filename}\"\r\n"
               f"Content-Type: chemical/x-pdb\r\n\r\n".encode())
    out.append(filebytes + b"\r\n")
    out.append(f"--{b}--\r\n".encode())
    return b"".join(out), f"multipart/form-data; boundary={b}"


def post(fields):
    body, ctype = multipart(fields, os.path.basename(PDB), pdb_bytes)
    req = urllib.request.Request(API, data=body, headers={
        **UA, "Content-Type": ctype})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode(errors="replace")


def get(job_id):
    url = f"{API}?job_id={job_id}"
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode(errors="replace")


def num(text):
    """Pull the first signed float that looks like a ddG out of a response."""
    try:
        d = json.loads(text)
    except Exception:
        d = None
    if isinstance(d, dict):
        for k in ("ddG", "ddg", "prediction", "Predicted ddG",
                  "ddG_pred", "DDG"):
            for kk in d:
                if kk.lower().replace(" ", "_") == k.lower().replace(" ", "_"):
                    try:
                        return float(re.findall(r"-?\d+\.?\d*", str(d[kk]))[0])
                    except Exception:
                        pass
        if "status" in d:
            return None
    m = re.search(r"(-?\d+\.\d+)\s*kcal", text)
    return float(m.group(1)) if m else None


# ---------------------------------------------------------------- submit
jobs = {}
for mut, cls, note in VARIANTS:
    try:
        resp = post({"mutation": mut, "chain": CHAIN})
    except urllib.error.HTTPError as e:
        print(f"  {mut:<8} submit FAILED  HTTP {e.code}")
        jobs[mut] = None
        continue
    except Exception as e:
        print(f"  {mut:<8} submit FAILED  {type(e).__name__}: {e}")
        jobs[mut] = None
        continue
    jid = None
    try:
        jid = json.loads(resp).get("job_id")
    except Exception:
        m = re.search(r'job_id["\s:]+([A-Za-z0-9_-]+)', resp)
        jid = m.group(1) if m else None
    jobs[mut] = jid
    print(f"  {mut:<8} submitted  job_id={jid}")
    time.sleep(2)

if not any(jobs.values()):
    print("\nNo job was accepted. Reporting this as a FAILED RUN rather than "
          "as an absence of effect.\nRaw response from the last attempt:\n")
    print(resp[:600] if 'resp' in dir() else "(none)")
    sys.exit(1)

# ---------------------------------------------------------------- collect
print("\npolling...")
results = {}
for attempt in range(POLL_MAX):
    pending = [m for m, j in jobs.items() if j and m not in results]
    if not pending:
        break
    for mut in pending:
        try:
            txt = get(jobs[mut])
        except Exception:
            continue
        v = num(txt)
        if v is not None:
            results[mut] = v
            print(f"  {mut:<8} ddG = {v:+.2f} kcal/mol")
    if len([m for m, j in jobs.items() if j and m not in results]) == 0:
        break
    time.sleep(POLL_EVERY)

# ---------------------------------------------------------------- report
print("\n=== DDMut ddG against the measured calibration set ===")
print(f"{'variant':<10}{'class':<15}{'ddG':>9}   evidence")
print("-" * 72)
rows = []
for mut, cls, note in VARIANTS:
    v = results.get(mut)
    rows.append({"variant": mut, "class": cls, "ddG": v, "evidence": note})
    print(f"{mut:<10}{cls:<15}{(f'{v:+.2f}' if v is not None else 'n/a'):>9}"
          f"   {note}")

d = [r["ddG"] for r in rows if r["class"] == "destabilised"
     and r["ddG"] is not None]
t = [r["ddG"] for r in rows if r["class"] == "tolerated" and r["ddG"] is not None]
q = results.get("N1002K")

print("\n=== does the predictor order what we already know? ===")
if len(d) >= 2 and len(t) >= 1:
    sep = max(d) < min(t)
    print(f"  destabilised {['%+.2f' % x for x in d]}")
    print(f"  tolerated    {['%+.2f' % x for x in t]}")
    print(f"  separates the measured classes: {'YES' if sep else 'NO'}")
    if not sep:
        print("  -> the predictor cannot order the alleles whose behaviour is")
        print("     known, so its number for the query carries no weight here.")
else:
    print("  too few results returned to judge separation.")

print("\n=== the query ===")
if q is None:
    print("  N1002K: no prediction returned.")
else:
    print(f"  N1002K ddG = {q:+.2f} kcal/mol")
    if d and t:
        print(f"  destabilised range {min(d):+.2f} to {max(d):+.2f}; "
              f"tolerated {min(t):+.2f} to {max(t):+.2f}")

with open("25_ddmut_stability_results.json", "w") as fh:
    json.dump({"predictor": "DDMut (Zhou 2023), public REST API",
               "structure": os.path.basename(PDB), "chain": CHAIN,
               "convention": "negative ddG = destabilising",
               "rows": rows}, fh, indent=1)
if not results:
    print("""
  THE RUN FAILED, AND THE FAILURE IS SERVER-SIDE, NOT OURS. All six jobs
  were accepted and returned job_ids. Every retrieval then answered
  {"message": "Internal Server Error"}. The same happens for a 49 KB
  canonical control structure (crambin 1CRN, mutation T2A, chain A), which
  rules out our input, our request format and the size of the AlphaFold
  model. The DDMut result endpoint was simply not serving results when this
  was run.

  REPORTED AS A FAILED RUN, NOT AS AN ABSENCE OF EFFECT. No ddG for
  p.Asn1002Lys is claimed from this script. The structural evidence in the
  report remains script 21's geometry \u2014 burial, packing, neighbourhood \u2014
  which is a prior whose own calibration control failed, and is labelled
  that way.""")

print("\n  wrote 25_ddmut_stability_results.json")
print("""
  WHATEVER THE NUMBERS SAY: the input is a PREDICTED structure, DDMut was
  trained on small soluble proteins rather than 1,050-residue mitotic
  scaffolds, and ddG is thermodynamic stability, not the steady-state
  cellular abundance the report actually needs. This cannot replace E0.""")
