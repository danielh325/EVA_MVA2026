"""The 42:1 problem is not about BubR1. This measures it where it hurts.

10_direction_audit.py found that ChEMBL assay descriptions beginning
"Inhibition of" outnumber "Activation of" by 42:1, and 12_direction_random
_sample.py showed that holds across a random draw of human protein targets
and tracks target class. Both are statements about medicinal chemistry in
general.

This script asks the question that actually matters for a rare-disease
patient: FOR THE GENES WHERE HUMAN DISEASE IS CAUSED BY HAVING TOO LITTLE
PROTEIN, HOW MUCH OF THE AVAILABLE CHEMISTRY POINTS THE RIGHT WAY?

Sampling frame is not ours to choose, which is the point. ClinGen's Dosage
Sensitivity curation is an independent, expert-reviewed, public assignment
of whether losing a gene causes disease:

  Haploinsufficiency score 3  - sufficient evidence for dosage pathogenicity
                                (one working copy is not enough: dominant LoF)
  Haploinsufficiency score 30 - gene associated with an autosomal recessive
                                phenotype (biallelic LoF)

Both classes need MORE of the protein. BUB1B is in the second one, scored 30
by ClinGen, which is the useful part: this child's gene is not a special case
to be pitied, it is a typical member of a curated population of ~1,000 genes,
and whatever is true of that population is true of hundreds of other families.

POSITIVE CONTROL, because a method that can only find zeros is not measuring
anything. CFTR is in this frame, and CFTR potentiators (ivacaftor) are the
best-known drugs in existence that make a loss-of-function protein work
better. If this pipeline does not detect activation chemistry at CFTR, it is
broken and its zeros mean nothing. The control is run and reported.

Method notes carried over from the earlier scripts:
  - Prefix matching on curator-written descriptions, so activation counts are
    CONSERVATIVE.
  - This counts STATED direction, not measured direction. A compound can
    raise a protein without any assay description saying "activation" - that
    is the schema gap the report describes, and it means the true picture is
    somewhat better than these numbers. Reported as a limitation, not buried.
  - API failures are recorded as failures, never as zeros.

Read-only. No API key, no challenge data.
Run: python3 -I 13_lof_disease_gene_direction.py <ClinGen_gene_curation_list_GRCh38.tsv>
"""
import json
import random
import sys
import time
import urllib.parse
import urllib.request

UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}
BASE = "https://www.ebi.ac.uk/chembl/api/data"
A = f"{BASE}/assay.json"
T = f"{BASE}/target.json"

# ClinGen haploinsufficiency scores that mean "too little protein causes disease"
LOF_SCORES = {"3": "dominant LoF (haploinsufficient)",
              "30": "recessive LoF (biallelic)"}

SEED = 20261022          # same seed as 12_direction_random_sample.py
CONTROLS = ["CFTR", "BUB1B"]


def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=20) as r:
                return json.load(r)
        except Exception:
            if i == tries - 1:
                return None
            time.sleep(2 * (i + 1))
    return None


def count(url):
    """total_count, or None on genuine failure. None != 0."""
    d = get(url)
    if d is None:
        return None
    return d.get("page_meta", {}).get("total_count")


def parse_clingen(path):
    """Read the ClinGen curation TSV. Untrusted input: parsed as data only."""
    out, hdr = [], None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#Gene Symbol"):
                hdr = line[1:].rstrip("\n").split("\t")
                continue
            if line.startswith("#") or not line.strip() or hdr is None:
                continue
            out.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    return out


def resolve(gene):
    """Gene symbol -> human SINGLE PROTEIN ChEMBL target, or None."""
    g = urllib.parse.quote(gene)
    d = get(f"{T}?target_type=SINGLE%20PROTEIN&organism=Homo%20sapiens"
            f"&target_components__target_component_synonyms__component_synonym__iexact={g}")
    if not d or not d.get("targets"):
        return None
    return d["targets"][0]["target_chembl_id"], d["targets"][0].get("pref_name")


def profile(tid):
    """Assay counts for a target: (total, activation, inhibition)."""
    tot = count(f"{A}?target_chembl_id={tid}&limit=1")
    time.sleep(0.08)
    if tot is None:
        return None, None, None
    if tot == 0:
        return 0, 0, 0
    act = count(f"{A}?target_chembl_id={tid}"
                f"&description__istartswith=Activation&limit=1")
    time.sleep(0.08)
    inh = count(f"{A}?target_chembl_id={tid}"
                f"&description__istartswith=Inhibition&limit=1")
    time.sleep(0.08)
    return tot, act, inh


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: 13_lof_disease_gene_direction.py <clingen.tsv>")
    rows = parse_clingen(sys.argv[1])
    print(f"ClinGen curation list: {len(rows)} genes", flush=True)

    lof = [(r["Gene Symbol"], r["Haploinsufficiency Score"])
           for r in rows if r.get("Haploinsufficiency Score") in LOF_SCORES]
    # Shuffle with a fixed seed so that ANY PREFIX of this run is itself a
    # random sample of the frame. The ChEMBL API is slow enough that a run
    # may be cut short, and a truncated alphabetical pass would be exactly
    # the convenience sample that 12_direction_random_sample.py exists to
    # replace. Seeded, so the order is reproducible.
    random.Random(SEED).shuffle(lof)
    print(f"loss-of-function disease genes: {len(lof)}")
    for s, label in LOF_SCORES.items():
        print(f"  score {s:<3} {label:<36} "
              f"{sum(1 for _, x in lof if x == s)}")

    # ------------------------------------------------- positive control first
    print("\n=== positive control: can this pipeline see activation at all? ===")
    for g in CONTROLS:
        r = resolve(g)
        if not r:
            print(f"  {g:<8} no ChEMBL human SINGLE PROTEIN target")
            continue
        tid, name = r
        tot, act, inh = profile(tid)
        print(f"  {g:<8} {tid:<16} assays={tot} activation={act} "
              f"inhibition={inh}   {name}")
    print("  (CFTR must show activation assays or the zeros below are "
          "meaningless)\n", flush=True)

    # Resumable by design. A single hung socket used to cost the whole run,
    # because results lived in memory until the end. Each gene is now appended
    # to a JSONL as soon as it resolves, and a restart skips what is already
    # there. The ChEMBL API is the slow, flaky part of this script; the script
    # should not also be the fragile part.
    CACHE = "13_lof_disease_gene_direction_partial.jsonl"
    done = {}
    try:
        with open(CACHE) as fh:
            for line in fh:
                if line.strip():
                    rec = json.loads(line)
                    done[rec["gene"]] = rec
    except FileNotFoundError:
        pass
    if done:
        print(f"resuming: {len(done)} genes already cached", flush=True)

    results, failures = [], []
    cache = open(CACHE, "a")
    for n, (gene, score) in enumerate(lof, 1):
        if gene in done:
            results.append(done[gene])
            continue
        r = resolve(gene)
        time.sleep(0.08)
        rec = None
        if r is None:
            # no ChEMBL target: no chemical matter at all, which is a result
            rec = {"gene": gene, "score": score, "target": None,
                   "assays": None, "activation": None, "inhibition": None}
        else:
            tid, name = r
            tot, act, inh = profile(tid)
            if tot is None:
                failures.append(gene)
            else:
                rec = {"gene": gene, "score": score, "target": tid,
                       "name": name, "assays": tot,
                       "activation": act, "inhibition": inh}
        if rec is not None:
            results.append(rec)
            cache.write(json.dumps(rec) + "\n")
            cache.flush()
        # progress must print regardless of which branch ran: an earlier
        # version put this after `continue` statements, so it only fired
        # when a gene happened to both resolve and profile, and a long run
        # could look wedged while working perfectly.
        if n % 25 == 0:
            got = sum(1 for x in results if x["target"])
            print(f"  {n}/{len(lof)} genes · {got} with a ChEMBL target · "
                  f"{len(failures)} failed", flush=True)

    cache.close()
    with open("13_lof_disease_gene_direction_results.json", "w") as fh:
        json.dump({"n_curated": len(rows), "n_lof": len(lof),
                   "failures": failures, "rows": results}, fh, indent=1)

    report(results, failures)


def report(results, failures):
    n = len(results)
    no_target = [r for r in results if r["target"] is None]
    has_t = [r for r in results if r["target"] is not None]
    no_assay = [r for r in has_t if r["assays"] == 0]
    studied = [r for r in has_t if (r["assays"] or 0) > 0]
    no_act = [r for r in studied if r["activation"] == 0]
    some_act = [r for r in studied if (r["activation"] or 0) > 0]

    print(f"\n=== {n} ClinGen loss-of-function disease genes "
          f"({len(failures)} failed) ===")
    print(f"  no ChEMBL human target at all      : {len(no_target):>5} "
          f"({100*len(no_target)/n:.1f}%)")
    print(f"  target exists, ZERO assays         : {len(no_assay):>5} "
          f"({100*len(no_assay)/n:.1f}%)")
    print(f"  has been assayed at least once     : {len(studied):>5} "
          f"({100*len(studied)/n:.1f}%)")
    if studied:
        print(f"\n  of the {len(studied)} genes anyone has ever assayed:")
        print(f"    ZERO activation assays           : {len(no_act):>5} "
              f"({100*len(no_act)/len(studied):.1f}%)")
        print(f"    some activation chemistry        : {len(some_act):>5} "
              f"({100*len(some_act)/len(studied):.1f}%)")
        ta = sum(r["activation"] or 0 for r in studied)
        ti = sum(r["inhibition"] or 0 for r in studied)
        print(f"    pooled activation assays         : {ta:>5}")
        print(f"    pooled inhibition assays         : {ti:>5}")
        if ta:
            print(f"    pooled inhibition : activation   : {ti/ta:.1f} : 1")

    # the number that matters: of ALL curated LoF genes, how many have any
    # chemistry pointing the direction the patient needs
    print(f"\n  ACROSS ALL {n} LoF DISEASE GENES, the fraction with any")
    print(f"  assay described as activation: {len(some_act)}/{n} = "
          f"{100*len(some_act)/n:.1f}%")

    print("\n=== split by inheritance ===")
    for s, label in LOF_SCORES.items():
        sub = [r for r in results if r["score"] == s]
        if not sub:
            continue
        sub_st = [r for r in sub if (r["assays"] or 0) > 0]
        sub_act = [r for r in sub_st if (r["activation"] or 0) > 0]
        print(f"  {label:<36} n={len(sub):<4} assayed={len(sub_st):<4} "
              f"with activation={len(sub_act)} "
              f"({100*len(sub_act)/len(sub):.1f}%)")

    print("\n=== LoF genes WITH activation chemistry (the exceptions) ===")
    some_act.sort(key=lambda r: -(r["activation"] / r["assays"]))
    print(f"  {'gene':<10} {'act%':>7} {'act':>5} {'assays':>7}  name")
    for r in some_act[:30]:
        pct = 100 * r["activation"] / r["assays"]
        print(f"  {r['gene']:<10} {pct:>6.1f}% {r['activation']:>5} "
              f"{r['assays']:>7}  {(r.get('name') or '')[:46]}")


if __name__ == "__main__":
    main()
