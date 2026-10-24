"""The 42-to-1 problem: measure the directional asymmetry of medicinal chemistry.

Section 4.5 of the report establishes that BubR1 has no usable chemical starting
point in the direction this genotype needs. Section 4.6 asks whether that is
bad luck or structural, and this script answers it against ChEMBL.

Three measurements, each independently checkable:

  1. GLOBAL. Count assays by the direction their description states. Across
     ~2M ChEMBL assays, "Inhibition of..." outnumbers "Activation of..." by
     roughly 42:1. The receptor vocabulary is the control: "Agonist..." vs
     "Antagonist..." is near parity and actually favours the agonist, because
     mimicking a ligand at a receptor is an ordinary thing to do while making
     an enzyme go faster is not.

  2. PANEL. The same split for the lesion, its checkpoint neighbours, and the
     mechanism-adjacent proteins where this report needs an INCREASE.

  3. BUB1B SPECIFICALLY. Both knowledge graphs report zero drug edges for
     BUB1B. ChEMBL has a BUB1B target with 24 activity records. Every one of
     them inhibits or destabilises; none raises it. That discrepancy is the
     coverage gap at the most important node in the report.

A note on method. An earlier version matched `description__icontains=activation`
and miscounted phrases like "Inhibition of Aurora B activation" as activation,
inflating AURKB from 2 to 3 and HDAC1 from 3 to 6. This uses
`description__istartswith`, so every activation count here is conservative and
the true asymmetry is if anything larger. ChEMBL descriptions are curator-written
free text, so this counts STATED direction, not measured direction - stated as a
limitation in 4.6 rather than worked around.

Read-only. No API key, no challenge data. Run: python3 -I this_script.py
"""
import json
import time
import urllib.parse
import urllib.request
import datetime

UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}
DATE = datetime.date.today().isoformat()
ASSAY = "https://www.ebi.ac.uk/chembl/api/data/assay.json"
ACT = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
TSEARCH = "https://www.ebi.ac.uk/chembl/api/data/target/search.json"


def _get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            if i == tries - 1:
                print(f"    ! {type(e).__name__} on {url[:90]}")
                return None
            time.sleep(2 * (i + 1))
    return None


def count(url):
    d = _get(url)
    time.sleep(0.15)
    return None if d is None else d.get("page_meta", {}).get("total_count")


def resolve(gene):
    """Gene symbol -> human SINGLE PROTEIN ChEMBL target, matched on synonym."""
    d = _get(f"{TSEARCH}?q={urllib.parse.quote(gene)}&limit=30")
    time.sleep(0.2)
    if not d:
        return None, "lookup failed"
    fallback = None
    for t in d.get("targets", []):
        if (t.get("organism") != "Homo sapiens"
                or t.get("target_type") != "SINGLE PROTEIN"):
            continue
        syns = {c.get("component_synonym", "").upper()
                for tc in t.get("target_components", [])
                for c in tc.get("target_component_synonyms", [])}
        if gene.upper() in syns:
            return t["target_chembl_id"], t.get("pref_name")
        fallback = fallback or (t["target_chembl_id"], t.get("pref_name"))
    return fallback or (None, "no human SINGLE PROTEIN match")


PANEL = [
    ("BUB1B", "the lesion"), ("BUB1", "checkpoint"), ("BUB3", "checkpoint"),
    ("CEP57", "MVA gene"), ("TRIP13", "MVA gene"),
    ("AURKB", "druggable / contraindicated"),
    ("PLK1", "druggable / contraindicated"),
    ("TTK", "druggable / contraindicated"),
    ("CENPE", "druggable / contraindicated"),
    ("SIRT2", "we need an INCREASE"), ("SIRT1", "comparator"),
    ("NAMPT", "we need an INCREASE"),
    ("CREBBP", "inhibit the writer"), ("EP300", "inhibit the writer"),
    ("HDAC1", "comparator"), ("HSP90AA1", "contraindicated"),
]

DIRECTIONS = ["Inhibition", "Activation", "Agonist", "Antagonist",
              "Stimulation", "Degradation"]


def main():
    print(f"=== ChEMBL direction-of-effect audit | retrieved {DATE} ===\n")

    # ---- 1. global ----
    print("--- 1. GLOBAL: all ChEMBL assays by stated direction ---")
    total = count(f"{ASSAY}?limit=1")
    print(f"  {str(total):>10}  all assays")
    glob = {}
    for word in DIRECTIONS:
        n = count(f"{ASSAY}?limit=1&description__istartswith={word}")
        glob[word] = n
        print(f"  {str(n):>10}  description starts \"{word}\"")
    if glob.get("Inhibition") and glob.get("Activation"):
        print(f"\n  enzyme language   inhibition : activation = "
              f"{glob['Inhibition'] / glob['Activation']:.1f} : 1")
    if glob.get("Agonist") and glob.get("Antagonist"):
        print(f"  receptor language antagonist : agonist   = "
              f"{glob['Antagonist'] / glob['Agonist']:.2f} : 1  "
              f"(the control - near parity)")

    # ---- 2. panel ----
    print(f"\n--- 2. PANEL: per-target split ---")
    print(f"{'gene':<10}{'target':<16}{'role':<30}"
          f"{'assays':>8}{'activ':>7}{'inhib':>7}{'act%':>7}")
    print("-" * 85)
    rows = []
    for gene, role in PANEL:
        tid, name = resolve(gene)
        if not tid:
            print(f"{gene:<10}{'-':<16}{role:<30}{'NO TARGET':>29}")
            continue
        tot = count(f"{ASSAY}?target_chembl_id={tid}&limit=1")
        act = count(f"{ASSAY}?target_chembl_id={tid}"
                    f"&description__istartswith=Activation&limit=1")
        inh = count(f"{ASSAY}?target_chembl_id={tid}"
                    f"&description__istartswith=Inhibition&limit=1")
        pct = (100.0 * act / tot) if (tot and act is not None) else 0.0
        print(f"{gene:<10}{tid:<16}{role:<30}"
              f"{str(tot):>8}{str(act):>7}{str(inh):>7}{pct:>6.1f}%")
        rows.append({"gene": gene, "target": tid, "role": role,
                     "assays": tot, "activation": act, "inhibition": inh})

    # ---- 3. BUB1B in detail ----
    print(f"\n--- 3. BUB1B: what the graphs say is zero ---")
    tid, _ = resolve("BUB1B")
    if tid:
        n = count(f"{ACT}?target_chembl_id={tid}&limit=1")
        print(f"  ChEMBL target {tid}: {n} activity records "
              f"(both knowledge graphs report 0 drug edges)")
        d = _get(f"{ASSAY}?target_chembl_id={tid}&limit=25")
        time.sleep(0.2)
        for a in (d or {}).get("assays", []):
            print(f"    - {(a.get('description') or '')[:118]}")

    with open("10_direction_audit_results.json", "w") as fh:
        json.dump({"date": DATE, "global": glob, "all_assays": total,
                   "panel": rows}, fh, indent=1)
    print("\nwrote 10_direction_audit_results.json")


if __name__ == "__main__":
    main()
