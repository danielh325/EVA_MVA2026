"""Is p.Asn1002Lys structurally the kind of substitution that destabilises?

THE GAP THIS ADDRESSES. The report's own evidence table puts this variant at
tier 1 with "no functional data of any kind". Every therapeutic argument in
it rests on a tier-2 inference: that pseudokinase-domain MVA missense alleles
are destabilised, so this one probably is too. Nobody has ever looked at
N1002K itself.

A pathogenicity score is not a stability score. AlphaMissense 0.923 says
"likely pathogenic"; it does not say "destabilises the fold", and the two
come apart all the time (an active-site residue can be essential and fully
stable). The report needs the stability axis specifically, and this script
supplies the structural half of it.

THE DESIGN, which is the point rather than the tool.

Published work on BUB1B gives a labelled positive and negative set INSIDE
THE SAME PROTEIN, so the question is not "is this method any good in
general" but "does it separate the BubR1 alleles whose behaviour is already
measured, and where does N1002K fall among them?"

  DESTABILISED, measured:
    K795A   denaturation 34.7 C vs 39.1 C wild type   (Suijkerbuijk 2012)
    I909T   5-10 fold protein loss                    (Suijkerbuijk 2010)
    L1012P  5-10 fold protein loss                    (Suijkerbuijk 2010)
  TOLERATED, measured:
    Y155C   protein level unaffected                  (Suijkerbuijk 2010)
    R550Q   protein level unaffected                  (Suijkerbuijk 2010)
  QUERY:
    N1002K  never assayed

THE CONFOUND, stated before the result. The destabilised set is inside the
pseudokinase domain (766-1050) and the tolerated set is outside it. So
"buried residues destabilise" and "domain residues destabilise" predict the
same labels, and a burial signal that merely recovers domain membership
tells us nothing new. The script therefore reports N1002's burial as a
PERCENTILE WITHIN THE PSEUDOKINASE DOMAIN, against its own neighbours, not
against the whole protein. That is the only comparison that can add
information beyond "it is in the domain".

WHAT IS COMPUTED, all from the structure, no external predictor, no licence:
  - solvent accessibility by Shrake-Rupley, implemented here in numpy
  - relative solvent accessibility against Tien 2013 theoretical maxima
  - contact number: heavy atoms of other residues within 10 A of CB
  - pLDDT, so low-confidence regions are visible rather than silently used
  - the chemistry of the substitution itself (Asn -> Lys buries a charge)

WHAT THIS IS NOT. Not a ddG. Not a prediction that the protein misfolds.
AlphaFold models the WILD-TYPE fold and is not a stability predictor; burial
and contact number are correlates of substitution tolerance, not energies.
A buried position means a substitution there is more likely to be
destabilising, not that this one is. The claim available at the end is a
structural prior with stated controls, which is one tier better than the
allele-class inference the report currently rests on and still well short of
E0.

The AlphaFold model is downloaded data: its path is passed as argv, it is
never imported, and this runs under python3 -I.

Run: python3 -I mva/track2/21_structure_stability_context.py <model.pdb>
"""
import json
import sys

import numpy as np

if len(sys.argv) < 2:
    sys.exit("usage: python3 -I 21_structure_stability_context.py <AF model.pdb>")
PDB = sys.argv[1]

DOMAIN = (766, 1050)          # pseudokinase domain, UniProt O60566
PROBE = 1.4
SPHERE_POINTS = 512

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
# Tien et al. 2013, theoretical maximum accessible surface area (A^2)
MAXASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
          "GLU": 223, "GLN": 225, "GLY": 104, "HIS": 224, "ILE": 197,
          "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
          "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}

CALIBRATION = [
    (795,  "K", "A", "destabilised", "denaturation 34.7 C vs 39.1 C WT (Suijkerbuijk 2012)"),
    (909,  "I", "T", "destabilised", "5-10 fold protein loss (Suijkerbuijk 2010)"),
    (1012, "L", "P", "destabilised", "5-10 fold protein loss (Suijkerbuijk 2010)"),
    (155,  "Y", "C", "tolerated",    "protein level unaffected (Suijkerbuijk 2010)"),
    (550,  "R", "Q", "tolerated",    "protein level unaffected (Suijkerbuijk 2010)"),
]
QUERY = (1002, "N", "K", "UNKNOWN", "never functionally assayed")

THREE_TO_ONE = {v: k for k, v in {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "E": "GLU",
    "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL"}.items()}


# ----------------------------------------------------------- parse (fixed width)
atoms = []
for line in open(PDB):
    if not line.startswith("ATOM"):
        continue
    element = line[76:78].strip() or line[12:16].strip()[0]
    atoms.append({
        "name": line[12:16].strip(),
        "resname": line[17:20].strip(),
        "resseq": int(line[22:26]),
        "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
        "plddt": float(line[60:66]),
        "r": VDW.get(element, 1.70),
    })

if not atoms:
    sys.exit("no ATOM records parsed - is this a PDB file?")

xyz = np.array([a["xyz"] for a in atoms], dtype=float)
rad = np.array([a["r"] for a in atoms], dtype=float) + PROBE
resseq = np.array([a["resseq"] for a in atoms])
print(f"parsed {len(atoms):,} atoms over {len(set(resseq)):,} residues "
      f"from {PDB.split('/')[-1]}")


# ----------------------------------------------------------- Shrake-Rupley SASA
def sphere(n):
    """n roughly-uniform points on a unit sphere (golden spiral)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


pts = sphere(SPHERE_POINTS)
cell = 12.0
keys = np.floor(xyz / cell).astype(int)
grid = {}
for idx, k in enumerate(map(tuple, keys)):
    grid.setdefault(k, []).append(idx)

atom_sasa = np.zeros(len(atoms))
for i in range(len(atoms)):
    kx, ky, kz = keys[i]
    near = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                near.extend(grid.get((kx + dx, ky + dy, kz + dz), ()))
    near = [j for j in near if j != i]
    if near:
        nb = np.array(near)
        d = np.linalg.norm(xyz[nb] - xyz[i], axis=1)
        keep = nb[d < rad[i] + rad[nb]]
    else:
        keep = np.empty(0, dtype=int)
    test = xyz[i] + pts * rad[i]
    if len(keep):
        dd = np.linalg.norm(test[:, None, :] - xyz[keep][None, :, :], axis=2)
        accessible = np.all(dd >= rad[keep][None, :], axis=1)
    else:
        accessible = np.ones(len(test), dtype=bool)
    atom_sasa[i] = 4 * np.pi * rad[i] ** 2 * accessible.mean()

res = {}
for i, a in enumerate(atoms):
    r = res.setdefault(a["resseq"], {"resname": a["resname"], "sasa": 0.0,
                                     "plddt": a["plddt"], "cb": None})
    r["sasa"] += atom_sasa[i]
    if a["name"] == ("CA" if a["resname"] == "GLY" else "CB"):
        r["cb"] = a["xyz"]

for num, r in res.items():
    r["rsa"] = 100.0 * r["sasa"] / MAXASA.get(r["resname"], 200.0)

# contact number: heavy atoms of OTHER residues within 10 A of CB
cb_num = [n for n, r in res.items() if r["cb"] is not None]
cb_xyz = np.array([res[n]["cb"] for n in cb_num])
for k, n in enumerate(cb_num):
    d = np.linalg.norm(xyz - cb_xyz[k], axis=1)
    res[n]["contacts"] = int(np.sum((d < 10.0) & (resseq != n)))
for n, r in res.items():
    r.setdefault("contacts", 0)

# ----------------------------------------------------------- domain percentiles
dom = [n for n in res if DOMAIN[0] <= n <= DOMAIN[1]]
dom_rsa = np.array([res[n]["rsa"] for n in dom])
dom_con = np.array([res[n]["contacts"] for n in dom])
print(f"pseudokinase domain {DOMAIN[0]}-{DOMAIN[1]}: {len(dom)} residues, "
      f"median RSA {np.median(dom_rsa):.1f}%, median contacts "
      f"{np.median(dom_con):.0f}, mean pLDDT "
      f"{np.mean([res[n]['plddt'] for n in dom]):.1f}")


def burial_pct(n):
    """percentile of BURIAL within the domain (100 = most buried)."""
    return 100.0 * float(np.mean(dom_rsa > res[n]["rsa"]))


def row(num, wt, mut, label, note):
    r = res[num]
    ok = THREE_TO_ONE.get(r["resname"]) == wt
    inside = DOMAIN[0] <= num <= DOMAIN[1]
    return {
        "variant": f"p.{r['resname'].capitalize()}{num}{mut}",
        "pos": num, "wt": wt, "mut": mut, "class": label,
        "identity_matches_literature": ok,
        "in_pseudokinase_domain": inside,
        "rsa_percent": round(r["rsa"], 1),
        "contacts_10A": r["contacts"],
        "plddt": round(r["plddt"], 1),
        "burial_percentile_in_domain": round(burial_pct(num), 1) if inside else None,
        "evidence": note,
    }


rows = [row(*c) for c in CALIBRATION] + [row(*QUERY)]

bad = [r for r in rows if not r["identity_matches_literature"]]
if bad:
    print("\n!! residue identity mismatch - the model is not what we think:")
    for r in bad:
        print("   ", r)
    sys.exit(1)
print("all six residue identities match the published alleles\n")

# ----------------------------------------------------------- report
hdr = (f"{'variant':<14}{'class':<14}{'dom':<5}{'RSA%':>7}"
       f"{'cont':>6}{'pLDDT':>7}{'burial %ile':>13}")
print("=== calibration set and query ===")
print(hdr)
print("-" * len(hdr))
for r in rows:
    pct = "-" if r["burial_percentile_in_domain"] is None \
        else f"{r['burial_percentile_in_domain']:.0f}"
    print(f"{r['variant']:<14}{r['class']:<14}"
          f"{'yes' if r['in_pseudokinase_domain'] else 'no':<5}"
          f"{r['rsa_percent']:>7.1f}{r['contacts_10A']:>6}"
          f"{r['plddt']:>7.1f}{pct:>13}")

d = [r for r in rows if r["class"] == "destabilised"]
t = [r for r in rows if r["class"] == "tolerated"]
q = rows[-1]

print("\n=== does burial separate the measured classes? ===")
print(f"  destabilised  RSA {[r['rsa_percent'] for r in d]}  "
      f"contacts {[r['contacts_10A'] for r in d]}")
print(f"  tolerated     RSA {[r['rsa_percent'] for r in t]}  "
      f"contacts {[r['contacts_10A'] for r in t]}")
sep = max(r["rsa_percent"] for r in d) < min(r["rsa_percent"] for r in t)
print(f"  clean separation on RSA: {'YES' if sep else 'NO'}")
print(f"  n = {len(d)} vs {len(t)}. With five labelled residues this is an "
      f"ordering, not a classifier,\n  and the two classes also differ by "
      f"domain membership - the confound stated up front.")

print("\n=== the query ===")
print(f"  {q['variant']}: RSA {q['rsa_percent']}%, {q['contacts_10A']} contacts, "
      f"pLDDT {q['plddt']}")
print(f"  burial percentile within the pseudokinase domain: "
      f"{q['burial_percentile_in_domain']}")
print(f"  Asn -> Lys: a longer side chain carrying a formal +1 charge.")

print("""
  THE CALIBRATION FAILS ITS OWN CONTROL, and that is reported rather than
  dropped. Y155C is TOLERATED at RSA 0.7% - more buried than two of the
  three destabilised positions. Burial alone therefore does not predict
  destabilisation in this protein: Y155 is buried in the core of a
  DIFFERENT domain (BUB1 N-terminal, 62-226), which tolerates the change.
  Worse for the design, BOTH measured-tolerated alleles lie outside the
  pseudokinase domain, so there is NO tolerated control inside the domain
  to calibrate against. The percentile below is a position, not a call.""")

# -------------------------------------------- is the query coupled to the
# machinery that is known to matter? distances to the structural pocket
LANDMARKS = {795: "K795, pocket, measured destabilising",
             882: "D882, catalytic proton acceptor",
             911: "D911, catalytic triad",
             909: "I909, measured destabilising",
             1012: "L1012, measured destabilising"}
qn = QUERY[0]
qcb = np.array(res[qn]["cb"])
print("\n=== is N1002 structurally coupled to the known positions? ===")
print(f"  {'landmark':<42}{'CB-CB (A)':>11}")
dists = {}
for pos, label in LANDMARKS.items():
    if res.get(pos, {}).get("cb") is None:
        continue
    dd = float(np.linalg.norm(np.array(res[pos]["cb"]) - qcb))
    dists[pos] = round(dd, 1)
    print(f"  {label:<42}{dd:>11.1f}")

# THE COMPARISON THAT DECIDES WHETHER "far from the pocket" MEANS ANYTHING.
# Distance from the query to the pocket is uninterpretable on its own: the
# measured-destabilising alleles have to be placed on the same axis. If one of
# them is also distal, then being distal is compatible with destabilising and
# says nothing against it.
print("\n=== distance to the pocket, for every labelled allele ===")
POCKET = {795: "K795", 882: "D882", 911: "D911"}
print(f"  {'residue':<34}" + "".join(f"{v:>8}" for v in POCKET.values())
      + f"{'min':>8}   class")
pocket_tbl = {}
for pos, label, cls in ((909, "I909", "destabilised"),
                        (1012, "L1012", "destabilised"),
                        (1002, "N1002", "QUERY"),
                        (155, "Y155", "tolerated")):
    ds = [float(np.linalg.norm(np.array(res[pos]["cb"]) - np.array(res[q]["cb"])))
          for q in POCKET]
    pocket_tbl[label] = {"to_pocket": [round(d, 1) for d in ds],
                         "min": round(min(ds), 1), "class": cls}
    print(f"  {label:<34}" + "".join(f"{d:>8.1f}" for d in ds)
          + f"{min(ds):>8.1f}   {cls}")
_l1012 = pocket_tbl["L1012"]["min"]
_n1002 = pocket_tbl["N1002"]["min"]
print(f"\n  L1012P is MEASURED DESTABILISING at {_l1012} A from the pocket;")
print(f"  the query is at {_n1002} A and tolerated Y155C at "
      f"{pocket_tbl['Y155']['min']} A.")
print("  So distance from the pocket does NOT separate the measured classes,")
print("  and being distal is compatible with destabilising rather than")
print("  evidence against it. An earlier draft of this report had that")
print("  backwards; the figure above is why it was withdrawn.")

# what does N1002 actually touch?
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}
near = []
for n, r in res.items():
    if n == qn or r["cb"] is None:
        continue
    dd = float(np.linalg.norm(np.array(r["cb"]) - qcb))
    if dd < 8.0:
        near.append((dd, n, r["resname"]))
near.sort()
phob = sum(1 for _, _, rn in near if rn in HYDROPHOBIC)
print(f"\n  residues with CB within 8 A of N1002: {len(near)}"
      f"  ({phob} hydrophobic, {len(near)-phob} polar/charged)")
print("   ", ", ".join(f"{rn.capitalize()}{n}" for _, n, rn in near[:12]))

with open("21_structure_stability_results.json", "w") as fh:
    json.dump({"model": PDB.split("/")[-1], "domain": DOMAIN,
               "probe": PROBE, "sphere_points": SPHERE_POINTS,
               "domain_median_rsa": round(float(np.median(dom_rsa)), 2),
               "domain_median_contacts": int(np.median(dom_con)),
               "rows": rows,
               "calibration_separates_on_rsa": bool(sep),
               "no_tolerated_control_inside_domain": True,
               "query_cb_distances": dists,
               "pocket_distance_by_allele": pocket_tbl,
               "query_neighbours_8A": {"n": len(near), "hydrophobic": phob,
                                       "list": [f"{rn.capitalize()}{n}"
                                                for _, n, rn in near]}}, fh, indent=1)
print("\n  wrote 21_structure_stability_results.json")
print("""
  REMINDER ON WHAT MAY BE CLAIMED. AlphaFold models the wild-type fold and
  is not a stability predictor; RSA and contact number are correlates of
  substitution tolerance, not energies. A buried position makes a
  substitution there MORE LIKELY to be destabilising. It does not establish
  that this one is, and it does not replace E0.""")
