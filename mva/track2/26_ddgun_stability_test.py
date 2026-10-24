"""Predict ddG for p.Asn1002Lys with DDGun, against the measured controls.

WHY DDGUN AND NOT A TRAINED MODEL. DDGun (Montanucci, Capriotti, Birolo,
Casadio et al., BMC Bioinformatics 2019) is an UNTRAINED, knowledge-based
ddG predictor: it combines substitution scores, hydrophobicity difference
and local sequence-profile conservation, with no parameters fitted to any
experimental ddG set. That matters here. A trained predictor could have
seen K795A, I909T or L1012P during training, which would make the
calibration circular. An untrained method cannot have memorised them, so
the controls mean what they appear to mean.

THE DESIGN, unchanged from scripts 21 and 25. Five BUB1B alleles whose
behaviour is MEASURED are predicted alongside the query, so the question
is not "what number comes out for N1002K" but "does the method order the
alleles we already know, and where does the query fall among them?" A
predictor that cannot separate the knowns says nothing about the unknown,
and that outcome is reported rather than buried.

  DESTABILISED, measured:  K795A, I909T, L1012P
  TOLERATED, measured:     Y155C, R550Q
  QUERY:                   N1002K

CONVENTION: DDGun reports ddG in kcal/mol with NEGATIVE = destabilising.

THE PIPELINE, and its one external dependency. DDGun needs a sequence
profile, which needs a multiple alignment of homologues. This script
fetches BUB1B homologues from UniProt, aligns them with EBI Clustal Omega
over REST, puts the human sequence first so that profile positions are
human residue numbers, and builds the profile from that. No personal data
and no challenge data is sent anywhere: the query is a public UniProt
accession and the contact address given to EBI is a non-personal
placeholder, never the submitter's address.

WHAT A NUMBER FROM THIS CANNOT DO. ddG is thermodynamic folding stability.
The report needs steady-state abundance in a dividing human cell, which
also depends on chaperone availability and on the proteasome - neither of
which this calculation models. Sequence-profile methods also degrade on
long multi-domain proteins with disordered regions, which BubR1 is. A
prediction here raises or lowers a prior. It does not replace E0.

Run: python3 -I 26_ddgun_stability_test.py
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import warnings

warnings.filterwarnings("ignore")

# DDGun ships against an older Biopython; newer releases removed
# three_to_one. Shim it before importing, rather than pinning Biopython down.
import Bio.PDB.Polypeptide as _pp
from Bio.Data.IUPACData import protein_letters_3to1 as _M
if not hasattr(_pp, "three_to_one"):
    _pp.three_to_one = lambda r: _M[(r if isinstance(r, str)
                                     else r.get_resname()).capitalize()]

from ddgun import ddgun_seq                      # noqa: E402
from ddgun.aa import Substitution                # noqa: E402
from ddgun.profile import Profile                # noqa: E402

UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}
ACC = "O60566"
CLUSTALO = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
CONTACT = "noreply@example.org"          # placeholder; never the submitter's
MSA_PATH = "26_bub1b_msa.tsv"
MAX_HOMOLOGS = 60

VARIANTS = [
    ("K795A",  "destabilised", "denaturation 34.7 C vs 39.1 C WT"),
    ("I909T",  "destabilised", "5-10 fold protein loss"),
    ("L1012P", "destabilised", "5-10 fold protein loss"),
    ("Y155C",  "tolerated",    "protein level unaffected"),
    ("R550Q",  "tolerated",    "protein level unaffected"),
    ("N1002K", "QUERY",        "never assayed"),
]


def fetch(url, data=None, timeout=120):
    req = urllib.request.Request(url, data=data, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode(errors="replace")


def read_fasta(text):
    out, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


# ------------------------------------------------------- 1. the sequences
print("fetching human BubR1 and homologues from UniProt ...")
human = read_fasta(fetch(f"https://rest.uniprot.org/uniprotkb/{ACC}.fasta"))[0]
print(f"  {ACC}  {len(human[1])} aa")
if len(human[1]) != 1050:
    sys.exit(f"unexpected length {len(human[1])}, expected 1050 - aborting")

q = urllib.parse.quote('gene:BUB1B AND length:[800 TO 1300]')
txt = fetch(f"https://rest.uniprot.org/uniprotkb/search?query={q}"
            f"&format=fasta&size={MAX_HOMOLOGS}")
homologs = [(n, s) for n, s in read_fasta(txt) if ACC not in n]
print(f"  {len(homologs)} homologues retrieved")
if len(homologs) < 10:
    sys.exit("too few homologues for a usable profile - aborting rather than "
             "building a profile that cannot support a conservation signal")

# human MUST be first: Profile.from_msa indexes on the first sequence, which
# is what makes profile positions equal human residue numbers.
fasta = f">{ACC}\n{human[1]}\n" + "".join(
    f">{n.replace('|','_')}\n{s}\n" for n, s in homologs)

# ------------------------------------------------------- 2. the alignment
# Done LOCALLY with Biopython's pairwise aligner, not through a web service.
# EBI's alignment REST endpoints require a working contact address, and the
# submitter's address is not ours to hand to an unrelated service. Aligning
# here also makes the step reproducible with no network and no account.
#
# Each homologue is aligned to human BubR1 and projected onto HUMAN
# coordinates, so column i of every row is human residue i+1. That is a
# star alignment against the query rather than a true progressive MSA: it
# is adequate for a position-specific frequency profile over orthologues,
# and it cannot drift the position mapping, which is the property that
# matters most here.
print("aligning homologues to human BubR1 locally (Biopython) ...")
from Bio import Align                                     # noqa: E402
from Bio.Align import substitution_matrices               # noqa: E402

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0

hs = human[1]
rows_msa = [(ACC, hs)]
for n, s_ in homologs:
    try:
        a = aligner.align(hs, s_)[0]
        gh, gq = str(a[0]), str(a[1])
    except Exception as e:
        print(f"  ! {n}: {type(e).__name__}, skipped")
        continue
    proj = [c for h, c in zip(gh, gq) if h != "-"]
    if len(proj) != len(hs):
        print(f"  ! {n}: projection length {len(proj)}, skipped")
        continue
    rows_msa.append((n.replace("|", "_"), "".join(proj)))
print(f"  {len(rows_msa)} sequences in human coordinates "
      f"({len(rows_msa)-1} homologues aligned)")
if len(rows_msa) < 11:
    sys.exit("too few aligned homologues - aborting rather than building a "
             "profile with no conservation signal")

with open(MSA_PATH, "w") as fh:
    for n, s_ in rows_msa:
        fh.write(f"{n} {s_}\n")
aln = rows_msa

# ------------------------------------------------------- 3. the profile
prof = Profile.from_msa(MSA_PATH)
n_pos = len(prof.profile)
print(f"  profile built: {n_pos} positions")
if n_pos != 1050:
    sys.exit(f"profile has {n_pos} positions, expected 1050 - the position "
             f"mapping would be wrong, so no prediction is made")

seqcol = "".join(prof.profile["SEQ"].tolist())
bad = [v for v, _, _ in VARIANTS if seqcol[int(re.findall(r"\d+", v)[0]) - 1]
       != v[0]]
if bad:
    sys.exit(f"wild-type residue mismatch for {bad} - aborting")
print("  all six wild-type residues match the profile sequence")

# ------------------------------------------------------- 4. predict
rows = []
for v, cls, note in VARIANTS:
    try:
        ddg = float(ddgun_seq(Substitution.parse(v), prof))
    except Exception as e:
        print(f"  {v}: FAILED ({type(e).__name__}: {e})")
        ddg = None
    rows.append({"variant": v, "class": cls, "ddG": ddg, "evidence": note})

print("\n=== DDGun (untrained) ddG against the measured calibration set ===")
print(f"{'variant':<10}{'class':<15}{'ddG':>9}   evidence")
print("-" * 74)
for r in rows:
    d = f"{r['ddG']:+.2f}" if r["ddG"] is not None else "n/a"
    print(f"{r['variant']:<10}{r['class']:<15}{d:>9}   {r['evidence']}")

D = [r["ddG"] for r in rows if r["class"] == "destabilised" and r["ddG"] is not None]
T = [r["ddG"] for r in rows if r["class"] == "tolerated" and r["ddG"] is not None]
Q = next((r["ddG"] for r in rows if r["class"] == "QUERY"), None)

print("\n=== does it order what we already know? ===")
sep = None
if len(D) >= 2 and len(T) >= 1:
    sep = max(D) < min(T)
    print(f"  destabilised  {['%+.2f' % x for x in D]}")
    print(f"  tolerated     {['%+.2f' % x for x in T]}")
    print(f"  separates the measured classes: {'YES' if sep else 'NO'}")
    if not sep:
        print("  -> it cannot order the alleles whose behaviour is known, so")
        print("     its number for the query carries no weight. Reported.")
else:
    print("  too few predictions returned to judge separation.")

print("\n=== the query ===")
if Q is None:
    print("  N1002K: no prediction.")
else:
    print(f"  N1002K  ddG = {Q:+.2f} kcal/mol")
    if D and T:
        inside = min(D) <= Q <= max(D)
        print(f"  destabilised span {min(D):+.2f} to {max(D):+.2f}; "
              f"tolerated span {min(T):+.2f} to {max(T):+.2f}")
        print(f"  the query falls {'INSIDE' if inside else 'OUTSIDE'} the "
              f"measured-destabilising span"
              + ("" if sep else " - but the separation test FAILED, so this "
                 "placement is not evidence"))

with open("26_ddgun_stability_results.json", "w") as fh:
    json.dump({"predictor": "DDGun seq (Montanucci 2019), untrained",
               "convention": "negative ddG = destabilising",
               "uniprot": ACC, "msa_sequences": len(aln),
               "profile_positions": n_pos,
               "separates_measured_classes": sep, "rows": rows}, fh, indent=1)
print(f"\n  wrote 26_ddgun_stability_results.json and {MSA_PATH}")
print("""
  STANDING LIMITS. ddG is thermodynamic folding stability, not the
  steady-state abundance in a dividing cell that the report actually needs;
  chaperone availability and proteasome flux are not modelled. Profile
  methods degrade on long multi-domain proteins with disordered regions,
  which BubR1 is. This moves a prior. It does not replace E0.""")
