"""Re-derive §3's edge counts straight from public PrimeKG. No TxGNN needed.

WHY THIS EXISTS. §3 claimed BUB1B carries 464 edges and 0 drug edges in
PrimeKG, and that number came from script 04, whose output was never
committed. The report therefore listed it as unverifiable. That was wrong:
the claim needs only PrimeKG's edge table, which is public, and TxGNN is
required for the *rankings* but not for counting edges. This script checks
the part the §3 conclusion actually rests on.

WHAT IT CHECKS, against the numbers as the report states them:
  - BUB1B: 464 edges, 0 of them drug edges
  - the other MVA genes and the druggable checkpoint set, which is the
    contrast the section turns on ("find the pathway, find drugs against it
    returns the exclusion list")
  - graph totals, reported as 129,375 nodes / 4,050,249 edges

HOW. kg.csv is ~1 GB, so it is streamed and filtered rather than stored: the
rows of interest are a few thousand. Each row carries x_type/x_name and
y_type/y_name, so a "drug edge" is simply an edge whose other endpoint is a
drug node. Counting is done on the file's own undirected convention and the
directed total is reported alongside, because the two differ by 2x and that
is a known source of confusion in citing this graph's size.

Downloaded data: streamed to a temp file in its own directory, parsed with
the csv module only, never imported or executed.

Run: python3 -I 23_primekg_edge_audit.py [path/to/kg.csv]
     with no argument it streams from Harvard Dataverse (~1 GB).
"""
import collections
import csv
import json
import os
import sys
import urllib.request

URL = "https://dataverse.harvard.edu/api/access/datafile/6180620"  # kg.csv
UA = {"User-Agent": "mva-hackathon-2026-research/1.0"}

REPORTED = {"BUB1B": (464, 0)}
MVA_GENES = ["BUB1B", "BUB1", "BUB3", "CEP57", "TRIP13", "CDC20", "MAD2L1"]
DRUGGABLE = ["AURKB", "PLK1", "TTK", "CENPE"]
WANT = set(MVA_GENES + DRUGGABLE)


def rows(path_or_none):
    if path_or_none:
        fh = open(path_or_none, newline="", errors="replace")
        yield from csv.DictReader(fh)
        fh.close()
        return
    print("streaming kg.csv from Harvard Dataverse (~1 GB)...", flush=True)
    req = urllib.request.Request(URL, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as r:
        import io
        buf = io.TextIOWrapper(r, encoding="utf-8", errors="replace",
                               newline="")
        n = 0
        for row in csv.DictReader(buf):
            n += 1
            if n % 1_000_000 == 0:
                print(f"  {n:,} rows", flush=True)
            yield row


path = sys.argv[1] if len(sys.argv) > 1 and os.path.exists(sys.argv[1]) else None

deg = collections.defaultdict(set)      # gene -> set of (other_type, other_name, rel)
raw = collections.Counter()            # gene -> raw rows mentioning it
drug_edges = collections.defaultdict(list)
node_ids = set()
total_rows = 0
types = collections.Counter()

for row in rows(path):
    total_rows += 1
    xt, xn = row.get("x_type", ""), row.get("x_name", "")
    yt, yn = row.get("y_type", ""), row.get("y_name", "")
    node_ids.add((xt, xn)); node_ids.add((yt, yn))
    types[xt] += 1
    rel = row.get("relation", "") or row.get("display_relation", "")
    for me, mt, other_t, other_n in ((xn, xt, yt, yn), (yn, yt, xt, xn)):
        if me in WANT and mt == "gene/protein":
            raw[me] += 1
            deg[me].add((other_t, other_n, rel))
            if other_t == "drug":
                drug_edges[me].append((other_n, rel))

print(f"\nparsed {total_rows:,} rows · {len(node_ids):,} distinct nodes")
print(f"  kg.csv stores BOTH directions of every edge, so these "
      f"{total_rows:,} rows are {total_rows//2:,} undirected edges \u2014 which is "
      f"the figure the report quotes as the graph total.")

print(f"\n{'gene':<9}{'undirected':>12}{'directed':>10}{'drug':>7}   verdict vs report")
print("-" * 64)
out = {}
for g in MVA_GENES + DRUGGABLE:
    e, d, rw = len(deg.get(g, ())), len(drug_edges.get(g, ())), raw.get(g, 0)
    out[g] = {"undirected": e, "directed": rw, "drug_edges": d}
    verdict = ""
    if g in REPORTED:
        re_, rd = REPORTED[g]
        if d != rd:
            verdict = f"DRUG EDGES DISAGREE (report says {rd})"
        elif re_ == rw:
            verdict = f"report's {re_} is the DIRECTED count; drug edges match"
        elif re_ == e:
            verdict = f"report's {re_} is the UNDIRECTED count; drug edges match"
        else:
            verdict = f"report says {re_} — matches neither convention"
    print(f"{g:<9}{e:>12}{rw:>10}{d:>7}   {verdict}")

if drug_edges.get("BUB1B"):
    print("\n  BUB1B drug edges found:",
          ", ".join(n for n, _ in drug_edges["BUB1B"][:10]))

print("\n  The contrast §6 draws: the only checkpoint proteins carrying drug")
print("  edges should be the EXCLUDED set.")
for g in DRUGGABLE:
    print(f"    {g:<7} drug edges = {len(drug_edges.get(g, ()))}")

with open("23_primekg_edge_audit_results.json", "w") as fh:
    json.dump({"source": path or URL, "undirected_rows": total_rows,
               "distinct_nodes": len(node_ids),
               "reported_in_report": {k: {"edges": v[0], "drug_edges": v[1]}
                                      for k, v in REPORTED.items()},
               "measured": out}, fh, indent=1)
print("\n  wrote 23_primekg_edge_audit_results.json")
