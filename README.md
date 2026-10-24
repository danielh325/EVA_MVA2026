# MVA Hackathon 2026 — Track 2: mechanism and repurposing for *BUB1B* MVA1

Daniel Hwang's Track 2 submission to [Rare Disease, Real Kid: MVA Hackathon 2026](https://sagebio-rare-disease-real-kid-mva-hackathon-2026.hf.space/) (Sage Bionetworks, the MVA Society, Hugging Face and BEACON).

It builds on team **`texdata`**'s earlier submission in this repository, which identified the causal genotype in Track 1. Their work is credited below and kept apart from ours.

> **Hypotheses for follow-up, not evidence that a medicine works.** Nothing here has been tested in anyone with MVA1, and nothing substitutes for the child's oncology team.

## Read this first

| | |
|---|---|
| **Report** | [`reports/track2-mechanism-and-repurposing.md`](reports/track2-mechanism-and-repurposing.md) (HTML build: [`reports/track2-report.html`](reports/track2-report.html)) |
| **Pitch script** | [`reports/track2-video-script.md`](reports/track2-video-script.md) |
| **Structure viewer** | [`reports/track2-structure-viewer.html`](reports/track2-structure-viewer.html) |
| **Code and outputs** | [`mva/track2/`](mva/track2/); every script's printed output is committed next to it as `NN_*_output.txt` |
| **Corrections** | Appendix C of the report lists every claim we withdrew or narrowed, and why |

## What is ours

Everything under `mva/track2/`, the Track 2 report and its HTML, figures named `t2-*`, the structure viewer and PDB, the pitch script, and `journals/2026-10-20.md` and `journals/2026-10-22.md`.

| Scripts | What they do |
|---|---|
| `01`–`03` | TxGNN zero-shot scores for MVA, a node-degree check on the ranking, and GraphMask path explanations (exact search; corrected mode keeps only drug→*target* edges) |
| `04`–`07` | OptimusKG coverage of the MVA genes, drugs two hops from them in PrimeKG and OptimusKG, and the labelled direction of ligands at alternative targets |
| `09`–`16` | ChEMBL direction-of-effect audits: SIRT2 activators, record- vs title-level direction, a seeded random sample, 1,018 ClinGen genes, BUB1B CETSA records |
| `18`–`19` | Decision-impact benchmark (exploratory + held-out) and its adjudication |
| `21`–`26` | Structural context on AlphaFold, PrimeKG edge audit, and two stability predictors, neither of which gave a usable ΔΔG (DDGun failed its controls; the DDMut server returned errors) |
| `17`, `20`, `22`, `24` | Figures and the HTML build |
| `mvakg.py` | Shared PrimeKG / OptimusKG loaders used by `02`–`07` |

`01`–`07` follow analyses first written by team `texdata`. We rewrote them as our own implementation on 24 October 2026 and re-ran them. Where the new runs disagreed with the old text, the report now follows the new runs (report §3 and Appendix C).

## What is team `texdata`'s

Team `texdata` (Tadej Fius) wrote the Track 1 analysis, its reports, their planning notes and journals, and Track 2 script `08`, between 28 August and 10 September 2026. **None of it is part of this submission, and none of it is in this tree.** It remains in the git history and in their repository, [github.com/tfius/mva-hackathon-2026](https://github.com/tfius/mva-hackathon-2026).

The report cites two of their results, each labelled as theirs: the Track 1 diagnosis (report §1) and the pharmacogenomic screen from script `08` (report Appendix A). Script `08` ran on the challenge VCF, which has since been deleted, so we could not re-run or check it.

The full authorship of every file is in `git log`.

## Running Track 2

```bash
pip install -r mva/track2/requirements.txt
cd mva/track2

# 09-26: public APIs only, no local data
python3 -I 10_direction_audit.py

# 01-07: need local graph data under $MVA_DATA (default ./data, git-ignored)
#   $MVA_DATA/primekg/kg.csv                    PrimeKG release
#   $MVA_DATA/txgnn-data/                       TxGNN's PrimeKG folder
#   $MVA_DATA/txgnn-ckpt/TxGNNExplorer/         pretrained checkpoint + GraphMask gates
#   $MVA_DATA/work/                             01 writes raw scores here; 02 reads them
# OptimusKG is fetched by the optimuskg package on first use.
MVA_DATA=/path/to/data python3 05_mechanism_anchored.py
```

Run `01`–`07` without `-I`, because they import `mvakg.py` from the script directory. Two workarounds are built in. TxGNN calls `DataFrame.append`, which pandas 2 removed, so `01` restores it on top of `pd.concat`. Harvard Dataverse rejects the default `python-requests` User-Agent, so `mvakg.py` sends curl's.

There are no unit tests. Each script is an analysis, and its committed output is the result; re-running it is the check. Report §10 records which scripts were re-run cold and what changed.

## Data handling

No challenge sequence data is in this repository. The VCF, BAM and FASTQ files stayed outside it and are excluded by `.gitignore`. The report contains derived findings (variant coordinates, a pharmacogenomic profile and clinical history), as the hackathon terms provide for. Report §10 states exactly what entered AI model context.

## AI disclosure

The Track 2 analysis, code and report were written and run with Anthropic's Claude through Claude Code: Claude Opus 5, and Claude Opus 5.5 for the 24 October rewrite of `01`–`07`. Direction, scope, judgement and the decision to submit are the submitter's. The models whose outputs are cited as evidence, and how much weight each gets, are listed in report §10.

## Licences

Code and reports are CC BY 4.0, matching the hackathon terms; team `texdata`'s material is also CC BY 4.0, and this README is the attribution. PrimeKG, TxGNN, OptimusKG, ChEMBL, ClinGen, UniProt and AlphaFold DB are used under their own licences.

## Acknowledgement

To the child and family who shared a genome with strangers; to the MVA Society, Sage Bionetworks, Hugging Face and BEACON; and to team `texdata`, whose Track 1 work made Track 2 possible.
