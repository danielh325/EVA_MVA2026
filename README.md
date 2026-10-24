# Rare Disease, Real Kid — MVA Hackathon 2026

Submission for the [Rare Disease, Real Kid: MVA Hackathon 2026](https://sagebio-rare-disease-real-kid-mva-hackathon-2026.hf.space/) — Sage Bionetworks, the MVA Society, Hugging Face and BEACON.

---

## Track 1 — solved

**Biallelic *BUB1B*, mosaic variegated aneuploidy syndrome 1** (MIM 257300, ORPHA:1052).

| | Allele A | Allele B |
|---|---|---|
| GRCh38 | `chr15:40209701 T>G` | `chr15:40220612 T>G` |
| HGVS | `ENST00000287598.11:c.2210T>G` | `ENST00000287598.11:c.3006T>G` |
| Protein | `p.Leu737Ter` | `p.Asn1002Lys` |
| Frequency | gnomAD genomes 3.3 × 10⁻⁵ | **1 allele in 1,461,878** gnomAD v4 exomes |
| Support | ClinVar Pathogenic/Likely_pathogenic for MVA1 | AlphaMissense 0.923, kinase domain |

> **Scored 100.0/100 rank points, F-max 1.000, full match at rank 1** against the clinically confirmed answer key. One null allele and one hypomorphic missense — the configuration MVA1 requires, since biallelic truncating *BUB1B* is not viable.

Allele B is in no clinical database. It surfaced only because the locus was worked exhaustively — 175 non-reference calls across *BUB1B* ±50 kb, classified by exon position — rather than filtered by predicted consequence. A conventional rare-variant chain finds the nonsense and stops.

An **Exomiser run given only the eight HPO terms and the sex**, with no gene panel, ranked *BUB1B* first of 363 genes and the compound-heterozygous pair second, matched to *Mosaic variegated aneuploidy syndrome*.

## Track 2 — mechanism and repurposing

Seven candidates reaching one target — **more working BubR1** — from independent directions: supply the cofactor (NAD⁺/SIRT2), inhibit the writer of the destabilising acetyl mark (CBP), rescue the null allele (PTC readthrough), correct it genomically, or **stabilise the mutant protein directly**. The first four need no drug that binds BubR1. The last one does, and it exists.

**The premise is not ours, and it was measured.** Pseudokinase-domain MVA missense alleles pushed back to wild-type levels show *fully restored* checkpoint function — they *"do not impose significant constraints on BUBR1 function other than affecting overall BUBR1 protein abundance"* ([Suijkerbuijk 2010, PMID 20516114](https://pubmed.ncbi.nlm.nih.gov/20516114/)). The lesion in this allele class is a **destabilised fold cleared by the proteasome**, not a broken protein.

Full analysis, candidate table, contraindications and validation plan in [`reports/track2-mechanism-and-repurposing.md`](reports/track2-mechanism-and-repurposing.md).

## Findings

**BubR1 is not undrugged — and we got the direction wrong twice before getting it right.** *BUB1B* has 464 edges in PrimeKG and **zero drug edges** — as do BUB1, BUB3, CEP57 and TRIP13, every gene the MVA node touches. OptimusKG, four times denser and drawn from 65 independent sources, reproduces that zero exactly. We first read the agreement as proof that no chemical matter exists. **It isn't.** ChEMBL `CHEMBL4295998` carries **12 assays and 24 activity records**, including approved **ceritinib at IC50 13 µM**. We then wrote that every record pushes the protein *down*. **That was wrong too.** Four cell-based records belong to a single CETSA assay *titled* **"Destabilization of BUB1B…"**, and **two of them report the opposite** — TG-101209 at fold-change **2.19** and BAY 1816032 at **1.33**, each carrying the curator's annotation `Increase In Heat Resistance Observed`. The source paper states it plainly: compound 8h destabilises BUB1B *"whereas TG-101209 stabilizes it"*.

**We missed them because of our own method.** The 42:1 audit counts direction by reading assay *descriptions*, and an assay titled for one direction can contain both. Direction is a property of the **record** — `standard_type`, the sign of the value, `activity_comment` — not of the title. That correction is the fourth blind spot in the report, and the only one that is about us rather than about a graph. [`14_bub1b_cetsa_audit.py`](mva/track2/14_bub1b_cetsa_audit.py)

**The 42-to-1 problem — the part that is not about BubR1.** Across all 1,970,438 ChEMBL assays, descriptions beginning "Inhibition of…" outnumber "Activation of…" by **42 to 1** (376,061 vs 8,952). At receptors the ratio *inverts* — agonist 24,287 against antagonist 17,788 — because you can mimic a ligand at a receptor and cannot easily make an enzyme go faster. Repeated as a **seeded random draw over all 5,869 human protein targets** — not the convenience sample an earlier draft used, which flattered us — **the median target has zero assays describing activation**, and so does the 75th percentile; **77.3 %** have none at all, 86.7 % are under 1 %, and the pooled ratio is **129.6 : 1**, three times the global figure. **Stratifying by ChEMBL's own protein classification tests the explanation rather than restating it:** transcription factors (which contain the nuclear receptors) top the table at 2.31 %, ion channels follow, and **enzymes sit at the bottom with 72 of 81 at exactly zero** — proteins with a native on-switch carry activation chemistry, proteins without one do not. A checkpoint scaffold has no switch, which is why our best candidates reach BubR1 indirectly. Every loss-of-function disease starts on the wrong side of this, and it is checkable in a handful of API calls before any model is run. [`10_direction_audit.py`](mva/track2/10_direction_audit.py) · [`12_direction_random_sample.py`](mva/track2/12_direction_random_sample.py).

**The druggable checkpoint proteins are the contraindicated ones.** AURKB, PLK1, TTK and CENPE carry drug edges; BUB1B carries none. A pipeline reasoning "find the pathway, find drugs against it" returns the contraindication list with a clean subgraph attached. (One ligand the graph attaches to AURKB is reversine — whose aneuploidy-inducing activity is in fact MPS1/TTK inhibition. We repeated that annotation uncritically at first, in a report arguing that annotations need auditing.)

**The domain we called "the kinase domain" for ten drafts is a pseudokinase domain, and the word was hiding a mechanism.** BubR1's catalytic activity is vestigial; what the C-terminal domain actually does is scaffold **PP2A-B56**, which opposes Aurora B at the kinetochore ([PMID 33207204](https://pubmed.ncbi.nlm.nih.gov/33207204/)). So `p.Asn1002Lys` plausibly causes a *phosphatase-recruitment* defect — a failure to **silence** the checkpoint and align chromosomes — on top of the dose deficit from the null allele. That prediction has an uncomfortable corollary: if Aurora B is left relatively unopposed, partial Aurora B inhibition might *help*, and our contraindication list says the opposite. We did not flip it — one cell-line paper, no measurement at this substitution, and a child who is also short of BubR1 outright — but we named the tension and added the experiment (E7) that settles it. Inference, not result.

**The model's top recommendations are contraindicated by mechanism.** TxGNN ranks paclitaxel 7th of 1,801 on *indications*, vinblastine 17th. The graph knows MVA is a cancer-predisposition syndrome and retrieves sarcoma chemotherapy; it cannot know the checkpoint is already at half dose, because that lives in the variant.

**Three distinct blind spots, and only one is about size.** A *schema* gap — nicotinamide riboside does not bind SIRT2, it raises the NAD⁺ SIRT2 consumes, and no graph encodes "increases availability of a cofactor", so substrate-pool interventions are invisible at any scale. A *semantics* gap — CREBBP and EP300 are already in PrimeKG as BUB1B interactors, and the graph still cannot say that *inhibiting* them should raise BubR1. And a *coverage* gap, found by auditing our own weakest claim: we wrote "no SIRT2 activator exists" on the strength of a graph lookup, and ChEMBL **CHEMBL4462** turns out to hold 4,236 activity records including a real activation set (EC150 15–25 µM, max activation 206–257 %) that neither graph carries. The practical conclusion survived — every one of those compounds is an unnamed preclinical chemical and there is no approved or clinical-stage SIRT2 activator — but the stated basis did not.

**A pathogenicity score does not find the second allele, and the score that looks like a second opinion is not one.** The full AlphaGenome Atlas was queried against this genome — 3.98 million called SNVs, plus every possible change across *BUB1B*. Allele B ranks **306th of 3.98 million** on the Atlas variant-impact score, and that score is **60 % AlphaMissense re-entered as a feature**, which the Atlas's own SHAP release makes measurable rather than arguable. Against all 2,382 ClinVar SNVs in *BUB1B* the classified sets separate completely — benign ≤ 23.78, pathogenic ≥ 27.12 — and **allele B falls in the gap**, which is what a hypomorph should do. What the Atlas did add is a splice negative calibrated against the known splice-pathogenic variants of the same gene, and 35 named deep-intronic positions for the next family with a missing second allele. [`reports/track1-alphagenome-atlas.md`](reports/track1-alphagenome-atlas.md).

**A standard secondary-findings filter asserts negatives it has not earned.** Filtering ClinVar to `Pathogenic|Likely_pathogenic` dropped **137 non-reference calls** on this genome, twelve reviewed by expert panel. It hid Factor V Leiden, and it let us write "no fluoropyrimidine toxicity risk allele" when the child carries an expert-panel `drug_response` *DPYD* variant.

**No mosaic aneuploidy is detectable in blood, and that is consistent with the diagnosis.** Bounded below f ≈ 0.054 by B-allele frequency and f ≈ 0.097 by read depth. *Variegated* aneuploidy puts a different random chromosome in each cell, so a 30% aneuploid population across 22 autosomes leaves each chromosome below any per-chromosome test. Both naive versions of that analysis produced false positives on the GC-rich chromosomes; the filters are the result.

## Independent review

Three independent agents reviewed the dossier cold — clinical management, adversarial refutation, and experimental feasibility. They found more than a dozen self-review passes had, and the corrections are documented in the reports and [`journals/`](journals/) rather than quietly applied:

- **A hard factual error, stated three times** — "K668 is deleted on the null allele" is false; `p.Leu737Ter` stops after residue 736 and 668 < 737.
- **A safety-level misreading risk** — vincristine is a backbone agent of every standard rhabdomyosarcoma regimen and appeared in our contraindication list. §6 now opens with an explicit guard.
- **A published surveillance guideline the dossier lacked** — renal ultrasound q3m from birth to age 7, radiation avoided, HPV vaccination.
- **A finding stated backwards** — an incidental *LZTR1* variant should not prompt surveillance imaging.
- **ACMG arithmetic that did not reach its stated conclusion** — 4 points is VUS, not Likely pathogenic.

## Layout

```
mva/
  env/          setup; clear_execstack.py (PF_X on PT_GNU_STACK, for PyTorch 1.10 on kernel 7)
  src/mva/      gene panels, Ensembl helpers, VEP REST client
  track1/       00 reference · 01 normalise · 02 ClinVar · 03 locus deep-dive
                04 mosaic aneuploidy (BAF) · 05 realign · 06 coverage+SV
                07 phasing attempt · 08 aneuploidy from depth · 09 SV panel intersect
                10-17 AlphaGenome Atlas: allele lookup · VCF scan · ClinVar benchmark
                      BUB1B saturation map · figures · base editing · incidental
                      findings · SHAP double-counting audit
                exomiser/  unbiased HPO-only run
  track2/       01 TxGNN zero-shot · 02 degree-controlled ranking · 03 GraphMask paths
                04 OptimusKG coverage · 05 mechanism-anchored reachability
                06 same on OptimusKG · 07 alternative targets · 08 pharmacogenomics
                09 SIRT2 activator audit · 10 direction-of-effect audit (the 42:1 result)
                11 direction breadth — 200 targets, median zero activation assays
  results/      submission CSV, both reports, calendar
reports/        Track 1 and Track 2 submissions, AlphaGenome companion, pitch script
  figures/      BUB1B calibration and saturation-map figures
journals/       daily working log, including what went wrong
external/       third-party checkouts (not committed)
```

## Reproducing

Everything runs from the challenge VCF and FASTQs plus public reference data. `mva/env/download_refs.sh` and `download_refs2.sh` fetch GRCh38, ClinVar, AlphaMissense, Exomiser 15.1.0 with the 2512 hg38 bundle, and the Ensembl VEP 116 cache. The AlphaGenome Atlas tables are a separate ~393 GB download from `deepmind.google.com/science/alphagenome/downloads`; point `ATLAS_DIR` at them.

Three workarounds worth knowing: pandas 2 removed `DataFrame.append`, which TxGNN needs (pin 1.5.3); Harvard Dataverse returns 403 to the default `python-requests` User-Agent while serving curl the identical URL; and PyTorch 1.10's executable-stack flag blocks loading on Linux 6+.

## Identifiers

`WGS_EX2312012` and `HGWCNDSX7` are the organisers' own library and flow-cell names, already present in the challenge dataset's filenames and required for reproducibility. They are not patient identifiers.

## Data handling

The challenge **sequence data** stays outside this repository, under `/mnt/data/`, excluded by `.gitignore`, and will be deleted at the close of the hackathon with notification to `MVAHackathon2026@synapse.org`. No attempt was made to re-identify or contact the family.

**This repository is nonetheless patient-derived, and says so rather than blurring it.** The reports contain variant coordinates, sex, gestational age, birth weight, a malignancy and a pharmacogenomic profile. The hackathon's terms provide for exactly this — submissions, code and reports are published CC BY 4.0 and the research output *is* the deliverable — but the honest framing is "the raw data stays private and the derived findings are published, as the terms require", not "nothing patient-derived leaves this machine".

## AI and machine learning

*Disclosed per the hackathon rules as updated 28 August 2026, and extended past the minimum: several models produced **evidence** here rather than assisting with it, and a reader should be able to see which conclusion rests on which.*

**Assistant:** Anthropic **Claude Opus 5** via **Claude Code** (agentic CLI, OAuth), **Claude Max** plan. It wrote and ran the whole pipeline and drafted the reports — every quantitative claim is produced by code in this repository and reproducible from it. Independent review was also AI-driven: three agents reviewed the dossier cold and found a hard factual error, a safety-level misreading risk, a missing surveillance guideline, a reversed finding and an ACMG miscalculation. Those corrections are documented in the reports and [`journals/`](journals/) rather than silently applied. *Data-handling setting is a per-account privacy preference and is recorded in the reports by the submitter.* No challenge data was pasted into a chat interface — the analysis ran locally against files on disk.

**Models whose output is cited as evidence:**

| Model | Type | Role | How it is treated |
|---|---|---|---|
| **AlphaMissense** | deep learning | 0.923 on `p.Asn1002Lys` | a primary pathogenicity argument (CC BY-NC-SA, non-commercial terms) |
| **REVEL** | ensemble ML | 0.472, same variant | reported as **disagreeing**; in ClinGen's indeterminate band, so PP3 is not claimed |
| **MVP** | deep learning | 0.852 | supporting only |
| **SpliceAI** | deep learning | via Exomiser | **not run standalone**; the standalone splice check is AlphaGenome's |
| **AlphaGenome Atlas** splicing | deep learning, precomputed | 0.087 / 0.048 on the two alleles, against 2.54–3.40 for the known splice-pathogenic records in the same gene; genome-wide scan clean | the **splice negative**, SNVs only (non-commercial terms) |
| **AlphaGenome Atlas** AVI | linear ensemble incl. AlphaMissense | PHRED 33.8 / 25.6; allele B ranks 306th of 4 M calls | **given no weight** — its own SHAP makes allele B's score 60 % AlphaMissense, so citing both would double count; a score does not find allele B, the recessive pairing does |
| **Exomiser** HiPhive + ACMG engine | semantic similarity; rule automation | *BUB1B* 1st of 363 from HPO alone | independence **qualified** (it reads ClinVar); ACMG output **corrected twice** |
| **TxGNN** + GraphMask | relational GNN; graph XAI | zero-shot repurposing and explanations | degree artefact corrected (ρ −0.51 → +0.054); top recommendations judged **wrong** on mechanism; raw explanation paths shown to be hub artefacts |

**Knowledge graphs** — the substrate of the Track 2 analysis, and the subject of several of its findings:

| Graph | Scale | Role |
|---|---|---|
| **PrimeKG** | 129,375 nodes / 4,050,249 edges | TxGNN's training substrate; the graph in which *BUB1B* has 464 edges and **zero drug edges** |
| **OptimusKG** | 190,531 nodes / 21,813,816 edges, 65 resources, 18 ontologies | **independent replication** — four times denser, same zero drug edges at BUB1B, which is what makes the gap a property of the graphs rather than of PrimeKG; and the test separating a schema gap from a size gap |

Neither graph was modified; both were queried as distributed. Scripts in [`mva/track2/`](mva/track2/).

## Licences

Code and reports CC BY 4.0, matching the hackathon terms. AlphaMissense is CC BY-NC-SA 4.0 and used under its non-commercial terms. AlphaGenome Atlas AVI scores are permissively licensed; the Atlas splicing and feature-importance sets are non-commercial only and used under those terms. ClinVar, gnomAD, Ensembl, UniProt, PrimeKG, OptimusKG, TxGNN and Exomiser under their respective licences.

## Acknowledgement

To the child and family who chose to share a genome with strangers, and to the MVA Society, Sage Bionetworks, Hugging Face and BEACON for making it possible to do anything useful with it.
