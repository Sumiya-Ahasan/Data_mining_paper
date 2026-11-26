# Method plan (v3, 2026-11-26) — single source of truth

v3: decoys failed G2 (v1 uniform and v2 homolog-structured; decision_log 2026-11-26) → calibration nulls are
labelled absent targets in simulated communities with fully known truth. Decoys stay only as a comparison arm.

**Title (working):** Calibrated, FDR-controlled detection of low-abundance ARGs in wastewater metagenomes under background shift.

## RQs
- **RQ1:** At fixed q, does conformal selection on a learned evidence score give a lower *valid* LoD50/LoD90 for rare ARGs than tuned cutoffs, existing tools, and target-decoy competition (TDC), with FDR ≤ q?
- **RQ2:** Under geographic / platform / matrix shift, which calibration (unweighted, weighted WCS / e-BH, within-sample decoys) keeps FDR ≤ q at least power cost; which features drive the shift?

## 1. Hypotheses & labels
- Unit = (sample s, ARG family j), m = |card_90| (1318). H0: family j absent from the sample.
- Labels (label v2, `src/window_id.py`), from BLAST of every genome in the sample: Y=1 strict (≥90% id HSPs cover ≥80% of a
  cluster member); excluded = fragment (a 150-bp window ≥90% id, not strict); Y=0 = everything else, distant homologs included.
- Simulated communities (`src/community.py`): every genome is known → nothing is masked; absent targets with homolog reads stay nulls.
- Real backgrounds (spike-in check only): background content unknown → evaluate only spiked targets / sequins there.

## 2. Decoys (comparison arm only — failed as null generator)
- G2 FAIL for uniform 75–85% (v1) and near/mosaic (v2) decoys: a decoy loses homolog reads to its parent target
  (decision_log 2026-11-26). Kept as the decoy-calibrated / TDC comparison and reported as a negative result.
- `src/make_decoys.py`: 3 decoys per target, 75–85% identity, substitution-only, 3rd-codon biased, no in-frame stops.
- `src/build_reference.sh`: drop decoys ≥90% id to any target (cd-hit-est-2d), dedupe decoys at 90%; index targets+decoys (KMA, Bowtie2); also targets-only index for the read-competition check.
- Key assumption (Gate G2): decoy nulls exchangeable with true target nulls.

## 2b. Simulated communities (main calibration and evaluation)
- Template per real background: Kraken2 standard-8 (2025-04-02 index) + Bracken species profile; species with a RefSeq
  assembly (≤40 per template, ≥0.05% of read pairs); ≤3 assemblies per species: pool A (train/cal/test_same), pool B
  (test_shift = unseen strains). Raphenya isolates excluded from backgrounds.
- Sample = background genomes at template share × log-normal noise (σ 0.5) of 5M nominal pairs + 3 spiked isolates (ladder).
- ART HS25 reads only from alignable regions (150-bp tiles step 50 that align to card90d, ±500 bp), at max(cov,10)× then
  thinned to the target coverage. Unclassified / non-RefSeq reads are not simulated (realism gap → gate R).
- Pilot: train 20 (DNK), cal 20 (PAK), test_same 10 (PAK), test_shift 10 (SEN). Notebook: notebooks/simcomm_kaggle.md.

## 3. Data (verified 2026-11-24, see data_manifest/)
| Role | Runs | Verified facts |
|---|---|---|
| Train bg | ERR1713349 (DNK), ERR1713397 (USA) | HiSeq 3000, PE ~2×151, 52.8M / 44.6M pairs |
| Cal bg | ERR1713384 (PAK), ERR1713346 (COL) | HiSeq 3000, 53.3M / 28.2M pairs |
| Test geo | ERR1713386 (SEN), ERR1713347 (CZE) | HiSeq 3000, 64.6M / 37.9M pairs |
| Test platform (same matrix) | SRR28520292 (CB-INF-1, PRJNA1095031) | NovaSeq 6000, US (VA) influent, 267.5M pairs |
| Test platform + matrix | SRR28520290 (CB-AS-1), SRR28520279 (CB-SE-1) | NovaSeq 6000, activated sludge / secondary effluent (labels verified in BioSample isolation_source) |
| ~~Dropped~~ | ~~SRR12326541, SRR13735366 (PRJNA648659)~~ | AS/E codes undefined in metadata and accessible text → matrix unverifiable |
| Isolates | Raphenya 2022, Zenodo 10.5281/zenodo.6543963 (AMR-Hackathon-2021/benchmarking_datasets v1.1.0) | genus-level split; target ~60/20/20 |
| External | PRJNA1095031 (Davis et al. 2025 ES&T) | 30 runs: 10 INF, 10 AS, 10 SE; NovaSeq; 242–426M pairs (73–129 Gb) each |

- Mixes: 3 isolates per spiked sample at different coverage rungs (0.1–10×); ~120 train / 60 cal / 120 test samples.
- Reads: InSilicoSeq, platform-matched error model (verify read length; else `iss model`).

## 4. Score
15 features (no is_decoy feature). LightGBM μ̂(X) on the train role, labels Y (target nulls + positives).

## 5. Calibration & selection
p_j = (1 + #{i∈C0: V_i ≥ V_j}) / (n+1); C0 = all labelled cal target nulls (zero-read ones count); zero-read test unit → p=1.
| Arm | Calibration nulls | Procedure |
|---|---|---|
| M1 | pooled labelled cal target nulls | BH per sample (`src/conformal_eval.py`) |
| M1s | same, per stratum (CARD mechanism efflux vs other; Mondrian) | BH per sample |
| M2 | M1 + w(x) | WCS (Jin & Candès 2023) — NOT plain BH (to write after the pilot) |
| M3 | M1 + w(x) | weighted conformal e-values + e-BH (verify construction) |
| M4 | cal decoys | BH — comparison arm (expected to fail) |
Weights: L2 logistic domain classifier on hits (cal vs test, unlabelled features only), clipped at 95th pct; n_eff = (Σw)²/Σw².
Power rule: need n_eff + 1 ≥ m/(q·k).

## 6. Baselines
rgi bwt; KMA default; KARGA; cutoff tuned on cal to FDP ≤ q; LightGBM+Platt with tuned threshold; **TDC (+1 correction) on same decoys**.

## 7. Evaluation
FDR (mean FDP, cluster bootstrap over mixes) vs q ∈ {0.01,0.05,0.10,0.20}; logistic TPR-vs-log(coverage) → LoD50/LoD90 with CI, reported only when FDR CI ≤ q; ΔFDR_shift; SHAP + domain-classifier coefficients.
External (PRJNA1095031): Mix A sequins = positives; **NTC sequins (in reference, absent from Mix A) = real nulls** (paper reports one NTC, S1067_MR_008_A, with 0–981 erroneous reads/sample); plus decoys. Sequins have no RefSeq homology → this validates FDR calibration of the *detection* procedure, not ARG biology. Subsample deeper (e.g. 20M pairs) and only a few runs (each is ~75–130 Gb). Use runs DISJOINT from the three test backgrounds (CB-*-1); choose them after reading the SI spike-level table so the validation set spans spike levels.
- Background subsampling for PRJNA1095031: stream first 50M pairs, sample 5M at random (`subsample_backgrounds.sh ... 5000000 50000000`); report the tile-order bias. Pilot check: sequin reads should not map to card_90 or decoys.

## 8. Pilot gates (weeks 2–4)
| Gate | Test | Pass | Fail → |
|---|---|---|---|
| G1 | n cal nulls | ≥ m/(q·k) at q=0.1 | more cal communities |
| ~~G2~~ | decoy vs target-null tails | FAILED twice → decoys dropped as calibration | (done) |
| V1 | FDR on test_same (no shift), M1/M1s | mean FDP ≤ q (q=0.05, 0.10) | check score / strata; Mondrian by target |
| V2 | FDR on test_shift, M1 | report; > q → weighted arms needed (RQ2) | — |
| R | realism: sim vs real background (units hit, reads on targets/M, read identity, efflux share) | same order of magnitude | add genomes / unclassified proxy; state limitation |

## 8b. Untested assumptions (check before claiming)
1. Simulated communities produce homolog reads like real wastewater (unclassified reads are not simulated) → gate R + real-data spike-in/sequin check.
2. Pooled calibration across targets: targets differ (efflux genes get far more null hits) → M1s Mondrian arm; report FDP per stratum.
3. Region-restricted simulation loses spurious alignments from the rest of the genome → main run: one community simulated full-genome vs restricted.
4. Units within a sample are dependent → FDR is estimated over samples (bootstrap), not assumed.
5. Kraken misassignment changes only realism, not truth (truth = genomes actually simulated).

## 9. Related work that MUST be cited/positioned
- Jin & Candès 2025, Biometrika (doi:10.1093/biomet/asaf066; arXiv 2307.09291) weighted conformalized selection.
- Bates et al. 2023 conformal p-values for outlier testing (PRDS).
- Davis, Vikesland & Pruden 2025 ES&T 59(12):6192 — LoD/LoQ via sequins (defines LoD as lowest conc detected in all 3 replicates; no FDR).
- Conformal for AMR *phenotype* (not metagenomic ARG detection): Corvelo Benz et al. RECOMB 2025 (MALDI-TOF); Inda-Díaz et al. mBio 2026; phylogeny-weighted conformal AST "four failures" (NeurIPS 2026) — shows non-exchangeability risk.
- Target-decoy FDR in proteomics (Elias & Gygi 2007) — methodological ancestor of decoys.
- Novelty searches (2026-11-24, ~10 queries, not systematic): no conformal or decoy-FDR method found for metagenomic ARG detection.
