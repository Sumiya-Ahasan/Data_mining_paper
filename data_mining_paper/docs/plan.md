# Method plan (v2, 2026-11-24) — single source of truth

**Title (working):** Calibrated, FDR-controlled detection of low-abundance ARGs in wastewater metagenomes under background shift.

## RQs
- **RQ1:** At fixed q, does conformal selection on a learned evidence score give a lower *valid* LoD50/LoD90 for rare ARGs than tuned cutoffs, existing tools, and target-decoy competition (TDC), with FDR ≤ q?
- **RQ2:** Under geographic / platform / matrix shift, which calibration (unweighted, weighted WCS / e-BH, within-sample decoys) keeps FDR ≤ q at least power cost; which features drive the shift?

## 1. Hypotheses & labels
- Unit = (sample s, ARG family j), m = |card_90| (~1318). H0: family j absent from spiked content.
- Y=1: in a spiked isolate's ARG set (RGI Perfect/Strict on genome). Y=0: not in any spiked isolate AND 0 reads in the unspiked background. Ambiguous (≥1 background read): tested, excluded from FDP numerator and denominator.
- Spiked sample reuses the identical 5M-pair background subsample → exact mask.
- Caveat to report: LoD applies to ARGs rare in that background.

## 2. Decoys (null generator)
- `src/make_decoys.py`: 3 decoys per target, 75–85% identity, substitution-only, 3rd-codon biased, no in-frame stops.
- `src/build_reference.sh`: drop decoys ≥90% id to any target (cd-hit-est-2d), dedupe decoys at 90%; index targets+decoys (KMA, Bowtie2); also targets-only index for the read-competition check.
- Key assumption (Gate G2): decoy nulls exchangeable with true target nulls.

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
14 features (targets and decoys identical, no is_decoy feature). LightGBM μ̂(X), train split only, decoys labelled Y=0.

## 5. Calibration & selection
p_j = (1 + #{i∈C0: V_i ≥ V_j}) / (n+1); zero-read → p=1.
| Arm | Calibration nulls | Procedure |
|---|---|---|
| M1 | pooled cal nulls | BH (valid: unweighted conformal p are PRDS) |
| M2 | same + w(x) | WCS (Jin & Candès 2023) — NOT plain BH |
| M3 | same + w(x) | weighted conformal e-values + e-BH (verify construction) |
| M4 | test sample's own decoys | BH |
Weights: L2 logistic domain classifier on **decoy hits only** (cal vs test), clipped at 95th pct; report n_eff = (Σw)²/Σw².
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
| G1 | n_eff of cal nulls | ≥ m/(q·k) at q=0.1 | d=5 decoys, more cal mixes, or report q≥0.1 |
| G2 | C2ST decoy-null vs target-null | AUC ≤ 0.6 | target nulls only; decoys only for weights; drop M4 |
| G3 | domain AUC cal vs test (nulls) | > 0.6 | RQ2 = "no correction needed" (still reportable) |

## 9. Related work that MUST be cited/positioned
- Jin & Candès 2025, Biometrika (doi:10.1093/biomet/asaf066; arXiv 2307.09291) weighted conformalized selection.
- Bates et al. 2023 conformal p-values for outlier testing (PRDS).
- Davis, Vikesland & Pruden 2025 ES&T 59(12):6192 — LoD/LoQ via sequins (defines LoD as lowest conc detected in all 3 replicates; no FDR).
- Conformal for AMR *phenotype* (not metagenomic ARG detection): Corvelo Benz et al. RECOMB 2025 (MALDI-TOF); Inda-Díaz et al. mBio 2026; phylogeny-weighted conformal AST "four failures" (NeurIPS 2026) — shows non-exchangeability risk.
- Target-decoy FDR in proteomics (Elias & Gygi 2007) — methodological ancestor of decoys.
- Novelty searches (2026-11-24, ~10 queries, not systematic): no conformal or decoy-FDR method found for metagenomic ARG detection.
