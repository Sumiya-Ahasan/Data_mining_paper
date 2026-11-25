# Literature table — FDR-controlled low-abundance ARG detection

Built 2026-11-25. **Every entry below was checked against Crossref/doi.org (DOI), arXiv, or Semantic Scholar for title, authors, year and venue.** BibTeX keys = `docs/literature/references.bib` (40 entries).
"Role" describes how *we* use the paper; it is based on title/abstract-level knowledge and must be confirmed when you read it. **Read** column: tick when you have actually read it.

Priority: **A** = core, must read fully · **B** = cite, skim methods/results · **C** = cite for a single fact/tool.

## 1. ARG databases, detection tools, and their fixed thresholds (the problem)
| Key | Reference | Venue | Role in our paper | Section | Pri | Read |
|---|---|---|---|---|---|---|
| alcock2023card | Alcock et al. 2023 (CARD 2023) | Nucleic Acids Res | Reference database (CARD 4.0.2 used); RGI baseline (rgi bwt) | Methods, Intro | A | [ ] |
| clausen2018kma | Clausen et al. 2018 (KMA) | BMC Bioinformatics | Aligner + ConClave score feature; KMA default baseline | Methods | A | [ ] |
| prosperi2021karga | Prosperi et al. 2021 (KARGA) | IEEE BHI 2021 | k-mer ARG baseline with its own threshold | Methods, Related | B | [ ] |
| arango2018deeparg | Arango-Argoty et al. 2018 (DeepARG) | Microbiome | Learned ARG caller; precision/recall trade-off of best-hit vs ML | Related | B | [ ] |
| marini2022amrmeta | Marini et al. 2022 (AMR-meta) | GigaScience | k-mer + metafeature classifier for short reads | Related | B | [ ] |
| lakin2019metamarc | Lakin et al. 2019 (Meta-MARC) | Commun Biol | HMM-based detection of divergent ARG sequences | Related | C | [ ] |
| bonin2023megares | Bonin et al. 2023 (MEGARes/AMR++ v3) | Nucleic Acids Res | Alternative DB/pipeline using coverage cut-offs | Related | C | [ ] |
| yin2018argsoap | Yin et al. 2018 (ARGs-OAP v2) | Bioinformatics | Widely used environmental ARG pipeline (fixed identity/length cut-offs) | Related | C | [ ] |
| bortolaia2020resfinder | Bortolaia et al. 2020 (ResFinder 4.0) | J Antimicrob Chemother | Genome-level ARG detection with identity/coverage thresholds | Related | C | [ ] |

## 2. Wastewater metagenomics, benchmarks, internal standards, detection limits
| Key | Reference | Venue | Role in our paper | Section | Pri | Read |
|---|---|---|---|---|---|---|
| hendriksen2019sewage | Hendriksen et al. 2019 | Nat Commun | Source of PRJEB13831 global sewage backgrounds (train/cal/geo-test) | Data | A | [ ] |
| davis2025qmg | Davis, Vikesland & Pruden 2025 | Environ Sci Technol | PRJNA1095031 backgrounds + sequin external check; LoD = "detected in all 3 replicates" (no error control) → our gap | Intro, Data, Discussion | A | [ ] |
| raphenya2022bench | Raphenya et al. 2022 | Sci Data | 174 isolate genomes for spike-ins | Data | A | [ ] |
| hardwick2018sequins | Hardwick et al. 2018 (meta sequins) | Nat Commun | Synthetic internal standards used by Davis 2025 | Data | B | [ ] |
| crossette2021istd | Crossette et al. 2021 | mBio | Internal-standard quantification in metagenomics | Related | B | [ ] |
| langenfeld2025qmg | Langenfeld et al. 2025 | Nucleic Acids Res | Quantitative limits in metagenomics (viruses) | Related | B | [ ] |
| ebinger2021lod | Ebinger et al. 2021 | Comput Struct Biotechnol J | Sample-specific LoD theory for clinical metagenomics — closest prior to "sample-specific LoD" | Related, Discussion | A | [ ] |

## 3. Multiple testing and FDR (the statistical target)
| Key | Reference | Venue | Role in our paper | Section | Pri | Read |
|---|---|---|---|---|---|---|
| bh1995 | Benjamini & Hochberg 1995 | J R Stat Soc B | BH procedure (arm M1, M4) | Methods | A | [ ] |
| by2001 | Benjamini & Yekutieli 2001 | Ann Stat | FDR under dependence (PRDS) — why plain BH on weighted p fails | Methods | B | [ ] |
| wang2022ebh | Wang & Ramdas 2022 | J R Stat Soc B | e-BH under arbitrary dependence (arm M3) | Methods | A | [ ] |
| barber2015knockoff | Barber & Candès 2015 | Ann Stat | Knockoffs: FDR via synthetic null copies — conceptual cousin of decoys | Related | B | [ ] |
| elias2007tdc | Elias & Gygi 2007 | Nat Methods | Target–decoy FDR (proteomics) — origin of decoy idea; TDC baseline | Related, Methods | A | [ ] |
| keich2015tdc | Keich, Kertesz-Farkas & Noble 2015 | J Proteome Res | Improved TDC FDR estimation (+1 correction) — TDC baseline details | Methods | B | [ ] |

## 4. Conformal inference, selection, covariate shift (the method)
| Key | Reference | Venue | Role in our paper | Section | Pri | Read |
|---|---|---|---|---|---|---|
| vovk2005alrw | Vovk, Gammerman & Shafer 2005 | Springer (book) | Foundation of conformal prediction | Intro | C | [ ] |
| angelopoulos2023gentle | Angelopoulos & Bates 2023 | Found Trends ML | Accessible conformal tutorial; terminology | Methods | B | [ ] |
| bates2023outliers | Bates et al. 2023 | Ann Stat | Null-only conformal p-values; PRDS → BH valid (arm M1) | Methods | A | [ ] |
| jin2023selection | Jin & Candès 2023 | JMLR 24 | Conformal selection with FDR control | Methods | A | [ ] |
| jin2025wcs | Jin & Candès 2025 | Biometrika (arXiv 2307.09291) | Weighted conformal selection (WCS) under covariate shift; weighted p not PRDS (arm M2) | Methods | A | [ ] |
| tibshirani2019shift | Tibshirani et al. 2019 | NeurIPS | Weighted conformal under covariate shift (likelihood-ratio weights) | Methods | A | [ ] |
| shimodaira2000 | Shimodaira 2000 | J Stat Plan Inference | Importance weighting under covariate shift | Methods | C | [ ] |
| bashari2025robust | Bashari, Sesia & Romano 2025 | ICML | Contaminated reference data in conformal outlier detection — relates to weight contamination | Methods, Discussion | B | [ ] |

## 5. Conformal methods in AMR / biology (closest applications; positions novelty)
| Key | Reference | Venue | Role in our paper | Section | Pri | Read |
|---|---|---|---|---|---|---|
| benz2025maldi | Corvelo Benz et al. 2025 | RECOMB 2025 (LNCS) | Conformal AMR detection from MALDI-TOF spectra (phenotype, not metagenomic ARG detection) | Related | A | [ ] |
| inda2026mbio | Inda-Díaz et al. 2026 | mBio | Confidence-based (conformal) patient-level resistance prediction | Related | B | [ ] |
| lopezdecastro2025scrna | López-De-Castro et al. 2025 | Bioinformatics | Conformal inference for scRNA-seq annotation (bioinformatics precedent) | Related | C | [ ] |

## 6. Tools and ML components (cite once in Methods)
| Key | Reference | Venue | Role in our paper | Pri | Read |
|---|---|---|---|---|---|
| fu2012cdhit | Fu et al. 2012 (CD-HIT) | Bioinformatics | card_90 clustering; decoy filtering | C | [ ] |
| langmead2012bowtie2 | Langmead & Salzberg 2012 (Bowtie 2) | Nat Methods | Second aligner | C | [ ] |
| gourle2019iss | Gourlé et al. 2019 (InSilicoSeq) | Bioinformatics | Platform-matched read simulation | C | [ ] |
| brown2022mobileog | Brown et al. 2022 (mobileOG-db) | Appl Environ Microbiol | Mobile-element linkage feature | C | [ ] |
| ke2017lightgbm | Ke et al. 2017 (LightGBM) | NeurIPS | Score model μ̂(x) | C | [ ] |
| lundberg2017shap | Lundberg & Lee 2017 (SHAP) | NeurIPS | Feature attribution (RQ2) | C | [ ] |
| lopezpaz2017c2st | Lopez-Paz & Oquab 2017 | ICLR | Classifier two-sample test (gates G2/G3) | C | [ ] |

## 7. Found but NOT yet verifiable — do not cite until checked
| Item | Status |
|---|---|
| "Phylogeny-weighted conformal prediction for genomic antibiotic susceptibility testing: four failures of a natural idea" (NeurIPS 2026, neurips.cc/virtual/2026/169434) | Title + venue seen in search; authors/proceedings record not yet retrieved (API rate limit). Highly relevant: shows conformal coverage breaking under clade shift. |

## 8. Gaps this table still has (search next)
- A benchmark comparing ARG-tool false-positive rates at low coverage (to quantify the "fixed thresholds" problem with numbers).
- Any use of decoy/null references in metagenomic read classification (none found in 2026-11-24 searches; record search terms in the paper).
- Recent (2024–2026) conformal FDR work on genomics/metagenomics beyond Section 5.

## Gap statement supported by this table (draft, 3 sentences)
ARG detection tools (§1) call genes with fixed identity/coverage/score cut-offs that carry no stated error rate, and wastewater detection limits (§2, e.g. Davis et al. 2025) are defined by replicate detection without false-discovery control. Conformal selection provides finite-sample FDR control (§3–4) and has reached AMR only for phenotype prediction (§5), not for metagenomic gene detection under background shift. We close this gap with decoy-calibrated, shift-weighted conformal selection and an FDR-coupled LoD.
