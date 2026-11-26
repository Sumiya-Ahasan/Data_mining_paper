# Status (update after every work session)

Workflow + timeline: see `docs/workflow.md`.

Legend: [x] done · [ ] todo · (YOU) = manual step for the researcher · (C) = Claude can do

## Week 1 — foundations
- [x] (C) Verify 8 background runs on ENA → `data_manifest/bg_hiseq_prjeb13831.tsv` (+ dropped_prjna648659.tsv)
- [x] (C) Verify PRJNA1095031 (30 runs) and PRJNA648659 (24 runs) → `data_manifest/`
- [x] (C) Davis 2025 key facts (Mix A, NTC sequins, LoD definition) → plan §7
- [x] (C) Zenodo isolate record exists (10.5281/zenodo.6543963)
- [x] (C) Narrow novelty check → plan §9
- [x] (C) `src/make_decoys.py` (tested locally on synthetic seqs: identity 0.76–0.85, 0 in-frame stops)
- [x] (C) `src/build_reference.sh`, `src/subsample_backgrounds.sh` (written, NOT yet run)
- [x] PRJNA648659 dropped (codes unverifiable); matrix test from PRJNA1095031 CB-INF-1 / CB-AS-1 / CB-SE-1
- [x] (C) Raphenya 2022 isolate table (174, GCF accessions) → `data_manifest/raphenya2022_strain_metadata.tsv`
- [x] Isolate split approved → data_manifest/isolate_split.tsv
- [x] Davis 2025 main PDF in data_manifest/davis2025 (methods extracted)
- [ ] (YOU) Download SI files es4c08284_si_001.pdf + es4c08284_si_002.xlsx (browser; scripts get 403)
- [ ] (YOU) Get meta sequin reference FASTA + NTC list (sequins.bio or Anaquin resources)
- [x] Kaggle build_reference.sh done (CARD 4.0.2, 1318 clusters, 3954 decoys)
- [ ] (YOU) Save reference output as Kaggle Dataset `card90-decoys-ref`
- [x] HiSeq backgrounds subsampled + validated (6 x 5M pairs, all checks pass) → results/backgrounds/bg_hiseq_check.json
- [ ] (YOU) Kaggle: `subsample_backgrounds.sh data_manifest/bg_test_prjna1095031.tsv bg 5000000 50000000` (3 runs)
- [ ] (YOU) Save both outputs as Kaggle Dataset `ww-backgrounds-5M`
- [x] GitHub repo created (files sit in subfolder data_mining_paper/; local push still rejected — tidy later)

- [x] (C) Literature table (40 verified refs) → `docs/literature/literature_table.md`
- [ ] (YOU) Read priority-A papers in the table; tick the Read column

## Week 2 — pilot (next for Claude to write)
- [x] (C) Pilot design + genomes + labels + ART simulation + alignment: src/pilot/pilot.py (data_manifest/pilot_design.tsv)
- [x] (C) Feature extraction: src/pilot/features.py (15 features + labels)
- [x] (C) Gates G1/G2/G3: src/pilot/gates.py
- [x] (YOU) Save HiSeq bg/ as private Kaggle Dataset ww-bg-5m (needed by pilot)
- [x] (YOU) Run pilot on Kaggle (dm_main v2, 2026-11-26, ~46 min, PILOT FINISHED OK)
- [x] (C) Check pilot: G1 pass (after fix), G2 inconclusive, G3 pass → results/pilot/gates_v1_relabelled.json
- [x] (C) Fix labels (spiked = positive even if background has reads) + G1 count + one-sided G2 + n_eff
- [x] (C) Write G2 extension: 174 isolates alone at 10x, no background (notebooks/g2ext_kaggle.md, data_manifest/g2ext_design.tsv)
- [x] (YOU) Run G2 extension on Kaggle (dm_g2ext, 174/174 genomes, finished OK)
- [x] (C) Check G2 extension: G2 FAIL (well powered) → label gaps (fragments, split HSPs) + real homologs with conserved blocks; results/g2ext/
- [x] (C) Decoy v2 (U/N/M types, 150-bp window cap) + label v2 (merged HSPs, fragment class, all cluster members) + per-type G2; tested locally (src/make_decoys_v2.py, filter_decoys.py, window_id.py, build_reference_v2.sh)
- [x] (YOU) Run G2 test v2 on Kaggle (dm_g2v2, 174/174, G2V2 FINISHED OK)
- [x] (C) Check G2 test v2: G2 FAIL for every decoy type (U, N, M, N+M), also with target nulls in score training and when restricted to non-efflux or inactivation-only targets → results/g2v2/
- [ ] (YOU, + supervisor) Decide the calibration-null source (decoys fail; see decision_log 2026-11-26 G2 v2)
- [ ] (YOU) Security: rotate the GitHub token (it was stored in repo/.git/config inside Kaggle outputs of all 3 notebooks; fixed in the notebook code); delete C:/Users/User/Downloads/g2v2_results.zip after Claude's check
