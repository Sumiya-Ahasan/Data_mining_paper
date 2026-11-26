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
- [ ] (C) Write G2 extension: isolate-only (no background) runs on many isolates → enough target-null tail events
- [ ] (YOU) Run G2 extension on Kaggle when ready
