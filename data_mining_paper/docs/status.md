# Status (update after every work session)

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
- [ ] (YOU) Approve the isolate split proposal in decision_log (genus-skew problem)
- [x] Davis 2025 main PDF in data_manifest/davis2025 (methods extracted)
- [ ] (YOU) Download SI files es4c08284_si_001.pdf + es4c08284_si_002.xlsx (browser; scripts get 403)
- [ ] (YOU) Get meta sequin reference FASTA + NTC list (sequins.bio or Anaquin resources)
- [ ] (YOU) Kaggle: run `build_reference.sh` → record CARD version + cluster/decoy counts in decision_log
- [ ] (YOU) Kaggle: `subsample_backgrounds.sh data_manifest/bg_hiseq_prjeb13831.tsv bg 5000000` (6 runs)
- [ ] (YOU) Kaggle: `subsample_backgrounds.sh data_manifest/bg_test_prjna1095031.tsv bg 5000000 50000000` (3 runs)
- [ ] (YOU) Save both outputs as Kaggle Dataset `ww-backgrounds-5M`
- [ ] (YOU) Create private GitHub repo, push this folder

## Week 2 — pilot (next for Claude to write)
- [ ] (C) mix-design + InSilicoSeq simulation script
- [ ] (C) alignment + 14-feature extraction script
- [ ] (C) gate script G1/G2/G3
- [ ] (YOU) run pilot on Kaggle (~20 samples), paste outputs to `results/pilot/`
