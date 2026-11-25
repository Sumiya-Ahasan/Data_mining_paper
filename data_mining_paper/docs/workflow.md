# Your paper workflow (Nov 24, 2026 → mid-March 2027)

The one rule that matters: **research and writing run in parallel, and nothing enters the paper unless it traces to a file in this repo.**

## 1. Every work session (5 minutes of overhead)
1. **Start:** open `docs/status.md`. Pick 1–3 unchecked items. Tell Claude "continue from status".
2. **Work:** Claude writes scripts or text; you run Kaggle, download things behind logins, decide.
3. **End:** paste outputs or errors to Claude → Claude ticks `status.md` and adds a dated line to `docs/decision_log.md`. Then commit:
   `git add . ; git commit -m "<what changed>" ; git push`

## 2. Who does what
| You | Claude |
|---|---|
| Run Kaggle notebooks, save outputs as Kaggle Datasets | Write and fix all scripts (`src/`) |
| Download files that need a browser/login (ACS SI, paywalled papers) | Verify accessions, numbers, citations (DOI/PubMed) |
| Make decisions (splits, journal, what to cut) | Propose options with trade-offs; log decisions |
| **Read the ~10 core papers yourself** (list in §5) — you must defend them to reviewers | Summarise papers, draft sections, simulated peer review |
| Talk to a supervisor/faculty co-author | Draft emails, cover letter, response to reviewers |

## 3. Timeline (16 weeks)
| Weeks | Dates | Research track | Writing track |
|---|---|---|---|
| 1–2 | Nov 24 – Dec 7 | ✅ reference built; subsample backgrounds; Claude writes pilot scripts | **Literature search + Related Work draft** |
| 3–4 | Dec 8 – Dec 21 | Pilot (~20 samples) → gates G1/G2/G3 → **go / adjust** | **Methods draft** (data, decoys, features, conformal) |
| 5–8 | Dec 22 – Jan 18 | Full simulation, alignment, features, baselines | Methods updated after full runs; figure plan |
| 9–11 | Jan 19 – Feb 8 | Conformal arms M1–M4, FDR/TPR/LoD, bootstrap | **Results** written per analysis block as soon as numbers exist |
| 12–13 | Feb 9 – Feb 22 | External sequin check (PRJNA1095031), SHAP | Results finished; **pick target journal**; read its author guidelines |
| 14 | Feb 23 – Mar 1 | Freeze code; Zenodo DOI for code | **Introduction, Discussion, Conclusions, Abstract** (written last) |
| 15 | Mar 2 – Mar 8 | — | Simulated peer review → revise; citation check; supervisor read |
| 16 | Mar 9 – Mar 15 | — | Format for journal, declarations, cover letter → **submit** (optional bioRxiv preprint) |

**Go/no-go points:** end of week 4 (pilot gates) and end of week 11 (does the method beat the TDC and tuned-cutoff baselines?). If the method fails at week 11, the paper becomes an honest "when does conformal control help / not help" paper — still publishable at Q2 if rigorous.

## 4. Paper rules (what reviewers will check)
1. **Every number** in the paper comes from a file in `results/` produced by a script in `src/` — no typing numbers by hand.
2. **Every figure** is made by a script in `src/figures/` (re-runnable).
3. **Every citation** is verified (DOI/PubMed) before it goes in; no citation from memory.
4. **Every claim of "first"/"novel"** states how the search was done (databases, date).
5. Fixed seeds everywhere; record tool versions (CARD 4.0.2, KMA, Bowtie2, InSilicoSeq) in the decision log.
6. Report failures honestly (e.g. a gate that failed and what you did).
7. Declarations required by most journals: data & code availability, AI-tool use, author contributions, funding, competing interests.

## 5. Papers you must read yourself (not just summaries)
Jin & Candès 2023 (weighted conformal selection) · Bates et al. 2023 (conformal p-values for outliers) · Benjamini & Hochberg 1995 · Wang & Ramdas 2022 (e-BH) · Elias & Gygi 2007 (target–decoy) · Davis et al. 2025 ES&T · Raphenya et al. 2022 Sci Data · Hendriksen et al. 2019 (global sewage, PRJEB13831) · one ARG-tool paper each for RGI/CARD, KMA, KARGA.
(Claude verifies every one of these before it is cited; this list is a reading list, not yet a verified bibliography.)

## 6. Tools
- **Zotero** (free) for your own reading library; Claude keeps `docs/literature/references.bib` in sync.
- **GitHub (private)** = backup + Kaggle access. **Kaggle Datasets (private)** = big data.
- Writing: Markdown in `docs/paper/` → LaTeX/Word at week 16 for the journal template.

## 7. Strongly recommended
- **Find a faculty supervisor/co-author at UIU now.** Journals and reviewers take a single-undergraduate-author methods paper less seriously; a supervisor also catches problems early.
- Weekly 15-minute review (e.g. every Sunday): what moved, which risk grew, is the timeline still real?
