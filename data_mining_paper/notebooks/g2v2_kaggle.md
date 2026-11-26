# G2 test v2 on Kaggle (decoy v2 + label v2) — notebook cells

**Why:** the G2 extension failed: uniform 75–85% decoys do not behave like absent ARGs that have real
homologs (conserved blocks) in the sample, and some partial genes were labelled "absent"
(decision_log 2026-11-26). This run uses the **same 174 genomes and the same ART seeds** with:
- **reference v2** (`src/build_reference_v2.sh`): same CARD 4.0.2 card_90 targets; per target 1 uniform (U),
  2 near (N, 86–89%) and 2 mosaic (M, conserved + variable blocks) decoys; no decoy has a 150-bp window ≥ 90%
  identical to any CARD gene;
- **label v2** (`pilot.py label`): strict = ≥90% HSPs on one contig cover ≥80%; fragment (any 150-bp window
  ≥ 90%) = excluded; every other target is a null (distant homologs included).
`gates.json` → `G2_decoy_validity.by_decoy_type` gives the G2 verdict for U, N, M and N+M.

**Settings:** new notebook (e.g. `dm_g2v2`) · Internet ON · Accelerator None · **Add-ons → Secrets: tick
GITHUB_TOKEN** · **+ Add Input** (two inputs):
1. Dataset `card90-decoys-ref`
2. **Notebook output** of `dm_g2ext` (Add Input → Your Work / Notebooks → dm_g2ext) — gives the 174 genomes,
   so nothing is downloaded again.

Run as **Save Version → Save & Run All (Commit)**. Expected runtime: similar to dm_g2ext (labels + simulation).

### Cell 1 — code + tools
```python
# Cell 1 — code + tools
%cd /kaggle/working
from kaggle_secrets import UserSecretsClient
tok = UserSecretsClient().get_secret("GITHUB_TOKEN")
!rm -rf /kaggle/working/repo && git clone -q https://{tok}@github.com/Sumiya-Ahasan/Data_mining_paper.git /kaggle/working/repo && git -C /kaggle/working/repo remote set-url origin https://github.com/Sumiya-Ahasan/Data_mining_paper.git && echo "CLONE OK"
%cd /kaggle/working/repo/data_mining_paper
!apt-get -qq update > /dev/null && apt-get -qq install -y bowtie2 samtools ncbi-blast+ cd-hit art-nextgen-simulation-tools > /dev/null
!pip -q install pysam
!which bowtie2 samtools blastn cd-hit-est art_illumina && ls src/build_reference_v2.sh notebooks/g2v2_kaggle.md
```

### Cell 2 — inputs + KMA
```python
# Cell 2 — inputs + KMA (reuse from the reference dataset if it runs, else build)
import os, glob
REF = "/kaggle/input/datasets/sumiyaahasan/card90-decoys-ref/ref"
assert os.path.exists(f"{REF}/card90.fasta.clstr"), "REF path wrong: run  !find /kaggle/input -maxdepth 6 -name card90.fasta.clstr"
G = [d for d in glob.glob("/kaggle/input/**/g2ext/genomes", recursive=True)]
assert G, "dm_g2ext output not found: run  !find /kaggle/input -maxdepth 6 -type d -name genomes"
GEN = G[0]; print("genomes:", GEN, len(glob.glob(f"{GEN}/*.fna")))
k = glob.glob("/kaggle/input/**/kma/kma", recursive=True)
ok = False
if k:  # dataset files can lose the executable bit -> copy + chmod, then test
    os.system(f"cp {k[0]} /kaggle/working/kma_bin && chmod +x /kaggle/working/kma_bin")
    ok = os.system("/kaggle/working/kma_bin -v > /dev/null 2>&1") == 0
    k = ["/kaggle/working/kma_bin"]
if not ok:
    os.system("git clone -q https://bitbucket.org/genomicepidemiology/kma.git /kaggle/working/kma && make -s -C /kaggle/working/kma")
    k = ["/kaggle/working/kma/kma"]
os.environ["KMA"] = k[0]; print("KMA:", k[0]); os.system(f"{k[0]} -v")
```
The genome count printed must be **174**.

### Cell 3 — build reference v2 (~10–30 min)
```python
# Cell 3 — build reference v2 (same targets, decoy v2)
R2 = "/kaggle/working/ref2"
!bash src/build_reference_v2.sh {REF} {R2} 2>&1 | tee /kaggle/working/ref2_build_log.txt | tail -n 12
```
Must end with `REF V2 DONE`. The lines above it show how many decoys of each type (U/N/M) were kept.

### Cell 4 — run G2 test v2 (stops at the first failing step; log in g2v2/run_log.txt)
```python
# Cell 4 — run G2 test v2 (label v2 + simulate + features + gates on reference v2)
W = "/kaggle/working/g2v2"
D = "data_manifest/g2ext_design.tsv"
A = f"--ref {R2} --work {W} --design {D} --threads 4"
!mkdir -p {W} && ln -sfn {GEN} {W}/genomes && ( python src/pilot/pilot.py label {A} && python src/pilot/pilot.py simulate {A} && python src/pilot/features.py --ref {R2} --work {W} --design {D} --out {W}/features.tsv.gz && python src/pilot/gates.py --features {W}/features.tsv.gz --out {W}/gates.json && echo "G2V2 FINISHED OK" ) 2>&1 | tee {W}/run_log.txt
```

### Cell 5 — pack ONLY the small results for Claude
```python
# Cell 5 — pack ONLY the small results for Claude
!cd /kaggle/working && tar czf /kaggle/working/g2v2_results.tgz ref2_build_log.txt ref2/decoys_raw.fasta.tsv ref2/cdhit_2d.log ref2/cdhit_decoys.log g2v2/run_log.txt g2v2/gates.json g2v2/features.tsv.gz g2v2/null_totals.tsv g2v2/labels/isolate_args.tsv 2>/dev/null; ls -la /kaggle/working/g2v2_results.tgz
```
Download only **`g2v2_results.tgz`** (Output tab → that file's ⋮ → Download), **not** "Download all".
Keep this notebook's output: if decoy v2 passes, `ref2/` becomes the reference for the main run.

### If something fails
Send Claude the last ~30 lines of the failing cell's output. A new version starts from an empty working folder
(keep both inputs attached).
