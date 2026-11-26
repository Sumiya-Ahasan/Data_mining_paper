# G2 extension on Kaggle — notebook cells

**Why:** the pilot could not decide G2 (are decoys a safe stand-in for absent ARGs?) because the top scores held only
3 target-null and 4 decoy hits. This run simulates each of the 174 isolates **alone at 10x, no background**,
so every absent ARG is known for sure and both target nulls and decoys get reads from the same source.
Design: `data_manifest/g2ext_design.tsv` (174 samples; made by `pilot.py g2design`).

**Settings:** new notebook (e.g. `dm_g2ext`) · Internet ON · Accelerator None · **+ Add Input**: only
`card90-decoys-ref` (no backgrounds needed). Run as **Save Version → Save & Run All (Commit)**.
Expected runtime ≈ 3–6 h (174 genome downloads + BLAST labels + simulation); Kaggle's limit is 12 h.

### Cell 1 — code + tools
```python
# Cell 1 — code + tools
%cd /kaggle/working
from kaggle_secrets import UserSecretsClient
tok = UserSecretsClient().get_secret("GITHUB_TOKEN")
!rm -rf /kaggle/working/repo && git clone -q https://{tok}@github.com/Sumiya-Ahasan/Data_mining_paper.git /kaggle/working/repo && echo "CLONE OK"
%cd /kaggle/working/repo/data_mining_paper
!apt-get -qq update > /dev/null && apt-get -qq install -y bowtie2 samtools ncbi-blast+ art-nextgen-simulation-tools > /dev/null
!pip -q install pysam
!which bowtie2 samtools blastn art_illumina && ls data_manifest/g2ext_design.tsv
```

### Cell 2 — paths + KMA (reuse from the reference dataset if it runs, else build)
```python
# Cell 2 — paths + KMA (reuse from the reference dataset if it runs, else build)
import os, glob
REF = "/kaggle/input/datasets/sumiyaahasan/card90-decoys-ref/ref"
assert os.path.exists(f"{REF}/card90_decoys.fasta"), "REF path wrong: run  !find /kaggle/input -maxdepth 6 -name card90_decoys.fasta"
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

### Cell 3 — run the G2 extension (stops at the first failing step; log in g2ext/run_log.txt)
```python
# Cell 3 — run the G2 extension (stops at the first failing step; log in g2ext/run_log.txt)
W = "/kaggle/working/g2ext"
D = "data_manifest/g2ext_design.tsv"
A = f"--ref {REF} --work {W} --design {D} --threads 4"
!mkdir -p {W} && ( python src/pilot/pilot.py genomes {A} && python src/pilot/pilot.py label {A} && python src/pilot/pilot.py simulate {A} && python src/pilot/features.py --ref {REF} --work {W} --design {D} --out {W}/features.tsv.gz && python src/pilot/gates.py --features {W}/features.tsv.gz --out {W}/gates.json && echo "G2EXT FINISHED OK" ) 2>&1 | tee {W}/run_log.txt
```

### Cell 4 — pack ONLY the small results for Claude
```python
# Cell 4 — pack ONLY the small results for Claude
!cd /kaggle/working/g2ext && tar czf /kaggle/working/g2ext_results.tgz run_log.txt gates.json features.tsv.gz null_totals.tsv labels/isolate_args.tsv samples/*/spike_pairs.txt samples/*/sample.res samples/*/spike.bam.log 2>/dev/null; ls -la /kaggle/working/g2ext_results.tgz
```
Download only **`g2ext_results.tgz`** (Output tab → that file's ⋮ → Download), not "Download all".

**What to read in gates.json:** only `G2_decoy_validity` matters here (G1/G3 are pilot gates). `verdict` is
`pass` / `fail` / `inconclusive`; `top_target_nulls` lists the highest-scoring absent ARGs (Claude checks them for
labelling misses).

### If something fails
Send Claude the last ~30 lines of the Cell 3 output. A genome that fails to download is listed as FAILED and
skipped; the rest still run.
