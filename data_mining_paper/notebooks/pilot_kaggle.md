# Pilot on Kaggle — notebook cells

**Settings:** Internet ON · Accelerator None · attach your two private datasets (**+ Add Input**):
`card90-decoys-ref` and `ww-backgrounds-5m`. Run as **Save Version → Save & Run All (Commit)** so it survives closing the browser. Expected runtime ≈ 2–4 h (20 samples).

### Cell 1 — find dataset paths (check once, then edit REF/BG below)
```python
!ls /kaggle/input/
!find /kaggle/input -maxdepth 3 -name "card90_decoys.fasta" -o -maxdepth 3 -name "ERR1713384_1.fastq.gz"
```

### Cell 2 — code + tools
```python
from kaggle_secrets import UserSecretsClient
tok = UserSecretsClient().get_secret("GITHUB_TOKEN")
!rm -rf /kaggle/working/repo && git clone -q https://{tok}@github.com/Sumiya-Ahasan/Data_mining_paper.git /kaggle/working/repo
%cd /kaggle/working/repo/data_mining_paper
!apt-get -qq update > /dev/null && apt-get -qq install -y bowtie2 samtools ncbi-blast+ art-nextgen-simulation-tools > /dev/null
!pip -q install pysam
!which bowtie2 samtools blastn art_illumina
```

### Cell 3 — KMA (reuse from the reference dataset if present, else build)
```python
import os, glob
REF = "/kaggle/input/card90-decoys-ref/ref"        # <- edit if Cell 1 shows a different path
BG  = "/kaggle/input/ww-backgrounds-5m/bg"          # <- edit if Cell 1 shows a different path
k = glob.glob("/kaggle/input/**/kma/kma", recursive=True)
if not k:
    os.system("git clone -q https://bitbucket.org/genomicepidemiology/kma.git /kaggle/working/kma && make -s -C /kaggle/working/kma")
    k = ["/kaggle/working/kma/kma"]
os.environ["KMA"] = k[0]; print("KMA:", k[0]); os.system(f"{k[0]} -v")
```

### Cell 4 — run the pilot (each step skips finished work, so re-running is safe)
```python
A = f"--ref {REF} --bg {BG} --work /kaggle/working/pilot"
!python src/pilot/pilot.py genomes {A}
!python src/pilot/pilot.py label {A}
!python src/pilot/pilot.py bgprep {A}
!python src/pilot/pilot.py simulate {A}
!python src/pilot/features.py --ref {REF} --work /kaggle/working/pilot --out /kaggle/working/pilot/features.tsv.gz
!python src/pilot/gates.py --features /kaggle/working/pilot/features.tsv.gz --out /kaggle/working/pilot/gates.json
```

### Cell 5 — pack the small results for Claude
```python
!cd /kaggle/working/pilot && tar czf /kaggle/working/pilot_results.tgz gates.json features.tsv.gz labels/isolate_args.tsv bg/*.log bg/*.res samples/*/spike_pairs.txt samples/*/sample.res samples/*/spike.bam.log
!ls -la /kaggle/working/pilot_results.tgz
```
Download **`pilot_results.tgz`** (Output tab) and send it to Claude. Large files (BAMs, FASTQs) stay on Kaggle.

### If something fails
Copy the last ~30 lines of the failing cell's output to Claude. Common fixes: wrong REF/BG path (Cell 1), `art_illumina` missing (apt package name), KMA flag differences (`kma -h`).
