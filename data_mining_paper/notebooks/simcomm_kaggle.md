# Simulated-community pilot on Kaggle (labelled calibration nulls) — notebook cells

**Why:** decoys failed G2 for every type (decision_log 2026-11-26), so calibration nulls are now absent targets
in simulated communities where truth is known for every gene. Each real background (Kraken2 profile) gives a
community template; samples = RefSeq genomes of its species at abundances around the template + 3 spiked ARG
isolates. Nothing is masked: absent targets with homolog reads stay nulls. Roles: train 20, cal 20 (PAK
template), test_same 10 (PAK, no shift), test_shift 10 (SEN template, strains unseen in train/cal).

**Gates of this pilot** (`eval.json`, `realism.tsv`):
- **V1 validity:** test_same, target arms: mean FDP ≤ q at q = 0.05 and 0.10.
- **V2 shift:** test_shift FDP (unweighted). Above q means weighted conformal is needed (RQ2), not a failure of the pilot.
- **R realism:** simulated vs real background: target units hit, reads on targets per M, read identity, efflux share.
  Large gaps mean the simulation misses real homolog sources; report and improve.
- **Decoy arm:** expected to fail (it is the comparison).

**Settings:** new notebook `dm_simcomm` · Internet ON · Accelerator None · **Secrets: GITHUB_TOKEN** ·
**+ Add Input** (three inputs):
1. **Notebook output** of `dm_g2v2` (reference v2 + isolate labels)
2. **Notebook output** of `dm_g2ext` (174 isolate genomes)
3. Dataset `ww-bg-5m` (HiSeq backgrounds)

Run as **Save Version → Save & Run All (Commit)**. Expected runtime: about 3–5 h (not measured yet).

### Cell 1 — code + tools
```python
# Cell 1 — code + tools
%cd /kaggle/working
from kaggle_secrets import UserSecretsClient
tok = UserSecretsClient().get_secret("GITHUB_TOKEN")
!rm -rf /kaggle/working/repo && git clone -q https://{tok}@github.com/Sumiya-Ahasan/Data_mining_paper.git /kaggle/working/repo && git -C /kaggle/working/repo remote set-url origin https://github.com/Sumiya-Ahasan/Data_mining_paper.git && echo "CLONE OK"
%cd /kaggle/working/repo/data_mining_paper
!apt-get -qq update > /dev/null && apt-get -qq install -y bowtie2 samtools ncbi-blast+ art-nextgen-simulation-tools > /dev/null
!apt-get -qq install -y kraken2 > /dev/null 2>&1; apt-get -qq install -y bracken > /dev/null 2>&1; true
!pip -q install pysam lightgbm
!which kraken2 || (git clone -q https://github.com/DerrickWood/kraken2.git /kaggle/working/k2src && cd /kaggle/working/k2src && ./install_kraken2.sh /usr/local/bin > /dev/null)
!which bowtie2 samtools blastn art_illumina kraken2; which bracken || echo "no bracken -> Kraken species counts are used"
```

### Cell 2 — inputs + KMA
```python
# Cell 2 — inputs + KMA
import os, glob
R2 = [d for d in glob.glob("/kaggle/input/**/ref2", recursive=True) if os.path.exists(f"{d}/idx/card90d.1.bt2")]
assert R2, "dm_g2v2 output (ref2) not found: run  !find /kaggle/input -maxdepth 5 -type d -name ref2"
REF = R2[0]
LAB = glob.glob("/kaggle/input/**/g2v2/labels", recursive=True)
G = glob.glob("/kaggle/input/**/g2ext/genomes", recursive=True)
BG = [os.path.dirname(p) for p in glob.glob("/kaggle/input/**/ERR1713384_1.f*q*", recursive=True)]
assert LAB and G and BG, ("labels", LAB, "genomes", G, "bg", BG)
GEN = G[0]; BGD = BG[0]
print("REF", REF, "| labels", LAB[0], "| genomes", GEN, len(glob.glob(f"{GEN}/*.fna")), "| bg", BGD)
k = glob.glob("/kaggle/input/**/kma_bin", recursive=True) + glob.glob("/kaggle/input/**/kma/kma", recursive=True)
ok = False
if k:
    os.system(f"cp {k[0]} /kaggle/working/kma_bin && chmod +x /kaggle/working/kma_bin")
    ok = os.system("/kaggle/working/kma_bin -v > /dev/null 2>&1") == 0
if not ok:
    os.system("git clone -q https://bitbucket.org/genomicepidemiology/kma.git /kaggle/working/kma && make -s -C /kaggle/working/kma && cp /kaggle/working/kma/kma /kaggle/working/kma_bin")
os.environ["KMA"] = "/kaggle/working/kma_bin"; os.system("/kaggle/working/kma_bin -v")
```
The genome count printed must be **174**.

### Cell 3 — Kraken2 database (8 GB, to /tmp, not saved)
```python
# Cell 3 — Kraken2 standard-8 database (streamed into /tmp/k2)
import re, urllib.request
page = urllib.request.urlopen("https://benlangmead.github.io/aws-indexes/k2").read().decode()
urls = sorted(set(re.findall(r"https://genome-idx\.s3\.amazonaws\.com/kraken/k2_standard_08gb_\d{8}\.tar\.gz", page)))
K2URL = urls[-1] if urls else "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz"
print(K2URL)
!mkdir -p /tmp/k2 && curl -sL {K2URL} | tar xz -C /tmp/k2 && ls /tmp/k2 && df -h /tmp | tail -1
```

### Cell 4 — backgrounds, profiles, genomes, design
```python
# Cell 4 — real backgrounds -> Kraken profiles -> RefSeq genomes -> design
W = "/kaggle/working/simcomm"
C = f"--work {W} --ref {REF} --bg {BGD} --k2db /tmp/k2 --threads 4"
!mkdir -p {W}/genomes {W}/labels && for f in {GEN}/*.fna; do ln -sfn $f {W}/genomes/; done && cp -n {LAB[0]}/*.hsp2.tsv {W}/labels/
!( python src/community.py realbg {C} && python src/community.py profile {C} && python src/community.py select {C} && python src/community.py genomes {C} && python src/community.py design {C} && echo "STEP4 OK" ) 2>&1 | tee {W}/log4.txt | grep -v "^+ " | tail -n 25
!rm -f {W}/bg/*.150.fq.gz
```

### Cell 5 — labels, regions, simulation (longest step)
```python
# Cell 5 — labels for all genomes, alignable regions, community simulation
!( python src/pilot/pilot.py label --all --work {W} --ref {REF} --threads 4 && python src/community.py regions {C} && python src/community.py simulate {C} && echo "STEP5 OK" ) 2>&1 | tee {W}/log5.txt | grep -v "^+ " | tail -n 15
```

### Cell 6 — features, evaluation, realism
```python
# Cell 6 — features + conformal evaluation + realism check
!( python src/pilot/features.py --ref {REF} --work {W} --design {W}/design_samples.tsv --members {W}/design_members.tsv --out {W}/features.tsv.gz && python src/conformal_eval.py --features {W}/features.tsv.gz --ref {REF} --out {W}/eval.json && python src/community.py realism {C} && echo "SIMCOMM FINISHED OK" ) 2>&1 | tee {W}/log6.txt | grep -v "units with reads" | tail -n 20
```

### Cell 7 — pack ONLY the small results for Claude
```python
# Cell 7 — pack ONLY the small results for Claude
!cd /kaggle/working && tar czf simcomm_results.tgz simcomm/log4.txt simcomm/log5.txt simcomm/log6.txt simcomm/profiles simcomm/community_genomes.tsv simcomm/design_samples.tsv simcomm/design_members.tsv simcomm/regions/regions.tsv simcomm/labels/isolate_args.tsv simcomm/features.tsv.gz simcomm/truth.tsv.gz simcomm/null_totals.tsv simcomm/eval.json simcomm/realism.tsv simcomm/realism_units.tsv 2>/dev/null; ls -la simcomm_results.tgz
```
Download only **`simcomm_results.tgz`** (Output tab → that file's ⋮ → Download), **not** "Download all".

### If something fails
Send Claude the last ~30 lines of the failing cell. Each step skips finished work, but a new version starts
from an empty working folder (keep all three inputs attached).
