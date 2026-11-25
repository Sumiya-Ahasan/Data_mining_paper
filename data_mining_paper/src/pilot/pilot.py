#!/usr/bin/env python3
"""Pilot driver: isolate split, mix design, genome download, ARG labels, background prep,
spike-in simulation and alignment. Feature extraction: features.py. Gates: gates.py.

Subcommands (run in order; each skips work already done):
  split     -> data_manifest/isolate_split.tsv           (local or Kaggle)
  design    -> data_manifest/pilot_design.tsv            (local or Kaggle)
  genomes   -> $WORK/genomes/<biosample>.fna             (Kaggle, internet ON)
  label     -> $WORK/labels/isolate_args.tsv             (Kaggle; BLAST of CARD homologs vs genomes)
  bgprep    -> $WORK/bg/<run>_{1,2}.150.fq.gz, bg.bam, kma .res, mask   (Kaggle)
  simulate  -> $WORK/samples/<sample>/spike.bam, sample.res              (Kaggle)

Paths come from flags (defaults suit Kaggle; see notebooks/pilot_kaggle.md).
"""
import argparse, csv, gzip, os, random, re, shutil, subprocess, sys, urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LADDER = [0.1, 0.25, 0.5, 1, 2, 5, 10]
# Approved 2026-11-25 (decision_log): genus-level split
SPLIT = {
    "train": {"Klebsiella", "Staphylococcus", "Salmonella"},
    "cal": {"Escherichia", "Citrobacter", "Aeromonas", "Kluyvera"},
    "test": {"Enterobacter", "Pseudomonas", "Acinetobacter", "Enterococcus", "Providencia"},
}
# Pilot: one calibration background, one geographic-test background (both HiSeq 3000)
PILOT_BG = {"cal": "ERR1713384", "test": "ERR1713386"}
N_MIX = 10
READ_LEN = 150  # backgrounds trimmed to 150 so spike (ART HS25 150 bp) and background match


def sh(cmd, **kw):
    print("+", cmd, flush=True)
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", **kw)


def tsv_write(path, rows, cols):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t"); w.writeheader(); w.writerows(rows)


def tsv_read(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- split / design
def cmd_split(a):
    rows = []
    for r in tsv_read(REPO / "data_manifest/raphenya2022_strain_metadata.tsv"):
        genus = r["Organism"].split()[0]
        split = next(k for k, v in SPLIT.items() if genus in v)
        rows.append(dict(biosample=r["biosample"], organism=r["Organism"], genus=genus, split=split,
                         assembly=r["assembly_filename"]))
    tsv_write(REPO / "data_manifest/isolate_split.tsv", rows, list(rows[0]))
    for s in SPLIT:
        print(s, sum(r["split"] == s for r in rows))


def cmd_design(a):
    iso = tsv_read(REPO / "data_manifest/isolate_split.tsv")
    rng = random.Random(a.seed)
    rows = []
    for role, bg in PILOT_BG.items():
        pool = [r["biosample"] for r in iso if r["split"] == role]
        for i in range(N_MIX):
            picks = rng.sample(pool, 3)                      # 3 distinct isolates per mix
            rungs = [LADDER[(3 * i + j) % len(LADDER)] for j in range(3)]  # balanced, distinct rungs
            row = dict(sample=f"pilot_{role}_{i:02d}", background=bg, role=role, art_seed=a.seed * 1000 + len(rows))
            for j in range(3):
                row[f"iso{j+1}"] = picks[j]; row[f"cov{j+1}"] = rungs[j]
            rows.append(row)
    cols = ["sample", "background", "role", "iso1", "cov1", "iso2", "cov2", "iso3", "cov3", "art_seed"]
    tsv_write(REPO / "data_manifest/pilot_design.tsv", rows, cols)
    print(f"{len(rows)} pilot samples -> data_manifest/pilot_design.tsv")


def design_isolates():
    d = tsv_read(REPO / "data_manifest/pilot_design.tsv")
    return sorted({r[f"iso{j}"] for r in d for j in (1, 2, 3)})


# ---------------------------------------------------------------- genomes
def ncbi_url(asm):
    acc = asm.split("_")[0] + "_" + asm.split("_")[1]          # GCF_016903215.1
    digits = acc.split("_")[1].split(".")[0]                    # 016903215
    p = "/".join(digits[i:i + 3] for i in range(0, 9, 3))
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{acc[:3]}/{p}/{asm}/{asm}_genomic.fna.gz"


def cmd_genomes(a):
    iso = {r["biosample"]: r for r in tsv_read(REPO / "data_manifest/isolate_split.tsv")}
    want = sorted(iso) if a.all else design_isolates()
    out = Path(a.work) / "genomes"; out.mkdir(parents=True, exist_ok=True)
    for b in want:
        f = out / f"{b}.fna"
        if f.exists() and f.stat().st_size > 0:
            continue
        url = ncbi_url(iso[b]["assembly"])
        print("get", b, url, flush=True)
        with urllib.request.urlopen(url, timeout=120) as r, open(f, "wb") as w:
            w.write(gzip.decompress(r.read()))
    print(len(want), "genomes in", out)


# ---------------------------------------------------------------- ARG labels per isolate
def cluster_map(clstr):
    """member id (first token of CARD header) -> representative id."""
    m, members, rep = {}, [], None
    for line in open(clstr, encoding="utf-8"):
        if line.startswith(">Cluster"):
            for x in members: m[x] = rep
            members, rep = [], None
            continue
        name = re.search(r">(.+?)\.\.\.", line).group(1)
        members.append(name)
        if line.rstrip().endswith("*"):
            rep = name
    for x in members: m[x] = rep
    return m


def cmd_label(a):
    """Y=1 ('strict'): a CARD homolog member hits the genome at >=90% id over >=80% of its length.
    'ambiguous': best hit >=70% id over >=50% length but not strict (excluded from scoring)."""
    ref = Path(a.ref); work = Path(a.work)
    cmap = cluster_map(ref / "card90.fasta.clstr")
    homolog = ref / "card/nucleotide_fasta_protein_homolog_model.fasta"
    out = work / "labels"; out.mkdir(parents=True, exist_ok=True)
    rows = []
    for b in design_isolates() if not a.all else sorted(p.stem for p in (work / "genomes").glob("*.fna")):
        g = work / "genomes" / f"{b}.fna"
        hits = out / f"{b}.blast.tsv"
        if not hits.exists():
            sh(f"makeblastdb -in {g} -dbtype nucl -out {out}/{b}.db > /dev/null")
            sh(f"blastn -task dc-megablast -query {homolog} -db {out}/{b}.db -evalue 1e-10 "
               f"-max_target_seqs 5 -num_threads {a.threads} "
               f"-outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' > {hits}")
        best = {}  # rep -> (status_rank, pident, qcov)
        for line in open(hits):
            q, s, pid, L, qlen, *_ = line.split("\t")
            pid, cov = float(pid), int(L) / int(qlen)
            rep = cmap.get(q)
            if rep is None:
                continue
            st = 2 if (pid >= 90 and cov >= 0.8) else (1 if (pid >= 70 and cov >= 0.5) else 0)
            if st and (rep not in best or (st, pid) > best[rep][:2]):
                best[rep] = (st, pid, round(cov, 3))
        for rep, (st, pid, cov) in best.items():
            rows.append(dict(biosample=b, unit=rep, status="strict" if st == 2 else "ambiguous", pident=pid, qcov=cov))
        print(b, "strict:", sum(v[0] == 2 for v in best.values()), "ambiguous:", sum(v[0] == 1 for v in best.values()), flush=True)
    tsv_write(out / "isolate_args.tsv", rows, ["biosample", "unit", "status", "pident", "qcov"])


# ---------------------------------------------------------------- background prep
def bt2(r1, r2, idx, bam, threads):
    sh(f"bowtie2 --very-sensitive-local -p {threads} -x {idx} -1 {r1} -2 {r2} --no-unal 2> {bam}.log "
       f"| samtools view -b -F 0x900 -o {bam} -")


def kma(r1, r2, db, prefix, threads):
    sh(f"{os.environ.get('KMA', 'kma')} -ipe {r1} {r2} -o {prefix} -t_db {db} -t {threads} -nf -na > {prefix}.log 2>&1")


def bg_file(bgdir, run, mate):
    """Kaggle unzips .fastq.gz on dataset upload, so accept .fastq.gz or .fastq."""
    for ext in (".fastq.gz", ".fastq", ".fq.gz", ".fq"):
        f = Path(bgdir) / f"{run}_{mate}{ext}"
        if f.exists():
            return f
    raise FileNotFoundError(f"no FASTQ for {run}_{mate} in {bgdir}")


def trim_cmd(src, dst):
    return f"zcat -f {src} | awk 'NR%2==0{{print substr($0,1,{READ_LEN});next}}{{print}}' | gzip -1 > {dst}"


def cmd_bgprep(a):
    ref, work, bgdir = Path(a.ref), Path(a.work), Path(a.bg)
    for run in sorted(set(PILOT_BG.values())):
        o = work / "bg"; o.mkdir(parents=True, exist_ok=True)
        r1, r2 = o / f"{run}_1.150.fq.gz", o / f"{run}_2.150.fq.gz"
        if not r2.exists():
            sh(trim_cmd(bg_file(bgdir, run, 1), r1)); sh(trim_cmd(bg_file(bgdir, run, 2), r2))
        if not (o / f"{run}.bam").exists():
            bt2(r1, r2, ref / "idx/card90d", o / f"{run}.bam", a.threads)
        if not (o / f"{run}.res").exists():
            kma(r1, r2, ref / "idx/card90d", o / run, a.threads)


# ---------------------------------------------------------------- spike simulation + alignment
def cmd_simulate(a):
    ref, work = Path(a.ref), Path(a.work)
    for r in tsv_read(REPO / "data_manifest/pilot_design.tsv"):
        if a.only and r["sample"] not in a.only.split(","):
            continue
        sd = work / "samples" / r["sample"]; sd.mkdir(parents=True, exist_ok=True)
        if (sd / "sample.res").exists() and (sd / "spike.bam").exists():
            print("skip", r["sample"]); continue
        for f in ("spike_1.fq.gz", "spike_2.fq.gz"):
            (sd / f).unlink(missing_ok=True)
        for j in (1, 2, 3):
            b, cov = r[f"iso{j}"], r[f"cov{j}"]
            pre = sd / f"art_{b}"
            sh(f"art_illumina -ss HS25 -i {work}/genomes/{b}.fna -p -l {READ_LEN} -f {cov} -m 300 -s 30 "
               f"-rs {int(r['art_seed']) + j} -na -o {pre} > {pre}.log 2>&1")
            for m in (1, 2):   # rename reads: SPK_<biosample>_<n>/<mate>; append to spike files
                sh(f"awk -v p=SPK_{b} -v m={m} 'NR%4==1{{print \"@\"p\"_\"(NR+3)/4\"/\"m; next}}{{print}}' "
                   f"{pre}{m}.fq | gzip -1 >> {sd}/spike_{m}.fq.gz")
                os.remove(f"{pre}{m}.fq")
        bt2(sd / "spike_1.fq.gz", sd / "spike_2.fq.gz", ref / "idx/card90d", sd / "spike.bam", a.threads)
        bg = work / "bg"; run = r["background"]
        for m in (1, 2):   # full sample = trimmed background + spike (gzip members concatenate)
            sh(f"cat {bg}/{run}_{m}.150.fq.gz {sd}/spike_{m}.fq.gz > {sd}/full_{m}.fq.gz")
        kma(sd / "full_1.fq.gz", sd / "full_2.fq.gz", ref / "idx/card90d", sd / "sample", a.threads)
        n = int(subprocess.run(f"zcat {sd}/spike_1.fq.gz | wc -l", shell=True, capture_output=True,
                               text=True, executable="/bin/bash").stdout) // 4
        (sd / "spike_pairs.txt").write_text(str(n))
        for m in (1, 2):
            (sd / f"full_{m}.fq.gz").unlink()
        print(r["sample"], "spike pairs:", n, flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["split", "design", "genomes", "label", "bgprep", "simulate"])
    p.add_argument("--ref", default="/kaggle/input/card90-decoys-ref/ref")
    p.add_argument("--bg", default="/kaggle/input/ww-bg-5m/bg")
    p.add_argument("--work", default="/kaggle/working/pilot")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--all", action="store_true", help="genomes/label: all 174 isolates, not only pilot ones")
    p.add_argument("--only", default="", help="simulate: comma-separated sample ids")
    a = p.parse_args()
    globals()["cmd_" + a.cmd](a)


if __name__ == "__main__":
    main()
