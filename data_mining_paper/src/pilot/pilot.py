#!/usr/bin/env python3
"""Pilot driver: isolate split, mix design, genome download, ARG labels, background prep,
spike-in simulation and alignment. Feature extraction: features.py. Gates: gates.py.

Subcommands (run in order; each skips work already done):
  split     -> data_manifest/isolate_split.tsv           (local or Kaggle)
  design    -> data_manifest/pilot_design.tsv            (local or Kaggle)
  genomes   -> $WORK/genomes/<biosample>.fna             (Kaggle, internet ON)
  label     -> $WORK/labels/isolate_args.tsv             (Kaggle; BLAST of CARD homologs vs genomes; label v2)
  bgprep    -> $WORK/bg/<run>_{1,2}.150.fq.gz, bg.bam, kma .res, mask   (Kaggle)
  simulate  -> $WORK/samples/<sample>/spike.bam, sample.res              (Kaggle)
  g2design  -> data_manifest/g2ext_design.tsv: every isolate alone, no background (G2 extension)

--design picks the design file (default pilot_design.tsv); genomes/label/simulate use its isolates.
A design row with background "none" = isolate reads only (truth fully known; used for G2).

Paths come from flags (defaults suit Kaggle; see notebooks/pilot_kaggle.md).
"""
import argparse, csv, gzip, os, random, re, shutil, subprocess, sys, urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from window_id import max_window_identity, classify  # noqa: E402
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


def cmd_g2design(a):
    """G2 extension: each of the 174 isolates simulated alone at --cov, no background."""
    iso = tsv_read(REPO / "data_manifest/isolate_split.tsv")
    rows = [dict(sample=f"g2x_{r['biosample']}", background="none", role=r["split"], iso1=r["biosample"],
                 cov1=a.cov, iso2="", cov2="", iso3="", cov3="", art_seed=a.seed * 1000 + 500 + i)
            for i, r in enumerate(iso)]
    cols = ["sample", "background", "role", "iso1", "cov1", "iso2", "cov2", "iso3", "cov3", "art_seed"]
    tsv_write(REPO / "data_manifest/g2ext_design.tsv", rows, cols)
    print(len(rows), "samples -> data_manifest/g2ext_design.tsv")


def design_isolates(design):
    d = tsv_read(Path(design))
    return sorted({r[f"iso{j}"] for r in d for j in (1, 2, 3) if r.get(f"iso{j}")})


# ---------------------------------------------------------------- genomes
def ncbi_url(asm):
    acc = asm.split("_")[0] + "_" + asm.split("_")[1]          # GCF_016903215.1
    digits = acc.split("_")[1].split(".")[0]                    # 016903215
    p = "/".join(digits[i:i + 3] for i in range(0, 9, 3))
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{acc[:3]}/{p}/{asm}/{asm}_genomic.fna.gz"


def cmd_genomes(a):
    iso = {r["biosample"]: r for r in tsv_read(REPO / "data_manifest/isolate_split.tsv")}
    want = sorted(iso) if a.all else design_isolates(a.design)
    out = Path(a.work) / "genomes"; out.mkdir(parents=True, exist_ok=True)
    failed = []
    for b in want:
        f = out / f"{b}.fna"
        if f.exists() and f.stat().st_size > 0:
            continue
        url = ncbi_url(iso[b]["assembly"])
        print("get", b, url, flush=True)
        for t in range(3):
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    data = gzip.decompress(r.read())
                f.write_bytes(data); break
            except Exception as e:
                print("  retry", t + 1, b, e, flush=True)
        else:
            failed.append(b)
    print(len(want) - len(failed), "genomes in", out, "| FAILED:", failed or "none")


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
    """Label v2 (2026-11-26, src/window_id.py): per CARD unit (card_90 representative) and genome, over ALL
    cluster members: 'strict' (y=1) = HSPs at >=90% id on one contig cover >=80% of the member;
    'fragment' (excluded) = not strict but some HSP has a 150-column window at >=90% identity.
    Everything else is a null, distant homologs included (v1's 70%/50% 'ambiguous' class hid those)."""
    ref = Path(a.ref); work = Path(a.work)
    cmap = cluster_map(ref / "card90.fasta.clstr")
    homolog = ref / "card/nucleotide_fasta_protein_homolog_model.fasta"
    out = work / "labels"; out.mkdir(parents=True, exist_ok=True)
    rows = []
    for b in design_isolates(a.design) if not a.all else sorted(p.stem for p in (work / "genomes").glob("*.fna")):
        g = work / "genomes" / f"{b}.fna"
        if not g.exists():
            print("no genome, skip", b); continue
        summ = out / f"{b}.hsp2.tsv"    # one line per HSP, alignment strings reduced to the max window identity
        if not summ.exists():
            raw = out / f"{b}.blast2.raw"
            sh(f"makeblastdb -in {g} -dbtype nucl -out {out}/{b}.db > /dev/null")
            sh(f"blastn -task dc-megablast -query {homolog} -db {out}/{b}.db -evalue 1e-10 "
               f"-max_target_seqs 5 -num_threads {a.threads} "
               f"-outfmt '6 qseqid sseqid pident length qlen qstart qend qseq sseq' > {raw}")
            with open(raw) as fi, open(str(summ) + ".tmp", "w") as fo:
                for line in fi:
                    q, s_, pid, L, qlen, qs, qe, qseq, sseq = line.rstrip("\n").split("\t")
                    mw = max_window_identity(qseq, sseq)
                    fo.write("\t".join([q, s_, pid, L, qlen, qs, qe, "" if mw is None else f"{mw:.4f}"]) + "\n")
            os.replace(str(summ) + ".tmp", summ); raw.unlink()
            for f in out.glob(f"{b}.db.*"): f.unlink()
        hs = {}
        for line in open(summ):
            q, s_, pid, L, qlen, qs, qe, mw = line.rstrip("\n").split("\t")
            hs.setdefault(q, []).append(dict(sseqid=s_, pident=float(pid), qlen=int(qlen), qstart=int(qs),
                                             qend=int(qe), maxwin=float(mw) if mw else None))
        best = {}  # rep -> (rank, cov90, maxwin, member)
        for q, h in hs.items():
            rep = cmap.get(q)
            if rep is None:
                continue
            st, cov, mw = classify(h)
            if st is None: continue
            rank = 2 if st == "strict" else 1
            if rep not in best or (rank, cov, mw) > best[rep][:3]:
                best[rep] = (rank, cov, mw, q)
        for rep, (rank, cov, mw, q) in best.items():
            rows.append(dict(biosample=b, unit=rep, status="strict" if rank == 2 else "fragment",
                             qcov90=cov, max_window_identity=mw, member=q))
        print(b, "strict:", sum(v[0] == 2 for v in best.values()), "fragment:", sum(v[0] == 1 for v in best.values()),
              flush=True)
    tsv_write(out / "isolate_args.tsv", rows, ["biosample", "unit", "status", "qcov90", "max_window_identity", "member"])


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
    for r in tsv_read(Path(a.design)):
        if a.only and r["sample"] not in a.only.split(","):
            continue
        isos = [(r[f"iso{j}"], r[f"cov{j}"]) for j in (1, 2, 3) if r.get(f"iso{j}")]
        if any(not (work / "genomes" / f"{b}.fna").exists() for b, _ in isos):
            print("missing genome, skip", r["sample"]); continue
        sd = work / "samples" / r["sample"]; sd.mkdir(parents=True, exist_ok=True)
        if (sd / "sample.res").exists() and (sd / "spike.bam").exists():
            print("skip", r["sample"]); continue
        for f in ("spike_1.fq.gz", "spike_2.fq.gz"):
            (sd / f).unlink(missing_ok=True)
        for j, (b, cov) in enumerate(isos, 1):
            pre = sd / f"art_{b}"
            sh(f"art_illumina -ss HS25 -i {work}/genomes/{b}.fna -p -l {READ_LEN} -f {cov} -m 300 -s 30 "
               f"-rs {int(r['art_seed']) + j} -na -o {pre} > {pre}.log 2>&1")
            for m in (1, 2):   # rename reads: SPK_<biosample>_<n>/<mate>; append to spike files
                sh(f"awk -v p=SPK_{b} -v m={m} 'NR%4==1{{print \"@\"p\"_\"(NR+3)/4\"/\"m; next}}{{print}}' "
                   f"{pre}{m}.fq | gzip -1 >> {sd}/spike_{m}.fq.gz")
                os.remove(f"{pre}{m}.fq")
        bt2(sd / "spike_1.fq.gz", sd / "spike_2.fq.gz", ref / "idx/card90d", sd / "spike.bam", a.threads)
        bg = work / "bg"; run = r["background"]
        if run == "none":   # isolate reads only
            kma(sd / "spike_1.fq.gz", sd / "spike_2.fq.gz", ref / "idx/card90d", sd / "sample", a.threads)
        else:
            for m in (1, 2):   # full sample = trimmed background + spike (gzip members concatenate)
                sh(f"cat {bg}/{run}_{m}.150.fq.gz {sd}/spike_{m}.fq.gz > {sd}/full_{m}.fq.gz")
            kma(sd / "full_1.fq.gz", sd / "full_2.fq.gz", ref / "idx/card90d", sd / "sample", a.threads)
        n = int(subprocess.run(f"zcat {sd}/spike_1.fq.gz | wc -l", shell=True, capture_output=True,
                               text=True, executable="/bin/bash").stdout) // 4
        (sd / "spike_pairs.txt").write_text(str(n))
        for f in ("full_1.fq.gz", "full_2.fq.gz", "spike_1.fq.gz", "spike_2.fq.gz"):   # BAM + .res are kept
            (sd / f).unlink(missing_ok=True)
        print(r["sample"], "spike pairs:", n, flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["split", "design", "genomes", "label", "bgprep", "simulate", "g2design"])
    p.add_argument("--design", default=str(REPO / "data_manifest/pilot_design.tsv"))
    p.add_argument("--cov", type=float, default=10, help="g2design: coverage per isolate")
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
