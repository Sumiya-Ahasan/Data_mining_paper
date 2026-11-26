#!/usr/bin/env python3
"""Simulated wastewater communities with fully known truth (decision_log 2026-11-26: decoys fail G2, so
calibration nulls are labelled absent targets from simulated communities).

Each real background run is profiled with Kraken2 (+Bracken); its species (those with a RefSeq genome) become
the template of a community. A simulated sample = background genomes at abundances drawn around the template
+ up to 3 spiked ARG isolates (ladder coverage). Truth for EVERY target comes from BLAST labels of every genome
in the sample (label v2, pilot.py label): nothing is masked, so absent targets with homolog reads stay nulls.

Speed-up: reads are simulated only from genome regions where an error-free 150-bp tile (step 50) aligns to the
reference with the same Bowtie2 settings (+-500 bp flanks); reads from elsewhere cannot align. ART runs at
max(cov, 10)x and pairs are then thinned to the target coverage (avoids ART's per-sequence rounding at low
coverage). Features use the nominal depth (--pairs) for RPKM.

Subcommands (run in order; each skips work already done):
  realbg   -> $WORK/bg/<run>.bam, trimmed reads           (template runs)
  profile  -> $WORK/profiles/<run>.species.tsv           (Kraken2 + Bracken, or Kraken report if no Bracken)
  select   -> $WORK/community_genomes.tsv                (RefSeq assemblies per species, strain pools A/B)
  genomes  -> $WORK/genomes/<GCF>.fna                    (background genomes; isolates are linked separately)
  design   -> $WORK/design_samples.tsv, design_members.tsv
  (labels: python src/pilot/pilot.py label --all --work $WORK --ref $REF)
  regions  -> $WORK/regions/<id>.fna + regions.tsv      (only genomes used in the design)
  simulate -> $WORK/samples/<sample>/spike.bam (all reads), sample.res, spike_pairs.txt (= nominal pairs)
  realism  -> $WORK/realism.tsv (real background vs its simulated communities)
"""
import argparse, csv, gzip, json, math, os, random, shutil, subprocess, sys, urllib.request
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src/pilot"))
import pilot  # noqa: E402  (sh, bt2, kma, tsv_read, tsv_write, bg_file, trim_cmd, ncbi_url, LADDER)

TEMPLATES = {"train": "ERR1713349", "cal": "ERR1713384", "test_same": "ERR1713384", "test_shift": "ERR1713386"}
N_SAMPLES = {"train": 20, "cal": 20, "test_same": 10, "test_shift": 10}
SPIKE_SPLIT = {"train": "train", "cal": "cal", "test_same": "cal", "test_shift": "test"}
POOL = {"train": "A", "cal": "A", "test_same": "A", "test_shift": "B"}   # B = strains unseen in train/cal
READ_LEN, TILE_STEP, FLANK, C_HI = 150, 50, 500, 10.0
LEVEL = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
ASM_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"


def read_fasta(path):
    seqs, name, buf = [], None, []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if name is not None: seqs.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip().upper())
    if name is not None: seqs.append((name, "".join(buf)))
    return seqs


# ---------------------------------------------------------------- real backgrounds
def runs_of(a):
    return sorted(set(a.runs.split(","))) if a.runs else sorted(set(TEMPLATES.values()))


def cmd_realbg(a):
    o = Path(a.work) / "bg"; o.mkdir(parents=True, exist_ok=True)
    for run in runs_of(a):
        r1, r2 = o / f"{run}_1.150.fq.gz", o / f"{run}_2.150.fq.gz"
        if not r2.exists():
            pilot.sh(pilot.trim_cmd(pilot.bg_file(a.bg, run, 1), r1)); pilot.sh(pilot.trim_cmd(pilot.bg_file(a.bg, run, 2), r2))
        if not (o / f"{run}.bam").exists():
            pilot.bt2(r1, r2, Path(a.ref) / "idx/card90d", o / f"{run}.bam", a.threads)


def parse_kraken_report(path):
    """-> (total reads, {species taxid: (name, clade reads)}). Kraken2 report: pct, clade, direct, rank, taxid, name."""
    total, sp = 0, {}
    for line in open(path, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6: continue
        rank, taxid, name = f[3].strip(), f[4].strip(), f[5].strip()
        if rank in ("U", "R"): total += int(f[1])
        if rank == "S": sp[taxid] = (name, int(f[1]))
    return total, sp


def cmd_profile(a):
    o = Path(a.work) / "profiles"; o.mkdir(parents=True, exist_ok=True)
    bg = Path(a.work) / "bg"
    for run in runs_of(a):
        out = o / f"{run}.species.tsv"
        if out.exists(): continue
        rep = o / f"{run}.k2report"
        if not rep.exists():
            pilot.sh(f"kraken2 --db {a.k2db} --threads {a.threads} --paired --gzip-compressed --report {rep} "
                     f"--output /dev/null {bg}/{run}_1.150.fq.gz {bg}/{run}_2.150.fq.gz")
        total, sp = parse_kraken_report(rep)
        reads, method = {t: n for t, (_, n) in sp.items()}, "kraken_clade"
        if shutil.which("bracken") and (Path(a.k2db) / f"database{READ_LEN}mers.kmer_distrib").exists():
            br = o / f"{run}.bracken"
            pilot.sh(f"bracken -d {a.k2db} -i {rep} -o {br} -r {READ_LEN} -l S > {br}.log 2>&1")
            reads = {r["taxonomy_id"]: int(r["new_est_reads"]) for r in pilot.tsv_read(br)}
            method = "bracken"
        rows = [dict(species_taxid=t, name=sp.get(t, ("?", 0))[0], reads=n, frac_total=round(n / total, 8))
                for t, n in sorted(reads.items(), key=lambda x: -x[1]) if n > 0]
        pilot.tsv_write(out, rows, ["species_taxid", "name", "reads", "frac_total"])
        json.dump(dict(run=run, total_pairs=total, method=method, classified_species_frac=round(sum(reads.values()) / total, 4)),
                  open(o / f"{run}.summary.json", "w"), indent=1)
        print(run, method, "pairs", total, "species", len(rows), flush=True)


# ---------------------------------------------------------------- genome selection
def raphenya_accessions():
    return {"_".join(r["assembly"].split("_")[:2]) for r in pilot.tsv_read(REPO / "data_manifest/isolate_split.tsv")}


def choose_assemblies(rows, k):
    """rows: assembly_summary dicts of one species -> best k (reference genome, then completeness, then newest)."""
    rows = sorted(rows, key=lambda r: r["seq_rel_date"], reverse=True)
    rows = sorted(rows, key=lambda r: (r["refseq_category"] != "reference genome", LEVEL.get(r["assembly_level"], 9)))
    return rows[:k]


def pools_for(n):
    """strain index -> pool. 1 strain: shared (AB); 2: A, B; 3: A, A, B."""
    return {1: ["AB"], 2: ["A", "B"], 3: ["A", "A", "B"]}[n]


def cmd_select(a):
    work = Path(a.work)
    want = {}   # taxid -> name, from every template profile (top max_species with frac >= min_frac)
    prof = {run: pilot.tsv_read(work / "profiles" / f"{run}.species.tsv") for run in runs_of(a)}
    asm_path = Path(a.asm_summary) if a.asm_summary else work / "assembly_summary_refseq.txt"
    if not asm_path.exists():
        print("download", ASM_URL, flush=True)
        urllib.request.urlretrieve(ASM_URL, asm_path)
    cols, by_sp = None, defaultdict(list)
    taxa = {r["species_taxid"] for rows in prof.values() for r in rows}
    excl = raphenya_accessions()
    with open(asm_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#assembly_accession"):
                cols = line[1:].rstrip("\n").split("\t"); continue
            if line.startswith("#") or cols is None: continue
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            if (r["species_taxid"] in taxa and r["group"] in ("bacteria", "archaea") and r["version_status"] == "latest"
                    and r["genome_rep"] == "Full" and r["excluded_from_refseq"] in ("na", "")
                    and r["assembly_accession"] not in excl and r["ftp_path"].startswith("http")):
                by_sp[r["species_taxid"]].append(r)
    out = []
    for run, rows in prof.items():
        kept = [r for r in rows if r["species_taxid"] in by_sp and float(r["frac_total"]) >= a.min_frac][:a.max_species]
        for r in kept: want[r["species_taxid"]] = r["name"]
        print(run, "species kept:", len(kept), "covering", round(sum(float(r["frac_total"]) for r in kept), 4),
              "of all read pairs", flush=True)
    for t, name in sorted(want.items()):
        ch = choose_assemblies(by_sp[t], a.strains)
        for r, pool in zip(ch, pools_for(len(ch))):
            out.append(dict(species_taxid=t, species=name, assembly=r["assembly_accession"], pool=pool,
                            genome_size=r["genome_size"], level=r["assembly_level"], ftp_path=r["ftp_path"]))
    pilot.tsv_write(work / "community_genomes.tsv", out, list(out[0]))
    print(len(want), "species,", len(out), "assemblies ->", work / "community_genomes.tsv")


def cmd_genomes(a):
    g = Path(a.work) / "genomes"; g.mkdir(parents=True, exist_ok=True)
    failed = []
    for r in pilot.tsv_read(Path(a.work) / "community_genomes.tsv"):
        f = g / f"{r['assembly']}.fna"
        if f.exists() and f.stat().st_size > 0: continue
        url = f"{r['ftp_path']}/{r['ftp_path'].rstrip('/').split('/')[-1]}_genomic.fna.gz"
        for t in range(3):
            try:
                with urllib.request.urlopen(url, timeout=180) as resp:
                    f.write_bytes(gzip.decompress(resp.read()))
                break
            except Exception as e:
                print("  retry", t + 1, r["assembly"], e, flush=True)
        else:
            failed.append(r["assembly"])
    print("background genomes in", g, "| FAILED:", failed or "none")


# ---------------------------------------------------------------- alignable regions
def merge_intervals(iv, flank, length):
    out = []
    for s, e in sorted((max(0, s - flank), min(length, e + flank)) for s, e in iv):
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out


def tiles(seqs, step=TILE_STEP, L=READ_LEN):
    for ci, (_, s) in enumerate(seqs):
        last = len(s) - L
        if last < 0: continue
        for p in list(range(0, last + 1, step)) + ([last] if last % step else []):
            yield f"{ci}|{p}", s[p:p + L]


def cmd_regions(a):
    work, ref = Path(a.work), Path(a.ref)
    rd = work / "regions"; rd.mkdir(parents=True, exist_ok=True)
    stats_path = rd / "regions.tsv"
    done = {r["genome"]: r for r in pilot.tsv_read(stats_path)} if stats_path.exists() else {}
    used = {r["genome"] for r in pilot.tsv_read(work / "design_members.tsv")} if (work / "design_members.tsv").exists() else None
    for g in sorted((work / "genomes").glob("*.fna")):
        gid = g.stem
        if used is not None and gid not in used: continue
        if gid in done and (rd / f"{gid}.fna").exists(): continue
        seqs = read_fasta(g)
        tf = rd / f"{gid}.tiles.fa"
        with open(tf, "w") as fh:
            for n, s in tiles(seqs): fh.write(f">{n}\n{s}\n")
        sam = subprocess.run(f"bowtie2 --very-sensitive-local -p {a.threads} -f -x {ref}/idx/card90d -U {tf} "
                             f"--no-unal --no-hd 2> /dev/null | cut -f1", shell=True, capture_output=True, text=True,
                             executable="/bin/bash", check=True).stdout.split()
        tf.unlink()
        hits = defaultdict(list)
        for n in set(sam):
            ci, p = n.split("|"); hits[int(ci)].append((int(p), int(p) + READ_LEN))
        nreg, bp = 0, 0
        with open(rd / f"{gid}.fna.tmp", "w") as fh:
            for ci, iv in sorted(hits.items()):
                name, s = seqs[ci]
                for st, en in merge_intervals(iv, FLANK, len(s)):
                    fh.write(f">{gid}_r{nreg} {name.split()[0]}:{st}-{en}\n{s[st:en]}\n"); nreg += 1; bp += en - st
        os.replace(rd / f"{gid}.fna.tmp", rd / f"{gid}.fna")
        done[gid] = dict(genome=gid, n_regions=nreg, region_bp=bp, genome_bp=sum(len(s) for _, s in seqs))
        pilot.tsv_write(stats_path, list(done.values()), ["genome", "n_regions", "region_bp", "genome_bp"])
        print(gid, "regions", nreg, "bp", bp, flush=True)


# ---------------------------------------------------------------- design
def cmd_design(a):
    work = Path(a.work); rng = random.Random(a.seed)
    gen = pilot.tsv_read(work / "community_genomes.tsv")
    strains = defaultdict(lambda: defaultdict(list))   # taxid -> pool -> [(assembly, size)]
    for r in gen:
        for p in ("A", "B"):
            if p in r["pool"]: strains[r["species_taxid"]][p].append((r["assembly"], int(r["genome_size"])))
    have = {p.stem for p in (work / "genomes").glob("*.fna")}
    iso = pilot.tsv_read(REPO / "data_manifest/isolate_split.tsv")
    samples, members = [], []
    for role, run in TEMPLATES.items():
        prof = [r for r in pilot.tsv_read(work / "profiles" / f"{run}.species.tsv") if r["species_taxid"] in strains]
        pool_iso = [r["biosample"] for r in iso if r["split"] == SPIKE_SPLIT[role]]
        for i in range(N_SAMPLES[role]):
            sid = f"sc_{role}_{i:02d}"
            samples.append(dict(sample=sid, role=role, template=run, background="none", n_pairs=a.pairs,
                                art_seed=a.seed * 1000 + len(samples)))
            sp = [(r, strains[r["species_taxid"]][POOL[role]]) for r in prof]
            sp = [(r, s) for r, s in sp if any(x[0] in have for x in s)]
            w = [float(r["frac_total"]) * math.exp(rng.gauss(0, a.sigma)) for r, _ in sp]
            scale = sum(float(r["frac_total"]) for r, _ in sp) / sum(w)
            for (r, s), wi in zip(sp, w):
                asm, size = rng.choice([x for x in s if x[0] in have])
                cov = a.pairs * wi * scale * 2 * READ_LEN / size
                members.append(dict(sample=sid, genome=asm, kind="bg", coverage=round(cov, 5)))
            for j, b in enumerate(rng.sample(pool_iso, 3)):
                members.append(dict(sample=sid, genome=b, kind="spike",
                                    coverage=pilot.LADDER[(3 * i + j) % len(pilot.LADDER)]))
    pilot.tsv_write(work / "design_samples.tsv", samples, list(samples[0]))
    pilot.tsv_write(work / "design_members.tsv", members, list(members[0]))
    print(len(samples), "samples,", len(members), "members ->", work)


# ---------------------------------------------------------------- simulation
def thin_rename(f1, f2, o1, o2, keep, prefix, rng):
    """Keep each read pair with probability keep; rename @<prefix>_<n>/<mate>. Returns pairs kept."""
    n = 0
    with open(f1) as a1, open(f2) as a2:
        while True:
            r1 = [a1.readline() for _ in range(4)]; r2 = [a2.readline() for _ in range(4)]
            if not r1[0]: break
            if keep < 1 and rng.random() >= keep: continue
            n += 1
            o1.write(f"@{prefix}_{n}/1\n{r1[1]}+\n{r1[3]}"); o2.write(f"@{prefix}_{n}/2\n{r2[1]}+\n{r2[3]}")
    return n


def cmd_simulate(a):
    work, ref = Path(a.work), Path(a.ref)
    mem = defaultdict(list)
    for r in pilot.tsv_read(work / "design_members.tsv"): mem[r["sample"]].append(r)
    for s in pilot.tsv_read(work / "design_samples.tsv"):
        if a.only and s["sample"] not in a.only.split(","): continue
        sd = work / "samples" / s["sample"]; sd.mkdir(parents=True, exist_ok=True)
        if (sd / "sample.res").exists() and (sd / "spike.bam").exists():
            print("skip", s["sample"]); continue
        rng = random.Random(int(s["art_seed"]))
        n_sim = 0
        with gzip.open(sd / "spike_1.fq.gz", "wt", compresslevel=1) as o1, gzip.open(sd / "spike_2.fq.gz", "wt", compresslevel=1) as o2:
            for j, m in enumerate(mem[s["sample"]]):
                reg = work / "regions" / f"{m['genome']}.fna"
                cov = float(m["coverage"])
                if cov <= 0 or not reg.exists() or reg.stat().st_size == 0: continue
                pre = sd / "art_"
                fc = max(cov, C_HI)
                pilot.sh(f"art_illumina -ss HS25 -i {reg} -p -l {READ_LEN} -f {fc} -m 300 -s 30 "
                         f"-rs {int(s['art_seed']) * 100 + j} -na -o {pre} > {pre}.log 2>&1")
                pfx = ("SPK_" if m["kind"] == "spike" else "BG_") + m["genome"]
                n_sim += thin_rename(f"{pre}1.fq", f"{pre}2.fq", o1, o2, cov / fc, pfx, rng)
                os.remove(f"{pre}1.fq"); os.remove(f"{pre}2.fq")
        pilot.bt2(sd / "spike_1.fq.gz", sd / "spike_2.fq.gz", ref / "idx/card90d", sd / "spike.bam", a.threads)
        pilot.kma(sd / "spike_1.fq.gz", sd / "spike_2.fq.gz", ref / "idx/card90d", sd / "sample", a.threads)
        (sd / "spike_pairs.txt").write_text(str(s["n_pairs"]))      # nominal depth for RPKM
        (sd / "simulated_pairs.txt").write_text(str(n_sim))
        for f in ("spike_1.fq.gz", "spike_2.fq.gz"): (sd / f).unlink(missing_ok=True)
        print(s["sample"], "simulated pairs (alignable regions only):", n_sim, flush=True)


# ---------------------------------------------------------------- realism check
def cmd_realism(a):
    """Per sample: target units with >=1 / >=10 reads, reads on targets per million pairs, median read identity of
    target hits, share of hit units in the efflux class; real background BAM vs its simulated communities."""
    import numpy as np
    sys.path.insert(0, str(REPO / "src/pilot")); from features import accumulate, read_fasta as rf
    work, ref = Path(a.work), Path(a.ref)
    targets = set(rf(ref / "card90.fasta"))
    aro = {r["ARO Accession"]: r["Resistance Mechanism"] for r in pilot.tsv_read(ref / "card/aro_index.tsv")}
    efflux = {t for t in targets if "efflux" in aro.get(t.split("|")[4], "")}

    def summ(bam, pairs):
        acc = defaultdict(lambda: defaultdict(float, starts=[], blocks=[]))
        accumulate(bam, acc)
        hit = {u: v for u, v in acc.items() if u in targets}
        idn = [v["ident"] / v["n"] for v in hit.values()]
        return dict(units_ge1=len(hit), units_ge10=sum(v["n"] >= 10 for v in hit.values()),
                    target_reads_per_M=round(sum(v["n"] for v in hit.values()) / (pairs / 1e6), 2),
                    median_identity=round(float(np.median(idn)), 4) if idn else "",
                    efflux_share=round(sum(u in efflux for u in hit) / max(1, len(hit)), 3))
    rows = []
    prof = {p.stem.split(".")[0]: json.load(open(p)) for p in (work / "profiles").glob("*.summary.json")}
    for run in runs_of(a):
        rows.append(dict(kind="real", sample=run, template=run, **summ(work / "bg" / f"{run}.bam", prof[run]["total_pairs"])))
    for s in pilot.tsv_read(work / "design_samples.tsv"):
        b = work / "samples" / s["sample"] / "spike.bam"
        if b.exists():
            rows.append(dict(kind="sim_" + s["role"], sample=s["sample"], template=s["template"],
                             **summ(b, int(s["n_pairs"]))))
    pilot.tsv_write(work / "realism.tsv", rows, list(rows[0]))
    print("realism ->", work / "realism.tsv")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["realbg", "profile", "select", "genomes", "regions", "design", "simulate", "realism"])
    p.add_argument("--work", default="/kaggle/working/simcomm")
    p.add_argument("--ref", default="/kaggle/working/ref2")
    p.add_argument("--bg", default="/kaggle/input/ww-bg-5m/bg")
    p.add_argument("--k2db", default="/tmp/k2")
    p.add_argument("--runs", default="", help="comma-separated runs (default: the template runs)")
    p.add_argument("--asm-summary", default="", help="local assembly_summary_refseq.txt (else downloaded)")
    p.add_argument("--max-species", type=int, default=40)
    p.add_argument("--min-frac", type=float, default=0.0005, help="min share of ALL read pairs per species")
    p.add_argument("--strains", type=int, default=3)
    p.add_argument("--pairs", type=int, default=5_000_000, help="nominal read pairs per sample")
    p.add_argument("--sigma", type=float, default=0.5, help="log-normal abundance noise around the template")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--seed", type=int, default=2027)
    p.add_argument("--only", default="")
    a = p.parse_args()
    globals()["cmd_" + a.cmd](a)


if __name__ == "__main__":
    main()
