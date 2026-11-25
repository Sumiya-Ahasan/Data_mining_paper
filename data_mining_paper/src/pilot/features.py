#!/usr/bin/env python3
"""Per-(sample, reference unit) evidence features + labels for the pilot.

Units = every reference in card90_decoys.fasta (1318 targets + decoys). A row is written for every unit
with >= 1 primary Bowtie2 alignment in (background BAM + spike BAM); all other units have p = 1 later.

Labels (column y): 1 = target in a spiked isolate (strict); 0 = target absent from spike and background,
or any decoy; NaN = ambiguous (background-present target, or ambiguous isolate hit) -> excluded from scoring.
Diagnostic columns (never model features): is_decoy, n_spike_reads, coverage, sample, role.

Usage: python features.py --ref REF --work WORK --out results/pilot/features.tsv.gz
"""
import argparse, csv, gzip, math, sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
FEATURES = ["n_reads", "rpkm", "breadth1", "breadth2", "evenness", "start_entropy", "gini_depth",
            "mean_mapq", "mean_identity", "indel_rate", "unique_frac", "softclip_frac", "gc_delta",
            "kma_score_per_kb", "kma_called"]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip().upper())
    if name: seqs[name] = "".join(buf)
    return seqs


def gc(s):
    s = s.upper(); n = sum(s.count(c) for c in "ACGT")
    return (s.count("G") + s.count("C")) / n if n else 0.0


def entropy_norm(starts, length, bin_bp=10):
    nb = max(1, math.ceil(length / bin_bp))
    if nb == 1 or not starts: return 0.0
    cnt = np.bincount(np.minimum(np.asarray(starts) // bin_bp, nb - 1), minlength=nb).astype(float)
    p = cnt[cnt > 0] / cnt.sum()
    return float(-(p * np.log(p)).sum() / math.log(nb))


def gini(x):
    x = np.sort(np.asarray(x, dtype=float))
    if x.sum() == 0: return 0.0
    n = len(x); i = np.arange(1, n + 1)
    return float((2 * (i * x).sum()) / (n * x.sum()) - (n + 1) / n)


def accumulate(bam_path, acc):
    import pysam  # imported here so gates.py can import FEATURES without pysam
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rd in bam:
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                continue
            u = acc[rd.reference_name]
            alen = rd.query_alignment_length or 1
            nm = rd.get_tag("NM") if rd.has_tag("NM") else 0
            indel = sum(l for op, l in rd.cigartuples if op in (1, 2))
            clip = sum(l for op, l in rd.cigartuples if op == 4)
            u["n"] += 1
            u["spike"] += rd.query_name.startswith("SPK_")
            u["mapq"] += rd.mapping_quality
            u["uniq"] += rd.mapping_quality >= 10
            u["ident"] += 1 - nm / alen
            u["indel"] += indel / alen
            u["clip"] += clip / (rd.query_length or 1)
            seq = rd.query_alignment_sequence or ""
            u["gc_n"] += sum(seq.count(c) for c in "GC"); u["acgt_n"] += sum(seq.count(c) for c in "ACGT")
            u["starts"].append(rd.reference_start)
            u["blocks"].extend(rd.get_blocks())


def kma_res(path):
    out = {}
    if not Path(path).exists(): return out
    for row in csv.reader(open(path, encoding="utf-8"), delimiter="\t"):
        if not row or row[0].startswith("#"): continue
        out[row[0].split()[0]] = float(row[1])   # Template -> Score
    return out


def tsv_read(path):
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref", default="/kaggle/input/card90-decoys-ref/ref")
    p.add_argument("--work", default="/kaggle/working/pilot")
    p.add_argument("--bg_pairs", type=int, default=5_000_000)
    p.add_argument("--out", default="/kaggle/working/pilot/features.tsv.gz")
    a = p.parse_args()
    ref, work = Path(a.ref), Path(a.work)

    seqs = read_fasta(ref / "card90_decoys.fasta")
    ref_gc = {k: gc(v) for k, v in seqs.items()}
    targets = set(read_fasta(ref / "card90.fasta"))

    lab = defaultdict(dict)   # biosample -> unit -> status
    for r in tsv_read(work / "labels/isolate_args.tsv"):
        lab[r["biosample"]][r["unit"]] = r["status"]

    mask_cache = {}
    out_rows = 0
    with gzip.open(a.out, "wt", newline="") as fh:
        cols = ["sample", "role", "background", "unit", "is_decoy", "y", "coverage", "n_spike_reads"] + FEATURES
        w = csv.DictWriter(fh, cols, delimiter="\t"); w.writeheader()
        for d in tsv_read(REPO / "data_manifest/pilot_design.tsv"):
            sd = work / "samples" / d["sample"]
            if not (sd / "spike.bam").exists():
                print("missing", d["sample"], file=sys.stderr); continue
            run = d["background"]
            if run not in mask_cache:   # background-present targets: >=1 read in unspiked background
                bacc = defaultdict(lambda: defaultdict(float, starts=[], blocks=[]))
                accumulate(work / "bg" / f"{run}.bam", bacc)
                mask_cache[run] = {u for u, v in bacc.items() if v["n"] >= 1 and u in targets}
            mask = mask_cache[run]
            acc = defaultdict(lambda: defaultdict(float, starts=[], blocks=[]))
            accumulate(work / "bg" / f"{run}.bam", acc)
            accumulate(sd / "spike.bam", acc)
            kscore = kma_res(sd / "sample.res")
            total_pairs = a.bg_pairs + int((sd / "spike_pairs.txt").read_text())
            # truth for this sample
            pos, amb = {}, set()
            for j in (1, 2, 3):
                b, cov = d[f"iso{j}"], float(d[f"cov{j}"])
                for u, st in lab[b].items():
                    if st == "strict": pos[u] = max(cov, pos.get(u, 0))
                    else: amb.add(u)
            for u, v in acc.items():
                L = len(seqs[u]); n = int(v["n"])
                depth = np.zeros(L, dtype=np.int32)
                for s, e in v["blocks"]:
                    depth[s:min(e, L)] += 1
                md = depth.mean()
                b1 = float((depth >= 1).mean()); b2 = float((depth >= 2).mean())
                expected = 1 - math.exp(-md) if md > 0 else 0
                is_decoy = u not in targets
                if is_decoy: y = 0
                elif u in mask: y = ""          # background-present: truth unknowable -> excluded
                elif u in pos: y = 1
                elif u in amb: y = ""
                else: y = 0
                w.writerow(dict(
                    sample=d["sample"], role=d["role"], background=run, unit=u, is_decoy=int(is_decoy), y=y,
                    coverage=pos.get(u, "") if not is_decoy else "", n_spike_reads=int(v["spike"]),
                    n_reads=n, rpkm=round(n / (L / 1000) / (2 * total_pairs / 1e6), 4),
                    breadth1=round(b1, 4), breadth2=round(b2, 4),
                    evenness=round(b1 / expected, 4) if expected > 0 else 0,
                    start_entropy=round(entropy_norm(v["starts"], L), 4), gini_depth=round(gini(depth), 4),
                    mean_mapq=round(v["mapq"] / n, 3), mean_identity=round(v["ident"] / n, 4),
                    indel_rate=round(v["indel"] / n, 5), unique_frac=round(v["uniq"] / n, 4),
                    softclip_frac=round(v["clip"] / n, 4),
                    gc_delta=round(abs((v["gc_n"] / v["acgt_n"] if v["acgt_n"] else 0) - ref_gc[u]), 4),
                    kma_score_per_kb=round(kscore.get(u, 0.0) / (L / 1000), 3), kma_called=int(u in kscore)))
                out_rows += 1
            # spiked positives with zero reads: kept so TPR has the right denominator (p = 1 later)
            for u, cov in pos.items():
                if u in acc or u in mask: continue
                row = {f: 0 for f in FEATURES}
                row.update(sample=d["sample"], role=d["role"], background=run, unit=u, is_decoy=0, y=1,
                           coverage=cov, n_spike_reads=0)
                w.writerow(row); out_rows += 1
            print(d["sample"], "units with reads:", len(acc), "positives present:",
                  sum(1 for u in pos if u in acc), "/", len(pos), flush=True)
    print("rows:", out_rows, "->", a.out)


if __name__ == "__main__":
    main()
