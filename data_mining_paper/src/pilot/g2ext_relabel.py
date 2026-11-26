"""Diagnostic (G2 extension, 2026-11-26): relabel target nulls of the G2 extension with fragment-aware rules (rep's own BLAST HSPs only):
 strict   : HSPs with pid>=90 summed per (query, subject) cover >=80% of the gene     -> y=1 (missed by single-HSP rule)
 fragment : any HSP pid>=90 and length>=150 bp (a read-length piece of the gene)       -> excluded
 homolog  : any HSP pid>=70 and length>=150 (local blocks can donate reads)             -> excluded only with --strict_homolog
args: features.tsv.gz null_totals.tsv blast_dir out_dir [--strict_homolog]"""
import sys, glob, os
from collections import defaultdict
import pandas as pd

feat, tots, bdir, out = sys.argv[1:5]; sh = "--strict_homolog" in sys.argv
os.makedirs(out, exist_ok=True)
cls = {}
for f in glob.glob(os.path.join(bdir, "*.blast.tsv")):
    b = os.path.basename(f).split(".")[0]
    cov90 = defaultdict(int); st = {}
    for line in open(f):
        q, s, pid, L, qlen, *_ = line.rstrip("\n").split("\t"); pid, L, qlen = float(pid), int(L), int(qlen)
        if pid >= 90: cov90[(q, s)] += L
        k = 1 if (pid >= 70 and L >= 150) else 0
        k = 2 if (pid >= 90 and L >= 150) else k
        if cov90[(q, s)] >= 0.8 * qlen: k = 3
        st[q] = max(st.get(q, 0), k)
    for q, k in st.items(): cls[(b, q)] = k
df = pd.read_csv(feat, sep="\t"); df["y"] = pd.to_numeric(df["y"], errors="coerce")
m = (df.y == 0) & (df.is_decoy == 0)
k = pd.Series([cls.get((s.replace("g2x_", ""), u), 0) for s, u in zip(df["sample"], df["unit"])], index=df.index)
newpos = m & (k == 3); excl = m & ((k == 2) | ((k == 1) & sh))
print("target-null rows:", int(m.sum()), "| -> strict:", int(newpos.sum()), "| -> excluded fragment/homolog:", int(excl.sum()))
df.loc[newpos, "y"] = 1; df.loc[excl, "y"] = float("nan")
# zero-read nulls that become excluded are not in features; count only hit rows (upper bound on change is small)
t = pd.read_csv(tots, sep="\t")
chg = (newpos | excl).groupby(df["sample"]).sum()
t["target_nulls"] = t["target_nulls"] - t["sample"].map(chg).fillna(0).astype(int)
df.to_csv(os.path.join(out, "features.tsv.gz"), sep="\t", index=False)
t.to_csv(os.path.join(out, "null_totals.tsv"), sep="\t", index=False)
