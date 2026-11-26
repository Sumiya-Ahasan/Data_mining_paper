"""Diagnostic: G2 tail test restricted to target subsets by CARD Resistance Mechanism.
Usage: python -P src/pilot/g2_by_mechanism.py {dec|tn} <extracted g2v2 results dir with ref2/ and g2v2/>
(dec = score trained on positives vs decoys as in gates.py; tn = also on target nulls)."""
import sys, csv, json, numpy as np, pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent)); from features import FEATURES, read_fasta
C = sys.argv[2]; train_tn = sys.argv[1] == "tn"
aro = pd.read_csv(f"{C}/ref2/card/aro_index.tsv", sep="\t")
mech = aro.groupby("ARO Accession")["Resistance Mechanism"].first().to_dict()
targets = list(read_fasta(f"{C}/ref2/card90.fasta"))
tm = {t: mech.get(t.split("|")[4], "NA") for t in targets}
print("unmapped targets:", sum(v == "NA" for v in tm.values()), "of", len(tm))
print(pd.Series(tm).str.contains("efflux").value_counts().to_dict())
f = pd.read_csv(f"{C}/g2v2/features.tsv.gz", sep="\t"); f["y"] = pd.to_numeric(f.y, errors="coerce")
f["parent"] = np.where(f.is_decoy == 1, f.unit.str.split("|", n=2).str[2], f.unit)
hit = f[f.n_reads >= 1].copy()
sc = pd.Series(np.nan, index=hit.index)
for te in hit.role.unique():
    T = hit[(hit.role != te) & (hit.y.notna() if train_tn else ((hit.y == 1) | (hit.is_decoy == 1)))]
    m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced")).fit(T[FEATURES].values, (T.y == 1).astype(int))
    E = hit[hit.role == te]; sc[E.index] = m.predict_proba(E[FEATURES].values)[:, 1]
hit["s"] = sc
# exact null counts per sample: labels + design
lab = defaultdict(dict)
for r in csv.DictReader(open(f"{C}/g2v2/labels/isolate_args.tsv"), delimiter="\t"): lab[r["biosample"]][r["unit"]] = r["status"]
des = list(csv.DictReader(open(Path(__file__).resolve().parents[2] / "data_manifest/g2ext_design.tsv"), delimiter="\t"))
decoys = [u for u in read_fasta(f"{C}/ref2/card90_decoys.fasta") if u.startswith("DECOY|")]
dpar = pd.Series([u.split("|", 2)[2] for u in decoys])
def run(name, keep):
    K = {t for t in targets if keep(tm[t])}
    TN = sum(len(K - set(lab[d["iso1"]])) for d in des)          # strict+fragment units removed
    D = len(des) * int(dpar.isin(K).sum())
    h = hit[hit.parent.isin(K)]
    tn = h[(h.y == 0) & (h.is_decoy == 0)].s; dn = h[h.is_decoy == 1].s
    out = []
    for q in (0.5, 0.9, 0.99, 0.999):
        t = dn.quantile(q); a = int((tn >= t).sum()); b = int((dn >= t).sum())
        out.append((q, a, b, round(b / D * TN, 2), fisher_exact([[a, TN - a], [b, D - b]], alternative="greater")[1]))
    print(f"{name}: targets {len(K)} TN {TN} DN {D} hits tn {len(tn)} dn {len(dn)}")
    for q, a, b, e, p in out: print(f"   q{q}: tn {a} dn {b} exp {e} p {p:.2g}")
run("all", lambda m: True)
run("no efflux", lambda m: "efflux" not in m)
run("inactivation only", lambda m: m == "antibiotic inactivation")
