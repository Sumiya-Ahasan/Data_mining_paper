#!/usr/bin/env python3
"""Conformal selection with labelled calibration nulls, evaluated on simulated communities (truth known).

Score: LightGBM (else logistic) P(positive) fit on role=train rows with reads, labels y (target nulls + positives).
Calibration nulls (role=cal), three arms:
  target_pooled     all cal target nulls (zero-read ones count, score -inf)
  target_mondrian   same, but per stratum (CARD mechanism: efflux vs other); p uses the test unit's stratum
  decoy_pooled      cal decoys (the method that failed G2; kept for comparison)
p = (1 + #{cal null scores >= s}) / (n_cal + 1); zero-read test units p = 1.
Per test sample: BH over every non-excluded target (m = |targets| - excluded) -> FDP, TPR (overall and by coverage).
Groups: test_same (calibration template, held-out communities, no shift) and test_shift (other template, unseen strains).
Gate V1 (validity, test_same, target arms): pass = mean FDP <= q; fail = bootstrap 95% lower bound > q; else inconclusive.

Usage: python conformal_eval.py --features F --ref REF --out eval.json   (truth.tsv.gz, null_totals.tsv next to F)
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "pilot"))
from features import FEATURES, read_fasta  # noqa: E402

QS = (0.01, 0.05, 0.10, 0.20)
COV_BINS = [0, 0.25, 1, 5, np.inf]


def fit_score(train):
    X, y = train[FEATURES].values, (train.y == 1).astype(int).values
    try:
        import lightgbm as lgb
        m = lgb.LGBMClassifier(n_estimators=300, learning_rate=0.05, num_leaves=15, min_child_samples=20,
                               class_weight="balanced", verbose=-1, random_state=0).fit(X, y)
        return m, "lightgbm"
    except ImportError:
        from sklearn.linear_model import LogisticRegression
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
        return make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000, class_weight="balanced")).fit(X, y), "logistic"


def pvals(test_scores, cal_scores, n_cal):
    """(1 + #{cal >= s}) / (n_cal + 1); cal_scores = scores of cal nulls WITH reads (others are -inf)."""
    c = np.sort(np.asarray(cal_scores))
    ge = len(c) - np.searchsorted(c, test_scores, side="left")
    return (1 + ge) / (n_cal + 1)


def bh(p, m, q):
    """indices selected by BH at level q among len(p) tested p-values plus (m - len(p)) p-values equal to 1."""
    if not len(p): return np.array([], dtype=int)
    o = np.argsort(p); ps = p[o]
    ok = np.nonzero(ps <= q * np.arange(1, len(ps) + 1) / m)[0]
    return o[:ok[-1] + 1] if len(ok) else np.array([], dtype=int)


def boot_ci(x, rng, B=2000):
    x = np.asarray(x, float)
    if len(x) < 2: return [float(x.mean()) if len(x) else None] * 2
    b = rng.choice(x, (B, len(x))).mean(1)
    return [round(float(np.quantile(b, .025)), 4), round(float(np.quantile(b, .975)), 4)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = Path(a.features).parent
    f = pd.read_csv(a.features, sep="\t"); f["y"] = pd.to_numeric(f.y, errors="coerce")
    truth = pd.read_csv(d / "truth.tsv.gz", sep="\t")
    tot = pd.read_csv(d / "null_totals.tsv", sep="\t")
    targets = list(read_fasta(Path(a.ref) / "card90.fasta"))
    aro = pd.read_csv(Path(a.ref) / "card/aro_index.tsv", sep="\t").groupby("ARO Accession")["Resistance Mechanism"].first()
    stratum = {t: "efflux" if "efflux" in str(aro.get(t.split("|")[4], "")) else "other" for t in targets}
    f["stratum"] = f.unit.map(stratum).fillna("decoy")
    hit = f[f.n_reads >= 1].copy()

    model, kind = fit_score(hit[(hit.role == "train") & (hit.is_decoy == 0) & hit.y.notna()])
    ev = hit[hit.role != "train"].copy()
    ev["s"] = model.predict_proba(ev[FEATURES].values)[:, 1]

    # calibration sets (counts include zero-read nulls)
    cal_samples = sorted(f[f.role == "cal"]["sample"].unique())
    ct = truth[truth["sample"].isin(cal_samples)]
    n_null_cal = {s_: sum(1 for t in targets if stratum[t] == s_) * len(cal_samples)
                  - int((ct.unit.map(stratum) == s_).sum()) for s_ in ("efflux", "other")}
    c = ev[ev.role == "cal"]
    cn = c[(c.is_decoy == 0) & (c.y == 0)]
    cal = {"target_pooled": {"all": (cn.s.values, sum(n_null_cal.values()))},
           "target_mondrian": {s_: (cn[cn.stratum == s_].s.values, n_null_cal[s_]) for s_ in n_null_cal},
           "decoy_pooled": {"all": (c[c.is_decoy == 1].s.values, int(tot[tot.role == "cal"].decoy_nulls.sum()))}}

    rng = np.random.default_rng(0)
    res = {"score_model": kind, "n_train_rows": int(((hit.role == "train") & (hit.is_decoy == 0) & hit.y.notna()).sum()),
           "calibration": {k: {s_: {"n_null": n, "n_with_reads": len(v)} for s_, (v, n) in arms.items()}
                           for k, arms in cal.items()}, "groups": {}}
    for grp in ("test_same", "test_shift"):
        g_samples = sorted(f[f.role == grp]["sample"].unique())
        if not g_samples: continue
        out = {}
        for arm, arms in cal.items():
            per_q = {}
            for q in QS:
                fdp, tpr, R, tp_bin, pos_bin = [], [], [], np.zeros(len(COV_BINS) - 1), np.zeros(len(COV_BINS) - 1)
                for smp in g_samples:
                    tr = truth[truth["sample"] == smp]
                    excl = set(tr[tr.status == "excluded"].unit)
                    pos = tr[tr.status == "positive"].set_index("unit").coverage.astype(float)
                    m = len(targets) - len(excl)
                    u = ev[(ev["sample"] == smp) & (ev.is_decoy == 0) & ~ev.unit.isin(excl)]
                    p = np.ones(len(u))
                    for s_, (cs, n) in arms.items():
                        mask = np.ones(len(u), bool) if s_ == "all" else (u.stratum == s_).values
                        if mask.any(): p[mask] = pvals(u.s.values[mask], cs, n)
                    sel = u.iloc[bh(p, m, q)]
                    nf = int((sel.y == 0).sum()); ntp = int((sel.y == 1).sum())
                    fdp.append(nf / max(1, len(sel))); R.append(len(sel)); tpr.append(ntp / max(1, len(pos)))
                    b_all = np.digitize(pos.values, COV_BINS[1:-1])
                    b_tp = np.digitize(pos[pos.index.isin(set(sel[sel.y == 1].unit))].values, COV_BINS[1:-1])
                    pos_bin += np.bincount(b_all, minlength=len(pos_bin)); tp_bin += np.bincount(b_tp, minlength=len(tp_bin))
                mf = float(np.mean(fdp)); ci = boot_ci(fdp, rng)
                per_q[str(q)] = {"mean_fdp": round(mf, 4), "fdp_ci95": ci, "mean_selected": round(float(np.mean(R)), 1),
                                 "mean_tpr": round(float(np.mean(tpr)), 4),
                                 "tpr_by_coverage": {f"{COV_BINS[i]}-{COV_BINS[i+1]}": round(float(tp_bin[i] / pos_bin[i]), 4)
                                                     if pos_bin[i] else None for i in range(len(pos_bin))},
                                 "verdict": "pass" if mf <= q else ("fail" if ci[0] is not None and ci[0] > q else "inconclusive")}
            out[arm] = per_q
        res["groups"][grp] = {"n_samples": len(g_samples), "arms": out}
    json.dump(res, open(a.out, "w"), indent=1)
    for grp, g in res["groups"].items():
        print(f"== {grp} ({g['n_samples']} samples)")
        for arm, per_q in g["arms"].items():
            print("  ", arm, " | ".join(f"q={q}: FDR {v['mean_fdp']} {v['fdp_ci95']} TPR {v['mean_tpr']} {v['verdict']}"
                                         for q, v in per_q.items()))
    print("->", a.out)


if __name__ == "__main__":
    main()
