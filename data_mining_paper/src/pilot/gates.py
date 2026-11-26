#!/usr/bin/env python3
"""Pilot gates (plan §8) from features.tsv.gz.

G1 power:  n calibration nulls (ALL cal nulls incl. zero-read ones, from null_totals.tsv; zero-read nulls
           tie at the minimum score but still count in the conformal denominator) vs m/(q*k) - 1.
G2 decoys: validity needs decoy null scores to be at least as heavy-tailed as target null scores (then
           p-values are conservative). One-sided test: tail rate per null (all nulls, incl. zero-read) above
           decoy-score quantiles, target nulls vs decoys, Fisher exact (greater). fail = any p < 0.05;
           pass = >= 10 target-null tail events at the median threshold and no fail; else inconclusive.
           Score = logistic P(positive) cross-fitted by role (train on cal, score test, and vice versa).
           The two-sided C2ST AUC (hits only) is kept as a diagnostic.
G3 shift:  classifier two-sample test, cal vs test decoy hits (pass: AUC > 0.60) + indicative n_eff.
Also prints label counts and where null reads come from (background vs spike).

Usage: python gates.py --features features.tsv.gz --out results/pilot/gates.json
"""
import argparse, json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from scipy.stats import fisher_exact
from pathlib import Path

from features import FEATURES

M = 1318


def c2st_auc(X0, X1, seed=0):
    """5-fold CV AUC of distinguishing two samples; LightGBM if available, else logistic regression.
    Returns (auc_gbm_or_None, auc_logistic, n0, n1)."""
    n0, n1 = len(X0), len(X1)
    if min(n0, n1) < 20:
        return None, None, n0, n1
    X = np.vstack([X0, X1]); y = np.r_[np.zeros(n0), np.ones(n1)]
    cv = StratifiedKFold(5, shuffle=True, random_state=seed)
    lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced"))
    auc_lr = roc_auc_score(y, cross_val_predict(lr, X, y, cv=cv, method="predict_proba")[:, 1])
    auc_gb = None
    try:
        import lightgbm as lgb
        gb = lgb.LGBMClassifier(n_estimators=200, learning_rate=0.05, num_leaves=15, min_child_samples=20,
                                class_weight="balanced", verbose=-1, random_state=seed)
        auc_gb = roc_auc_score(y, cross_val_predict(gb, X, y, cv=cv, method="predict_proba")[:, 1])
    except ImportError:
        pass
    return auc_gb, auc_lr, n0, n1


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--features", default="/kaggle/working/pilot/features.tsv.gz")
    p.add_argument("--out", default="/kaggle/working/pilot/gates.json")
    p.add_argument("--totals", default=None, help="null_totals.tsv (default: next to --features)")
    a = p.parse_args()
    df = pd.read_csv(a.features, sep="\t")
    tot = pd.read_csv(a.totals or Path(a.features).parent / "null_totals.tsv", sep="\t")
    if "bg_present" not in df: df["bg_present"] = 0
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    hit = df[df.n_reads >= 1]
    res = {}

    # label overview
    res["counts"] = {
        role: {"targets_pos_total": int(((g.y == 1) & (g.is_decoy == 0)).sum()),
               "targets_pos_with_reads": int(((g.y == 1) & (g.n_reads >= 1)).sum()),
               "targets_pos_bg_present": int(((g.y == 1) & (g.bg_present == 1)).sum()),
               "target_nulls_with_reads": int(((g.y == 0) & (g.is_decoy == 0) & (g.n_reads >= 1)).sum()),
               "decoy_nulls_with_reads": int(((g.is_decoy == 1) & (g.n_reads >= 1)).sum()),
               "excluded": int(g.y.isna().sum()),
               "target_nulls_total": int(tot[tot.role == role].target_nulls.sum()),
               "decoy_nulls_total": int(tot[tot.role == role].decoy_nulls.sum()),
               "samples": int(g["sample"].nunique())}
        for role, g in df.groupby("role")}
    nulls = hit[hit.y == 0]
    res["null_read_origin_spike_frac"] = {
        "target_nulls": round(float(nulls[nulls.is_decoy == 0].n_spike_reads.sum() / max(1, nulls[nulls.is_decoy == 0].n_reads.sum())), 3),
        "decoy_nulls": round(float(nulls[nulls.is_decoy == 1].n_spike_reads.sum() / max(1, nulls[nulls.is_decoy == 1].n_reads.sum())), 3)}

    # G1
    cal_nulls = nulls[nulls.role == "cal"]
    ct = tot[tot.role == "cal"]
    n_cal = int(ct.decoy_nulls.sum() + ct.target_nulls.sum())
    k_per_sample = df[(df.role == "test") & (df.y == 1) & (df.n_reads >= 1)].groupby("sample").size()
    k = float(k_per_sample.median()) if len(k_per_sample) else 0.0
    g1 = {"n_cal_nulls_all": n_cal, "n_cal_decoy_nulls_all": int(ct.decoy_nulls.sum()),
          "n_cal_nulls_with_reads": len(cal_nulls), "median_k_test": k}
    for q in (0.05, 0.10):
        need = M / (q * k) - 1 if k > 0 else float("inf")
        g1[f"needed_q{q}"] = round(need, 1); g1[f"pass_q{q}"] = bool(n_cal >= need)
    res["G1_power"] = g1

    # G2 (one-sided tail dominance; C2ST kept as diagnostic)
    gb, lr, n0, n1 = c2st_auc(nulls[nulls.is_decoy == 0][FEATURES].values, nulls[nulls.is_decoy == 1][FEATURES].values)
    g2 = {"c2st_auc_lightgbm": gb, "c2st_auc_logistic": lr, "n_target_null_hits": n0, "n_decoy_null_hits": n1}
    sc = pd.Series(np.nan, index=hit.index)
    for tr, te in (("cal", "test"), ("test", "cal")):
        T = hit[(hit.role == tr) & ((hit.y == 1) | (hit.is_decoy == 1))]
        if T.y.nunique() < 2: continue
        mdl = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced"))
        mdl.fit(T[FEATURES].values, (T.y == 1).astype(int).values)
        E = hit[hit.role == te]; sc[E.index] = mdl.predict_proba(E[FEATURES].values)[:, 1]
    tn_s = sc[(hit.y == 0) & (hit.is_decoy == 0)].dropna(); dn_s = sc[hit.is_decoy == 1].dropna()
    TN, DN = int(tot.target_nulls.sum()), int(tot.decoy_nulls.sum())
    tails, worst_p = [], 1.0
    if len(dn_s):
        for qq in (0.5, 0.9, 0.99):
            t = float(dn_s.quantile(qq)); a_ = int((tn_s >= t).sum()); b_ = int((dn_s >= t).sum())
            pv = float(fisher_exact([[a_, TN - a_], [b_, DN - b_]], alternative="greater")[1])
            worst_p = min(worst_p, pv)
            tails.append({"decoy_quantile": qq, "threshold": round(t, 4), "target_null_events": a_,
                          "decoy_events": b_, "target_null_rate_per10k": round(1e4 * a_ / TN, 2),
                          "decoy_rate_per10k": round(1e4 * b_ / DN, 2), "fisher_p_greater": round(pv, 4)})
    g2.update(tails=tails, n_target_nulls_all=TN, n_decoy_nulls_all=DN)
    g2["verdict"] = ("fail" if worst_p < 0.05 else
                     "pass" if tails and tails[0]["target_null_events"] >= 10 else "inconclusive")
    res["G2_decoy_validity"] = g2

    # G3
    dn = nulls[nulls.is_decoy == 1]
    gb, lr, n0, n1 = c2st_auc(dn[dn.role == "cal"][FEATURES].values, dn[dn.role == "test"][FEATURES].values)
    auc = gb if gb is not None else lr
    res["G3_shift"] = {"auc_lightgbm": gb, "auc_logistic": lr, "n_cal": n0, "n_test": n1,
                       "pass": None if auc is None else bool(auc > 0.60)}
    if min(n0, n1) >= 20:   # indicative n_eff of likelihood-ratio weights on cal decoy hits (logistic, cross-fitted)
        X = np.vstack([dn[dn.role == "cal"][FEATURES].values, dn[dn.role == "test"][FEATURES].values])
        yy = np.r_[np.zeros(n0), np.ones(n1)]
        pr = cross_val_predict(make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)), X, yy,
                               cv=StratifiedKFold(5, shuffle=True, random_state=0), method="predict_proba")[:n0, 1]
        wt = pr / (1 - pr) * (n0 / n1)
        res["G3_shift"]["n_eff_cal_hits"] = round(float(wt.sum() ** 2 / (wt ** 2).sum()), 1)

    # per-feature means for the write-up (null types and positives)
    grp = hit.assign(kind=np.where(hit.y == 1, "positive", np.where(hit.is_decoy == 1, "decoy_null",
                     np.where(hit.y == 0, "target_null", "excluded"))))
    res["feature_means"] = grp.groupby("kind")[FEATURES].mean().round(4).to_dict(orient="index")

    print(json.dumps(res, indent=1, default=float))
    json.dump(res, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
