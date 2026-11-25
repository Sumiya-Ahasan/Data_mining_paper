#!/usr/bin/env python3
"""Pilot gates (plan §8) from features.tsv.gz.

G1 power:  n calibration nulls (y=0, >=1 read, cal samples) vs required m/(q*k), q in {0.05, 0.10}
G2 decoys: classifier two-sample test, decoy nulls vs target nulls (pass: AUC <= 0.60)
G3 shift:  classifier two-sample test, cal nulls vs test nulls, decoy hits only (pass: AUC > 0.60)
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
    a = p.parse_args()
    df = pd.read_csv(a.features, sep="\t")
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    hit = df[df.n_reads >= 1]
    res = {}

    # label overview
    res["counts"] = {
        role: {"targets_pos_total": int(((g.y == 1) & (g.is_decoy == 0)).sum()),
               "targets_pos_with_reads": int(((g.y == 1) & (g.n_reads >= 1)).sum()),
               "target_nulls_with_reads": int(((g.y == 0) & (g.is_decoy == 0) & (g.n_reads >= 1)).sum()),
               "decoy_nulls_with_reads": int(((g.is_decoy == 1) & (g.n_reads >= 1)).sum()),
               "ambiguous": int(g.y.isna().sum()),
               "samples": int(g["sample"].nunique())}
        for role, g in df.groupby("role")}
    nulls = hit[hit.y == 0]
    res["null_read_origin_spike_frac"] = {
        "target_nulls": round(float(nulls[nulls.is_decoy == 0].n_spike_reads.sum() / max(1, nulls[nulls.is_decoy == 0].n_reads.sum())), 3),
        "decoy_nulls": round(float(nulls[nulls.is_decoy == 1].n_spike_reads.sum() / max(1, nulls[nulls.is_decoy == 1].n_reads.sum())), 3)}

    # G1
    cal_nulls = nulls[nulls.role == "cal"]
    n_cal = len(cal_nulls)
    k_per_sample = df[(df.role == "test") & (df.y == 1) & (df.n_reads >= 1)].groupby("sample").size()
    k = float(k_per_sample.median()) if len(k_per_sample) else 0.0
    g1 = {"n_cal_nulls": n_cal, "n_cal_target_nulls": int((cal_nulls.is_decoy == 0).sum()),
          "n_cal_decoy_nulls": int((cal_nulls.is_decoy == 1).sum()), "median_k_test": k}
    for q in (0.05, 0.10):
        need = M / (q * k) - 1 if k > 0 else float("inf")
        g1[f"needed_q{q}"] = round(need, 1); g1[f"pass_q{q}"] = bool(n_cal >= need)
    res["G1_power"] = g1

    # G2
    gb, lr, n0, n1 = c2st_auc(nulls[nulls.is_decoy == 0][FEATURES].values, nulls[nulls.is_decoy == 1][FEATURES].values)
    auc = gb if gb is not None else lr
    res["G2_decoy_exchangeability"] = {"auc_lightgbm": gb, "auc_logistic": lr, "n_target_nulls": n0,
                                       "n_decoy_nulls": n1, "pass": None if auc is None else bool(auc <= 0.60)}

    # G3
    dn = nulls[nulls.is_decoy == 1]
    gb, lr, n0, n1 = c2st_auc(dn[dn.role == "cal"][FEATURES].values, dn[dn.role == "test"][FEATURES].values)
    auc = gb if gb is not None else lr
    res["G3_shift"] = {"auc_lightgbm": gb, "auc_logistic": lr, "n_cal": n0, "n_test": n1,
                       "pass": None if auc is None else bool(auc > 0.60)}

    # per-feature means for the write-up (null types and positives)
    grp = hit.assign(kind=np.where(hit.y == 1, "positive", np.where(hit.is_decoy == 1, "decoy_null",
                     np.where(hit.y == 0, "target_null", "ambiguous"))))
    res["feature_means"] = grp.groupby("kind")[FEATURES].mean().round(4).to_dict(orient="index")

    print(json.dumps(res, indent=1, default=float))
    json.dump(res, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
