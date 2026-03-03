#!/usr/bin/env python3
"""
Theta and Psi Sign Prediction Using Machine Learning
======================================================
Predicts whether Pitzer mixing parameters θ (theta) and ψ (psi) are
positive or negative for electrolyte mixtures, using electrolyte
physicochemical features from baseline_numeric_only.csv.

Data sources:
  - data/pitzers_mixing_parameters.xlsx  : labeled mixing systems (θ, ψ)
  - data/baseline_numeric_only.csv       : per-electrolyte feature vectors

Pipeline:
  1. Load & parse the mixing-system table (split "HCl-NaCl" → [HCl, NaCl])
  2. Join each component to its baseline feature vector
  3. Build symmetric pairwise features (diffs, sums, absolutes)
  4. Binary-classify sign(θ) and sign(ψ)
  5. Train/validation/test split with stratification
  6. Train five classifiers, select best on validation, evaluate once on test
  7. Plot confusion matrices & feature importances
  8. SHAP analysis (beeswarm summary, mean-|SHAP| bar, waterfall per test sample)
"""

import os
import warnings
import pickle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from sklearn.ensemble import (
    RandomForestClassifier,
    GradientBoostingClassifier,
)
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    roc_auc_score,
    confusion_matrix,
    ConfusionMatrixDisplay,
    classification_report,
)

warnings.filterwarnings("ignore")

# ──────────────────────────────────────────────────────────────────────────────
# 0.  PATHS
# ──────────────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "data")

MIXING_PATH   = os.path.join(DATA_DIR, "pitzers_mixing_parameters.xlsx")
BASELINE_PATH = os.path.join(DATA_DIR, "baseline_numeric_only.csv")
OUT_DIR       = SCRIPT_DIR                 # figures saved alongside this script


# ──────────────────────────────────────────────────────────────────────────────
# 1.  LOAD DATA
# ──────────────────────────────────────────────────────────────────────────────
def load_mixing(path: str) -> pd.DataFrame:
    """Load and clean the Pitzer mixing-parameter table."""
    df = pd.read_excel(path, header=1)
    df.columns = ["system", "exptl", "theta", "psi", "sd_tp", "sd_0", "maxI", "ref"]
    df = df.dropna(subset=["system"])
    df["theta"] = pd.to_numeric(df["theta"], errors="coerce")
    df["psi"]   = pd.to_numeric(df["psi"],   errors="coerce")
    df = df.dropna(subset=["theta", "psi"])
    df = df.reset_index(drop=True)
    return df


def load_baseline(path: str) -> pd.DataFrame:
    """Load the electrolyte feature table and de-duplicate."""
    df = pd.read_csv(path)
    df = df.drop_duplicates(subset="electrolyte")
    df = df.set_index("electrolyte")
    return df


# ──────────────────────────────────────────────────────────────────────────────
# 2.  PARSE MIXING SYSTEMS  ("HCl-NaCl" → ["HCl", "NaCl"])
# ──────────────────────────────────────────────────────────────────────────────
def parse_system(system_str: str) -> list:
    """
    Split a mixing-system label on '-', respecting parentheses.
    E.g. "CaCl2-Ca(NO3)2" → ["CaCl2", "Ca(NO3)2"]
    """
    parts = system_str.split("-")
    elecs, buf = [], ""
    for p in parts:
        buf = f"{buf}-{p}" if buf else p
        if buf.count("(") == buf.count(")"):
            elecs.append(buf)
            buf = ""
    return elecs


# ──────────────────────────────────────────────────────────────────────────────
# 3.  FEATURE ENGINEERING
# ──────────────────────────────────────────────────────────────────────────────
BASELINE_FEATURES = [
    "r_M_angstrom", "r_X_angstrom",
    "anion_1_delta_G_hydration", "anion_1_diffusion_coefficient",
    "anion_1_molecular_weight", "anion_1_n_atoms",
    "anion_1_radius_hydrated", "anion_1_viscosity_jones_dole",
    "cation_1_delta_G_hydration", "cation_1_diffusion_coefficient",
    "cation_1_molecular_weight", "cation_1_n_atoms",
    "cation_1_radius_hydrated", "cation_1_viscosity_jones_dole",
    "B_MX_0_original", "B_MX_1_original", "B_MX_2_original",
    "C_MX_phi_original",
    "electrolyte_type_numeric", "cation_type_numeric", "anion_type_numeric",
    "max_molality_original",
]


def build_feature_matrix(mixing: pd.DataFrame, baseline: pd.DataFrame) -> pd.DataFrame:
    """
    For each mixing pair, join per-electrolyte features and compute
    symmetric pairwise features (sum, difference, abs-difference).
    Returns a feature DataFrame aligned with the mixing table index.
    """
    rows = []
    valid_idx = []

    # keep only features actually present in baseline
    feat_cols = [c for c in BASELINE_FEATURES if c in baseline.columns]

    for idx, row in mixing.iterrows():
        elecs = parse_system(row["system"])
        if len(elecs) != 2:
            continue
        e1, e2 = elecs
        if e1 not in baseline.index or e2 not in baseline.index:
            continue

        f1 = baseline.loc[e1, feat_cols].astype(float)
        f2 = baseline.loc[e2, feat_cols].astype(float)

        feat = {}
        for c in feat_cols:
            feat[f"e1_{c}"]    = f1[c]
            feat[f"e2_{c}"]    = f2[c]
            feat[f"sum_{c}"]   = f1[c] + f2[c]
            feat[f"diff_{c}"]  = f1[c] - f2[c]
            feat[f"absd_{c}"]  = abs(f1[c] - f2[c])

        rows.append(feat)
        valid_idx.append(idx)

    X = pd.DataFrame(rows, index=valid_idx)
    return X


# ──────────────────────────────────────────────────────────────────────────────
# 4.  BUILD LABELS  (sign: 1 = positive, 0 = negative/zero)
# ──────────────────────────────────────────────────────────────────────────────
def build_labels(mixing: pd.DataFrame, valid_idx, param: str) -> pd.Series:
    sub = mixing.loc[valid_idx, param]
    # For psi: the one exactly-zero value is dropped later via valid_idx filtering
    return (sub > 0).astype(int)


# ──────────────────────────────────────────────────────────────────────────────
# 5.  PREPARE DATASET FOR ONE PARAMETER
# ──────────────────────────────────────────────────────────────────────────────
def prepare_dataset(X_full: pd.DataFrame, mixing: pd.DataFrame, param: str):
    """
    For a given parameter ('theta' or 'psi'):
      - Drop rows where the value is exactly 0 (ambiguous sign)
      - Build binary labels
      - Return aligned X, y
    """
    valid_idx = X_full.index.tolist()
    param_vals = mixing.loc[valid_idx, param]

    # Drop zero values
    nonzero_mask = param_vals != 0
    valid_idx_nz = [i for i in valid_idx if nonzero_mask.loc[i]]

    X = X_full.loc[valid_idx_nz].copy()
    y = (mixing.loc[valid_idx_nz, param] > 0).astype(int)

    print(f"\n  {param.upper()} dataset: {len(X)} samples")
    print(f"    Positive (>0): {y.sum()}  |  Negative (<0): {(y==0).sum()}")
    return X, y


# ──────────────────────────────────────────────────────────────────────────────
# 6.  SPLIT & PREPROCESS
# ──────────────────────────────────────────────────────────────────────────────
def _safe_stratified_split(X, y, test_size, random_state):
    """
    Try a stratified split; fall back to random if any class has too few members.
    Prints a warning when the fallback is triggered.
    """
    min_class_count = y.value_counts().min()
    # Need at least 2 of the minority class to put ≥1 in each partition
    needed = int(np.ceil(1.0 / test_size))
    if min_class_count >= needed:
        return train_test_split(X, y, test_size=test_size,
                                stratify=y, random_state=random_state)
    else:
        print(f"    ⚠  Minority class has only {min_class_count} sample(s); "
              f"falling back to non-stratified split.")
        return train_test_split(X, y, test_size=test_size,
                                random_state=random_state)


def split_preprocess(X: pd.DataFrame, y: pd.Series,
                     train_frac=0.60, val_frac=0.20, test_frac=0.20,
                     random_state=42):
    """
    Stratified (where possible) 60 / 20 / 20 train / validation / test split.
    Falls back to random splitting when the minority class is too small for
    stratification. Imputation (mean) and StandardScaler fitted on train only.
    """
    assert abs(train_frac + val_frac + test_frac - 1.0) < 1e-9

    # First split: training vs. temp (val + test)
    X_tr, X_tmp, y_tr, y_tmp = _safe_stratified_split(
        X, y,
        test_size=(val_frac + test_frac),
        random_state=random_state,
    )

    # Second split: validation vs. test
    val_prop = val_frac / (val_frac + test_frac)
    X_val, X_te, y_val, y_te = _safe_stratified_split(
        X_tmp, y_tmp,
        test_size=(1.0 - val_prop),
        random_state=random_state,
    )

    # Impute then scale (fit on train only)
    imputer = SimpleImputer(strategy="mean")
    scaler  = StandardScaler()

    Xtr_s  = scaler.fit_transform(imputer.fit_transform(X_tr))
    Xval_s = scaler.transform(imputer.transform(X_val))
    Xte_s  = scaler.transform(imputer.transform(X_te))

    print(f"    Train {len(X_tr)} | Val {len(X_val)} | Test {len(X_te)}")
    print(f"    Train class dist:  {dict(y_tr.value_counts().sort_index())}")
    print(f"    Val   class dist:  {dict(y_val.value_counts().sort_index())}")
    print(f"    Test  class dist:  {dict(y_te.value_counts().sort_index())}")

    return (X_tr, X_val, X_te,
            Xtr_s, Xval_s, Xte_s,
            y_tr, y_val, y_te,
            imputer, scaler)


# ──────────────────────────────────────────────────────────────────────────────
# 7.  MODELS
# ──────────────────────────────────────────────────────────────────────────────
def build_models(random_state=42) -> dict:
    return {
        "Logistic Regression": LogisticRegression(
            max_iter=2000, class_weight="balanced", random_state=random_state
        ),
        "Random Forest": RandomForestClassifier(
            n_estimators=200, max_depth=6, min_samples_leaf=2,
            class_weight="balanced", random_state=random_state
        ),
        "Gradient Boosting": GradientBoostingClassifier(
            n_estimators=150, max_depth=3, learning_rate=0.05,
            random_state=random_state
        ),
        "SVM (RBF)": SVC(
            kernel="rbf", C=1.0, class_weight="balanced",
            probability=True, random_state=random_state
        ),
        "K-Nearest Neighbours": KNeighborsClassifier(n_neighbors=5),
    }


# ──────────────────────────────────────────────────────────────────────────────
# 8.  TRAIN & EVALUATE ON VALIDATION
# ──────────────────────────────────────────────────────────────────────────────
def train_evaluate(models: dict,
                   Xtr_s, y_tr,
                   Xval_s, y_val,
                   param_name: str) -> pd.DataFrame:
    print(f"\n{'='*70}")
    print(f"  Training models for sign({param_name})")
    print(f"{'='*70}")

    records = []
    for name, clf in models.items():
        clf.fit(Xtr_s, y_tr)

        y_tr_pred  = clf.predict(Xtr_s)
        y_val_pred = clf.predict(Xval_s)

        # AUC (needs predict_proba)
        try:
            val_auc = roc_auc_score(y_val, clf.predict_proba(Xval_s)[:, 1])
        except Exception:
            val_auc = float("nan")

        tr_acc  = accuracy_score(y_tr, y_tr_pred)
        val_acc = accuracy_score(y_val, y_val_pred)
        val_f1  = f1_score(y_val, y_val_pred, zero_division=0)

        print(f"  {name:28s}  "
              f"train_acc={tr_acc:.3f}  "
              f"val_acc={val_acc:.3f}  "
              f"val_F1={val_f1:.3f}  "
              f"val_AUC={val_auc:.3f}")

        records.append({
            "model": name,
            "train_acc": tr_acc,
            "val_acc":   val_acc,
            "val_f1":    val_f1,
            "val_auc":   val_auc,
        })

    results = pd.DataFrame(records).sort_values("val_f1", ascending=False)
    return results


# ──────────────────────────────────────────────────────────────────────────────
# 9.  CROSS-VALIDATE BEST MODEL
# ──────────────────────────────────────────────────────────────────────────────
def cross_validate(clf, Xtr_s, y_tr, cv=5) -> None:
    # Cap folds so every fold has at least 1 minority-class member
    min_class = int(pd.Series(y_tr).value_counts().min())
    cv = min(cv, min_class)
    if cv < 2:
        print(f"\n  ⚠  Only {min_class} minority sample(s) in training set — "
              "skipping cross-validation.")
        return
    skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=42)
    scores = cross_val_score(clf, Xtr_s, y_tr, cv=skf, scoring="f1")
    print(f"\n  {cv}-fold CV F1 scores: "
          + "  ".join(f"{s:.3f}" for s in scores)
          + f"  →  mean={scores.mean():.3f} ±{scores.std():.3f}")


# ──────────────────────────────────────────────────────────────────────────────
# 10.  FINAL TEST EVALUATION
# ──────────────────────────────────────────────────────────────────────────────
def evaluate_test(clf, Xte_s, y_te, model_name: str, param_name: str) -> dict:
    print(f"\n{'='*70}")
    print(f"  FINAL TEST EVALUATION  —  sign({param_name})  —  {model_name}")
    print(f"{'='*70}")

    y_pred = clf.predict(Xte_s)
    acc  = accuracy_score(y_te, y_pred)
    f1   = f1_score(y_te, y_pred, zero_division=0)

    unique_classes = np.unique(y_te)
    try:
        if len(unique_classes) < 2:
            auc = float("nan")
            print(f"  ⚠  Test set contains only class(es) {unique_classes.tolist()} "
                  "— AUC undefined.")
        else:
            auc = roc_auc_score(y_te, clf.predict_proba(Xte_s)[:, 1])
    except Exception:
        auc = float("nan")

    print(f"  Test accuracy : {acc:.4f}")
    print(f"  Test F1       : {f1:.4f}")
    print(f"  Test AUC      : {'n/a' if np.isnan(auc) else f'{auc:.4f}'}")
    print(f"\n{classification_report(y_te, y_pred, labels=[0, 1], target_names=['Negative','Positive'], zero_division=0)}")

    return {"acc": acc, "f1": f1, "auc": auc, "y_pred": y_pred}


# ──────────────────────────────────────────────────────────────────────────────
# 11.  PLOTTING
# ──────────────────────────────────────────────────────────────────────────────
def plot_confusion_matrices(all_results: dict, out_path: str) -> None:
    """
    One confusion-matrix panel per (parameter, split) combination.
    all_results: { param_name: { 'clf': clf, 'Xte_s': ..., 'y_te': ...,
                                  'Xval_s': ..., 'y_val': ..., 'best_name': ... } }
    """
    params  = list(all_results.keys())          # ["theta", "psi"]
    splits  = ["Validation", "Test"]
    n_rows, n_cols = len(splits), len(params)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
    fig.suptitle("Confusion Matrices — Sign Prediction (θ and ψ)",
                 fontsize=15, fontweight="bold", y=1.01)

    for col, param in enumerate(params):
        res  = all_results[param]
        clf  = res["clf"]
        name = res["best_name"]

        for row, split in enumerate(splits):
            ax = axes[row, col]
            Xs = res["Xval_s"] if split == "Validation" else res["Xte_s"]
            yt = res["y_val"]  if split == "Validation" else res["y_te"]

            y_pred = clf.predict(Xs)
            # Force both classes [0,1] in the confusion matrix
            cm = confusion_matrix(yt, y_pred, labels=[0, 1])
            disp = ConfusionMatrixDisplay(cm, display_labels=["Negative", "Positive"])
            disp.plot(ax=ax, colorbar=False, cmap="Blues")

            acc = accuracy_score(yt, y_pred)
            f1  = f1_score(yt, y_pred, zero_division=0, labels=[0, 1], average="binary")
            ax.set_title(
                f"sign({param})  |  {split}\n"
                f"{name}\nAcc={acc:.3f}  F1={f1:.3f}",
                fontsize=10,
            )

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\n  ✓  Confusion matrices saved → {out_path}")
    plt.close()


def plot_feature_importance(clf, feature_names: list,
                            param_name: str, out_path: str,
                            top_n=20, model_label: str = None) -> None:
    """Bar chart of feature importances (tree models) or |coefficients| (linear)."""
    if hasattr(clf, "feature_importances_"):
        importances = clf.feature_importances_
        kind = "Gini importance"
    elif hasattr(clf, "coef_"):
        importances = np.abs(clf.coef_[0])
        kind = "|coefficient|"
    else:
        print(f"  (feature importance not available for {type(clf).__name__})")
        return

    label = model_label or type(clf).__name__
    idx = np.argsort(importances)[::-1][:top_n]
    names  = [feature_names[i] for i in idx]
    values = importances[idx]

    fig, ax = plt.subplots(figsize=(9, max(4, top_n * 0.35)))
    ax.barh(range(len(names))[::-1], values, color="steelblue", edgecolor="k", linewidth=0.4)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names[::-1], fontsize=8)
    ax.set_xlabel(kind, fontsize=10)
    ax.set_title(
        f"Top-{top_n} Feature Importances — sign({param_name})\n({label})",
        fontsize=12, fontweight="bold"
    )
    ax.grid(axis="x", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  ✓  Feature importance saved → {out_path}")
    plt.close()


def plot_class_distributions(mixing: pd.DataFrame, valid_idx, out_path: str) -> None:
    """Bar charts showing class balance for θ and ψ sign."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    fig.suptitle("Class Balance — Sign of θ (theta) and ψ (psi)",
                 fontsize=13, fontweight="bold")

    for ax, param, color in zip(axes, ["theta", "psi"], ["steelblue", "salmon"]):
        vals = mixing.loc[valid_idx, param]
        vals = vals[vals != 0]
        counts = pd.Series(
            ["Positive" if v > 0 else "Negative" for v in vals]
        ).value_counts()

        ax.bar(counts.index, counts.values, color=color, edgecolor="k", linewidth=0.6)
        ax.set_title(f"sign({param})", fontsize=12)
        ax.set_ylabel("Count", fontsize=10)
        for i, (label, cnt) in enumerate(counts.items()):
            ax.text(i, cnt + 0.2, str(cnt), ha="center", fontsize=11, fontweight="bold")
        ax.set_ylim(0, counts.max() + 5)
        ax.grid(axis="y", alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  ✓  Class distribution plot saved → {out_path}")
    plt.close()


def plot_validation_comparison(results_theta: pd.DataFrame,
                               results_psi:   pd.DataFrame,
                               out_path: str) -> None:
    """Grouped bar chart comparing validation metrics across models for θ and ψ."""
    metrics = ["val_acc", "val_f1", "val_auc"]
    labels  = ["Accuracy", "F1", "AUC"]
    params  = [("theta", results_theta), ("psi", results_psi)]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
    fig.suptitle("Validation Metrics Comparison Across Models",
                 fontsize=14, fontweight="bold")

    bar_width = 0.22
    x = np.arange(len(results_theta))

    for ax, (param, results) in zip(axes, params):
        for j, (metric, label) in enumerate(zip(metrics, labels)):
            ax.bar(x + j * bar_width,
                   results[metric].values,
                   width=bar_width,
                   label=label,
                   edgecolor="k", linewidth=0.4)

        ax.set_xticks(x + bar_width)
        ax.set_xticklabels(results["model"].values, rotation=20, ha="right", fontsize=8)
        ax.set_ylim(0, 1.1)
        ax.set_title(f"sign({param})", fontsize=12)
        ax.set_ylabel("Score", fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(axis="y", alpha=0.3)

        # highlight best model
        best_idx = results["val_f1"].values.argmax()
        ax.axvspan(best_idx - 0.15, best_idx + 3 * bar_width + 0.05,
                   color="gold", alpha=0.25, zorder=0, label="Best (F1)")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  ✓  Validation comparison saved → {out_path}")
    plt.close()


# ──────────────────────────────────────────────────────────────────────────────
# 12.  TEST PERFORMANCE PLOT
# ──────────────────────────────────────────────────────────────────────────────
def plot_test_performance(all_results: dict,
                          mixing: pd.DataFrame,
                          out_path: str) -> None:
    """
    Two-row figure — one row per parameter (θ, ψ).
    Each row has three panels:
      Panel 1 — Predicted probability per test sample, sorted by probability,
                 coloured by true class; correctly vs. wrongly classified
                 are distinguished by marker shape.
      Panel 2 — ROC curve for all five trained models (where both classes are
                 present in the test set) or a simple accuracy/F1/AUC bar if
                 only one class is present.
      Panel 3 — Confusion matrix on the test set.
    """
    params      = list(all_results.keys())
    n_rows      = len(params)
    fig         = plt.figure(figsize=(18, 6 * n_rows))
    fig.suptitle("Best-Model Performance on Test Data — Sign of θ and ψ",
                 fontsize=15, fontweight="bold", y=1.01)

    COLORS  = {1: "#2ecc71", 0: "#e74c3c"}   # green=positive, red=negative
    MARKERS = {"correct": "o", "wrong": "X"}

    for row_idx, param in enumerate(params):
        res       = all_results[param]
        clf       = res["clf"]
        Xte_s     = res["Xte_s"]
        y_te      = res["y_te"]
        best_name = res["best_name"]
        X_te_df   = res["X_te"]
        models    = res["models"]
        tm        = res["test_metrics"]

        # system labels for the test set
        sys_labels = mixing.loc[X_te_df.index, "system"].tolist()

        # predicted probabilities
        try:
            proba = clf.predict_proba(Xte_s)[:, 1]
        except Exception:
            proba = clf.decision_function(Xte_s)
            proba = (proba - proba.min()) / (proba.ptp() + 1e-9)

        y_pred   = clf.predict(Xte_s)
        y_true   = y_te.values
        n_test   = len(y_true)
        correct  = (y_pred == y_true)

        # sort samples by predicted probability (ascending)
        order    = np.argsort(proba)
        proba_s  = proba[order]
        y_true_s = y_true[order]
        correct_s= correct[order]
        labels_s = [sys_labels[i] for i in order]

        # ── Panel 1: predicted probability per sample ─────────────────────
        ax1 = fig.add_subplot(n_rows, 3, row_idx * 3 + 1)

        for i in range(n_test):
            color  = COLORS[y_true_s[i]]
            marker = MARKERS["correct"] if correct_s[i] else MARKERS["wrong"]
            msize  = 110 if correct_s[i] else 140
            ax1.scatter(i, proba_s[i], c=color, marker=marker,
                        s=msize, zorder=3, edgecolors="k", linewidths=0.5)

        # decision threshold
        ax1.axhline(0.5, color="grey", linestyle="--", linewidth=1.2,
                    label="Threshold = 0.50")

        ax1.set_xticks(range(n_test))
        ax1.set_xticklabels(labels_s, rotation=45, ha="right", fontsize=7)
        ax1.set_ylim(-0.05, 1.05)
        ax1.set_ylabel("P(positive)", fontsize=10)
        ax1.set_title(
            f"sign({param}) — Predicted Probability per Test Sample\n"
            f"Best model: {best_name}",
            fontsize=10, fontweight="bold"
        )
        ax1.grid(axis="y", alpha=0.35)

        # custom legend
        from matplotlib.lines import Line2D
        legend_handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=COLORS[1],
                   markersize=9, markeredgecolor="k", label="True Positive"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=COLORS[0],
                   markersize=9, markeredgecolor="k", label="True Negative"),
            Line2D([0], [0], marker="X", color="w", markerfacecolor="grey",
                   markersize=9, markeredgecolor="k", label="Misclassified"),
            Line2D([0], [0], linestyle="--", color="grey", label="Threshold"),
        ]
        ax1.legend(handles=legend_handles, fontsize=7, loc="upper left",
                   framealpha=0.85)

        # ── Panel 2: ROC curves (all models) or metric bars ──────────────
        ax2 = fig.add_subplot(n_rows, 3, row_idx * 3 + 2)
        unique_classes = np.unique(y_true)

        if len(unique_classes) >= 2:
            from sklearn.metrics import roc_curve

            for mname, mclf in models.items():
                try:
                    mproba = mclf.predict_proba(Xte_s)[:, 1]
                except Exception:
                    continue
                fpr, tpr, _ = roc_curve(y_true, mproba)
                try:
                    auc_val = roc_auc_score(y_true, mproba)
                except Exception:
                    auc_val = float("nan")
                lw  = 2.5 if mname == best_name else 1.0
                lst = "-"  if mname == best_name else "--"
                ax2.plot(fpr, tpr, lw=lw, linestyle=lst,
                         label=f"{mname}  (AUC={auc_val:.2f})")

            ax2.plot([0, 1], [0, 1], "k:", lw=0.8, label="Random")
            ax2.set_xlabel("False Positive Rate", fontsize=10)
            ax2.set_ylabel("True Positive Rate", fontsize=10)
            ax2.set_title(f"ROC Curves on Test Set — sign({param})",
                          fontsize=10, fontweight="bold")
            ax2.legend(fontsize=7, loc="lower right")
            ax2.grid(alpha=0.3)
            ax2.set_xlim([0, 1])
            ax2.set_ylim([0, 1.05])

        else:
            # Only one class in test — show accuracy / F1 bar chart instead
            metric_names = ["Accuracy", "F1"]
            metric_vals  = [tm["acc"], tm["f1"]]
            bars = ax2.bar(metric_names, metric_vals,
                           color=["steelblue", "salmon"],
                           edgecolor="k", linewidth=0.7, width=0.4)
            for bar, val in zip(bars, metric_vals):
                ax2.text(bar.get_x() + bar.get_width() / 2,
                         val + 0.02, f"{val:.3f}",
                         ha="center", va="bottom", fontsize=11, fontweight="bold")
            ax2.set_ylim(0, 1.2)
            ax2.set_title(
                f"Test Metrics — sign({param})\n"
                f"(ROC n/a: single class in test set)",
                fontsize=10, fontweight="bold"
            )
            ax2.grid(axis="y", alpha=0.3)
            ax2.text(0.5, 0.5,
                     f"⚠  Only 1 class in test set\n"
                     f"(theta is heavily skewed\n"
                     f"toward positive values)",
                     transform=ax2.transAxes,
                     ha="center", va="center", fontsize=9,
                     color="grey", style="italic")

        # ── Panel 3: confusion matrix ─────────────────────────────────────
        ax3 = fig.add_subplot(n_rows, 3, row_idx * 3 + 3)
        cm   = confusion_matrix(y_true, y_pred, labels=[0, 1])
        disp = ConfusionMatrixDisplay(cm, display_labels=["Negative", "Positive"])
        disp.plot(ax=ax3, colorbar=False, cmap="Blues")

        acc = accuracy_score(y_true, y_pred)
        f1  = f1_score(y_true, y_pred, zero_division=0, labels=[0, 1],
                       average="binary")
        auc_str = ("n/a" if np.isnan(tm["auc"])
                   else f"{tm['auc']:.3f}")
        ax3.set_title(
            f"Test Confusion Matrix — sign({param})\n"
            f"{best_name}\n"
            f"Acc={acc:.3f}  F1={f1:.3f}  AUC={auc_str}",
            fontsize=10, fontweight="bold"
        )

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  ✓  Test performance plot saved → {out_path}")
    plt.close()


# ──────────────────────────────────────────────────────────────────────────────
# 13.  SHAP ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────
def _shap_values_pos_class(clf, Xtr_s: np.ndarray,
                            X_explain: np.ndarray) -> np.ndarray:
    """
    Return a 2-D SHAP array (n_samples × n_features) for the POSITIVE class,
    using the most appropriate explainer for the model type.
    """
    import shap

    if isinstance(clf, LogisticRegression):
        exp = shap.LinearExplainer(clf, Xtr_s)
        sv  = exp.shap_values(X_explain)
        # LinearExplainer returns 2-D (n, p) directly for binary logistic
        return np.array(sv)

    if isinstance(clf, (RandomForestClassifier, GradientBoostingClassifier)):
        exp = shap.TreeExplainer(clf)
        sv  = exp.shap_values(X_explain)
        # TreeExplainer may return list [neg, pos] or 3-D array (n, p, 2).
        if isinstance(sv, list):
            sv = np.array(sv[1])
        else:
            sv = np.array(sv)
        if sv.ndim == 3:
            sv = sv[:, :, 1]
        return sv

    # SVM / KNN / anything else → KernelExplainer with a kmeans background
    import shap
    n_bg  = min(10, len(Xtr_s))
    bg    = shap.kmeans(Xtr_s, n_bg)
    exp   = shap.KernelExplainer(clf.predict_proba, bg)
    sv    = exp.shap_values(X_explain, nsamples=256)
    if isinstance(sv, list):
        sv = np.array(sv[1])
    else:
        sv = np.array(sv)

    # Normalise to 2-D positive-class slice regardless of SHAP version
    if sv.ndim == 3:
        sv = sv[:, :, 1]
    return sv


# ── matplotlib helpers ────────────────────────────────────────────────────────
def _shorten(name: str) -> str:
    """Compact display label for a pairwise feature name."""
    return (name
            .replace("absd_", "Δ|")
            .replace("diff_", "Δ")
            .replace("sum_",  "Σ")
            .replace("e1_",   "A·")
            .replace("e2_",   "B·"))


def _beeswarm_ax(ax, sv: np.ndarray, fv: np.ndarray,
                 feat_names: list, title: str, top_n: int = 15) -> None:
    """
    Draw a SHAP beeswarm plot on *ax*.
    sv  : (n_samples, n_features) SHAP values
    fv  : (n_samples, n_features) raw feature values (for colouring)
    """
    n_feat = min(top_n, sv.shape[1])
    order  = np.argsort(np.abs(sv).mean(axis=0))[::-1][:n_feat]
    order  = order[::-1]                              # bottom → top
    cmap   = plt.cm.coolwarm
    rng    = np.random.RandomState(0)

    for row, fi in enumerate(order):
        col_sv = sv[:, fi]
        col_fv = fv[:, fi]
        vlo    = np.nanpercentile(col_fv, 5)
        vhi    = np.nanpercentile(col_fv, 95)
        normed = (np.clip(col_fv, vlo, vhi) - vlo) / max(vhi - vlo, 1e-9)
        jitter = rng.uniform(-0.22, 0.22, len(col_sv))
        ax.scatter(col_sv, row + jitter,
                   c=cmap(normed), s=28, alpha=0.75,
                   edgecolors="none", zorder=3)

    ax.axvline(0, color="k", lw=0.8, ls="--", alpha=0.5)
    ax.set_yticks(range(n_feat))
    ax.set_yticklabels([_shorten(feat_names[i]) for i in order], fontsize=7)
    ax.set_xlabel("SHAP value  →  positive class", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, shrink=0.55, pad=0.02)
    cb.set_label("Feature value\nlow ← → high", fontsize=7)
    cb.ax.tick_params(labelsize=7)


def _bar_ax(ax, sv: np.ndarray, feat_names: list,
            title: str, top_n: int = 15) -> None:
    """
    Draw a mean-|SHAP| bar chart on *ax*.
    """
    mean_abs = np.abs(sv).mean(axis=0)
    top_idx  = np.argsort(mean_abs)[::-1][:top_n][::-1]   # bottom → top
    ax.barh(range(top_n), mean_abs[top_idx],
            color="steelblue", edgecolor="k", linewidth=0.3, alpha=0.85)
    ax.set_yticks(range(top_n))
    ax.set_yticklabels([_shorten(feat_names[i]) for i in top_idx], fontsize=7)
    ax.set_xlabel("Mean |SHAP value|", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)


def _waterfall_ax(ax, sv_1d: np.ndarray, feat_names: list,
                  base_val: float, pred_prob: float,
                  title: str, top_n: int = 10) -> None:
    """
    Draw a SHAP waterfall chart on *ax*.
    Shows the top *top_n* features and a single "All others" bar.
    """
    n   = len(sv_1d)
    idx = np.argsort(np.abs(sv_1d))[::-1]

    # Top features + remainder
    top_i  = idx[:top_n]
    rest   = sv_1d[idx[top_n:]].sum() if n > top_n else 0.0

    items = [(feat_names[i], sv_1d[i]) for i in top_i]
    if n > top_n:
        items.append(("All others", rest))

    items = items[::-1]          # bottom → top

    running = base_val
    for row, (name, val) in enumerate(items):
        color = "#e74c3c" if val >= 0 else "#3498db"
        ax.barh(row, val, left=running,
                color=color, edgecolor="k", linewidth=0.35, height=0.6,
                alpha=0.85)
        running += val
        ax.text(running + (0.002 if val >= 0 else -0.002), row,
                f"{val:+.3f}",
                va="center",
                ha="left" if val >= 0 else "right",
                fontsize=6.5, color="k")

    ax.axvline(base_val, color="grey", ls="--", lw=0.9,
               label=f"E[f(X)] = {base_val:.3f}")
    ax.axvline(pred_prob, color="k", ls="-", lw=1.5,
               label=f"f(x) = {pred_prob:.3f}")

    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([_shorten(it[0]) for it in items], fontsize=7)
    ax.set_xlabel("Model output  P(positive)", fontsize=9)
    ax.set_title(title, fontsize=9, fontweight="bold")
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(axis="x", alpha=0.3)

    # Shade positive / negative regions
    ax.axvspan(0.5, ax.get_xlim()[1] if ax.get_xlim()[1] > 0.5 else 1,
               alpha=0.05, color="green")
    ax.axvspan(ax.get_xlim()[0] if ax.get_xlim()[0] < 0.5 else 0, 0.5,
               alpha=0.05, color="red")


# ── Main SHAP orchestration ───────────────────────────────────────────────────
def plot_shap_analysis(all_results: dict,
                       mixing: pd.DataFrame,
                       out_dir: str) -> None:
    """
    Generates two SHAP figures:

    1. ``theta_psi_shap_summary.png``
       2 × 2 grid: beeswarm (left) + mean-|SHAP| bar (right) for each parameter.
       SHAP values computed on the TRAINING set so the summary reflects the
       patterns the model actually learnt.

    2. ``theta_psi_shap_waterfall.png``
       One waterfall panel per TEST SAMPLE for each parameter, showing how
       individual features drive each prediction away from the base value.
    """
    import shap                                       # noqa – imported lazily

    print("\n" + "=" * 70)
    print("  SHAP ANALYSIS")
    print("=" * 70)

    params = list(all_results.keys())

    # ── Figure 1: Summary (beeswarm + bar) ───────────────────────────────
    n_rows = len(params)
    fig_s, axes_s = plt.subplots(n_rows, 2, figsize=(16, 7 * n_rows))
    if n_rows == 1:
        axes_s = axes_s[np.newaxis, :]
    fig_s.suptitle(
        "SHAP Feature Contributions — Sign of θ (theta) and ψ (psi)\n"
        "SHAP values computed on training set",
        fontsize=14, fontweight="bold", y=1.01,
    )

    # ── Figure 2: Waterfall per test sample ──────────────────────────────
    max_te  = max(len(all_results[p]["y_te"]) for p in params)
    fig_w   = plt.figure(figsize=(14, 5 * n_rows * max_te))
    gs_w    = gridspec.GridSpec(
        n_rows * max_te, 1, figure=fig_w, hspace=0.55
    )
    fig_w.suptitle(
        "SHAP Waterfall — Individual Test-Sample Predictions",
        fontsize=13, fontweight="bold", y=1.005,
    )
    wf_row = 0

    for row_idx, param in enumerate(params):
        res       = all_results[param]
        clf       = res["clf"]
        Xtr_s     = res["Xtr_s"]
        Xte_s     = res["Xte_s"]
        y_te      = res["y_te"]
        X_te_df   = res["X_te"]
        feat_names = res["feature_names"]
        best_name  = res["best_name"]

        # ── Compute SHAP on TRAIN (for summary) ──────────────────────────
        print(f"\n  sign({param})  [{best_name}]")
        print(f"    Computing SHAP on training set ({len(Xtr_s)} samples) …",
              end="", flush=True)
        sv_tr = _shap_values_pos_class(clf, Xtr_s, Xtr_s)
        print(" done.")

        # ── Compute SHAP on TEST (for waterfall) ─────────────────────────
        print(f"    Computing SHAP on test set    ({len(Xte_s)} samples) …",
              end="", flush=True)
        sv_te = _shap_values_pos_class(clf, Xtr_s, Xte_s)
        print(" done.")

        # Expected value (base value) from training predictions
        try:
            base_val = clf.predict_proba(Xtr_s)[:, 1].mean()
        except Exception:
            base_val = 0.5

        # ── Summary beeswarm (Panel left) ────────────────────────────────
        _beeswarm_ax(
            axes_s[row_idx, 0], sv_tr, Xtr_s, feat_names,
            title=f"SHAP Beeswarm — sign({param})\n{best_name}  (training set)",
        )

        # ── Mean |SHAP| bar (Panel right) ────────────────────────────────
        _bar_ax(
            axes_s[row_idx, 1], sv_tr, feat_names,
            title=f"Mean |SHAP| — sign({param})\n{best_name}  (training set)",
        )

        # ── Waterfall per test sample ─────────────────────────────────────
        sys_labels = mixing.loc[X_te_df.index, "system"].tolist()
        try:
            proba_te = clf.predict_proba(Xte_s)[:, 1]
        except Exception:
            proba_te = np.full(len(Xte_s), np.nan)

        for i, (sv_row, sys_name, prob, ytrue) in enumerate(
            zip(sv_te, sys_labels, proba_te, y_te.values)
        ):
            ax_wf = fig_w.add_subplot(gs_w[wf_row])
            wf_row += 1

            true_str = "Positive (+)" if ytrue == 1 else "Negative (−)"
            pred_str = "Positive (+)" if prob >= 0.5 else "Negative (−)"
            correct  = (ytrue == 1) == (prob >= 0.5)
            status   = "✓ correct" if correct else "✗ wrong"

            _waterfall_ax(
                ax_wf, sv_row, feat_names, base_val, prob,
                title=(f"sign({param})  |  system: {sys_name}  |  "
                       f"True: {true_str}   Pred: {pred_str}   {status}"),
            )

        # pad unused rows in the waterfall figure
        for _ in range(len(y_te), max_te):
            ax_wf = fig_w.add_subplot(gs_w[wf_row])
            wf_row += 1
            ax_wf.axis("off")

    # ── Save summary figure ───────────────────────────────────────────────
    fig_s.tight_layout()
    sum_path = os.path.join(out_dir, "theta_psi_shap_summary.png")
    fig_s.savefig(sum_path, dpi=150, bbox_inches="tight")
    print(f"\n  ✓  SHAP summary saved → {sum_path}")
    plt.close(fig_s)

    # ── Save waterfall figure ─────────────────────────────────────────────
    fig_w.tight_layout()
    wf_path = os.path.join(out_dir, "theta_psi_shap_waterfall.png")
    fig_w.savefig(wf_path, dpi=150, bbox_inches="tight")
    print(f"  ✓  SHAP waterfall saved → {wf_path}")
    plt.close(fig_w)


# ──────────────────────────────────────────────────────────────────────────────
# 14.  SAVE ARTIFACTS
# ──────────────────────────────────────────────────────────────────────────────
def save_model(clf, imputer, scaler, feature_names: list,
               param: str, out_dir: str) -> None:
    artifact = {
        "model":         clf,
        "imputer":       imputer,
        "scaler":        scaler,
        "feature_names": feature_names,
        "param":         param,
        "task":          f"sign({param}) classification",
    }
    path = os.path.join(out_dir, f"best_{param}_sign_model.pkl")
    with open(path, "wb") as f:
        pickle.dump(artifact, f)
    print(f"  ✓  Model artifact saved → {path}")


# ──────────────────────────────────────────────────────────────────────────────
# 15.  MAIN
# ──────────────────────────────────────────────────────────────────────────────
def run_pipeline(param: str,
                 X_full: pd.DataFrame,
                 mixing: pd.DataFrame,
                 all_results: dict,
                 val_results_store: dict) -> None:
    """
    Full classification pipeline for one parameter ('theta' or 'psi').
    Populates all_results[param] and val_results_store[param] in-place.
    """
    print(f"\n{'#'*70}")
    print(f"  PARAMETER: sign({param.upper()})")
    print(f"{'#'*70}")

    # 4 + 5. Prepare dataset
    X, y = prepare_dataset(X_full, mixing, param)

    if len(X) < 10:
        print(f"  ⚠ Too few samples ({len(X)}) — skipping {param}")
        return

    feature_names = X.columns.tolist()

    # 6. Split & preprocess
    (X_tr, X_val, X_te,
     Xtr_s, Xval_s, Xte_s,
     y_tr, y_val, y_te,
     imputer, scaler) = split_preprocess(X, y, random_state=42)

    # 7. Build models
    models = build_models(random_state=42)

    # 8. Train & evaluate on validation
    results_df = train_evaluate(models, Xtr_s, y_tr, Xval_s, y_val, param)
    val_results_store[param] = results_df

    # Best model by validation F1
    best_name = results_df.iloc[0]["model"]
    best_clf  = models[best_name]
    print(f"\n  🏆  Best model (val F1): {best_name}")

    # 9. Cross-validate best model on training data
    cross_validate(best_clf, Xtr_s, y_tr)

    # 10. Final test evaluation (ONCE)
    test_metrics = evaluate_test(best_clf, Xte_s, y_te, best_name, param)

    # Store for plotting
    all_results[param] = {
        "clf":        best_clf,
        "best_name":  best_name,
        "Xtr_s":      Xtr_s,
        "Xval_s":     Xval_s,
        "Xte_s":      Xte_s,
        "X_te":       X_te,          # original DataFrame, index = mixing row ids
        "y_tr":       y_tr,
        "y_val":      y_val,
        "y_te":       y_te,
        "test_metrics": test_metrics,
        "feature_names": feature_names,
        "imputer":    imputer,
        "scaler":     scaler,
        "models":     models,        # all trained models (for ROC comparison)
    }

    # Save model artifact
    save_model(best_clf, imputer, scaler, feature_names, param, OUT_DIR)

    # Feature importance plot — use best model if possible, else fallback
    fi_clf, fi_name = best_clf, best_name
    if not (hasattr(best_clf, "feature_importances_") or
            hasattr(best_clf, "coef_")):
        # Try to find a tree-based model from the trained dict for this param
        for alt_name, alt_clf in models.items():
            if hasattr(alt_clf, "feature_importances_"):
                fi_clf, fi_name = alt_clf, alt_name
                print(f"  (feature importance: falling back to {fi_name})")
                break
    fi_path = os.path.join(OUT_DIR, f"{param}_sign_feature_importance.png")
    plot_feature_importance(fi_clf, feature_names, param, fi_path,
                            model_label=fi_name)


def main():
    print("=" * 70)
    print("  THETA & PSI SIGN PREDICTION — Pitzer Mixing Parameters")
    print("=" * 70)

    # ── 1. Load data ──────────────────────────────────────────────────────
    print("\n[1] Loading data …")
    mixing   = load_mixing(MIXING_PATH)
    baseline = load_baseline(BASELINE_PATH)
    print(f"  Mixing table    : {len(mixing)} rows")
    print(f"  Baseline table  : {len(baseline)} electrolytes, "
          f"{len(baseline.columns)} columns")

    # ── 2 + 3. Build feature matrix ───────────────────────────────────────
    print("\n[2] Parsing systems and building pairwise features …")
    X_full = build_feature_matrix(mixing, baseline)
    valid_idx = X_full.index.tolist()

    n_dropped = len(mixing) - len(X_full)
    print(f"  Matched pairs   : {len(X_full)} "
          f"({n_dropped} dropped — electrolyte not in baseline)")
    print(f"  Feature columns : {X_full.shape[1]}")

    # ── Class-balance overview ────────────────────────────────────────────
    dist_path = os.path.join(OUT_DIR, "theta_psi_sign_class_distribution.png")
    plot_class_distributions(mixing, valid_idx, dist_path)

    # ── Run pipeline for both parameters ──────────────────────────────────
    all_results      = {}
    val_results_store = {}

    for param in ["theta", "psi"]:
        run_pipeline(param, X_full, mixing, all_results, val_results_store)

    if len(all_results) < 2:
        print("\n⚠  Insufficient data for one or both parameters — exiting.")
        return

    # ── Validation comparison plot ────────────────────────────────────────
    vc_path = os.path.join(OUT_DIR, "theta_psi_sign_validation_comparison.png")
    plot_validation_comparison(
        val_results_store["theta"],
        val_results_store["psi"],
        vc_path,
    )

    # ── Confusion matrices (val + test, both params) ──────────────────────
    cm_path = os.path.join(OUT_DIR, "theta_psi_sign_confusion_matrices.png")
    plot_confusion_matrices(all_results, cm_path)

    # ── Best-model test performance (prob strip + ROC + confusion matrix) ─
    tp_path = os.path.join(OUT_DIR, "theta_psi_sign_test_performance.png")
    plot_test_performance(all_results, mixing, tp_path)

    # ── SHAP analysis ─────────────────────────────────────────────────────
    plot_shap_analysis(all_results, mixing, OUT_DIR)

    # ── Summary ───────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    for param, res in all_results.items():
        tm = res["test_metrics"]
        print(f"\n  sign({param})  —  best model: {res['best_name']}")
        print(f"    Test accuracy : {tm['acc']:.4f}")
        print(f"    Test F1       : {tm['f1']:.4f}")
        auc_str = "n/a" if np.isnan(tm['auc']) else f"{tm['auc']:.4f}"
        print(f"    Test AUC      : {auc_str}")

    print("\n  Generated files:")
    files = [
        "theta_psi_sign_class_distribution.png",
        "theta_psi_sign_validation_comparison.png",
        "theta_psi_sign_confusion_matrices.png",
        "theta_psi_sign_test_performance.png",
        "theta_psi_shap_summary.png",
        "theta_psi_shap_waterfall.png",
        "theta_sign_feature_importance.png",
        "psi_sign_feature_importance.png",
        "best_theta_sign_model.pkl",
        "best_psi_sign_model.pkl",
    ]
    for f in files:
        print(f"    {os.path.join(OUT_DIR, f)}")

    print("\n" + "=" * 70)
    print("  DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
