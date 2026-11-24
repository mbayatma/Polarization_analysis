"""Random-forest feature ranking for W_L Z_L vs. W_T Z_T (SM events only).

Usage:
    python -m ml.rf_feature_importance --level LHE \
        --sm-path /eos/user/m/mbayatma/ml_pol_files_regen/polarization_SM.parquet \
        --output-dir ml_outputs_lhe_regen/rf_WLZL_vs_WTZT
    python -m ml.rf_feature_importance --level reco \
        --sm-path /eos/user/m/mbayatma/ml_pol_files_reco/polarization_SM.parquet \
        --output-dir ml_outputs_reco/rf_WLZL_vs_WTZT

Trains a random forest on the two pure polarization states (the mixed
W_L Z_T / W_T Z_L samples are dropped) and plots its built-in importance:
the mean decrease in Gini impurity each variable gives, averaged over trees,
with the spread across trees as the error bar. The held-out ROC AUC is
reported alongside so the ranking can be read against how well the forest
separates the two states at all, plus the confusion matrix of the held-out
events at a 0.5 score threshold.
"""

import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, roc_auc_score
from sklearn.model_selection import train_test_split

from ml import data as ml_data
from ml.feature_importance import label_for


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sm-path", required=True)
    parser.add_argument("--level", choices=["LHE", "reco"], required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--n-estimators", type=int, default=300)
    parser.add_argument("--max-depth", type=int, default=12)
    parser.add_argument("--min-samples-leaf", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def plot_confusion_matrix(cm, cm_frac, accuracy, args):
    """Row-normalized confusion matrix; each cell also shows the raw count."""
    names = [r"$W_T Z_T$", r"$W_L Z_L$"]
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(cm_frac, cmap="Blues", vmin=0, vmax=1)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{cm_frac[i, j]:.3f}\n({cm[i, j]})", ha="center", va="center",
                    color="white" if cm_frac[i, j] > 0.5 else "black", fontsize=12)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(names)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(names)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title(f"Random forest ({args.level}, SM test split)\naccuracy = {accuracy:.3f}, threshold 0.5")
    fig.colorbar(im, ax=ax, fraction=0.046, label="Fraction of true class")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.output_dir, f"rf_confusion_matrix.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    df = ml_data.add_labels(pd.read_parquet(args.sm_path))
    df = df[df["W_is_L"] == df["Z_is_L"]]          # keep only W_L Z_L and W_T Z_T
    y = (df["W_is_L"] == 1).astype(int).values      # 1 = W_L Z_L, 0 = W_T Z_T
    features = ml_data.get_feature_columns(df)
    print(f"{len(df)} SM events (W_L Z_L: {y.sum()}, W_T Z_T: {len(y) - y.sum()}), {len(features)} features")

    X_train_raw, X_test_raw, y_train, y_test = train_test_split(
        df[features], y, test_size=0.2, stratify=y, random_state=args.seed
    )
    X_train, medians = ml_data.clean_features(X_train_raw)
    X_test, _ = ml_data.clean_features(X_test_raw, medians)

    rf = RandomForestClassifier(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        min_samples_leaf=args.min_samples_leaf,
        n_jobs=-1,
        random_state=args.seed,
    )
    rf.fit(X_train, y_train)
    p_test = rf.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, p_test)
    print(f"Held-out ROC AUC (W_L Z_L vs W_T Z_T): {auc:.4f}")

    # Rows = true class, columns = predicted class, order (W_T Z_T, W_L Z_L).
    cm = confusion_matrix(y_test, (p_test > 0.5).astype(int), labels=[0, 1])
    cm_frac = cm / cm.sum(axis=1, keepdims=True)
    accuracy = np.trace(cm) / cm.sum()
    print(f"Confusion matrix (rows true, cols predicted; W_T Z_T, W_L Z_L):\n{cm}\naccuracy {accuracy:.4f}")
    plot_confusion_matrix(cm, cm_frac, accuracy, args)

    per_tree = np.array([t.feature_importances_ for t in rf.estimators_])
    imp = pd.DataFrame({
        "feature": features,
        "label": [label_for(f, args.level) for f in features],
        "importance": rf.feature_importances_,
        "std_over_trees": per_tree.std(axis=0),
    }).sort_values("importance", ascending=False).reset_index(drop=True)
    imp.to_csv(os.path.join(args.output_dir, "rf_feature_importance.csv"), index=False)
    print(imp[["feature", "importance"]].to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    fig, ax = plt.subplots(figsize=(8, 0.32 * len(imp) + 1.5))
    ax.barh(imp["label"][::-1], imp["importance"][::-1], xerr=imp["std_over_trees"][::-1],
            color="tab:blue", ecolor="gray")
    ax.set_xlabel("Feature importance (mean decrease in impurity)")
    ax.set_title(f"Random forest, $W_L Z_L$ vs $W_T Z_T$ ({args.level}, SM)\nheld-out AUC = {auc:.3f}")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.output_dir, f"rf_feature_importance.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    with open(os.path.join(args.output_dir, "rf_metrics.json"), "w") as f:
        json.dump({"auc_test": float(auc), "accuracy_test": float(accuracy),
                   "confusion_matrix": {"rows_true_cols_pred": ["W_T Z_T", "W_L Z_L"], "counts": cm.tolist()},
                   "n_train": int(len(y_train)), "n_test": int(len(y_test)),
                   "params": {k: getattr(args, k) for k in ("n_estimators", "max_depth", "min_samples_leaf", "seed")}},
                  f, indent=2)
    print(f"Saved {args.output_dir}/rf_feature_importance.(png|pdf|csv)")


if __name__ == "__main__":
    main()
