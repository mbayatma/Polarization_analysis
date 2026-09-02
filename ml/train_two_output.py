"""Train the SM-only two-output polarization DNN, then apply it to EFT events.

Usage:
    python -m ml.train_two_output --output-dir ml_outputs/two_output

Reproduces the methodology of the original two_outpu_training.ipynb notebook:
P(W_L) and P(Z_L) heads are trained on SM events ONLY (no EFT anywhere in
training), then the frozen model is run in inference-only mode on each EFT
operator's events. Any shift in the resulting scores relative to the SM
baseline is a side effect of applying an SM-trained model out of
distribution, not something the model was ever supervised to detect -- so a
per-operator z-score (shift / combined standard error) is used to judge
which shifts are large relative to the sampling noise, rather than reading
raw shift magnitude off a bar chart.

Contrast with ml/train.py's multi-task model, which trains on SM+EFT
together with an explicit SM-vs-EFT head; that answers "can a supervised
classifier tell these apart at all" (it can't -- AUC 0.502). This script
answers a different, narrower question: "does an SM-only polarization model
drift when shown each EFT operator, and by how much relative to noise".
"""

import argparse
import json
import os

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.metrics import auc, roc_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from ml import data as ml_data
from ml.model import build_two_output_dnn

DEFAULT_SM_PATH = "/eos/user/m/mbayatma/ml_pol_files/polarization_SM.parquet"
DEFAULT_EFT_PATH = "/eos/user/m/mbayatma/ml_pol_files/polarization_EFT.parquet"

SCORE_COLS = ["P_WL", "P_ZL", "P_LL", "P_TT", "P_mixed"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sm-path", default=DEFAULT_SM_PATH)
    parser.add_argument("--eft-path", default=DEFAULT_EFT_PATH)
    parser.add_argument("--output-dir", default="ml_outputs/two_output")
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--sample-frac",
        type=float,
        default=None,
        help="Subsample SM and EFT to this fraction before training/inference (for a fast smoke test).",
    )
    return parser.parse_args()


def add_derived_scores(df, P_WL, P_ZL):
    df = df.copy()
    df["P_WL"] = P_WL
    df["P_ZL"] = P_ZL
    df["P_LL"] = P_WL * P_ZL
    df["P_TT"] = (1 - P_WL) * (1 - P_ZL)
    df["P_mixed"] = P_WL * (1 - P_ZL) + (1 - P_WL) * P_ZL
    return df


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    plots_dir = os.path.join(args.output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    tf.random.set_seed(args.seed)
    np.random.seed(args.seed)

    print(f"Loading SM only from {args.sm_path} (EFT is not used for training)")
    df_sm = pd.read_parquet(args.sm_path)
    if args.sample_frac is not None:
        df_sm = df_sm.sample(frac=args.sample_frac, random_state=args.seed)
    df_sm = ml_data.add_labels(df_sm)
    feature_columns = ml_data.get_feature_columns(df_sm)
    print(f"SM (post-filter): {df_sm.shape}, {len(feature_columns)} features")

    X = df_sm[feature_columns]
    y = df_sm[["W_is_L", "Z_is_L"]].values.astype("float32")
    stratify_label = ml_data.stratify_key(df_sm)

    X_train_raw, X_test_raw, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=stratify_label, random_state=args.seed
    )

    X_train_clean, feature_medians = ml_data.clean_features(X_train_raw)
    X_test_clean, _ = ml_data.clean_features(X_test_raw, feature_medians)

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train_clean).astype("float32")
    X_test = scaler.transform(X_test_clean).astype("float32")

    model = build_two_output_dnn(input_dim=X_train.shape[1])
    model.summary()

    early_stop = tf.keras.callbacks.EarlyStopping(
        monitor="val_loss", patience=args.patience, restore_best_weights=True
    )
    history = model.fit(
        X_train, y_train,
        validation_data=(X_test, y_test),
        epochs=args.epochs,
        batch_size=args.batch_size,
        callbacks=[early_stop],
        verbose=2,
    )
    print(f"Trained {len(history.history['loss'])} epochs (early stopping, patience={args.patience})")

    model.save(os.path.join(args.output_dir, "two_output_polarization_DNN.keras"))
    joblib.dump(scaler, os.path.join(args.output_dir, "scaler.joblib"))
    joblib.dump(feature_medians, os.path.join(args.output_dir, "feature_medians.joblib"))
    with open(os.path.join(args.output_dir, "history.json"), "w") as f:
        json.dump({k: [float(v) for v in vals] for k, vals in history.history.items()}, f, indent=2)

    # --- Evaluate on the held-out SM test split ---
    pred_test = model.predict(X_test, batch_size=4096)
    P_WL_test, P_ZL_test = pred_test[:, 0], pred_test[:, 1]

    aucs = {}
    plt.figure(figsize=(7, 7))
    for i, label in enumerate(["W_L", "Z_L"]):
        fpr, tpr, _ = roc_curve(y_test[:, i], pred_test[:, i])
        roc_auc = auc(fpr, tpr)
        aucs[label] = roc_auc
        plt.plot(fpr, tpr, label=f"{label} (AUC = {roc_auc:.3f})")
    plt.plot([0, 1], [0, 1], "k--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("SM-only two-output DNN: polarization ROC (SM test split)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "roc_polarization.png"), dpi=150)
    plt.close()
    print(f"SM test AUCs: {aucs}")

    sm_test_scored = add_derived_scores(
        pd.DataFrame(index=range(len(P_WL_test))), P_WL_test, P_ZL_test
    )
    sm_baseline_mean = sm_test_scored[SCORE_COLS].mean()
    sm_baseline_std = sm_test_scored[SCORE_COLS].std()
    n_sm_test = len(sm_test_scored)

    # --- Apply the frozen, SM-only-trained model to each EFT operator ---
    print(f"Loading EFT from {args.eft_path} for inference only (no training on it)")
    df_eft = pd.read_parquet(args.eft_path)
    if args.sample_frac is not None:
        df_eft = df_eft.sample(frac=args.sample_frac, random_state=args.seed)
    df_eft = ml_data.add_labels(df_eft)
    X_eft = df_eft[feature_columns]
    X_eft_clean, _ = ml_data.clean_features(X_eft, feature_medians)
    X_eft_scaled = scaler.transform(X_eft_clean).astype("float32")

    eft_pred = model.predict(X_eft_scaled, batch_size=4096)
    df_eft_scored = add_derived_scores(df_eft, eft_pred[:, 0], eft_pred[:, 1])
    df_eft_scored.to_parquet(os.path.join(args.output_dir, "eft_scores.parquet"), index=False)

    operator_mean = df_eft_scored.groupby("operator")[SCORE_COLS].mean()
    operator_std = df_eft_scored.groupby("operator")[SCORE_COLS].std()
    operator_n = df_eft_scored.groupby("operator")[SCORE_COLS].count()

    shift = operator_mean.subtract(sm_baseline_mean, axis=1)
    # Combined standard error of (operator mean - SM baseline mean), treating
    # both as independent samples -- this is what turns "the bar is above
    # zero" into an actual claim about whether it's above *noise*.
    sm_stderr = sm_baseline_std / np.sqrt(n_sm_test)
    eft_stderr = operator_std.divide(np.sqrt(operator_n))
    combined_stderr = np.sqrt(eft_stderr**2 + sm_stderr**2)
    significance = shift / combined_stderr

    shift.to_csv(os.path.join(plots_dir, "operator_shifts.csv"))
    significance.to_csv(os.path.join(plots_dir, "operator_significance.csv"))
    print("\nOperator shifts relative to SM baseline (SM-only-trained model):")
    print(shift)
    print("\nSignificance (shift / combined stderr):")
    print(significance)

    shift.plot(kind="bar", figsize=(12, 6))
    plt.axhline(0, color="black", linewidth=1)
    plt.ylabel("Mean score shift relative to SM")
    plt.title("EFT-induced shifts: SM-only-trained two-output DNN")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "operator_shifts.png"), dpi=150)
    plt.close()

    significance.plot(kind="bar", figsize=(12, 6))
    plt.axhline(0, color="black", linewidth=1)
    for line_y in (2, -2, 3, -3):
        plt.axhline(line_y, color="gray", linewidth=0.7, linestyle=":")
    plt.ylabel("Significance (shift / combined std. error)")
    plt.title("EFT operator sensitivity: SM-only-trained two-output DNN")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "operator_significance.png"), dpi=150)
    plt.close()

    metrics = {
        "n_sm_train": int(len(X_train)),
        "n_sm_test": int(n_sm_test),
        "epochs_trained": len(history.history["loss"]),
        "aucs": aucs,
        "sm_baseline_mean": sm_baseline_mean.to_dict(),
        "max_abs_significance_by_operator": significance.abs().max(axis=1).to_dict(),
    }
    with open(os.path.join(args.output_dir, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"\nDone. Artifacts written to {args.output_dir}")


if __name__ == "__main__":
    main()
