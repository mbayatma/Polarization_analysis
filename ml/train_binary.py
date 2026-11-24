"""Train the SM-only binary W_L Z_L vs. rest classifier, then apply it to EFT events.

Usage:
    python -m ml.train_binary --output-dir ml_outputs_reco/binary \
        --sm-path  /eos/user/m/mbayatma/ml_pol_files_reco/polarization_SM.parquet \
        --eft-path /eos/user/m/mbayatma/ml_pol_files_reco/polarization_EFT.parquet

This is step 1 of the README's progression: a single sigmoid asking "is this
event fully longitudinal (W_L Z_L)". It is trained on SM events only and then
run, frozen, over each EFT operator -- the same protocol as
ml/train_two_output.py, so the two models' per-operator significances are
directly comparable and differ only in what the network was asked to predict.

Training on SM alone matters for the interpretation: the model never sees an
EFT event, so any shift in its scores on EFT samples is the response of a
fixed function to a changed input distribution, not something it was
supervised to find. That is why the per-operator number reported here is a
z-score against the SM baseline's own sampling spread rather than a raw shift.
"""

import argparse
import json
import os

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.metrics import auc, roc_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from ml import data as ml_data
from ml.model import build_binary_dnn

DEFAULT_SM_PATH = "/eos/user/m/mbayatma/ml_pol_files_reco/polarization_SM.parquet"
DEFAULT_EFT_PATH = "/eos/user/m/mbayatma/ml_pol_files_reco/polarization_EFT.parquet"

SCORE_COL = "P_LL_binary"


BACKGROUND_LABELS = {"rest": "other polarizations", "WTZT": "$W_T Z_T$"}


def select_background(df, background):
    """Keep the events that take part in the chosen W_L Z_L vs. X comparison.

    For "WTZT" only the two pure states survive (W and Z both longitudinal, or
    both transverse), so the mixed W_L Z_T / W_T Z_L samples are dropped.
    """
    if background == "rest":
        return df
    return df[df["W_is_L"] == df["Z_is_L"]].copy()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sm-path", default=DEFAULT_SM_PATH)
    parser.add_argument("--eft-path", default=DEFAULT_EFT_PATH)
    parser.add_argument("--output-dir", default="ml_outputs_reco/binary")
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--background",
        choices=["rest", "WTZT"],
        default="rest",
        help="What W_L Z_L is separated from: every other polarization (rest), "
             "or only the fully transverse W_T Z_T state (the other two states "
             "are then dropped from training and from EFT scoring).",
    )
    parser.add_argument("--sample-frac", type=float, default=None,
                        help="Subsample before training/inference (fast smoke test).")
    return parser.parse_args()


def main():
    args = parse_args()
    plots_dir = os.path.join(args.output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    tf.random.set_seed(args.seed)
    np.random.seed(args.seed)

    print(f"Loading SM only from {args.sm_path} (EFT is not used for training)")
    df_sm = pd.read_parquet(args.sm_path)
    if args.sample_frac is not None:
        df_sm = df_sm.sample(frac=args.sample_frac, random_state=args.seed)
    df_sm = select_background(ml_data.add_labels(df_sm), args.background)
    bkg_label = BACKGROUND_LABELS[args.background]

    feature_columns = ml_data.get_feature_columns(df_sm)
    print(f"SM (post-filter): {df_sm.shape}, {len(feature_columns)} features")
    print(f"features: {feature_columns}")

    # The target: fully longitudinal (W_L AND Z_L) vs. the chosen background.
    y = ((df_sm["W_is_L"] == 1) & (df_sm["Z_is_L"] == 1)).astype("float32").values
    print(f"Class balance: W_L Z_L = {int(y.sum())} / {len(y)} ({y.mean():.3f})")

    X = df_sm[feature_columns]
    X_train_raw, X_test_raw, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=ml_data.stratify_key(df_sm), random_state=args.seed
    )

    X_train_clean, feature_medians = ml_data.clean_features(X_train_raw)
    X_test_clean, _ = ml_data.clean_features(X_test_raw, feature_medians)

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train_clean).astype("float32")
    X_test = scaler.transform(X_test_clean).astype("float32")

    model = build_binary_dnn(input_dim=X_train.shape[1])
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
    print(f"Trained {len(history.history['loss'])} epochs (patience={args.patience})")

    model.save(os.path.join(args.output_dir, "binary_WLZL_DNN.keras"))
    joblib.dump(scaler, os.path.join(args.output_dir, "scaler.joblib"))
    joblib.dump(feature_medians, os.path.join(args.output_dir, "feature_medians.joblib"))
    with open(os.path.join(args.output_dir, "feature_columns.json"), "w") as f:
        json.dump(feature_columns, f, indent=2)
    with open(os.path.join(args.output_dir, "history.json"), "w") as f:
        json.dump({k: [float(v) for v in vals] for k, vals in history.history.items()}, f, indent=2)

    # --- Held-out SM test split ---
    p_test = model.predict(X_test, batch_size=4096).ravel()
    fpr, tpr, _ = roc_curve(y_test, p_test)
    roc_auc = auc(fpr, tpr)
    print(f"SM test AUC (W_L Z_L vs {args.background}): {roc_auc:.4f}")

    plt.figure(figsize=(6.5, 6.5))
    plt.plot(fpr, tpr, label=f"$W_L Z_L$ vs {bkg_label} (AUC = {roc_auc:.3f})")
    plt.plot([0, 1], [0, 1], "k--", lw=1)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title(f"SM-only binary DNN: $W_L Z_L$ vs {bkg_label}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "roc_binary.png"), dpi=150)
    plt.close()

    plt.figure(figsize=(7, 5))
    bins = np.linspace(0, 1, 51)
    plt.hist(p_test[y_test == 1], bins=bins, histtype="step", density=True,
             label="$W_L Z_L$ (signal)", lw=1.6)
    plt.hist(p_test[y_test == 0], bins=bins, histtype="step", density=True,
             label=bkg_label, lw=1.6)
    plt.xlabel("DNN score $P(W_L Z_L)$")
    plt.ylabel("normalized events")
    plt.title("SM test split: binary classifier output")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "score_distribution.png"), dpi=150)
    plt.close()

    sm_mean, sm_std, n_sm = p_test.mean(), p_test.std(), len(p_test)

    # --- Apply the frozen model to each EFT operator ---
    print(f"Loading EFT from {args.eft_path} for inference only")
    df_eft = pd.read_parquet(args.eft_path)
    if args.sample_frac is not None:
        df_eft = df_eft.sample(frac=args.sample_frac, random_state=args.seed)
    df_eft = select_background(ml_data.add_labels(df_eft), args.background)

    X_eft_clean, _ = ml_data.clean_features(df_eft[feature_columns], feature_medians)
    p_eft = model.predict(scaler.transform(X_eft_clean).astype("float32"), batch_size=4096).ravel()
    df_eft = df_eft.copy()
    df_eft[SCORE_COL] = p_eft
    df_eft[["polarization", "theory", "operator", SCORE_COL]].to_parquet(
        os.path.join(args.output_dir, "eft_scores.parquet"), index=False
    )

    grouped = df_eft.groupby("operator")[SCORE_COL]
    op_mean, op_std, op_n = grouped.mean(), grouped.std(), grouped.count()

    shift = op_mean - sm_mean
    # Combined standard error of (operator mean - SM baseline mean), treating
    # the two as independent samples. This is what turns "the bar is above
    # zero" into a statement about whether it is above *noise*.
    combined_stderr = np.sqrt(op_std**2 / op_n + sm_std**2 / n_sm)
    significance = shift / combined_stderr

    summary = pd.DataFrame({
        "mean_score": op_mean,
        "shift_vs_SM": shift,
        "combined_stderr": combined_stderr,
        "significance": significance,
        "n_events": op_n,
    }).sort_values("significance", key=np.abs, ascending=False)
    summary.to_csv(os.path.join(plots_dir, "operator_significance.csv"))

    print(f"\nSM baseline mean score: {sm_mean:.6f} (std {sm_std:.6f}, n = {n_sm})")
    print("\nPer-operator shift of the SM-trained binary score:")
    print(summary.to_string(float_format=lambda v: f"{v:.4f}"))

    ordered = summary.sort_values("significance")
    plt.figure(figsize=(10, 5.5))
    plt.bar(ordered.index, ordered["significance"],
            color=["tab:red" if abs(v) >= 2 else "tab:blue" for v in ordered["significance"]])
    plt.axhline(0, color="black", lw=1)
    for line_y in (2, -2):
        plt.axhline(line_y, color="gray", lw=0.8, ls=":")
    plt.ylabel("significance (shift / combined std. error)")
    plt.title(f"EFT operator sensitivity: SM-only binary $W_L Z_L$ vs {bkg_label} classifier")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "operator_significance.png"), dpi=150)
    plt.close()

    metrics = {
        "model": "binary_WLZL",
        "background": args.background,
        "sm_path": args.sm_path,
        "auc_sm_test": float(roc_auc),
        "n_sm_train": int(len(X_train)),
        "n_sm_test": int(n_sm),
        "epochs_trained": len(history.history["loss"]),
        "sm_baseline_mean": float(sm_mean),
        "max_abs_significance": float(significance.abs().max()),
        "significance_by_operator": {k: float(v) for k, v in significance.items()},
    }
    with open(os.path.join(args.output_dir, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"\nDone. Artifacts written to {args.output_dir}")


if __name__ == "__main__":
    main()
