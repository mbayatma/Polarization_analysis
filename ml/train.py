"""Train the multi-task polarization + SM-vs-EFT DNN.

Usage:
    python -m ml.train --output-dir ml_outputs
    python -m ml.train --output-dir ml_outputs/smoke --sample-frac 0.01 --epochs 2

Loads the SM and EFT parquet exports directly (not the combined
"polarization_full.parquet", whose contents depend on which --mode
export_ml.py was last run with), builds W/Z polarization + SM-vs-EFT
labels, trains the multi-task network in ml/model.py, and writes the
trained model, preprocessing artifacts, and evaluation plots/tables to
--output-dir.
"""

import argparse
import json
import os
import time

import joblib
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_class_weight

from ml import data as ml_data
from ml import evaluate as ml_evaluate
from ml.model import build_multitask_dnn

DEFAULT_SM_PATH = "/eos/user/m/mbayatma/ml_pol_files/polarization_SM.parquet"
DEFAULT_EFT_PATH = "/eos/user/m/mbayatma/ml_pol_files/polarization_EFT.parquet"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sm-path", default=DEFAULT_SM_PATH)
    parser.add_argument("--eft-path", default=DEFAULT_EFT_PATH)
    parser.add_argument("--output-dir", default="ml_outputs")
    parser.add_argument("--epochs", type=int, default=60)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--test-size", type=float, default=0.2)
    parser.add_argument("--patience", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--sample-frac",
        type=float,
        default=None,
        help="Subsample the dataset to this fraction before training (for a fast smoke test).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    plots_dir = os.path.join(args.output_dir, "plots")

    tf.random.set_seed(args.seed)
    np.random.seed(args.seed)

    print(f"Loading SM from {args.sm_path}")
    print(f"Loading EFT from {args.eft_path}")
    df = ml_data.load_dataset(args.sm_path, args.eft_path)
    print(f"Combined shape: {df.shape}, theory counts:\n{df['theory'].value_counts()}")

    if args.sample_frac is not None:
        df = df.sample(frac=args.sample_frac, random_state=args.seed)
        print(f"Subsampled to {df.shape} (sample_frac={args.sample_frac})")

    df = ml_data.add_labels(df)
    print(f"After label filtering: {df.shape}")
    print(f"Theory label counts:\n{df['is_EFT'].value_counts()}")
    print(f"W/Z label counts:\n{df[['W_is_L', 'Z_is_L']].value_counts()}")

    feature_columns = ml_data.get_feature_columns(df)
    print(f"{len(feature_columns)} feature columns: {feature_columns}")

    X = df[feature_columns]
    y_pol = df[["W_is_L", "Z_is_L"]].values.astype("float32")
    y_theory = df["is_EFT"].values.astype("float32")
    stratify_label = ml_data.stratify_key(df)

    (
        X_train_raw, X_test_raw,
        y_pol_train, y_pol_test,
        y_theory_train, y_theory_test,
        meta_train, meta_test,
    ) = train_test_split(
        X, y_pol, y_theory, df[["theory", "operator", "polarization"]],
        test_size=args.test_size,
        stratify=stratify_label,
        random_state=args.seed,
    )

    X_train_clean, feature_medians = ml_data.clean_features(X_train_raw)
    X_test_clean, _ = ml_data.clean_features(X_test_raw, feature_medians)

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train_clean).astype("float32")
    X_test = scaler.transform(X_test_clean).astype("float32")

    # The theory label is heavily imbalanced (SM is a small minority of the
    # combined SM+EFT sample), so give the theory head balanced sample weights.
    # The polarization labels are already ~50/50 by construction.
    theory_classes = np.unique(y_theory_train)
    theory_class_weights = compute_class_weight(
        class_weight="balanced", classes=theory_classes, y=y_theory_train
    )
    theory_weight_map = dict(zip(theory_classes, theory_class_weights))
    theory_sample_weight = np.array([theory_weight_map[v] for v in y_theory_train], dtype="float32")
    print(f"Theory class weights: {theory_weight_map}")

    model = build_multitask_dnn(input_dim=X_train.shape[1])
    model.summary()

    early_stop = tf.keras.callbacks.EarlyStopping(
        monitor="val_loss", patience=args.patience, restore_best_weights=True
    )

    start = time.time()
    history = model.fit(
        X_train,
        {"pol_output": y_pol_train, "theory_output": y_theory_train},
        sample_weight={
            "pol_output": np.ones(len(y_pol_train), dtype="float32"),
            "theory_output": theory_sample_weight,
        },
        validation_data=(
            X_test,
            {"pol_output": y_pol_test, "theory_output": y_theory_test},
        ),
        epochs=args.epochs,
        batch_size=args.batch_size,
        callbacks=[early_stop],
        verbose=2,
    )
    train_time = time.time() - start
    print(f"Training took {train_time:.1f}s over {len(history.history['loss'])} epochs")

    model.save(os.path.join(args.output_dir, "multitask_polarization_EFT_DNN.keras"))
    joblib.dump(scaler, os.path.join(args.output_dir, "scaler.joblib"))
    joblib.dump(feature_medians, os.path.join(args.output_dir, "feature_medians.joblib"))
    with open(os.path.join(args.output_dir, "feature_columns.json"), "w") as f:
        json.dump(feature_columns, f, indent=2)
    with open(os.path.join(args.output_dir, "history.json"), "w") as f:
        json.dump({k: [float(v) for v in vals] for k, vals in history.history.items()}, f, indent=2)

    ml_evaluate.plot_training_curves(history, plots_dir)

    pred = model.predict(X_test, batch_size=4096)
    pol_pred = pred["pol_output"]
    theory_pred = pred["theory_output"].ravel()

    aucs = ml_evaluate.plot_roc_curves(y_pol_test, pol_pred, y_theory_test, theory_pred, plots_dir)
    ml_evaluate.plot_theory_confusion_matrix(y_theory_test, theory_pred, plots_dir)
    print(f"Test AUCs: {aucs}")

    P_WL, P_ZL = pol_pred[:, 0], pol_pred[:, 1]
    df_scored = meta_test.reset_index(drop=True).copy()
    df_scored["P_WL"] = P_WL
    df_scored["P_ZL"] = P_ZL
    df_scored["P_LL"] = P_WL * P_ZL
    df_scored["P_TT"] = (1 - P_WL) * (1 - P_ZL)
    df_scored["P_mixed"] = P_WL * (1 - P_ZL) + (1 - P_WL) * P_ZL
    df_scored["P_EFT"] = theory_pred

    score_cols = ["P_WL", "P_ZL", "P_LL", "P_TT", "P_mixed", "P_EFT"]
    df_scored.to_parquet(os.path.join(args.output_dir, "test_scores.parquet"), index=False)

    operator_summary, operator_shifts = ml_evaluate.compute_operator_shifts(
        df_scored, score_cols, plots_dir
    )
    print("Operator shifts relative to SM baseline:")
    print(operator_shifts)

    metrics = {
        "train_time_seconds": train_time,
        "epochs_trained": len(history.history["loss"]),
        "n_train": int(len(X_train)),
        "n_test": int(len(X_test)),
        "aucs": aucs,
        "final_val_loss": float(history.history["val_loss"][-1]),
    }
    with open(os.path.join(args.output_dir, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"\nDone. Artifacts written to {args.output_dir}")


if __name__ == "__main__":
    main()
