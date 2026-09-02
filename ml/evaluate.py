"""Evaluation plots and EFT-sensitivity summaries for the multi-task DNN."""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc, confusion_matrix, ConfusionMatrixDisplay, roc_curve


def plot_training_curves(history, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    curves = [
        ("loss", "val_loss", "Total loss", "total_loss.png"),
        ("pol_output_loss", "val_pol_output_loss", "Polarization head loss", "pol_loss.png"),
        ("theory_output_loss", "val_theory_output_loss", "Theory head loss", "theory_loss.png"),
        ("pol_output_auc", "val_pol_output_auc", "Polarization head AUC", "pol_auc.png"),
        ("theory_output_auc", "val_theory_output_auc", "Theory head AUC", "theory_auc.png"),
    ]

    for train_key, val_key, title, filename in curves:
        if train_key not in history.history:
            continue
        plt.figure(figsize=(7, 5))
        plt.plot(history.history[train_key], label="Train")
        plt.plot(history.history[val_key], label="Validation")
        plt.xlabel("Epoch")
        plt.ylabel(title)
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, filename), dpi=150)
        plt.close()


def plot_roc_curves(y_pol_test, pol_pred, y_theory_test, theory_pred, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    aucs = {}

    plt.figure(figsize=(7, 7))
    for i, label in enumerate(["W_L", "Z_L"]):
        fpr, tpr, _ = roc_curve(y_pol_test[:, i], pol_pred[:, i])
        roc_auc = auc(fpr, tpr)
        aucs[label] = roc_auc
        plt.plot(fpr, tpr, label=f"{label} (AUC = {roc_auc:.3f})")
    plt.plot([0, 1], [0, 1], "k--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Polarization head ROC")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "roc_polarization.png"), dpi=150)
    plt.close()

    fpr, tpr, _ = roc_curve(y_theory_test, theory_pred)
    roc_auc = auc(fpr, tpr)
    aucs["EFT_vs_SM"] = roc_auc
    plt.figure(figsize=(7, 7))
    plt.plot(fpr, tpr, label=f"SM vs EFT (AUC = {roc_auc:.3f})")
    plt.plot([0, 1], [0, 1], "k--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Theory head ROC")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "roc_theory.png"), dpi=150)
    plt.close()

    return aucs


def plot_theory_confusion_matrix(y_theory_test, theory_pred, output_dir, threshold=0.5):
    os.makedirs(output_dir, exist_ok=True)
    y_pred = (theory_pred > threshold).astype(int)
    cm = confusion_matrix(y_theory_test, y_pred, labels=[0, 1], normalize="true")
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=["SM", "EFT"])
    disp.plot(cmap="Blues", values_format=".2f")
    plt.title("SM vs EFT Confusion Matrix (normalized)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "confusion_theory.png"), dpi=150)
    plt.close()


def compute_operator_shifts(df_scored, score_cols, output_dir):
    """Mean score per EFT operator vs. the SM baseline, plus a bar-chart plot.

    This is the physics payoff of the model: does the SM-vs-EFT / polarization
    score systematically shift under each EFT operator, relative to pure SM?
    """
    os.makedirs(output_dir, exist_ok=True)

    sm_baseline = df_scored.loc[df_scored["theory"] == "SM", score_cols].mean()

    eft_rows = df_scored[df_scored["theory"] == "EFT"]
    operator_summary = eft_rows.groupby("operator")[score_cols].mean()
    operator_shifts = operator_summary.subtract(sm_baseline, axis=1)

    operator_shifts.plot(kind="bar", figsize=(12, 6))
    plt.axhline(0, color="black", linewidth=1)
    plt.ylabel("Mean score shift relative to SM")
    plt.title("EFT operator shifts in multi-task DNN scores")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "operator_shifts.png"), dpi=150)
    plt.close()

    summary = pd.concat(
        [sm_baseline.to_frame("SM_baseline").T, operator_summary],
        axis=0,
    )
    summary.to_csv(os.path.join(output_dir, "operator_summary.csv"))
    operator_shifts.to_csv(os.path.join(output_dir, "operator_shifts.csv"))

    return operator_summary, operator_shifts
