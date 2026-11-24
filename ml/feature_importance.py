"""Permutation feature importance for a trained SM-only DNN (two-output or binary).

Usage:
    python -m ml.feature_importance --model-dir ml_outputs_lhe_regen/two_output \
        --sm-path /eos/user/m/mbayatma/ml_pol_files_regen/polarization_SM.parquet --level LHE
    python -m ml.feature_importance --model-dir ml_outputs_reco/two_output \
        --sm-path /eos/user/m/mbayatma/ml_pol_files_reco/polarization_SM.parquet --level reco
    python -m ml.feature_importance --model-type binary --background WTZT \
        --model-dir ml_outputs_reco/binary_WLZL_vs_WTZT \
        --sm-path /eos/user/m/mbayatma/ml_pol_files_reco/polarization_SM.parquet --level reco

Loads the model, scaler and feature medians written by ml/train_two_output.py
or ml/train_binary.py, rebuilds the identical SM test split (same seed,
event selection and stratification), then shuffles one feature at a time and
records how much each output's ROC AUC drops. A large drop means the network leans on that variable; ~0 means it
ignores it. Permutation importance is model-agnostic and is measured on
held-out events, so it reflects what the trained network actually uses rather
than anything about the training procedure.
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
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

from ml import data as ml_data
from ml.train_binary import BACKGROUND_LABELS, select_background

# Axis labels following the analysis' observable table (Table 2): masses and
# transverse momenta as M_{..} / p_{T,..}, rapidities and separations as
# y / |Δy| / Δφ / ΔR with comma-separated subscripts, angles as cos θ*.
# The table is written for the e+ μ+ μ- channel; the LHE variables are named
# the same way, although the samples contain every 3-lepton flavour channel.
FEATURE_LABELS = {
    # --- LHE level (plotter/physics_variables.py) ---
    "M_mumu": r"$M_{\mu^+\mu^-}$ [GeV]",
    "M_epmup": r"$M_{e^+\mu^+}$ [GeV]",
    "M_event_rep": r"$M_{\ell^+\ell^-}$ [GeV]",
    "M_ee_Z": r"$M_{e^+e^-}$ [GeV]",
    "dphi_mumu": r"$\Delta\phi_{\mu^+,\mu^-}$ [$^\circ$]",
    "dphi_epmup": r"$\Delta\phi_{e^+,\mu^+}$ [$^\circ$]",
    "DR_mupmum": r"$\Delta R_{\mu^+,\mu^-}$",
    "DR_epmup": r"$\Delta R_{e^+,\mu^+}$",
    "yep": r"$y_{e^+}$",
    "ymm": r"$y_{\mu^-}$",
    "dyuu": r"$|\Delta y_{\mu^+,\mu^-}|$",
    "dypp": r"$|\Delta y_{e^+,\mu^+}|$",
    "pt_leading_electron": r"$p_{T,e^+}$ [GeV]",
    "pt_leading_muon": r"$p_{T,\mu^-}$ [GeV]",
    "eta_leading_electron": r"$\eta_{e^+}$",
    "eta_leading_muon": r"$\eta_{\mu^-}$",
    "phi_leading_electron": r"$\phi_{e^+}$",
    "phi_leading_muon": r"$\phi_{\mu^-}$",
    "MET": r"$p_{T,\nu_e}$ [GeV]",
    "MET_phi": r"$\phi_{\nu_e}$",
    "dphi_ee_event": r"$\Delta\phi_{e^+,e^-}$ [$^\circ$]",
    "DR_ee_event": r"$\Delta R_{e^+,e^-}$",
    "dphi_mumu_event": r"$\Delta\phi_{\mu^+,\mu^-}$ [$^\circ$] (Z cand.)",
    "DR_mumu_event": r"$\Delta R_{\mu^+,\mu^-}$ (Z cand.)",
    "dphi_emu_event": r"$\Delta\phi_{e,\mu}$ [$^\circ$]",
    "DR_emu_event": r"$\Delta R_{e,\mu}$",
    "cthep": r"$\cos\theta^{*}_{e^+}$",
    "cthmm": r"$\cos\theta^{*}_{\mu^+}$",
    "cscatt": r"$|\cos\theta_{W}|$",
    "cosuu": r"$\cos\Delta\theta_{\mu^+,\mu^-}$",
    "cospp": r"$\cos\Delta\theta_{e^+,\mu^+}$",
    # --- reco level (plotter/reco.py); "rec" marks quantities that need the
    # neutrino p_z from the W-mass constraint, as in the table ---
    "M_Z": r"$M_{\ell^+\ell^-}$ [GeV]",
    "pt_Z": r"$p_{T,\ell^+\ell^-}$ [GeV]",
    "y_Z": r"$y_{\ell^+\ell^-}$",
    "cos_theta_star_Z": r"$\cos\theta^{*}_{\ell^-}$",
    "dphi_Z_leptons": r"$\Delta\phi_{\ell^+,\ell^-}$",
    "dR_Z_leptons": r"$\Delta R_{\ell^+,\ell^-}$",
    "pt_Z_lep1": r"$p_{T,\ell_1^Z}$ [GeV]",
    "pt_Z_lep2": r"$p_{T,\ell_2^Z}$ [GeV]",
    "mT_W": r"$m_{T,W}$ [GeV]",
    "pt_W": r"$p^{\mathrm{rec}}_{T,W}$ [GeV]",
    "y_W": r"$y^{\mathrm{rec}}_{W}$",
    "cos_theta_star_W": r"$\cos\theta^{*\,\mathrm{rec}}_{\ell_W}$",
    "pt_lep_W": r"$p_{T,\ell_W}$ [GeV]",
    "eta_lep_W": r"$\eta_{\ell_W}$",
    "dphi_lepW_MET": r"$\Delta\phi_{\ell_W,\nu}$",
    "M_3l": r"$M_{3\ell}$ [GeV]",
    "pt_3l": r"$p_{T,3\ell}$ [GeV]",
    "M_WZ": r"$M^{\mathrm{rec}}_{4\ell}$ [GeV]",
    "pt_WZ": r"$p^{\mathrm{rec}}_{T,4\ell}$ [GeV]",
    "dphi_Z_W": r"$\Delta\phi_{Z,W}$",
    "dy_Z_W": r"$|\Delta y^{\mathrm{rec}}_{Z,W}|$",
    "dR_lepW_Z": r"$\Delta R_{\ell_W,Z}$",
    "HT": r"$H_T$ [GeV]",
    "n_jet": r"$N_{\mathrm{jets}}$",
    "pt_lead_jet": r"$p_{T,j_1}$ [GeV]",
    "pt_lead_lep": r"$p_{T,\ell_1}$ [GeV]",
    "eta_lead_lep": r"$\eta_{\ell_1}$",
}
# Reco MET is the measured missing p_T, not the true neutrino.
RECO_LABEL_OVERRIDES = {"MET": r"$p_{T}^{\mathrm{miss}}$ [GeV]"}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--sm-path", required=True)
    parser.add_argument("--level", choices=["LHE", "reco"], required=True)
    parser.add_argument("--model-type", choices=["two_output", "binary"], default="two_output")
    parser.add_argument("--background", choices=["rest", "WTZT"], default="rest",
                        help="Binary models only: must match the --background used in training.")
    parser.add_argument("--n-repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42, help="Must match the training seed.")
    return parser.parse_args()


def label_for(feature, level):
    if level == "reco" and feature in RECO_LABEL_OVERRIDES:
        return RECO_LABEL_OVERRIDES[feature]
    return FEATURE_LABELS.get(feature, feature)


def main():
    args = parse_args()
    tf.random.set_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    plots_dir = os.path.join(args.model_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    model_file = ("two_output_polarization_DNN.keras" if args.model_type == "two_output"
                  else "binary_WLZL_DNN.keras")
    model = tf.keras.models.load_model(os.path.join(args.model_dir, model_file))
    scaler = joblib.load(os.path.join(args.model_dir, "scaler.joblib"))
    feature_medians = joblib.load(os.path.join(args.model_dir, "feature_medians.joblib"))

    # Same selection and split as the training script, so these are events the
    # model never trained on.
    df_sm = ml_data.add_labels(pd.read_parquet(args.sm_path))
    if args.model_type == "binary":
        df_sm = select_background(df_sm, args.background)
        y = ((df_sm["W_is_L"] == 1) & (df_sm["Z_is_L"] == 1)).values.astype("float32")[:, None]
        heads = [r"$W_L Z_L$ vs " + BACKGROUND_LABELS[args.background]]
        model_desc = f"binary $W_L Z_L$ vs {BACKGROUND_LABELS[args.background]}"
    else:
        y = df_sm[["W_is_L", "Z_is_L"]].values.astype("float32")
        heads = [r"$W_L$ head", r"$Z_L$ head"]
        model_desc = "two-output"
    feature_columns = ml_data.get_feature_columns(df_sm)
    if list(feature_medians.index) != feature_columns:
        raise SystemExit("Feature columns in the parquet do not match the trained model.")
    _, X_test_raw, _, y_test = train_test_split(
        df_sm[feature_columns], y, test_size=0.2,
        stratify=ml_data.stratify_key(df_sm), random_state=args.seed,
    )
    X_test = scaler.transform(ml_data.clean_features(X_test_raw, feature_medians)[0]).astype("float32")

    def aucs(X):
        pred = model.predict(X, batch_size=8192, verbose=0).reshape(len(X), -1)
        return np.array([roc_auc_score(y_test[:, i], pred[:, i]) for i in range(len(heads))])

    base = aucs(X_test)
    print("Baseline SM test AUC: " + ", ".join(f"{h} {a:.4f}" for h, a in zip(heads, base)))

    rows = []
    for j, feat in enumerate(feature_columns):
        drops = []
        for _ in range(args.n_repeats):
            X_perm = X_test.copy()
            X_perm[:, j] = rng.permutation(X_perm[:, j])
            drops.append(base - aucs(X_perm))
        drops = np.array(drops)
        row = {"feature": feat, "label": label_for(feat, args.level)}
        for i in range(len(heads)):
            row[f"drop_auc_{i}"] = drops[:, i].mean()
            row[f"drop_auc_{i}_std"] = drops[:, i].std()
        rows.append(row)
        print(f"  {feat:22s} dAUC " + "  ".join(f"{drops[:, i].mean():+.4f}" for i in range(len(heads))))

    imp = pd.DataFrame(rows)
    imp["drop_auc_mean"] = imp[[f"drop_auc_{i}" for i in range(len(heads))]].mean(axis=1)
    imp = imp.sort_values("drop_auc_mean", ascending=False).reset_index(drop=True)
    imp.to_csv(os.path.join(plots_dir, "feature_importance.csv"), index=False)

    x = np.arange(len(imp))
    width = 0.8 / len(heads)
    fig, ax = plt.subplots(figsize=(max(10, 0.45 * len(imp)), 6))
    for i, head in enumerate(heads):
        offset = (i - (len(heads) - 1) / 2) * width
        ax.bar(x + offset, imp[f"drop_auc_{i}"], width, yerr=imp[f"drop_auc_{i}_std"], label=head)
    ax.axhline(0, color="black", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(imp["label"], rotation=70, ha="right")
    ax.set_ylabel("Drop in ROC AUC when the variable is shuffled")
    ax.set_title(f"Permutation feature importance, SM-only {model_desc} DNN ({args.level}, SM test split)")
    ax.legend()
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(plots_dir, f"feature_importance.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    with open(os.path.join(plots_dir, "feature_importance_baseline.json"), "w") as f:
        json.dump({"heads": heads, "auc": [float(a) for a in base], "n_test": int(len(y_test)),
                   "n_repeats": args.n_repeats}, f, indent=2)
    print(f"Saved {plots_dir}/feature_importance.(png|pdf|csv)")


if __name__ == "__main__":
    main()
