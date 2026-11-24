"""Loading and labeling for the polarization/EFT multi-task dataset.

The exported parquet datasets (see export_ml.py) store one row per event with
kinematic variables plus three string metadata columns: `polarization`,
`theory` (SM or EFT) and `operator` (SM or the EFT operator name).
"""

import numpy as np
import pandas as pd

# Columns that are metadata/labels, never model features.
#
# `weight` and `xsec_weight` come from the reco-level export and must stay out
# of the feature set: xsec_weight takes only one magnitude per sample, so it
# identifies the sample an event came from exactly, and a model given it would
# score the label rather than the kinematics.
META_COLUMNS = [
    "polarization", "theory", "operator",
    "W_is_L", "Z_is_L", "is_EFT",
    "weight", "xsec_weight",
]


def load_dataset(sm_path, eft_path):
    """Load and concatenate the SM and EFT parquet datasets.

    Reads the theory-pure SM/EFT files directly rather than a combined
    "full" export, since export_ml.py's combined file only reflects
    whichever --mode was last used to generate it (see its `mode` guard).
    """
    df_sm = pd.read_parquet(sm_path)
    df_eft = pd.read_parquet(eft_path)
    return pd.concat([df_sm, df_eft], ignore_index=True)


def get_wz_binary_labels(polarization):
    """Map a polarization label to (W_is_L, Z_is_L) in {0, 1}, or (nan, nan).

    The Madspin-showered SM sample isn't a distinct polarization state and
    has no clean W/Z decomposition, so it's excluded (nan).
    """
    if "Madspin" in polarization:
        return np.nan, np.nan
    if "W_L Z_L" in polarization:
        return 1, 1
    if "W_L Z_T" in polarization:
        return 1, 0
    if "W_T Z_L" in polarization:
        return 0, 1
    if "W_T Z_T" in polarization:
        return 0, 0
    return np.nan, np.nan


def add_labels(df):
    """Attach W_is_L, Z_is_L and is_EFT label columns, dropping unlabelable rows.

    Only a handful of distinct polarization strings ever appear, so
    get_wz_binary_labels is evaluated once per distinct value and broadcast
    with a vectorized dict lookup, rather than once per row (a row-wise
    .apply constructing a pd.Series per call is 100x+ slower at ~5M rows).
    """
    df = df.copy()
    labels = {pol: get_wz_binary_labels(pol) for pol in df["polarization"].unique()}
    w_map = {pol: w for pol, (w, _) in labels.items()}
    z_map = {pol: z for pol, (_, z) in labels.items()}

    df["W_is_L"] = df["polarization"].map(w_map)
    df["Z_is_L"] = df["polarization"].map(z_map)
    df = df.dropna(subset=["W_is_L", "Z_is_L"]).copy()
    df["W_is_L"] = df["W_is_L"].astype(int)
    df["Z_is_L"] = df["Z_is_L"].astype(int)
    df["is_EFT"] = (df["theory"] == "EFT").astype(int)
    return df


def get_feature_columns(df):
    return [c for c in df.columns if c not in META_COLUMNS]


def clean_features(X, feature_medians=None):
    """Replace inf with nan, fill missing values, and clip to a safe range.

    Pass `feature_medians` (fit on the training split) when cleaning any
    other split, so held-out data is never used to compute the fill values.
    """
    X = X.replace([np.inf, -np.inf], np.nan)
    if feature_medians is None:
        feature_medians = X.median()
    X = X.fillna(feature_medians)
    X = X.clip(-1e6, 1e6)
    return X.astype(np.float32), feature_medians


def stratify_key(df):
    return (
        df["theory"].astype(str)
        + "_"
        + df["W_is_L"].astype(str)
        + "_"
        + df["Z_is_L"].astype(str)
    )
