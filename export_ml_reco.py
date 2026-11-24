"""Export the reco-level (Pythia8 + Delphes) WZ samples to Parquet for ML.

The LHE-level counterpart is export_ml.py. This writes the same three-column
label schema (`polarization`, `theory`, `operator`) that ml/data.py expects,
so the existing training scripts can read these files unchanged -- but the
features are the reco observables from plotter/reco.py rather than the
truth-level ones.

Two extra columns are carried through and excluded from the feature set by
ml/data.py's META_COLUMNS: `xsec_weight` (pb, for cross-section-weighted
work) and `weight` (the raw generator weight).

Usage:
    python3 export_ml_reco.py --output-dir /eos/user/m/mbayatma/ml_pol_files_reco
"""

import argparse
import os

import numpy as np
import pandas as pd

from plotter.config import base_path_RECO, eft_operators, polarizations, pol_to_short, reco_file_template
from plotter import reco


def sample_to_frame(pol_label, sample, theory, operator):
    """Select one sample and return it as a labelled DataFrame."""
    path = os.path.join(
        base_path_RECO,
        reco_file_template.format(pol=pol_to_short[pol_label], sample=sample),
    )
    if not os.path.exists(path):
        print(f"[WARN] missing: {path}")
        return None

    selected = reco.select_wz_3l(reco.load_arrays(path))
    if selected is None:
        print(f"[WARN] no events survive selection: {path}")
        return None

    columns = {var: selected[var] for var in reco.RECO_VARIABLES}
    columns["xsec_weight"] = selected["xsec_weight"]
    columns["weight"] = selected["weight"]

    df = pd.DataFrame(columns)
    n = len(df)
    df["polarization"] = pol_label
    df["theory"] = theory
    df["operator"] = operator

    print(f"  {pol_to_short[pol_label]:<6} {operator:<7} {n:6d} events, "
          f"sigma = {selected['xsec_weight'].sum():.6f} pb")
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="ml_datasets_reco")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print("=== SM ===")
    sm_frames = [
        f for f in (sample_to_frame(pol, "SM", "SM", "SM") for pol in polarizations)
        if f is not None
    ]

    print("\n=== EFT ===")
    eft_frames = []
    for op in eft_operators:
        for pol in polarizations:
            frame = sample_to_frame(pol, op, "EFT", op)
            if frame is not None:
                eft_frames.append(frame)

    if not sm_frames or not eft_frames:
        raise SystemExit("Nothing exported: missing SM or EFT samples.")

    df_sm = pd.concat(sm_frames, ignore_index=True)
    df_eft = pd.concat(eft_frames, ignore_index=True)

    sm_path = os.path.join(args.output_dir, "polarization_SM.parquet")
    eft_path = os.path.join(args.output_dir, "polarization_EFT.parquet")
    df_sm.to_parquet(sm_path, index=False)
    df_eft.to_parquet(eft_path, index=False)

    print(f"\nSM : {df_sm.shape} -> {sm_path}")
    print(f"EFT: {df_eft.shape} -> {eft_path}")
    print("\nEvents per operator:")
    print(df_eft.groupby("operator").size().to_string())
    print("\nEvents per polarization (SM):")
    print(df_sm.groupby("polarization").size().to_string())

    features = [c for c in df_sm.columns
                if c not in ("polarization", "theory", "operator", "xsec_weight", "weight")]
    print(f"\n{len(features)} feature columns: {features}")
    nan_frac = df_sm[features].isna().mean()
    if (nan_frac > 0).any():
        print("\nNaN fraction (SM, non-zero only):")
        print(nan_frac[nan_frac > 0].to_string())
    else:
        print("\nNo NaNs in SM features.")


if __name__ == "__main__":
    main()
