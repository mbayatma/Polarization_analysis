"""Combine a SM export and per-operator partial_eft_<op>.parquet files
(written by `export_ml.py --operator <op>`) into the standard
polarization_SM / polarization_EFT / polarization_full parquet files.

Usage:
    python3 merge_partial_exports.py --dir ml_datasets_regen
"""
import argparse
import glob
import os

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", required=True, help="Directory with sm.parquet and partial_eft_*.parquet")
    args = parser.parse_args()

    sm_path = os.path.join(args.dir, "polarization_SM.parquet")
    if not os.path.exists(sm_path):
        raise SystemExit(f"Missing {sm_path} (run export_ml.py --mode SM first)")
    df_sm = pd.read_parquet(sm_path)

    parts = sorted(glob.glob(os.path.join(args.dir, "partial_eft_*.parquet")))
    if not parts:
        raise SystemExit(f"No partial_eft_*.parquet files found in {args.dir}")

    df_eft = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    eft_path = os.path.join(args.dir, "polarization_EFT.parquet")
    df_eft.to_parquet(eft_path, index=False)
    print(f"EFT: {df_eft.shape} <- {len(parts)} parts -> {eft_path}")
    print(df_eft.groupby("operator").size().to_string())

    full_path = os.path.join(args.dir, "polarization_full.parquet")
    pd.concat([df_sm, df_eft], ignore_index=True).to_parquet(full_path, index=False)
    print(f"\nfull: -> {full_path}")


if __name__ == "__main__":
    main()
