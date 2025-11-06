# export_ml.py
import os
import argparse
import numpy as np
import pandas as pd
from plotter.io import load_SM_and_EFT_trees
from plotter.physics_variables import dispatch_map, build_vectors

def export_variables(trees_SM, trees_EFT, output_dir):
    all_dfs = []

    for tag, trees in [("SM", trees_SM), ("EFT", trees_EFT)]:
        for pol, tree in trees.items():
            print(f"➡️ Processing {tag}-{pol} with {tree.num_entries} events")

            # Load all 4-vectors once
            arrays = build_vectors(tree)
            n_events = len(arrays["px"])

            row_dict = {}
            for var_name, func in dispatch_map.items():
                try:
                    values = func(tree)  # compute variable for the whole tree
                    values = np.asarray(values)

                    if len(values) < n_events:
                        # pad with NaNs to match event count
                        padded = np.full(n_events, np.nan)
                        padded[:len(values)] = values
                        row_dict[var_name] = padded
                        print(f"⚠️ {var_name} shorter ({len(values)}) than {n_events}, padded with NaNs")
                    elif len(values) > n_events:
                        # truncate if something went wrong
                        row_dict[var_name] = values[:n_events]
                        print(f"⚠️ {var_name} longer ({len(values)}) than {n_events}, truncated")
                    else:
                        row_dict[var_name] = values
                except Exception as e:
                    print(f"⚠️ Failed {var_name} for {tag}-{pol}: {e}")
                    row_dict[var_name] = np.full(n_events, np.nan)

            # Add labels
            row_dict["polarization"] = [pol] * n_events
            row_dict["theory"] = [tag] * n_events

            df = pd.DataFrame(row_dict)
            all_dfs.append(df)

    full_df = pd.concat(all_dfs, ignore_index=True)

    os.makedirs(output_dir, exist_ok=True)

    # Save SM-only dataset
    df_sm = full_df[full_df["theory"] == "SM"]
    sm_file = os.path.join(output_dir, "polarization_SM.parquet")
    df_sm.to_parquet(sm_file, index=False)
    print(f"✅ Saved SM dataset: {sm_file}")

    # Save full dataset
    full_file = os.path.join(output_dir, "polarization_full.parquet")
    full_df.to_parquet(full_file, index=False)
    print(f"✅ Saved full dataset: {full_file}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", default="ml_datasets",
                        help="Directory to save Parquet files")
    args = parser.parse_args()

    trees_SM, trees_EFT = load_SM_and_EFT_trees(source="spin")
    export_variables(trees_SM, trees_EFT, args.output_dir)


if __name__ == "__main__":
    main()

