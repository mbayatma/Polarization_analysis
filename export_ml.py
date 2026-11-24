import os
import argparse
import numpy as np
import pandas as pd

from plotter.io import load_SM_and_EFT_trees
from plotter.physics_variables import dispatch_map, build_vectors


def process_tree_group(trees, theory, operator, all_dfs):
    """
    Process one group of trees.

    trees format:
        {
            polarization_label: uproot_tree
        }

    Example:
        trees["$W_L Z_L$"]
        trees["$W_T Z_T$"]
    """

    for pol, tree in trees.items():

        print(
            f"➡️ Processing "
            f"{theory} | operator={operator} | "
            f"pol={pol} | events={tree.num_entries}"
        )

        # Build vectors once
        arrays = build_vectors(tree)
        n_events = len(arrays["px"])

        row_dict = {}

        for var_name, func in dispatch_map.items():

            try:
                values = func(tree)
                values = np.asarray(values)

                # Variable shorter than event count
                if len(values) < n_events:

                    padded = np.full(n_events, np.nan)
                    padded[:len(values)] = values

                    row_dict[var_name] = padded

                    print(
                        f"⚠️ {var_name} shorter "
                        f"({len(values)}) than {n_events}, padded with NaNs"
                    )

                # Variable longer than event count
                elif len(values) > n_events:

                    row_dict[var_name] = values[:n_events]

                    print(
                        f"⚠️ {var_name} longer "
                        f"({len(values)}) than {n_events}, truncated"
                    )

                else:
                    row_dict[var_name] = values

            except Exception as e:

                print(
                    f"⚠️ Failed {var_name} "
                    f"for {theory}-{operator}-{pol}: {e}"
                )

                row_dict[var_name] = np.full(n_events, np.nan)

        # Labels
        row_dict["polarization"] = [pol] * n_events
        row_dict["theory"] = [theory] * n_events
        row_dict["operator"] = [operator] * n_events

        df = pd.DataFrame(row_dict)

        all_dfs.append(df)


def export_variables(
    trees_SM,
    trees_EFT,
    output_dir,
    mode="ALL"
):

    all_dfs = []

    print("\nVariables in dispatch_map:\n")

    for v in dispatch_map.keys():
        print("  ", v)

    print("\nNumber of variables:", len(dispatch_map))

    # =========================================================
    # Process SM
    # =========================================================
    if mode in ["SM", "ALL"]:

        print("\n==============================")
        print("Processing SM samples")
        print("==============================\n")

        process_tree_group(
            trees=trees_SM,
            theory="SM",
            operator="SM",
            all_dfs=all_dfs
        )

    # =========================================================
    # Process EFT
    # =========================================================
    if mode in ["EFT", "ALL"]:

        print("\n==============================")
        print("Processing EFT samples")
        print("==============================\n")

        if trees_EFT is None or len(trees_EFT) == 0:

            print("⚠️ No EFT trees found.")

        else:

            first_value = next(iter(trees_EFT.values()))

            # Nested structure:
            # trees_EFT[operator][polarization]
            if isinstance(first_value, dict):

                for operator, trees_op in trees_EFT.items():

                    process_tree_group(
                        trees=trees_op,
                        theory="EFT",
                        operator=operator,
                        all_dfs=all_dfs
                    )

            # Flat structure fallback
            else:

                print(
                    "⚠️ trees_EFT is flat. "
                    "Setting operator='EFT'."
                )

                process_tree_group(
                    trees=trees_EFT,
                    theory="EFT",
                    operator="EFT",
                    all_dfs=all_dfs
                )

    # =========================================================
    # Combine all dataframes
    # =========================================================
    if len(all_dfs) == 0:

        print("⚠️ No dataframes produced.")
        return

    full_df = pd.concat(all_dfs, ignore_index=True)

    print("\n===================================")
    print("Final dataframe summary")
    print("===================================\n")

    print("Shape:", full_df.shape)

    print("\nColumns:")
    print(full_df.columns.tolist())

    print("\nCounts by theory/operator/polarization:\n")

    print(
        full_df.groupby(
            ["theory", "operator", "polarization"]
        ).size().reset_index(name="n_events")
    )

    # =========================================================
    # Save outputs
    # =========================================================
    os.makedirs(output_dir, exist_ok=True)

    # -------------------------
    # Save SM only
    # -------------------------
    if "SM" in full_df["theory"].unique():

        df_sm = full_df[
            full_df["theory"] == "SM"
        ].copy()

        sm_file = os.path.join(
            output_dir,
            "polarization_SM.parquet"
        )

        df_sm.to_parquet(sm_file, index=False)

        print(f"\n✅ Saved SM dataset:\n{sm_file}")

    # -------------------------
    # Save EFT only
    # -------------------------
    if "EFT" in full_df["theory"].unique():

        df_eft = full_df[
            full_df["theory"] == "EFT"
        ].copy()

        eft_file = os.path.join(
            output_dir,
            "polarization_EFT.parquet"
        )

        df_eft.to_parquet(eft_file, index=False)

        print(f"\n✅ Saved EFT dataset:\n{eft_file}")

    # -------------------------
    # Save full dataset
    # -------------------------
    # Only write polarization_full.parquet when it actually contains both
    # theories. Otherwise a `--mode SM` or `--mode EFT` run would silently
    # overwrite the combined file with a theory-only one; downstream code
    # that expects "full" to mean SM+EFT would then break quietly. To
    # rebuild the true full dataset, concatenate the SM and EFT parquet
    # files directly (see ml/data.py:load_dataset).
    if {"SM", "EFT"}.issubset(set(full_df["theory"].unique())):

        full_file = os.path.join(
            output_dir,
            "polarization_full.parquet"
        )

        full_df.to_parquet(full_file, index=False)

        print(f"\n✅ Saved full dataset:\n{full_file}")

    else:
        print(
            "\n⚠️ Skipped polarization_full.parquet: this run only "
            f"produced theory={sorted(full_df['theory'].unique())}, not "
            "both SM and EFT. Run with --mode ALL to regenerate the "
            "true combined file."
        )


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--output_dir",
        default="ml_datasets",
        help="Directory to save parquet files"
    )

    parser.add_argument(
        "--mode",
        choices=["SM", "EFT", "ALL"],
        default="ALL",
        help="Choose which datasets to export"
    )

    parser.add_argument(
        "--operator",
        default=None,
        help=(
            "Export only this EFT operator (skips SM; writes "
            "<output_dir>/partial_eft_<operator>.parquet instead of the "
            "combined files). Lets the 12 operators be exported as separate "
            "parallel processes -- see scratchpad/export_lhe_logs/ for the "
            "driver used to do that and merge_partial_exports.py to combine "
            "the pieces afterwards -- since each operator's per-event "
            "variable computation is the slow part, not the file I/O."
        ),
    )

    args = parser.parse_args()

    if args.operator is not None:
        # Single-operator EFT-only path for parallel export.
        _, trees_EFT = load_SM_and_EFT_trees(source="spin", operator=args.operator)
        if not trees_EFT:
            raise SystemExit(f"No EFT trees found for operator {args.operator}")

        all_dfs = []
        process_tree_group(
            trees=trees_EFT, theory="EFT", operator=args.operator, all_dfs=all_dfs
        )
        df = pd.concat(all_dfs, ignore_index=True)
        os.makedirs(args.output_dir, exist_ok=True)
        out_path = os.path.join(args.output_dir, f"partial_eft_{args.operator}.parquet")
        df.to_parquet(out_path, index=False)
        print(f"\n✅ Saved partial EFT dataset ({args.operator}): {df.shape} -> {out_path}")
        return

    # Load ROOT trees
    trees_SM, trees_EFT = load_SM_and_EFT_trees(
        source="spin"
    )

    # Export parquet datasets
    export_variables(
        trees_SM=trees_SM,
        trees_EFT=trees_EFT,
        output_dir=args.output_dir,
        mode=args.mode
    )


if __name__ == "__main__":
    main()