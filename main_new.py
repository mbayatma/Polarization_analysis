# main_new.py

import os
import argparse
import numpy as np
from plotter.io import load_SM_and_EFT_trees
from plotter.plot_utils import compute_histogram, compute_ratio, plot_comparison
from plotter.config import polarizations, xsecs_SM, xsecs_EFT
from plotter.physics_variables import dispatch_map, label_map  # New variable map and labels

def process_variable(var_name, trees_SM, trees_EFT):
    all_data = []
    hist_SM_all = {}
    hist_EFT_all = {}
    ratio_uncertainties = {}

    all_keys = set(trees_SM.keys()).union(set(trees_EFT.keys()))

    if var_name not in dispatch_map:
        print(f"❌ Unknown variable: {var_name}")
        return

    compute_func = dispatch_map[var_name]

    for pol in all_keys:
        in_SM = pol in trees_SM
        in_EFT = pol in trees_EFT

        data_SM = compute_func(trees_SM[pol]) if in_SM else []
        data_EFT = compute_func(trees_EFT[pol]) if in_EFT else []

        if len(data_SM) == 0 and len(data_EFT) == 0:
            print(f"⚠️ Skipping {pol} — no valid data for {var_name}")
            continue

        if len(data_SM) > 0:
            hist_SM_all[pol] = data_SM
            all_data.extend(data_SM)
        if len(data_EFT) > 0:
            hist_EFT_all[pol] = data_EFT
            all_data.extend(data_EFT)

    if len(all_data) == 0:
        print(f"❌ No valid data found for variable: {var_name}")
        return

    # Define binning
    if var_name.startswith("M_"):
        bins = np.linspace(50, 130, 60)
    elif "pt" in var_name:
        bins = np.append(np.arange(0, 960, 40), [960, 1040])
        #bins = np.arange(0, 1000, 40)  # Coarser: 50 GeV bins
    elif "dphi" in var_name:
        bins = np.linspace(0, 180, 36)
    elif "DR" in var_name:
        bins = np.linspace(0, 6, 30)
    elif "y" in var_name:
        bins = np.linspace(-3, 3, 40) if "d" not in var_name else np.linspace(0, 5, 40)
    elif "pt" in var_name:
        bins = np.concatenate([
            np.linspace(0, 500, 26),
            np.linspace(500, 1000, 10),
            np.linspace(1000, 1400, 5)
        ])
    else:
        bin_min = np.floor(np.min(all_data) / 10) * 10
        bin_max = np.ceil(np.max(all_data) / 10) * 10
        bins = np.linspace(bin_min, bin_max, 50)

    # Compute histograms and ratios
    hists_SM, hists_EFT, ratios = {}, {}, {}

    for pol in hist_SM_all:
        tree_SM = trees_SM[pol]
        h_SM = compute_histogram(hist_SM_all[pol], xsecs_SM.get(pol, 1), tree_SM.num_entries, bins)
        hists_SM[pol] = h_SM

        if pol in hist_EFT_all:
            tree_EFT = trees_EFT[pol]
            h_EFT = compute_histogram(hist_EFT_all[pol], xsecs_EFT.get(pol, 1), tree_EFT.num_entries, bins)
            hists_EFT[pol] = h_EFT
            ratio, ratio_unc = compute_ratio(h_EFT, h_SM)
            ratios[pol] = (ratio, ratio_unc)

    # Save plot
    os.makedirs("polarization_plots_new", exist_ok=True)
    output_path = os.path.join("polarization_plots", f"{var_name}_comparison.png")

    plot_comparison(
        variable_name=var_name,
        hist_sm_dict=hists_SM,
        hist_eft_dict=hists_EFT,
        ratio_dict=ratios,
        bins=bins,
        xlabel=label_map.get(var_name, var_name),
        output_path=output_path,
        ratio_uncertainties=ratio_uncertainties
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", help="Physics variable to plot (e.g. M_ee, DR_epmup)")
    parser.add_argument("--all", action="store_true", help="Plot all physics variables")
    args = parser.parse_args()

    trees_SM, trees_EFT = load_SM_and_EFT_trees(source="spin")

    print("✅ SM trees loaded:")
    for pol, tree in trees_SM.items():
        print(f"  {pol}: {list(tree.keys())[:5]} ...")

    print("\n✅ EFT trees loaded:")
    for pol, tree in trees_EFT.items():
        print(f"  {pol}: {list(tree.keys())[:5]} ...")

    from plotter.physics_variables import dispatch_map

    if args.all:
        print("\n📋 Variables to be plotted:")
        for v in dispatch_map:
            print(" -", v)
        for var_name in dispatch_map:
            print(f"\n📊 Plotting: {var_name}")
            process_variable(var_name, trees_SM, trees_EFT)
    elif args.var:
        process_variable(args.var, trees_SM, trees_EFT)
    else:
        print("❌ Please specify --var <variable> or use --all")

if __name__ == "__main__":
    main()

