# main.py

import os
import argparse
import numpy as np
from plotter.io import load_SM_and_EFT_trees
from plotter.plot_utils import compute_histogram, compute_ratio, plot_comparison
from plotter.config import polarizations, xsecs_SM, xsecs_EFT
from plotter.mass_utils import compute_mass

#list of special variables

mass_variables = ["M_ee", "M_mumu", "M_epmup"]
# Axis labels (with LaTeX formatting)
label_map = {
    "M_ee": r"$M_{e^+e^-}$ [GeV]",
    "M_mumu": r"$M_{\mu^+\mu^-}$ [GeV]",
    "M_epmup": r"$M_{e^+\mu^+}$ [GeV]",
}

def extract_leading(tree, branch):
    array = tree[branch].arrays(library="np")[branch]
    if isinstance(array[0], (np.ndarray, list)):
        return np.array([x[0] for x in array if len(x) > 0])
    else:
        return array[np.isfinite(array)]

#Helper to get valid branch names
def get_valid_branches(tree):
    skip_prefixes = ("GenParticle_", "ph_", "photon")
    skip_exact = ("event", "eventWeight", "el_charge", "mu_charge", "nPhoton", "charge")
    return [
        b for b in tree.keys()
        if not any(b.startswith(p) for p in skip_prefixes) and b not in skip_exact
    ]

def process_variable(var_name, trees_SM, trees_EFT):
    all_data = []
    hist_SM_all = {}
    hist_EFT_all = {}

    for pol in polarizations:
        if var_name in mass_variables:
            data_SM = compute_mass(trees_SM[pol], var_name)
            data_EFT = compute_mass(trees_EFT[pol], var_name)
        else:  
            data_SM = extract_leading(trees_SM[pol], var_name)
            data_EFT = extract_leading(trees_EFT[pol], var_name)

        if len(data_SM) == 0 or len(data_EFT) == 0:
            print(f"⚠️ Skipping {pol} — no valid data for {var_name}")
            return

        hist_SM_all[pol] = data_SM
        hist_EFT_all[pol] = data_EFT
        all_data.extend(data_SM)
        all_data.extend(data_EFT)
    
    if var_name in mass_variables:
       bins = np.linspace(50, 130, 60)
    else:
       bin_min = np.floor(np.min(all_data) / 10) * 10
       bin_max = np.ceil(np.max(all_data) / 10) * 10
       bins = np.linspace(bin_min, bin_max, 50)

    hists_SM, hists_EFT, ratios = {}, {}, {}
    for pol in polarizations:
        n_SM = len(hist_SM_all[pol])
        n_EFT = len(hist_EFT_all[pol])
        tree_SM = trees_SM[pol]
        tree_EFT = trees_EFT[pol]

        h_SM = compute_histogram(hist_SM_all[pol], xsecs_SM[pol], tree_SM.num_entries, bins)
        h_EFT = compute_histogram(hist_EFT_all[pol], xsecs_EFT[pol], tree_EFT.num_entries, bins)

        hists_SM[pol] = h_SM
        hists_EFT[pol] = h_EFT
        ratios[pol] = compute_ratio(h_EFT, h_SM)

    # Save plot
    os.makedirs("polarization_plots", exist_ok=True)
    output_path = os.path.join("polarization_plots", f"{var_name}_comparison.png")

    plot_comparison(
        variable_name=var_name,
        hist_sm_dict=hists_SM,
        hist_eft_dict=hists_EFT,
        ratio_dict=ratios,
        bins=bins,
        xlabel=var_name,
        output_path=output_path
    )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", help="Branch name to plot (e.g. el_pt, MET)")
    parser.add_argument("--all", action="store_true", help="Plot all branches")
    args = parser.parse_args()

    trees_SM, trees_EFT = load_SM_and_EFT_trees()

    if args.all:
        branches = get_valid_branches(trees_SM[r"$W_L Z_L$"])
        for var_name in branches:
            print(f"\n📊 Plotting: {var_name}")
            process_variable(var_name, trees_SM, trees_EFT)
    elif args.var:
        process_variable(args.var, trees_SM, trees_EFT)
    else:
        print("❌ Please specify --var <branch> or use --all")

if __name__ == "__main__":
    main()
