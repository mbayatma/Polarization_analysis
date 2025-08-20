# main.py

import os
import argparse
import numpy as np
from plotter.io import load_SM_and_EFT_trees
from plotter.plot_utils import compute_histogram, compute_ratio, plot_comparison
from plotter.config import polarizations, xsecs_SM, xsecs_EFT
from plotter.mass_utils import compute_mass
from plotter.delta_phi_utils import delta_phi_mumu, delta_phi_epmup
from plotter.deltar_utils import compute_delta_r
from plotter.rapidity_utils import compute_rapidity
from plotter.config import rapidity_variables
from plotter.fourvector_utils import compute_pt, compute_eta, compute_phi, compute_mass  as compute_mass_fourvec
from plotter.fourvector_utils import compute_leading_pt, compute_leading_eta, compute_leading_phi

fourvec_map = {
    "pt": compute_pt,
    "eta": compute_eta,
    "phi": compute_phi,
    "mass": compute_mass_fourvec,
}



# Computed variables from four-vectors
computed_kinematics = ["pt", "eta", "phi", "mass", "pt_leading_electron", "pt_leading_muon", "eta_leading_electron", "eta_leading_muon", "phi_leading_electron", "phi_leading_muon"]


# Variable categories
mass_variables = ["M_ee", "M_mumu", "M_epmup"]
delta_phi_variables = ["dphi_mumu", "dphi_epmup"]
dr_variables = ["DR_mupmum", "DR_epmup"]
special_variables = mass_variables + delta_phi_variables + dr_variables + rapidity_variables + computed_kinematics

label_map = {
    "M_ee": r"$M_{e^+e^-}$ [GeV]",
    "M_mumu": r"$M_{\mu^+\mu^-}$ [GeV]",
    "M_epmup": r"$M_{e^+\mu^+}$ [GeV]",
    "dphi_mumu": r"$\Delta\phi(\mu^+, \mu^-)$ [degrees]",
    "dphi_epmup": r"$\Delta\phi(e^+, \mu^+)$ [degrees]",
    "DR_mupmum": r"$\Delta R(\mu^+, \mu^-)$",
    "DR_epmup": r"$\Delta R(e^+, \mu^+)$",
    "yep": r"$y_{e^+}$",
    "ymm": r"$y_{\mu^-}$",
    "yuu": r"$y_{\mu^+\mu^-}$",
    "dyuu": r"$|\Delta y(\mu^+, \mu^-)|$",
    "dypp": r"$|\Delta y(e^+, \mu^+)|$",
    "pt_leading_electron": r"$p_T(e^{+})$ [GeV]",
    "pt_leading_muon": r"$p_T(\mu^{+})$ [GeV]",
}

delta_phi_map = {
    "dphi_mumu": delta_phi_mumu,
    "dphi_epmup": delta_phi_epmup,
}

def extract_leading(tree, branch):
    array = tree[branch].arrays(library="np")[branch]
    if isinstance(array[0], (np.ndarray, list)):
        return np.array([x[0] for x in array if len(x) > 0])
    else:
        return array[np.isfinite(array)]

def get_valid_branches(tree):
    skip_prefixes = ("GenParticle_", "ph_", "photon")
    skip_exact = ("event", "eventWeight", "el_charge", "mu_charge", "nPhotons", "charge", "nElectrons", "nMuons", "nJets", "jet_btag")
    return [b for b in tree.keys() if not any(b.startswith(p) for p in skip_prefixes) and b not in skip_exact]



def process_variable(var_name, trees_SM, trees_EFT):
    all_data = []
    hist_SM_all = {}
    hist_EFT_all = {}
    ratio_uncertainties={}
    
    all_keys = set(trees_SM.keys()).union(set(trees_EFT.keys()))

    for pol in all_keys:
        in_SM = pol in trees_SM
        in_EFT = pol in trees_EFT

        if var_name in mass_variables:
            data_SM = compute_mass(trees_SM[pol], var_name) if in_SM else []
            data_EFT = compute_mass(trees_EFT[pol], var_name) if in_EFT else []
        elif var_name in delta_phi_variables:
            func = delta_phi_map[var_name]
            data_SM = func(trees_SM[pol]) if in_SM else []
            data_EFT = func(trees_EFT[pol]) if in_EFT else []
        elif var_name in dr_variables:
            data_SM = compute_delta_r(trees_SM[pol], var_name) if in_SM else []
            data_EFT = compute_delta_r(trees_EFT[pol], var_name) if in_EFT else []
        elif var_name in rapidity_variables:
            data_SM = compute_rapidity(trees_SM[pol], var_name) if in_SM else []
            data_EFT = compute_rapidity(trees_EFT[pol], var_name) if in_EFT else []
        elif var_name in fourvec_map:
             func = fourvec_map[var_name]
             data_SM = func(trees_SM[pol]) if in_SM else []
             data_EFT = func(trees_EFT[pol]) if in_EFT else []
        elif var_name == "pt_leading_electron":
             data_SM = compute_leading_pt(trees_SM[pol])[0] if in_SM else []
             data_EFT = compute_leading_pt(trees_EFT[pol])[0] if in_EFT else []

        elif var_name == "pt_leading_muon":
             data_SM = compute_leading_pt(trees_SM[pol])[1] if in_SM else []
             data_EFT = compute_leading_pt(trees_EFT[pol])[1] if in_EFT else []
        elif var_name == "eta_leading_electron":
             data_SM = compute_leading_eta(trees_SM[pol])[0] if in_SM else []
             data_EFT = compute_leading_eta(trees_EFT[pol])[0] if in_EFT else []

        elif var_name == "eta_leading_muon":
             data_SM = compute_leading_eta(trees_SM[pol])[1] if in_SM else []
             data_EFT = compute_leading_eta(trees_EFT[pol])[1] if in_EFT else []
        
        elif var_name == "phi_leading_electron":
             data_SM = compute_leading_phi(trees_SM[pol])[0] if in_SM else []
             data_EFT = compute_leading_phi(trees_EFT[pol])[0] if in_EFT else []

        elif var_name == "phi_leading_muon":
            data_SM = compute_leading_phi(trees_SM[pol])[1] if in_SM else []
            data_EFT = compute_leading_phi(trees_EFT[pol])[1] if in_EFT else []
        
        else:
            data_SM = extract_leading(trees_SM[pol], var_name) if in_SM else []
            data_EFT = extract_leading(trees_EFT[pol], var_name) if in_EFT else []
        
        if len(data_SM)==0 and len(data_EFT)==0:
            print(f"⚠️ Skipping {pol} — no valid data for {var_name}")
            continue

        if len(data_SM)>0:
            hist_SM_all[pol] = data_SM
            all_data.extend(data_SM)
        if len(data_EFT)>0:
            hist_EFT_all[pol] = data_EFT
            all_data.extend(data_EFT)

    if len(all_data) == 0:
        print(f"❌ No valid data found for variable: {var_name}")
        return

    # Define binning
    if var_name in mass_variables:
        bins = np.linspace(50, 130, 60)
    elif var_name in delta_phi_variables:
        bins = np.linspace(0, 180, 36)
    elif var_name in dr_variables:
        bins = np.linspace(0, 6, 30)
    elif var_name in rapidity_variables:
        bins = np.linspace(-3, 3, 40) if "y" in var_name and "d" not in var_name else np.linspace(0, 5, 40)
    elif var_name == "pt":
         bins = np.concatenate([
         np.linspace(0, 500, 26),         # fine bins for low pt
         np.linspace(500, 1000, 10),       # medium
         np.linspace(1000, 1400, 5)       # coarse in the tail
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
            ratios[pol]=(ratio, ratio_unc)
            #ratios[pol] = ratio
            #ratio_uncertainties[pol] = ratio_unc
            #ratios[pol] = compute_ratio(h_EFT, h_SM)
            #ratios[pol], ratio_unc[pol]= compute_ratio(h_EFT, h_SM)
    # Save plot
    os.makedirs("polarization_plots", exist_ok=True)
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
    
    
# Variables you want to exclude when plotting with --all
skip_variables = [
    "event", "eventWeight", "eventScale", "alphaEM", "alphaS",
    "mother1", "mother2", "color1", "color2",
    "pdgID", "pdgStatus", "px", "py", "pz", "E",
    "numParticles",  # add more if needed
]    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", help="Branch name to plot (e.g. el_pt, MET)")
    parser.add_argument("--all", action="store_true", help="Plot all branches")
    args = parser.parse_args()

    trees_SM, trees_EFT = load_SM_and_EFT_trees(source="spin")

    print("✅ SM trees loaded:")
    for pol, tree in trees_SM.items():
        print(f"  {pol}: {list(tree.keys())[:5]} ...")

    print("\n✅ EFT trees loaded:")
    for pol, tree in trees_EFT.items():
        print(f"  {pol}: {list(tree.keys())[:5]} ...")

    if args.all:
        ref_pol = polarizations[0] if polarizations[0] in trees_SM else list(trees_SM.keys())[0]
        branches = get_valid_branches(trees_SM[ref_pol])
        #all_variables = branches + special_variables
        all_variables = special_variables
        print("\n📋 Variables to be plotted:")
        for v in all_variables:
            print(" -", v)
        for var_name in all_variables:
            print(f"\n📊 Plotting: {var_name}")
            process_variable(var_name, trees_SM, trees_EFT)
    elif args.var:
        process_variable(args.var, trees_SM, trees_EFT)
    else:
        print("❌ Please specify --var <branch> or use --all")

if __name__ == "__main__":
    main()


