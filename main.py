# main.py

import os
import argparse
import numpy as np
from plotter.physics_variables import build_vectors
from plotter.io import load_SM_and_EFT_trees
from plotter.plot_utils import compute_histogram, compute_ratio, plot_comparison
from plotter.config import (polarizations, xsecs_SM, xsecs_EFT_all, eft_operators)
from plotter.physics_variables import (
    dispatch_map,
    label_map,
    count_event_categories,   # NEW: event counters
)

# Set from --diagnostics; gates the slow per-event category counting.
SHOW_COUNTS = False

# Variables that are naturally event-level (useful for coverage debug prints)
EVENT_LEVEL_VARS = {
    "M_event_rep",
    "pt_leading_electron", "pt_leading_muon",
    "eta_leading_electron", "eta_leading_muon",
    "phi_leading_electron", "phi_leading_muon",
    "MET", "MET_phi",
}

def process_variable(var_name, trees_SM, trees_EFT, xsecs_EFT, operator):
    all_data = []
    hist_SM_all = {}
    hist_EFT_all = {}
    ratio_uncertainties = {}

    all_keys = set(trees_SM.keys()).union(set(trees_EFT.keys()))

    if var_name not in dispatch_map:
        print(f"Unknown variable: {var_name}")
        return

    compute_func = dispatch_map[var_name]

    for pol in sorted(all_keys):
        in_SM = pol in trees_SM
        in_EFT = pol in trees_EFT

        data_SM = compute_func(trees_SM[pol]) if in_SM else []
        data_EFT = compute_func(trees_EFT[pol]) if in_EFT else []
        
                # --- NEW DEBUG: count SFOS and DF events actually contributing to histograms ---
        def flavour_category(var_name):
            """Return 'SFOS' if same-flavour variable, 'DF' if different-flavour."""
            if any(x in var_name for x in ["ee", "mumu"]):
                return "SFOS"
            elif any(x in var_name for x in ["epmup", "em", "mue"]):
                return "DF"
            else:
                return "All"

        category = flavour_category(var_name)

        if len(data_SM) > 0:
            n_SM_entries = len(data_SM)
            total_SM_events = getattr(trees_SM[pol], "num_entries", len(data_SM))
            print(f"[FILL][{pol}][SM]  {category} events filled = "
                  f"{n_SM_entries:>7}/{total_SM_events:<7}  "
                  f"({n_SM_entries/total_SM_events*100:5.2f}%)")

        if len(data_EFT) > 0:
            n_EFT_entries = len(data_EFT)
            total_EFT_events = getattr(trees_EFT[pol], "num_entries", len(data_EFT))
            print(f"[FILL][{pol}][EFT] {category} events filled = "
                  f"{n_EFT_entries:>7}/{total_EFT_events:<7}  "
                  f"({n_EFT_entries/total_EFT_events*100:5.2f}%)")

       

        # Debug: event categories & coverage (only print if we actually have this pol).
        # count_event_categories walks every event, and process_variable is
        # called once per variable, so this is gated to avoid repeating the
        # same per-polarization count ~30 times per operator.
        if in_SM and SHOW_COUNTS:
            counts = count_event_categories(trees_SM[pol])  # unweighted by default
            cov_msg = ""
            if var_name in EVENT_LEVEL_VARS and hasattr(trees_SM[pol], "num_entries"):
                cov = (len(data_SM) / max(1, trees_SM[pol].num_entries)) * 100.0
                cov_msg = f" | coverage(SM) ≈ {cov:.1f}%"
            print(f"[COUNTS][{pol}][SM]  SFOS={counts['unweighted']['N_SFOS']}, "
                  f"DF={counts['unweighted']['N_DF']}, DFOS={counts['unweighted']['N_DFOS']}{cov_msg}")

        if in_EFT and SHOW_COUNTS:
            counts = count_event_categories(trees_EFT[pol])
            cov_msg = ""
            if var_name in EVENT_LEVEL_VARS and hasattr(trees_EFT[pol], "num_entries"):
                cov = (len(data_EFT) / max(1, trees_EFT[pol].num_entries)) * 100.0
                cov_msg = f" | coverage(EFT) ≈ {cov:.1f}%"
            print(f"[COUNTS][{pol}][EFT] SFOS={counts['unweighted']['N_SFOS']}, "
                  f"DF={counts['unweighted']['N_DF']}, DFOS={counts['unweighted']['N_DFOS']}{cov_msg}")

        if len(data_SM) == 0 and len(data_EFT) == 0:
            print(f"⚠️  Skipping {pol} — no valid data for {var_name}")
            continue

        if len(data_SM) > 0:
            hist_SM_all[pol] = np.asarray(data_SM, dtype=float)
            all_data.extend(hist_SM_all[pol])
        if len(data_EFT) > 0:
            hist_EFT_all[pol] = np.asarray(data_EFT, dtype=float)
            all_data.extend(hist_EFT_all[pol])

    if len(all_data) == 0:
        print(f"❌ No valid data found for variable: {var_name}")
        return

    # -----------------
    # Define binning
    # -----------------
    if var_name.startswith("M_"):  # includes M_ee, M_mumu, M_epmup, and M_event_rep
        bins = np.linspace(50, 130, 60)
    elif "pt_leading_electron" in var_name:
        bins = np.concatenate([
            np.arange(0, 600, 40),
            [600, 800, 1000, 1200],
        ])
    elif "pt_leading_muon" in var_name:
        bins = np.concatenate([
            np.arange(0, 400, 40),
            np.arange(400, 800, 100),
            [1000, 1200],
        ])
    elif "dphi" in var_name:
        bins = np.linspace(0, 180, 36)
    elif "DR" in var_name:
        bins = np.linspace(0, 6, 30)
    elif "y" in var_name:
        bins = np.linspace(-3, 3, 40) if "d" not in var_name else np.linspace(0, 5, 40)
    elif var_name in ("cthep", "cthmm", "cosuu", "cospp"):
        bins = np.linspace(-1.0, 1.0, 41)
    elif var_name in ("cscatt", "crecscatt"):
        bins = np.linspace(0,1.0, 41)
    else:
        bin_min = np.floor(np.min(all_data) / 10) * 10
        bin_max = np.ceil(np.max(all_data) / 10) * 10
        if bin_min == bin_max:
            bin_min -= 1
            bin_max += 1
        bins = np.linspace(bin_min, bin_max, 50)

    # -----------------
    # Compute histograms and ratios
    # -----------------
    hists_SM, hists_EFT, ratios = {}, {}, {}

    for pol in hist_SM_all:
        tree_SM = trees_SM[pol]
        h_SM = compute_histogram(hist_SM_all[pol], xsecs_SM.get(pol, 1), tree_SM.num_entries, bins)
        hists_SM[pol] = h_SM

        # Debug print for SM
        xsec_SM = xsecs_SM.get(pol, 1)
        print(f"[DEBUG] {pol} (SM): xsec = {xsec_SM}, "
              f"N_events = {tree_SM.num_entries}, "
              f"xsec/event = {xsec_SM / tree_SM.num_entries:.3e}")
              
        #AREA CONSISTENCY CHECK for SM
        # ============================================================
        nfill_SM = len(hist_SM_all[pol])
        fill_frac_SM = nfill_SM / max(1, tree_SM.num_entries)
        area_SM = float(h_SM.sum())
        expected_area_SM = xsec_SM * fill_frac_SM
        rel_err_SM = (area_SM - expected_area_SM) / expected_area_SM if expected_area_SM > 0 else 0.0
        print(f"[AREA ABS][{pol}][SM] area={area_SM:.6e}  expected={expected_area_SM:.6e}  "
          f"fill_frac={fill_frac_SM:.4f}  rel.err={rel_err_SM:.2%}")


        if pol in hist_EFT_all:
            tree_EFT = trees_EFT[pol]
            xsec_EFT = xsecs_EFT.get(pol, 1)
            h_EFT = compute_histogram(hist_EFT_all[pol], xsec_EFT, tree_EFT.num_entries, bins)
            hists_EFT[pol] = h_EFT

            ratio, ratio_unc = compute_ratio(h_EFT, h_SM)
            ratios[pol] = (ratio, ratio_unc)

            # Debug print for EFT
            print(f"[DEBUG] {pol} (EFT): xsec = {xsec_EFT}, "
                  f"N_events = {tree_EFT.num_entries}, "
                  f"xsec/event = {xsec_EFT / tree_EFT.num_entries:.3e}")
            # AREA CONSISTENCY CHECK for EFT
            # ============================================================
            nfill_EFT = len(hist_EFT_all[pol])
            fill_frac_EFT = nfill_EFT / max(1, tree_EFT.num_entries)
            area_EFT = float(h_EFT.sum())
            expected_area_EFT = xsec_EFT * fill_frac_EFT
            rel_err_EFT = (area_EFT - expected_area_EFT) / expected_area_EFT if expected_area_EFT > 0 else 0.0
            print(f"[AREA ABS][{pol}][EFT] area={area_EFT:.6e}  expected={expected_area_EFT:.6e}  "
              f"fill_frac={fill_frac_EFT:.4f}  rel.err={rel_err_EFT:.2%}")
        

            # AREA CHECK
            sm_area = h_SM.sum()
            eft_area = h_EFT.sum()
            expected_ratio = xsec_EFT / xsec_SM if xsec_SM else np.nan
            area_ratio = (eft_area / sm_area) if sm_area > 0 else np.nan
            print(f"[AREA CHECK] {pol}: EFT/SM = {area_ratio:.3f}, expected ≈ {expected_ratio:.3f}")

    # -----------------
    # Save plot
    # -----------------
    operator_dir =os.path.join("polariztion_plots_new", operator)
    os.makedirs(operator_dir, exist_ok=True)
    output_path = os.path.join(operator_dir, f"{var_name}_comparison.png")

    plot_comparison(
        variable_name=var_name,
        hist_sm_dict=hists_SM,
        hist_eft_dict=hists_EFT,
        ratio_dict=ratios,
        bins=bins,
        xlabel=label_map.get(var_name, var_name),
        output_path=output_path,
        ratio_uncertainties=ratio_uncertainties,
        operator=operator, 
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", help="Physics variable to plot (e.g. M_event_rep, M_ee, DR_epmup)")
    parser.add_argument("--all", action="store_true", help="Plot all physics variables")
    parser.add_argument("--operator", help="EFT operator name (e.g. cHWB). If omitted with --all, runs all operators.")
    parser.add_argument(
        "--diagnostics",
        action="store_true",
        help="Print per-event selection/decay-composition diagnostics (slow: loops over every event).",
    )
    args = parser.parse_args()

    global SHOW_COUNTS
    SHOW_COUNTS = args.diagnostics

    # Decide which operators to run
    if args.operator:
        operators_to_run = [args.operator]
    else:
        # No operator specified: if --all, run all; otherwise we don't know what to do
        if not args.all:
            print("❌ Please specify --operator <OP> or use --all to run all operators.")
            return
        operators_to_run = eft_operators

    if not args.all and not args.var:
        print("❌ Please specify --var <variable> or use --all.")
        return

    # Loop over operators. Everything below is per-operator: loading the
    # trees, the optional diagnostics, and the plotting all have to stay
    # inside this loop, or only the last operator ends up being plotted.
    for op in operators_to_run:
        print(f"\n🔧 Using EFT operator: {op}")

        # Get cross sections for this operator
        if op not in xsecs_EFT_all:
            print(f"⚠️  No xsecs_EFT found for operator {op}, skipping.")
            continue
        xsecs_EFT = xsecs_EFT_all[op]

        # Load trees for this operator
        trees_SM, trees_EFT = load_SM_and_EFT_trees(source="spin", operator=op)

        print("✅ SM trees loaded:")
        for pol, tree in trees_SM.items():
            print(f"  {pol}: {list(tree.keys())[:5]} ...")

        print("\n✅ EFT trees loaded:")
        for pol, tree in trees_EFT.items():
            print(f"  {pol}: {list(tree.keys())[:5]} ...")

        if args.diagnostics:
            run_diagnostics(trees_SM, trees_EFT)

        variables = list(dispatch_map) if args.all else [args.var]
        if args.all:
            print(f"\n📋 Variables to be plotted for operator {op}:")
            for v in variables:
                print(" -", v)

        for var_name in variables:
            print(f"\n📊 Plotting: {var_name} (operator {op})")
            process_variable(var_name, trees_SM, trees_EFT, xsecs_EFT, op)


def run_diagnostics(trees_SM, trees_EFT):
    """Selection-efficiency and decay-composition printouts.

    Both loop in Python over every event of every polarization, so this is
    opt-in (--diagnostics) rather than something every plotting run pays for.
    """
    from plotter.physics_variables import _event_mask_3l_trigger

    print("\n✅ Applying 3-lepton + single-lepton trigger selection check:")
    for tag, trees in (("SM", trees_SM), ("EFT", trees_EFT)):
        for pol, t in trees.items():
            raw = t.arrays(["px", "py", "pz", "E", "pdgID"], library="np")
            mask = _event_mask_3l_trigger(raw)
            n_pass = sum(mask)
            n_total = len(raw["pdgID"])
            eff = 100.0 * n_pass / max(1, n_total)
            print(f"[SELECTION][{pol}][{tag}]  pass = {n_pass:6d} / {n_total:6d}  ({eff:5.1f}%)")
    print("-------------------------------------------------------------\n")

    def breakdown(tree):
        arr = build_vectors(tree)
        pdg = arr["pdgID"]
        z_ee = z_mm = w_e = w_mu = 0
        for ids in pdg:
            has_ep = any(pid == -11 for pid in ids)
            has_em = any(pid ==  11 for pid in ids)
            has_mup = any(pid == -13 for pid in ids)
            has_mum = any(pid ==  13 for pid in ids)
            if has_ep and has_em:
               z_ee += 1
            if has_mup and has_mum:
               z_mm += 1
            if any(abs(pid)==11 for pid in ids):
               w_e += 1
            if any(abs(pid)==13 for pid in ids):
              w_mu += 1
        return dict(Zee=z_ee, Zmm=z_mm, We=w_e, Wmu=w_mu, total=len(pdg))

    print("\n📊 Decay-channel composition check:")
    for pol, t in trees_SM.items():
        b = breakdown(t)
        print(f"[DECAYS][{pol}][SM]  Zee={b['Zee']:>6}  Zmm={b['Zmm']:>6}  W→e={b['We']:>6}  W→μ={b['Wmu']:>6}  total={b['total']:>6}")
    for pol, t in trees_EFT.items():
        b = breakdown(t)
        print(f"[DECAYS][{pol}][EFT] Zee={b['Zee']:>6}  Zmm={b['Zmm']:>6}  W→e={b['We']:>6}  W→μ={b['Wmu']:>6}  total={b['total']:>6}")
    print("-------------------------------------------------------------\n")

if __name__ == "__main__":
    main()

