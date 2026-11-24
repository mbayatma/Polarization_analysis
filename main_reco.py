"""SM vs. EFT comparison plots for the Pythia8-showered + Delphes samples.

The LHE-level counterpart is main.py. This script reads the flat reco
ntuples (see plotter/reco.py), applies the WZ 3-lepton selection, and
produces figures in exactly the same style as main.py's: it calls
plotter.plot_utils.plot_comparison directly, so the layout (step histograms,
dashed SM / solid EFT, log-scale top panel, a combined ratio panel plus one
ratio panel per polarization) is pixel-for-pixel the same code path, not a
lookalike.

The one thing that has to differ is how a histogram is normalized. The
LHE-level compute_histogram() assumes one flat per-event weight (xsec /
n_total), because every LHE event carries the same weight. Reco events don't:
the shower/reco step leaves each event with its own `xsec_weight` (and a
small fraction are negative), so compute_histogram_weighted() below sums those
directly instead. Everything downstream -- compute_ratio, plot_comparison --
is reused unchanged from plotter/plot_utils.py.

Usage:
    python3 main_reco.py --all                        # every operator, every variable
    python3 main_reco.py --operator cHq1              # one operator, every variable
    python3 main_reco.py --operator cHq1 --var M_Z    # one figure
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import numpy as np

from plotter.config import (
    base_path_RECO,
    eft_operators,
    polarizations,
    pol_to_short,
    reco_file_template,
)
from plotter.plot_utils import compute_ratio, plot_comparison
from plotter import reco

OUTPUT_DIR = "polarization_plots_reco"


def reco_path(pol_tag, sample):
    return os.path.join(base_path_RECO, reco_file_template.format(pol=pol_tag, sample=sample))


def load_selected(pol_tag, sample):
    """Read one sample and return the post-selection observable dict."""
    path = reco_path(pol_tag, sample)
    if not os.path.exists(path):
        print(f"[WARN] missing: {path}")
        return None
    return reco.select_wz_3l(reco.load_arrays(path))


def compute_histogram_weighted(data, weights, bins):
    """Differential-cross-section histogram from per-event weights.

    The reco-level analog of plotter.plot_utils.compute_histogram: that
    function multiplies a flat (xsec / n_total) weight onto every LHE event,
    since all LHE events in a sample carry the same weight. Reco events each
    carry their own xsec_weight (pb) from the shower/reco step instead, so
    the histogram is filled with those directly. Dividing by bin width keeps
    the same pb/bin units compute_histogram produces, so the two feed
    plot_comparison identically.
    """
    finite = np.isfinite(data)
    bin_width = bins[1] - bins[0]
    hist, _ = np.histogram(data[finite], bins=bins, weights=weights[finite])
    return hist / bin_width


def plot_variable(var_name, data_SM, data_EFT, operator, output_dir):
    """One figure, built by the same plot_comparison() main.py uses."""
    xlabel, bins = reco.VAR_SPECS[var_name]

    hist_sm_dict, hist_eft_dict, ratio_dict = {}, {}, {}
    # Sorted order matches main.py's `for pol in sorted(all_keys)`, which is
    # what fixes the legend/color ordering seen in polariztion_plots_new.
    for pol in sorted(polarizations):
        sm, eft = data_SM.get(pol), data_EFT.get(pol)
        if sm is not None:
            hist_sm_dict[pol] = compute_histogram_weighted(sm[var_name], sm["xsec_weight"], bins)
        if eft is not None:
            hist_eft_dict[pol] = compute_histogram_weighted(eft[var_name], eft["xsec_weight"], bins)
        if pol in hist_sm_dict and pol in hist_eft_dict:
            ratio_dict[pol] = compute_ratio(hist_eft_dict[pol], hist_sm_dict[pol])

    if not hist_sm_dict:
        print(f"[WARN] nothing to draw for {var_name} / {operator}")
        return None

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{var_name}_comparison.png")
    plot_comparison(
        variable_name=var_name,
        hist_sm_dict=hist_sm_dict,
        hist_eft_dict=hist_eft_dict,
        ratio_dict=ratio_dict,
        bins=bins,
        xlabel=xlabel,
        output_path=out_path,
        operator=operator,
    )
    return out_path


def report_rates(data_SM, data_EFT, operator):
    """Print the selected cross section per polarization, SM vs. EFT.

    This is the rate-level cross-check that the ratio panels should be
    consistent with: summing xsec_weight over selected events gives the
    fiducial cross section directly.
    """
    print(f"\n[RATES] selected fiducial cross section (pb), SM vs. {operator}")
    print(f"  {'polarization':<14} {'SM':>12} {operator:>12} {'ratio':>8}")
    for pol in polarizations:
        sm, eft = data_SM.get(pol), data_EFT.get(pol)
        if sm is None or eft is None:
            continue
        s_sm = sm["xsec_weight"].sum()
        s_eft = eft["xsec_weight"].sum()
        ratio = s_eft / s_sm if s_sm else np.nan
        print(f"  {pol_to_short[pol]:<14} {s_sm:12.6f} {s_eft:12.6f} {ratio:8.4f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--var", help="One observable (default: all)")
    parser.add_argument("--operator", help="One EFT operator (default: all)")
    parser.add_argument("--all", action="store_true", help="All operators and observables")
    parser.add_argument("--output-dir", default=OUTPUT_DIR)
    args = parser.parse_args()

    operators = [args.operator] if args.operator else list(eft_operators)
    variables = [args.var] if args.var else list(reco.RECO_VARIABLES)

    unknown = [v for v in variables if v not in reco.VAR_SPECS]
    if unknown:
        parser.error(f"unknown variable(s): {unknown}")

    # The SM samples are shared across every operator, so they are read once.
    print("Reading SM samples (shared across operators)")
    data_SM = {}
    for pol in polarizations:
        sel = load_selected(pol_to_short[pol], "SM")
        if sel is not None:
            data_SM[pol] = sel
            print(f"  {pol_to_short[pol]:<6} selected {len(sel['M_Z']):6d} events, "
                  f"sigma = {sel['xsec_weight'].sum():.6f} pb")

    for op in operators:
        print(f"\n=== {op} ===")
        data_EFT = {}
        for pol in polarizations:
            sel = load_selected(pol_to_short[pol], op)
            if sel is not None:
                data_EFT[pol] = sel
                print(f"  {pol_to_short[pol]:<6} selected {len(sel['M_Z']):6d} events, "
                      f"sigma = {sel['xsec_weight'].sum():.6f} pb")
        if not data_EFT:
            print(f"[WARN] no samples for {op}, skipping")
            continue

        report_rates(data_SM, data_EFT, op)

        out_dir = os.path.join(args.output_dir, op)
        written = 0
        for var_name in variables:
            if plot_variable(var_name, data_SM, data_EFT, op, out_dir):
                written += 1
        print(f"  wrote {written} figures to {out_dir}/")


if __name__ == "__main__":
    main()
