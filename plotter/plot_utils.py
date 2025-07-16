# plotter/plot_utils.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from .config import polarizations, xsecs_SM, xsecs_EFT, colors

def compute_histogram(data, xsec, n_total, bins):
    bin_width = bins[1] - bins[0]
    weight = (xsec / n_total) / bin_width
    weights = np.full_like(data, weight)
    hist, _ = np.histogram(data, bins=bins, weights=weights)
    return hist

def compute_ratio(hist_eft, hist_sm):
    return np.divide(hist_eft, hist_sm, out=np.ones_like(hist_sm, dtype=float), where=hist_sm != 0)

def plot_comparison(
    variable_name,
    hist_sm_dict,
    hist_eft_dict,
    ratio_dict,
    bins,
    xlabel,
    output_path=None
):
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    n_ratio_panels = len(polarizations)
    fig = plt.figure(figsize=(10, 5 + 1.2 * (1 + n_ratio_panels)))
    gs = GridSpec(n_ratio_panels + 2, 1, height_ratios=[5] + [1]*(n_ratio_panels + 1), hspace=0.15)

    # --- Top: histogram comparison ---
    ax0 = fig.add_subplot(gs[0])
    for pol in polarizations:
        ax0.step(bin_centers, hist_sm_dict[pol], where='mid', linestyle='--',
                 color=colors[pol], label=pol + " (SM)")
        ax0.step(bin_centers, hist_eft_dict[pol], where='mid', linestyle='-',
                 color=colors[pol], label=pol + " (EFT)")

    ax0.set_ylabel(r"$\frac{d\sigma}{dX}$ [pb/bin]")
    ax0.set_yscale("log")
    ax0.set_title(f"{variable_name} Distribution (Cross-Section Normalized)")
    ax0.grid(True)
    ax0.legend(ncol=2)
    ax0.tick_params(labelbottom=False)

    # --- Middle: combined ratio plot ---
    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    for pol in polarizations:
        valid = hist_sm_dict[pol] > 0
        ax1.plot(bin_centers[valid], ratio_dict[pol][valid], color=colors[pol], label=pol)
    ax1.axhline(1.0, color='black', linestyle='--', linewidth=1.2)
    ax1.set_ylabel("EFT / SM")
    ax1.set_ylim(0.0, 4)
    ax1.grid(True)
    ax1.tick_params(labelbottom=False)
    ax1.legend(ncol=2)

    # --- Lower: individual ratio plots ---
    for i, pol in enumerate(polarizations):
        ax = fig.add_subplot(gs[i + 2], sharex=ax0)
        valid = hist_sm_dict[pol] > 0
        ax.plot(bin_centers[valid], ratio_dict[pol][valid], color=colors[pol])
        ax.axhline(1.0, color='black', linestyle='--', linewidth=1.0)
        ax.set_ylabel(pol)
        ax.set_ylim(0.0, 1.5)
        ax.grid(True)
        if i == len(polarizations) - 1:
            ax.set_xlabel(xlabel)
        else:
            ax.tick_params(labelbottom=False)

    fig.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.show()  # Disabled for batch use
