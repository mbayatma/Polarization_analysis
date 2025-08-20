import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.signal import savgol_filter
from .config import polarizations, xsecs_SM, xsecs_EFT, colors


def compute_histogram(data, xsec, n_total, bins):
    bin_width = bins[1] - bins[0]
    weight = (xsec / n_total) / bin_width
    weights = np.full_like(data, weight)
    hist, _ = np.histogram(data, bins=bins, weights=weights)
    return hist


def compute_ratio(hist_eft, hist_sm):
    ratio = np.divide(hist_eft, hist_sm, out=np.ones_like(hist_sm, dtype=float), where=hist_sm != 0)
    err_sm = np.sqrt(hist_sm)
    err_eft = np.sqrt(hist_eft)
    ratio_unc = np.zeros_like(ratio)
    nonzero = hist_sm > 0
    ratio_unc[nonzero] = ratio[nonzero] * np.sqrt(
        (err_eft[nonzero] / hist_eft[nonzero])**2 +
        (err_sm[nonzero] / hist_sm[nonzero])**2
    )
    return ratio, ratio_unc


def smooth_histogram(hist, window=7, polyorder=2):
    if len(hist) < window:
        return hist
    return savgol_filter(hist, window_length=window, polyorder=polyorder)


def plot_comparison(
    variable_name,
    hist_sm_dict,
    hist_eft_dict,
    ratio_dict,
    bins,
    xlabel,
    output_path=None,
    ratio_uncertainties=None
):
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    ratio_keys = [pol for pol in polarizations if pol in hist_sm_dict and pol in ratio_dict]
    n_ratio_panels = len(ratio_keys)

    fig = plt.figure(figsize=(10, 5 + 1.2 * (1 + n_ratio_panels)))
    gs = GridSpec(n_ratio_panels + 2, 1, height_ratios=[5] + [1] * (n_ratio_panels + 1), hspace=0.2)

    # --- Top: SM vs EFT ---
    ax0 = fig.add_subplot(gs[0])
    for pol in hist_sm_dict:
        ax0.step(bin_centers, hist_sm_dict[pol], where='mid', linestyle='--',
                 color=colors.get(pol, 'black'), label=pol + " (SM)")
        if pol in hist_eft_dict:
            ax0.step(bin_centers, hist_eft_dict[pol], where='mid', linestyle='-',
                     color=colors.get(pol, 'black'), label=pol + " (EFT)")
    ax0.set_ylabel(r"$\frac{d\sigma}{dX}$ [pb/bin]")
    ax0.set_yscale("log")
    ax0.set_title(f"{variable_name} Distribution (Cross-Section Normalized)")
    ax0.grid(True)
    ax0.legend(ncol=2)
    ax0.tick_params(labelbottom=False)

    # --- Middle: Combined Ratio Panel ---
    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    for pol in ratio_keys:
        #valid = hist_sm_dict[pol] > 0
        ratio_vals, _ = ratio_dict[pol]
        valid = (hist_sm_dict[pol] > 0) & (hist_eft_dict[pol] > 0)
        #smoothed = smooth_histogram(ratio_vals[valid])
        ax1.plot(bin_centers[valid], ratio_vals[valid] ,color=colors.get(pol, 'black'), label=pol)
    ax1.axhline(1.0, color='black', linestyle='--', linewidth=1.0)
    ax1.set_ylabel("EFT / SM")
    ax1.set_ylim(0.0, 4)
    ax1.grid(True)
    ax1.tick_params(labelbottom=False)
    ax1.legend(ncol=2)

    # --- Lower: Individual Ratio Panels ---
    for i, pol in enumerate(ratio_keys):
        ax = fig.add_subplot(gs[i + 2], sharex=ax0)
        ratio_vals, _ = ratio_dict[pol]
        valid = (hist_sm_dict[pol] > 0) & (hist_eft_dict[pol] > 0)
        #valid = hist_sm_dict[pol] > 0
        #ratio_vals, ratio_unc = ratio_dict[pol]
        #smoothed = smooth_histogram(ratio_vals[valid])
        ax.plot(bin_centers[valid], ratio_vals[valid] ,color=colors.get(pol, 'black'))
        ax.axhline(1.0, color='black', linestyle='--', linewidth=1.0)
        ax.set_ylabel(pol)
        ax.set_ylim(0.0, 2.5)
        ax.grid(True)
        if i == len(ratio_keys) - 1:
            ax.set_xlabel(xlabel)
        else:
            ax.tick_params(labelbottom=False)

    fig.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.show()

