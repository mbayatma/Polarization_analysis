import csv
import os

# OLD SHOWERED FILES (kept for compatibility)
base_path = "/eos/user/m/mbayatma/polarization_files/"

# NEW FILES: SM and EFT (LHE → ROOT)
base_path_SM  = "/eos/user/m/mbayatma/polarization_files/SMRootFiles"
base_path_EFT = "/eos/user/m/mbayatma/polarization_files/BSMRootFiles"
base_path_BSM= base_path_EFT

# Regenerated EFT files (2026-09-11). 23 of the 44 original EFT files were
# byte-identical duplicates of each other and were regenerated; plotter/io.py
# prefers this directory and falls back to base_path_EFT for the rest.
base_path_EFT_REGEN = "/eos/user/m/mbayatma/polarization_files/ROOT_OUTPUT_regenerated_20260911"

# Pythia8-showered + Delphes samples, 0j and 1j merged (ATLAS slim card).
# A complete, self-consistent set: 4 polarizations x (SM + 12 operators).
# Flat reco ntuples ("<sample>_01j_flat.root", tree "events") are what the
# reco-level analysis reads; the "_01j_delphes.root" files next to them are
# the full Delphes output and are not needed downstream.
base_path_RECO = "/eos/user/m/mbayatma/polarization_files/WZ_01j_merged_ROOT"
reco_file_template = "WZ_3L_{pol}_{sample}_01j_flat.root"

polarizations = [
    r"$W_L Z_L$",
    r"$W_T Z_T$",
    r"$W_T Z_L$",
    r"$W_L Z_T$"
]

# mapping between fancy label → short name used in file names
pol_to_short = {
    r"$W_L Z_L$": "WL_ZL",
    r"$W_T Z_T$": "WT_ZT",
    r"$W_T Z_L$": "WT_ZL",
    r"$W_L Z_T$": "WL_ZT"
}

files_SM = {
    r"$W_L Z_L$": "Polarization_WZ_3L_WL_ZL_SM.root",
    r"$W_T Z_T$": "Polarization_WZ_3L_WT_ZT_SM.root",
    r"$W_T Z_L$": "Polarization_WZ_3L_WT_ZL_SM.root",
    r"$W_L Z_T$": "Polarization_WZ_3L_WL_ZT_SM.root"
}

files_SM_spin = {
    r"$W_L Z_L$": "WZ_3L_WL_ZL_SM.root",
    r"$W_T Z_T$": "WZ_3L_WT_ZT_SM.root",
    r"$W_T Z_L$": "WZ_3L_WT_ZL_SM.root",
    r"$W_L Z_T$": "WZ_3L_WL_ZT_SM.root",
}

xsecs_SM = {
    r"$W_L Z_L$": 0.0222885,
    r"$W_T Z_T$": 0.274366,
    r"$W_T Z_L$": 0.0312744,
    r"$W_L Z_T$": 0.03439
}



eft_operators = [
    "cW", "cHW", "cHWB", "cHbox", "cHDD",
    "cHq3", "cHq1", "cHu", "cHd", "cHl3",
    "cll1", "cWtil"
]

# ROOT file naming convention:
#   WZ_3L_<POL>_<OP>.root
#
# Example:
#   WZ_3L_WT_ZT_cHWB.root
#   WZ_3L_WL_ZL_cHu.root

files_EFT_all = {
    op: {
        pol: f"WZ_3L_{pol_to_short[pol]}_{op}.root"
        for pol in polarizations
    }
    for op in eft_operators
}

def _load_xsecs_eft_from_csv(csv_path):
    """Load per-operator, per-polarization EFT cross sections from xsections_all.csv.

    The CSV is the single source of truth for these numbers (previously they
    were hand-copied into a hardcoded dict here, operator by operator).
    """
    short_to_pol = {short: pol for pol, short in pol_to_short.items()}
    xsecs = {}
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            pol = short_to_pol.get(row["pol"])
            if pol is None:
                continue
            xsecs.setdefault(row["operator"], {})[pol] = float(row["xsec_pb"])
    return xsecs


_XSECTIONS_CSV = os.path.join(os.path.dirname(__file__), "..", "xsections_all.csv")

xsecs_EFT_all = {
    # cW isn't in xsections_all.csv (different production sample), so it
    # stays hardcoded here.
    "cW": {
        r"$W_L Z_L$": 0.0211174627,
        r"$W_T Z_T$": 0.2588394876,
        r"$W_T Z_L$": 0.0295997177,
        r"$W_L Z_T$": 0.0324668316
    },
}
xsecs_EFT_all.update(_load_xsecs_eft_from_csv(_XSECTIONS_CSV))

colors = {
    r"$W_L Z_L$": "red",
    r"$W_T Z_T$": "blue",
    r"$W_T Z_L$": "green",
    r"$W_L Z_T$": "goldenrod",
}


# ============================================
#   VARIABLE DEFINITIONS
# ============================================

mass_variables = {
    "M_ee": {
        "branches": ("el_pt", "el_eta", "el_phi", "el_charge"),
        "mass": 0.000511,
        "charges": (-1, +1),
        "label": r"$M_{e^+e^-}$"
    },
    "M_mumu": {
        "branches": ("mu_pt", "mu_eta", "mu_phi", "mu_charge"),
        "mass": 0.10566,
        "charges": (-1, +1),
        "label": r"$M_{\mu^+\mu^-}$"
    },
    "M_epmup": {
        "branches": (
            ("el_pt", "el_eta", "el_phi", "el_charge", 0.000511, +1),
            ("mu_pt", "mu_eta", "mu_phi", "mu_charge", 0.10566, +1)
        ),
        "label": r"$M_{e^+\mu^+}$"
    }
}

delta_phi_variables = {
    "dphi_mumu": {
        "function": "delta_phi_mumu",
        "xlabel": r"$\Delta\phi(\mu^+, \mu^-)$ [degrees]"
    },
    "dphi_epmup": {
        "function": "delta_phi_epmup",
        "xlabel": r"$\Delta\phi(e^+, \mu^+)$ [degrees]"
    }
}

rapidity_variables = ["yep", "ymm", "yuu", "dyuu", "dypp", "dyez"]

