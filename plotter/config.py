# OLD SHOWERED FILES (kept for compatibility)
base_path = "/eos/user/m/mbayatma/polarization_files/"

# NEW FILES: SM and EFT (LHE → ROOT)
base_path_SM  = "/eos/user/m/mbayatma/polarization_files/SMRootFiles"
base_path_EFT = "/eos/user/m/mbayatma/polarization_files/BSMRootFiles" 
base_path_BSM= base_path_EFT

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

xsecs_EFT_all = {

    # ---- existing operator cW ----
    "cW": {
        r"$W_L Z_L$": 0.0211174627,
        r"$W_T Z_T$": 0.2588394876,
        r"$W_T Z_L$": 0.0295997177,
        r"$W_L Z_T$": 0.0324668316
    },

    # ---- operator cHW (values you gave manually) ----
    "cHW": {
        r"$W_L Z_L$": 0.0211,
        r"$W_T Z_L$": 0.02957,
        r"$W_L Z_T$": 0.03251,
        r"$W_T Z_T$": 0.2592
    },
    

    "cHWB": {

        r"$W_T Z_T$": 0.26749807343959997,
        r"$W_L Z_L$": 0.021774346896949997,
        r"$W_T Z_L$": 0.03051100526077,
        r"$W_L Z_T$": 0.03354925146509,
     },

    "cHbox": {
        r"$W_T Z_T$": 0.25925791952620003,
        r"$W_L Z_L$": 0.021104686925369997,
        r"$W_T Z_L$": 0.029573201498610004,
        r"$W_L Z_T$": 0.03251220758310001,
    },

    "cHDD": {
        r"$W_T Z_T$": 0.2588608870859,
        r"$W_L Z_L$": 0.02107693044458999,
        r"$W_T Z_L$": 0.02953199588852,
        r"$W_L Z_T$": 0.03247489722058,
    },

    "cHq3": {
        r"$W_T Z_T$": 0.25925791952620003,
        r"$W_L Z_L$": 0.021104686925369997,
        r"$W_T Z_L$": 0.029573201498610004,
        r"$W_L Z_T$": 0.03251220758310001,
    },

    "cHq1": {
        r"$W_T Z_T$": 0.25925791952620003,
        r"$W_L Z_L$": 0.021104686925369997,
        r"$W_T Z_L$": 0.029573201498610004,
        r"$W_L Z_T$": 0.03251220758310001,
    },

    "cHu": {
        r"$W_T Z_T$": 0.25925791952620003,
        r"$W_L Z_L$": 0.021104686925369997,
        r"$W_T Z_L$": 0.029573201498610004,
        r"$W_L Z_T$": 0.03251220758310001,
    },

    "cHd": {
        r"$W_T Z_T$": 0.25925791952620003,
        r"$W_L Z_L$": 0.021104686925369997,
        r"$W_T Z_L$": 0.029573201498610004,
        r"$W_L Z_T$": 0.03251220758310001,
    },

    "cHl3": {
        r"$W_T Z_T$": 0.2629925953954001,
        r"$W_L Z_L$": 0.02141031754726999,
        r"$W_T Z_L$": 0.030001542141300006,
        r"$W_L Z_T$": 0.03300296391238,
    },

    "cll1": {
        r"$W_T Z_T$": 0.2752154066874,
        r"$W_L Z_L$": 0.02240371235962001,
        r"$W_T Z_L$": 0.03139321441688999,
        r"$W_L Z_T$": 0.03455161327375,
    },

    "cWtil": {
        r"$W_T Z_T$": 0.25925791952620003,
        r"$W_L Z_L$": 0.021104686925369997,
        r"$W_T Z_L$": 0.029573201498610004,
        r"$W_L Z_T$": 0.03251220758310001,
    },

}   
    # ---- others: auto-fill after reading xsections_all.csv ----
    # (you can add them manually or load them dynamically in main_new.py)
 

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

