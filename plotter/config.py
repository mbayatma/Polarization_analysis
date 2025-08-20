# plotter/config.py
#files with parton shower
base_path = "/eos/user/m/mbayatma/polarization_files/"

#new files without parton shower madspin and no madspin 

base_path_SM="/eos/user/m/mbayatma/polarization_files/SMRootFiles"
base_path_BSM="/eos/user/m/mbayatma/polarization_files/BSMRootFiles"

polarizations = [r"$W_L Z_L$", r"$W_T Z_T$", r"$W_T Z_L$", r"$W_L Z_T$", r"$W_T Z_T$ (Madspin)"]

files_SM_spin= { 
    r"$W_T Z_T$ (Madspin)": "MadSpin_WZ_3L_WT_ZT_SM.root", 
    r"$W_L Z_L$":  "WZ_3L_WL_ZL_SM.root" ,
    r"$W_T Z_T$": "WZ_3L_WT_ZT_SM.root",  
    r"$W_T Z_L$": "WZ_3L_WT_ZL_SM.root",  
    r"$W_L Z_T$": "WZ_3L_WL_ZT_SM.root"  
} 
files_BSM_spin= {
    r"$W_L Z_L$":  "WZ_3L_WL_ZL_NPcW.root" ,
    r"$W_T Z_T$": "WZ_3L_WT_ZT_NPcW.root",  
    r"$W_T Z_L$": "WZ_3L_WT_ZL_NPcW.root",  
    r"$W_L Z_T$": "WZ_3L_WL_ZT_NPcW.root"  
}

files_SM = {
    r"$W_L Z_L$": "Polarization_WZ_3L_WL_ZL_SM.root",
    r"$W_T Z_T$": "Polarization_WZ_3L_WT_ZT_SM.root",
    r"$W_T Z_L$": "Polarization_WZ_3L_WT_ZL_SM.root",
    r"$W_L Z_T$": "Polarization_WZ_3L_WL_ZT_SM.root"
}

files_EFT = {
    r"$W_L Z_L$": "Polarization_WZ_3L_WL_ZL_NPcW_pythia8_events.root",
    r"$W_T Z_T$": "Polarization_WZ_3L_WT_ZT_NPcW_pythia8_events.root",
    r"$W_T Z_L$": "Polarization_WZ_3L_WT_ZL_NPcW_pythia8_events.root",
    r"$W_L Z_T$": "Polarization_WZ_3L_WL_ZT_NPcW_pythia8_events.root"
}

xsecs_SM = {
    r"$W_L Z_L$": 0.0222885,
    r"$W_T Z_T$": 0.274366,
    r"$W_T Z_L$": 0.0312744,
    r"$W_L Z_T$": 0.0951278
}

xsecs_EFT = {
    r"$W_L Z_L$": 0.0211174627,
    r"$W_T Z_T$": 0.2588394876,
    r"$W_T Z_L$": 0.0295997177,
    r"$W_L Z_T$": 0.0324668316
}

colors = {
    r"$W_L Z_L$": "red",
    r"$W_T Z_T$": "blue",
    r"$W_T Z_L$": "green",
    r"$W_L Z_T$": "goldenrod",
    r"$W_T Z_T$ (Madspin)": "purple" 
}


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


