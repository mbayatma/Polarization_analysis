import numpy as np
import vector

# Define mass and charge logic per variable
mass_definitions = {
    "M_ee": {
        "pt_branch": "el_pt", "eta_branch": "el_eta", "phi_branch": "el_phi",
        "charge_branch": "el_charge", "mass": 0.000511, "q1": -1, "q2": +1
    },
    "M_mumu": {
        "pt_branch": "mu_pt", "eta_branch": "mu_eta", "phi_branch": "mu_phi",
        "charge_branch": "mu_charge", "mass": 0.10566, "q1": -1, "q2": +1
    },
    "M_epmup": {
        "e_pt": "el_pt", "e_eta": "el_eta", "e_phi": "el_phi", "e_charge": "el_charge",
        "mu_pt": "mu_pt", "mu_eta": "mu_eta", "mu_phi": "mu_phi", "mu_charge": "mu_charge"
    }
}


def compute_mass(tree, varname):
    mass_vals = []

    if varname == "M_epmup":
        # Special case: electron+ and muon+
        e_pt = tree["el_pt"].array(library="np")
        e_eta = tree["el_eta"].array(library="np")
        e_phi = tree["el_phi"].array(library="np")
        e_charge = tree["el_charge"].array(library="np")
        mu_pt = tree["mu_pt"].array(library="np")
        mu_eta = tree["mu_eta"].array(library="np")
        mu_phi = tree["mu_phi"].array(library="np")
        mu_charge = tree["mu_charge"].array(library="np")

        for i in range(len(e_charge)):
            try:
                # Find e+ and mu+
                e_idx = list(e_charge[i]).index(+1)
                mu_idx = list(mu_charge[i]).index(+1)

                v1 = vector.obj(pt=e_pt[i][e_idx], eta=e_eta[i][e_idx],
                                phi=e_phi[i][e_idx], mass=0.000511)
                v2 = vector.obj(pt=mu_pt[i][mu_idx], eta=mu_eta[i][mu_idx],
                                phi=mu_phi[i][mu_idx], mass=0.10566)

                mass_vals.append((v1 + v2).mass)
            except:
                continue

    elif varname in mass_definitions:
        cfg = mass_definitions[varname]
        pt = tree[cfg["pt_branch"]].array(library="np")
        eta = tree[cfg["eta_branch"]].array(library="np")
        phi = tree[cfg["phi_branch"]].array(library="np")
        charge = tree[cfg["charge_branch"]].array(library="np")

        for i in range(len(charge)):
            pos = [j for j in range(len(charge[i])) if charge[i][j] == cfg["q2"]]
            neg = [j for j in range(len(charge[i])) if charge[i][j] == cfg["q1"]]
            if not pos or not neg:
                continue
            try:
                i_pos, i_neg = pos[0], neg[0]
                v1 = vector.obj(pt=pt[i][i_pos], eta=eta[i][i_pos],
                                phi=phi[i][i_pos], mass=cfg["mass"])
                v2 = vector.obj(pt=pt[i][i_neg], eta=eta[i][i_neg],
                                phi=phi[i][i_neg], mass=cfg["mass"])
                mass_vals.append((v1 + v2).mass)
            except:
                continue
    else:
        raise ValueError(f"Unknown mass variable: {varname}")

    return np.array(mass_vals)

