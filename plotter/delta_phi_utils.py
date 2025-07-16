import numpy as np

def compute_delta_phi_deg(phi1, phi2):
    dphi = np.abs(phi1 - phi2)
    dphi = np.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
    return np.degrees(dphi)

def delta_phi_mumu(tree):
    mu_phi = tree["mu_phi"].array(library="np")
    mu_charge = tree["mu_charge"].array(library="np")
    dphi_list = []

    for phi, charge in zip(mu_phi, mu_charge):
        pos = [i for i, c in enumerate(charge) if c == +1]
        neg = [i for i, c in enumerate(charge) if c == -1]
        if pos and neg:
            dphi = compute_delta_phi_deg(phi[pos[0]], phi[neg[0]])
            dphi_list.append(dphi)
    return np.array(dphi_list)

def delta_phi_epmup(tree):
    el_phi = tree["el_phi"].array(library="np")
    el_charge = tree["el_charge"].array(library="np")
    mu_phi = tree["mu_phi"].array(library="np")
    mu_charge = tree["mu_charge"].array(library="np")
    dphi_list = []

    for eph, ech, mph, mch in zip(el_phi, el_charge, mu_phi, mu_charge):
        ep = [i for i, c in enumerate(ech) if c == +1]
        mup = [i for i, c in enumerate(mch) if c == +1]
        if ep and mup:
            dphi = compute_delta_phi_deg(eph[ep[0]], mph[mup[0]])
            dphi_list.append(dphi)
    return np.array(dphi_list)

