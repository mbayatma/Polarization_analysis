import numpy as np

def delta_r(eta1, phi1, eta2, phi2):
    dphi = np.abs(phi1 - phi2)
    dphi = np.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
    return np.sqrt((eta1 - eta2)**2 + dphi**2)

def compute_delta_r(tree, var_name):
    if var_name == "DR_mupmum":
        mu_eta = tree["mu_eta"].array(library="np")
        mu_phi = tree["mu_phi"].array(library="np")
        mu_charge = tree["mu_charge"].array(library="np")
        dr_list = []

        for eta, phi, charge in zip(mu_eta, mu_phi, mu_charge):
            pos = [i for i, c in enumerate(charge) if c == +1]
            neg = [i for i, c in enumerate(charge) if c == -1]
            if pos and neg:
                dr = delta_r(eta[pos[0]], phi[pos[0]], eta[neg[0]], phi[neg[0]])
                dr_list.append(dr)
        return np.array(dr_list)

    elif var_name == "DR_epmup":
        el_eta = tree["el_eta"].array(library="np")
        el_phi = tree["el_phi"].array(library="np")
        el_charge = tree["el_charge"].array(library="np")

        mu_eta = tree["mu_eta"].array(library="np")
        mu_phi = tree["mu_phi"].array(library="np")
        mu_charge = tree["mu_charge"].array(library="np")

        dr_list = []

        for eeta, ephi, ech, meta, mphi, mch in zip(el_eta, el_phi, el_charge, mu_eta, mu_phi, mu_charge):
            ep = [i for i, c in enumerate(ech) if c == +1]
            mup = [i for i, c in enumerate(mch) if c == +1]
            if ep and mup:
                dr = delta_r(eeta[ep[0]], ephi[ep[0]], meta[mup[0]], mphi[mup[0]])
                dr_list.append(dr)
        return np.array(dr_list)

    else:
        raise ValueError(f"Unknown delta R variable: {var_name}")

