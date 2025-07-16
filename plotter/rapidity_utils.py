import numpy as np
import vector

def compute_rapidity(tree, var_name):
    el_pt = tree["el_pt"].array(library="np")
    el_eta = tree["el_eta"].array(library="np")
    el_phi = tree["el_phi"].array(library="np")
    el_charge = tree["el_charge"].array(library="np")

    mu_pt = tree["mu_pt"].array(library="np")
    mu_eta = tree["mu_eta"].array(library="np")
    mu_phi = tree["mu_phi"].array(library="np")
    mu_charge = tree["mu_charge"].array(library="np")

    results = []

    for i in range(len(el_charge)):
        ep_idx = [j for j, c in enumerate(el_charge[i]) if c == +1]
        mm_idx = [j for j, c in enumerate(mu_charge[i]) if c == -1]
        mp_idx = [j for j, c in enumerate(mu_charge[i]) if c == +1]

        if not ep_idx or not mm_idx or not mp_idx:
            continue

        ep = vector.obj(pt=el_pt[i][ep_idx[0]], eta=el_eta[i][ep_idx[0]], phi=el_phi[i][ep_idx[0]], mass=0.000511)
        mm = vector.obj(pt=mu_pt[i][mm_idx[0]], eta=mu_eta[i][mm_idx[0]], phi=mu_phi[i][mm_idx[0]], mass=0.105)
        mp = vector.obj(pt=mu_pt[i][mp_idx[0]], eta=mu_eta[i][mp_idx[0]], phi=mu_phi[i][mp_idx[0]], mass=0.105)
        z = mp + mm

        if var_name == "yep":
            results.append(ep.rapidity)
        elif var_name == "ymm":
            results.append(mm.rapidity)
        elif var_name == "yuu":
            results.append(z.rapidity)
        elif var_name == "dyuu":
            results.append(np.abs(mp.rapidity - mm.rapidity))
        elif var_name == "dypp":
            results.append(np.abs(ep.rapidity - mp.rapidity))


    return np.array(results)
