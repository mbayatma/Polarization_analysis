import numpy as np
import vector



def build_vectors(event):
    """Construct Lorentz vectors from px, py, pz, E"""
    return vector.array({
        "px": np.array(event["px"]),
        "py": np.array(event["py"]),
        "pz": np.array(event["pz"]),
        "E":  np.array(event["E"])
    })

def compute_leading_pt(tree):
    arrays = tree.arrays(["px", "py", "pdgID"], library="np")
    electron_pts = []
    muon_pts = []

    for i in range(len(arrays["px"])):
        px = arrays["px"][i]
        py = arrays["py"][i]
        pdg = arrays["pdgID"][i]

        pt = np.hypot(px, py)

        e_pts = pt[(np.abs(pdg) == 11)]
        mu_pts = pt[(np.abs(pdg) == 13)]

        if len(e_pts) > 0:
            electron_pts.append(np.max(e_pts))
        if len(mu_pts) > 0:
            muon_pts.append(np.max(mu_pts))

    return np.array(electron_pts), np.array(muon_pts)

def compute_leading_eta(tree):
    arrays = tree.arrays(["px", "py", "pz", "E", "pdgID"], library="np")
    electron_etas = []
    muon_etas = []

    for i in range(len(arrays["px"])):
        px = arrays["px"][i]
        py = arrays["py"][i]
        pz = arrays["pz"][i]
        E = arrays["E"][i]
        pdg = arrays["pdgID"][i]

        eta_list = []
        for j in range(len(px)):
            vec = vector.obj(px=px[j], py=py[j], pz=pz[j], E=E[j])
            eta_list.append((vec.eta, pdg[j]))

        eta_e = [eta for eta, pid in eta_list if abs(pid) == 11]
        eta_mu = [eta for eta, pid in eta_list if abs(pid) == 13]

        if eta_e:
            electron_etas.append(eta_e[0])  # pick leading (first) for now
        if eta_mu:
            muon_etas.append(eta_mu[0])

    return np.array(electron_etas), np.array(muon_etas)

def compute_pt(tree):
    arrays = tree.arrays(["px", "py"], library="np")
    pts = []
    for i in range(len(arrays["px"])):
        px = np.array(arrays["px"][i])
        py = np.array(arrays["py"][i])
        pts.extend(np.hypot(px, py))
    return np.array(pts)

def compute_eta(tree):
    arrays = tree.arrays(["px", "py", "pz", "E"], library="np")
    etas = []
    for i in range(len(arrays["px"])):
        px = arrays["px"][i]
        py = arrays["py"][i]
        pz = arrays["pz"][i]
        E  = arrays["E"][i]
        for j in range(len(px)):
            vec = vector.obj(px=px[j], py=py[j], pz=pz[j], E=E[j])
            etas.append(vec.eta)
    return np.array(etas)

def compute_leading_phi(tree):
    arrays = tree.arrays(["px", "py", "pz", "E", "pdgID"], library="np")
    electron_phis = []
    muon_phis = []

    for i in range(len(arrays["px"])):
        px = arrays["px"][i]
        py = arrays["py"][i]
        pz = arrays["pz"][i]
        E = arrays["E"][i]
        pdg = arrays["pdgID"][i]

        phi_list = []
        for j in range(len(px)):
            vec = vector.obj(px=px[j], py=py[j], pz=pz[j], E=E[j])
            phi_list.append((vec.phi, pdg[j]))

        phi_e = [phi for phi, pid in phi_list if abs(pid) == 11]
        phi_mu = [phi for phi, pid in phi_list if abs(pid) == 13]

        if phi_e:
            electron_phis.append(phi_e[0])  # first occurrence for now
        if phi_mu:
            muon_phis.append(phi_mu[0])

    return np.array(electron_phis), np.array(muon_phis)

def compute_phi(tree):
    arrays = tree.arrays(["px", "py", "pz", "E"], library="np")
    phis = []
    for i in range(len(arrays["px"])):
        px = arrays["px"][i]
        py = arrays["py"][i]
        pz = arrays["pz"][i]
        E  = arrays["E"][i]
        for j in range(len(px)):
            vec = vector.obj(px=px[j], py=py[j], pz=pz[j], E=E[j])
            phis.append(np.degrees(vec.phi))
    return np.array(phis)

def compute_mass(tree):
    arrays = tree.arrays(["px", "py", "pz", "E"], library="np")
    masses = []
    for i in range(len(arrays["px"])):
        px = arrays["px"][i]
        py = arrays["py"][i]
        pz = arrays["pz"][i]
        E  = arrays["E"][i]
        if len(px) >= 2:
            vec1 = vector.obj(px=px[0], py=py[0], pz=pz[0], E=E[0])
            vec2 = vector.obj(px=px[1], py=py[1], pz=pz[1], E=E[1])
            masses.append((vec1 + vec2).mass)
    return np.array(masses)

