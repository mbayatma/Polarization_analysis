import numpy as np
import vector

# PDG IDs
ELECTRON_ID = 11
MUON_ID = 13

def build_vectors(tree):
    arrays = tree.arrays(["px", "py", "pz", "E", "pdgID"], library="np")
    return arrays

def get_particles(arrays, pdg_id):
    """Returns list of four-vectors for particles with given PDG ID"""
    px = arrays["px"]
    py = arrays["py"]
    pz = arrays["pz"]
    E  = arrays["E"]
    pid = arrays["pdgID"]

    selected = []
    for i in range(len(px)):
        vecs = []
        for j in range(len(px[i])):
            if abs(pid[i][j]) == pdg_id:
                vecs.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
        selected.append(vecs)
    return selected

# Leading pt / eta / phi

def compute_leading_pt_electron(tree):
    arrays = build_vectors(tree)
    electrons = get_particles(arrays, ELECTRON_ID)
    pts = [max([v.pt for v in e]) if e else np.nan for e in electrons]
    return np.array([p for p in pts if not np.isnan(p)])

def compute_leading_pt_muon(tree):
    arrays = build_vectors(tree)
    muons = get_particles(arrays, MUON_ID)
    pts = [max([v.pt for v in m]) if m else np.nan for m in muons]
    return np.array([p for p in pts if not np.isnan(p)])

def compute_leading_eta_electron(tree):
    arrays = build_vectors(tree)
    electrons = get_particles(arrays, ELECTRON_ID)
    etas = [e[0].eta if e else np.nan for e in electrons]
    return np.array([e for e in etas if not np.isnan(e)])

def compute_leading_eta_muon(tree):
    arrays = build_vectors(tree)
    muons = get_particles(arrays, MUON_ID)
    etas = [m[0].eta if m else np.nan for m in muons]
    return np.array([e for e in etas if not np.isnan(e)])

def compute_leading_phi_electron(tree):
    arrays = build_vectors(tree)
    electrons = get_particles(arrays, ELECTRON_ID)
    phis = [e[0].phi if e else np.nan for e in electrons]
    return np.array([p for p in phis if not np.isnan(p)])

def compute_leading_phi_muon(tree):
    arrays = build_vectors(tree)
    muons = get_particles(arrays, MUON_ID)
    phis = [m[0].phi if m else np.nan for m in muons]
    return np.array([p for p in phis if not np.isnan(p)])

# Invariant mass

def compute_mass_ee(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]
    py = arrays["py"]
    pz = arrays["pz"]
    E = arrays["E"]
    pdg = arrays["pdgID"]

    masses = []
    for i in range(len(pdg)):
        e_plus = []
        e_minus = []
        for j in range(len(pdg[i])):
            if pdg[i][j] == -11:  # e⁻
                e_minus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
            elif pdg[i][j] == 11:  # e⁺
                e_plus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
        if e_plus and e_minus:
            masses.append((e_plus[0] + e_minus[0]).mass)
    return np.array(masses)

def compute_mass_mumu(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]
    py = arrays["py"]
    pz = arrays["pz"]
    E = arrays["E"]
    pdg = arrays["pdgID"]

    masses = []
    for i in range(len(pdg)):
        mu_plus = []
        mu_minus = []
        for j in range(len(pdg[i])):
            if pdg[i][j] == -13:  # μ⁻
                mu_minus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
            elif pdg[i][j] == 13:  # μ⁺
                mu_plus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
        if mu_plus and mu_minus:
            masses.append((mu_plus[0] + mu_minus[0]).mass)
    return np.array(masses)

def compute_mass_epmup(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]
    py = arrays["py"]
    pz = arrays["pz"]
    E = arrays["E"]
    pdg = arrays["pdgID"]

    masses = []
    for i in range(len(pdg)):
        e_plus = []
        mu_plus = []
        for j in range(len(pdg[i])):
            if pdg[i][j] == 11:  # e⁺
                e_plus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
            elif pdg[i][j] == 13:  # μ⁺
                mu_plus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
        if e_plus and mu_plus:
            masses.append((e_plus[0] + mu_plus[0]).mass)
    return np.array(masses)
    
# ΔR and Δφ

def compute_delta_r(vecs1, vecs2):
    drs = []
    for v1, v2 in zip(vecs1, vecs2):
        if v1 and v2:
            drs.append(v1[0].deltaR(v2[0]))
    return np.array(drs)

def compute_delta_phi(vecs1, vecs2):
    dphis = []
    for v1, v2 in zip(vecs1, vecs2):
        if v1 and v2:
            delta_phi = v1[0].phi - v2[0].phi
            delta_phi = np.arctan2(np.sin(delta_phi), np.cos(delta_phi))  # wrap to [-π, π]
            dphis.append(np.degrees(np.abs(delta_phi)))
    return np.array(dphis)

def compute_dr_mumu(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]
    py = arrays["py"]
    pz = arrays["pz"]
    E = arrays["E"]
    pdg = arrays["pdgID"]

    dr_vals = []
    for i in range(len(pdg)):
        mu_plus = []
        mu_minus = []
        for j in range(len(pdg[i])):
            v = vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j])
            if pdg[i][j] == -13:  # μ⁻
                mu_minus.append(v)
            elif pdg[i][j] == 13:  # μ⁺
                mu_plus.append(v)

        if mu_plus and mu_minus:
            dr_vals.append(mu_plus[0].deltaR(mu_minus[0]))

    return np.array(dr_vals)
    
def compute_dr_epmup(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]
    py = arrays["py"]
    pz = arrays["pz"]
    E = arrays["E"]
    pdg = arrays["pdgID"]

    dr_vals = []
    for i in range(len(pdg)):
        e_plus = []
        mu_plus = []
        for j in range(len(pdg[i])):
            v = vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j])
            if pdg[i][j] == 11:  # e⁺
                e_plus.append(v)
            elif pdg[i][j] == 13:  # μ⁺
                mu_plus.append(v)

        if e_plus and mu_plus:
            dr_vals.append(e_plus[0].deltaR(mu_plus[0]))

    return np.array(dr_vals)

def compute_dphi_mumu(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]
    py = arrays["py"]
    pz = arrays["pz"]
    E = arrays["E"]
    pdg = arrays["pdgID"]

    dphi_list = []
    for i in range(len(pdg)):
        mu_plus = []
        mu_minus = []
        for j in range(len(pdg[i])):
            v = vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j])
            if pdg[i][j] == -13:  # μ⁻
                mu_minus.append(v)
            elif pdg[i][j] == 13:  # μ⁺
                mu_plus.append(v)

        if mu_plus and mu_minus:
            delta_phi = mu_plus[0].phi - mu_minus[0].phi
            delta_phi = np.arctan2(np.sin(delta_phi), np.cos(delta_phi))  # wrap to [-π, π]
            dphi_list.append(np.degrees(np.abs(delta_phi)))

    return np.array(dphi_list)


def compute_dphi_epmup(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]
    py = arrays["py"]
    pz = arrays["pz"]
    E = arrays["E"]
    pdg = arrays["pdgID"]

    dphi_list = []
    for i in range(len(pdg)):
        ep = []
        mup = []
        for j in range(len(pdg[i])):
            v = vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j])
            if pdg[i][j] == 11:  # e⁺
                ep.append(v)
            elif pdg[i][j] == 13:  # μ⁺
                mup.append(v)

        if ep and mup:
            delta_phi = ep[0].phi - mup[0].phi
            delta_phi = np.arctan2(np.sin(delta_phi), np.cos(delta_phi))
            dphi_list.append(np.degrees(np.abs(delta_phi)))

    return np.array(dphi_list)

# Rapidity and Δy

def compute_rapidity(vecs):
    y_vals = []
    for v in vecs:
        if v:
            y_vals.append(v[0].rapidity)
    return np.array(y_vals)

def compute_yep(tree):
    arrays = build_vectors(tree)
    electrons = get_particles(arrays, ELECTRON_ID)
    return compute_rapidity(electrons)

def compute_ymm(tree):
    arrays = build_vectors(tree)
    muons = get_particles(arrays, MUON_ID)
    return compute_rapidity(muons)

def compute_dyuu(tree):
    arrays = build_vectors(tree)
    muons = get_particles(arrays, MUON_ID)
    dy_vals = []
    for v in muons:
        if len(v) >= 2:
            dy_vals.append(np.abs(v[0].rapidity - v[1].rapidity))
    return np.array(dy_vals)

def compute_dypp(tree):
    arrays = build_vectors(tree)
    electrons = get_particles(arrays, ELECTRON_ID)
    muons = get_particles(arrays, MUON_ID)
    dy_vals = []
    for e, m in zip(electrons, muons):
        if e and m:
            dy_vals.append(np.abs(e[0].rapidity - m[0].rapidity))
    return np.array(dy_vals)
#compute MET
NEUTRINO_IDS = [12, -12, 14, -14, 16, -16]

def compute_met(tree):
    arrays = build_vectors(tree)
    pdg_ids = arrays["pdgID"]
    pxs = arrays["px"]
    pys = arrays["py"]

    met = []
    for ids, px, py in zip(pdg_ids, pxs, pys):
        met_px = sum(px[i] for i, pid in enumerate(ids) if pid in NEUTRINO_IDS)
        met_py = sum(py[i] for i, pid in enumerate(ids) if pid in NEUTRINO_IDS)
        met_val = np.sqrt(met_px**2 + met_py**2)
        met.append(met_val)

    return np.array(met)

def compute_met_phi(tree):
    arrays = build_vectors(tree)
    pdg_ids = arrays["pdgID"]
    pxs = arrays["px"]
    pys = arrays["py"]

    met_phi = []
    for ids, px, py in zip(pdg_ids, pxs, pys):
        met_px = sum(px[i] for i, pid in enumerate(ids) if pid in NEUTRINO_IDS)
        met_py = sum(py[i] for i, pid in enumerate(ids) if pid in NEUTRINO_IDS)
        if met_px == 0 and met_py == 0:
            continue  # Skip zero-MET events
        phi = np.arctan2(met_py, met_px)
        met_phi.append(phi)

    return np.array(met_phi)



    
dispatch_map = {
    "M_ee": compute_mass_ee,
    "M_mumu": compute_mass_mumu,
    "M_epmup": compute_mass_epmup,
    "dphi_mumu": compute_dphi_mumu,
    "dphi_epmup": compute_dphi_epmup,
    "DR_mupmum": compute_dr_mumu,
    "DR_epmup": compute_dr_epmup,
    "yep": compute_yep,
    "ymm": compute_ymm,
    "dyuu": compute_dyuu,
    "dypp": compute_dypp,
    "pt_leading_electron": compute_leading_pt_electron,
    "pt_leading_muon": compute_leading_pt_muon,
    "eta_leading_electron": compute_leading_eta_electron,
    "eta_leading_muon": compute_leading_eta_muon,
    "phi_leading_electron": compute_leading_phi_electron,
    "phi_leading_muon": compute_leading_phi_muon,
    "MET": compute_met,
    "MET_phi": compute_met_phi,
}


label_map = {
    "M_ee": r"$M_{e^+e^-}$ [GeV]",
    "M_mumu": r"$M_{\mu^+\mu^-}$ [GeV]",
    "M_epmup": r"$M_{e^+\mu^+}$ [GeV]",
    "dphi_mumu": r"$\Delta\phi(\mu^+, \mu^-)$ [deg]",
    "dphi_epmup": r"$\Delta\phi(e^+, \mu^+)$ [deg]",
    "DR_mupmum": r"$\Delta R(\mu^+, \mu^-)$",
    "DR_epmup": r"$\Delta R(e^+, \mu^+)$",
    "yep": r"$y_{e^+}$",
    "ymm": r"$y_{\mu^-}$",
    "dyuu": r"$|\Delta y(\mu^+, \mu^-)|$",
    "dypp": r"$|\Delta y(e^+, \mu^+)|$",
    "pt_leading_electron": r"$p_T(e^{+})$ [GeV]",
    "pt_leading_muon": r"$p_T(\mu^{+})$ [GeV]",
    "eta_leading_electron": r"$\eta(e^{+})$",
    "eta_leading_muon": r"$\eta(\mu^{+})$",
    "phi_leading_electron": r"$\phi(e^{+})$ [rad]",
    "phi_leading_muon": r"$\phi(\mu^{+})$ [rad]",
    "MET": r"MET [GeV]",
    "MET_phi": r"$\phi(\text{MET})$ [rad]",
}





