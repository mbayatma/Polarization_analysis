import numpy as np
import vector

# =========================
# Analysis / trigger config
# =========================
# Master switch: apply selection to every variable
APPLY_3L_TRIGGER = False

# "Exactly 3 charged leptons (e, μ)" selection
EXACT_N_LEPTONS = 3  # set to None to disable the multiplicity requirement

# Simple single-lepton trigger emulation (truth-level)
# ATLAS-like: e ~27 GeV, μ ~26 GeV (you can tune these)
PT_E_SINGLE = 27.0  # GeV
PT_MU_SINGLE = 26.0  # GeV

# Fiducial acceptance (very rough)
ETA_MAX = 2.5  # reject leptons outside |eta|<ETA_MAX for trigger logic


# =========================
# PDG IDs (sign convention)
#  11 = e−,  -11 = e+
#  13 = μ−,  -13 = μ+
# =========================

ELECTRON_ID = 11
MUON_ID = 13
NEUTRINO_IDS = [12, -12, 14, -14, 16, -16]


# -----------------------
# Core array extraction
# -----------------------
# build_vectors is called once per variable per tree, and every dispatch_map
# entry starts with it -- so plotting ~30 variables re-read and re-decompressed
# the same branches ~30 times per tree. The arrays depend only on the tree (and
# the trigger switch), so they are cached by file path + tree name. Caching the
# arrays rather than the derived variables keeps each compute_* function
# unchanged.
_VECTOR_CACHE = {}


def _tree_cache_key(tree):
    try:
        return (tree.file.file_path, tree.object_path, APPLY_3L_TRIGGER)
    except AttributeError:
        return None


def clear_vector_cache():
    """Drop cached event arrays (use when iterating over many samples)."""
    _VECTOR_CACHE.clear()


def build_vectors(tree):
    key = _tree_cache_key(tree)
    if key is not None and key in _VECTOR_CACHE:
        return _VECTOR_CACHE[key]

    arrays = tree.arrays(["px", "py", "pz", "E", "pdgID"], library="np")
    if APPLY_3L_TRIGGER:
       mask = _event_mask_3l_trigger(arrays)
       arrays = _apply_event_mask(arrays, mask)

    if key is not None:
        _VECTOR_CACHE[key] = arrays
    return arrays


def get_particles(arrays, pdg_id):
    """Returns list of four-vectors for particles with |pdgID| == pdg_id (per event)."""
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


# -----------------------
# Helpers for per-event logic
# -----------------------
def _split_leptons_with_charge(arrays):
    """Return per-event lists: e_plus, e_minus, mu_plus, mu_minus (as vector.objs)."""
    px, py, pz, E, pid = arrays["px"], arrays["py"], arrays["pz"], arrays["E"], arrays["pdgID"]
    e_plus, e_minus, mu_plus, mu_minus = [], [], [], []
    for ev in range(len(pid)):
        ep, em, mp, mm = [], [], [], []
        for j in range(len(pid[ev])):
            v = vector.obj(px=px[ev][j], py=py[ev][j], pz=pz[ev][j], E=E[ev][j])
            if pid[ev][j] == -11:     # e+
                ep.append(v)
            elif pid[ev][j] ==  11:   # e−
                em.append(v)
            elif pid[ev][j] == -13:   # μ+
                mp.append(v)
            elif pid[ev][j] ==  13:   # μ−
                mm.append(v)
        e_plus.append(ar_sort_pt(ep))
        e_minus.append(ar_sort_pt(em))
        mu_plus.append(ar_sort_pt(mp))
        mu_minus.append(ar_sort_pt(mm))
    return e_plus, e_minus, mu_plus, mu_minus


def _p3(v):
    return np.array([v.px, v.py, v.pz], dtype=float)

def _boost_p4(p4, beta):
    """
    Boost a four-vector p4=(E,px,py,pz) by velocity beta (3-vector).
    beta is the velocity of the *new frame* w.r.t. old frame.
    To go to the parent rest frame, use beta = parent_p / parent_E, and boost by -beta.
    """
    E, px, py, pz = p4
    b = np.array(beta, dtype=float)
    b2 = np.dot(b, b)
    if b2 >= 1.0 or b2 == 0.0:
        return np.array([E, px, py, pz], dtype=float)

    gamma = 1.0 / np.sqrt(1.0 - b2)
    p = np.array([px, py, pz], dtype=float)
    bp = np.dot(b, p)

    E_prime = gamma * (E - bp)
    factor = (gamma - 1.0) * bp / b2 - gamma * E
    p_prime = p + factor * b
    return np.array([E_prime, p_prime[0], p_prime[1], p_prime[2]], dtype=float)

def _to_rest_frame(p4_child, p4_parent):
    """
    Return child's four-vector in parent's rest frame.
    """
    Epar, pxpar, pypar, pzpar = p4_parent
    if Epar == 0:
        return np.array(p4_child, dtype=float)
    beta = np.array([pxpar, pypar, pzpar], dtype=float) / float(Epar)
    return _boost_p4(p4_child, -beta)

def _unit(v):
    n = np.linalg.norm(v)
    if n == 0:
        return None
    return v / n

def _cos_opening_angle_3d(v1, v2):
    p1 = _p3(v1); p2 = _p3(v2)
    n1 = np.linalg.norm(p1); n2 = np.linalg.norm(p2)
    if n1 == 0 or n2 == 0:
        return None
    return float(np.dot(p1, p2) / (n1 * n2))

def _sum_neutrinos(arrays):
    px, py, pz, E, pid = arrays["px"], arrays["py"], arrays["pz"], arrays["E"], arrays["pdgID"]
    out = []
    for ev in range(len(pid)):
        nu_px = 0.0; nu_py = 0.0; nu_pz = 0.0; nu_E = 0.0
        for j in range(len(pid[ev])):
            if pid[ev][j] in NEUTRINO_IDS:
                nu_px += px[ev][j]
                nu_py += py[ev][j]
                nu_pz += pz[ev][j]
                nu_E  += E[ev][j]
        out.append(vector.obj(px=nu_px, py=nu_py, pz=nu_pz, E=nu_E) if nu_E != 0 else None)
    return out

def _pick_best_z_pair(arrays, mZ=91.1876):
    """
    Pick best Z candidate among ee and mumu (closest to mZ).
    Returns per-event dict:
      {"flavor": "ee"/"mumu", "lp": v_plus, "lm": v_minus, "Z": vZ}
    or None if no SFOS pair.
    """
    ee_pairs = _pick_z_pair_vectors(arrays, lep_id=11, mZ=mZ)   # list of (e+,e-) or None
    mm_pairs = _pick_z_pair_vectors(arrays, lep_id=13, mZ=mZ)   # list of (mu+,mu-) or None

    out = []
    for pr_ee, pr_mm in zip(ee_pairs, mm_pairs):
        cand = []
        if pr_ee is not None:
            lp, lm = pr_ee
            Z = lp + lm
            cand.append(("ee", lp, lm, Z, abs(Z.mass - mZ)))
        if pr_mm is not None:
            lp, lm = pr_mm
            Z = lp + lm
            cand.append(("mumu", lp, lm, Z, abs(Z.mass - mZ)))

        if not cand:
            out.append(None)
        else:
            cand.sort(key=lambda x: x[4])
            fl, lp, lm, Z, _ = cand[0]
            out.append({"flavor": fl, "lp": lp, "lm": lm, "Z": Z})
    return out

def _pick_w_lepton(arrays, z_choice):
    """
    W-lepton = highest-pT charged lepton not used in the Z pair.
    Returns per-event vector.obj or None.
    """
    e_plus, e_minus, mu_plus, mu_minus = _split_leptons_with_charge(arrays)
    out = []
    for ev, zc in enumerate(z_choice):
        if zc is None:
            out.append(None); continue

        used = set()
        # identify used objects by (px,py,pz,E) tuple (good enough at truth level)
        def key(v): return (float(v.px), float(v.py), float(v.pz), float(v.E))
        used.add(key(zc["lp"]))
        used.add(key(zc["lm"]))

        cands = []
        for v in (e_plus[ev] + e_minus[ev] + mu_plus[ev] + mu_minus[ev]):
            if key(v) not in used:
                cands.append(v)

        out.append(max(cands, key=lambda v: v.pt) if cands else None)
    return out



def ar_sort_pt(vlist):
    """Sort a flat list of vector.objs by descending pT."""
    return sorted(vlist, key=lambda v: v.pt, reverse=True)


def _all_masses(vlistA, vlistB):
    """All pair masses between two event-level lists."""
    out = []
    for A, B in zip(vlistA, vlistB):
        masses = []
        for a in A:
            for b in B:
                masses.append((a + b).mass)
        out.append(masses)
    return out


def _pick_closest_to_target_per_event(candidates_per_event, target):
    """Pick the candidate closest to target (e.g., mZ) per event; None if empty."""
    picked = []
    for masses in candidates_per_event:
        if masses:
            arr = np.asarray(masses, dtype=float)
            picked.append(arr[np.argmin(np.abs(arr - target))])
        else:
            picked.append(None)
    return picked


def _pick_first_per_event(candidates_per_event):
    """Pick the first candidate per event; None if empty."""
    return [m[0] if m else None for m in candidates_per_event]


def _leading_vector_by_pt(arrays, pdg_id):
    """Pick the |pdg|==pdg_id object with max pT per event (or None)."""
    px, py, pz, E, pid = arrays["px"], arrays["py"], arrays["pz"], arrays["E"], arrays["pdgID"]
    lead = []
    for i in range(len(pid)):
        cand = []
        for j in range(len(pid[i])):
            if abs(pid[i][j]) == pdg_id:
                cand.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
        lead.append(max(cand, key=lambda v: v.pt) if cand else None)
    return lead

def _pt(px, py):
    return np.sqrt(px*px + py*py)

def _eta(px, py, pz, E):
    # safe eta using vector for consistency with the rest of the code
    v = vector.obj(px=px, py=py, pz=pz, E=E)
    return v.eta

def _event_mask_3l_trigger(arrays,
                           pt_e=PT_E_SINGLE,
                           pt_mu=PT_MU_SINGLE,
                           eta_max=ETA_MAX,
                           exact_n_leptons=EXACT_N_LEPTONS):
    """
    Returns a boolean list mask (per event):
      - exactly N charged leptons (e or mu), if exact_n_leptons is not None
      - AND passes single-lepton trigger emulation:
          any electron with pt>pt_e and |eta|<eta_max
          OR any muon     with pt>pt_mu and |eta|<eta_max
    """
    px, py, pz, E, pid = arrays["px"], arrays["py"], arrays["pz"], arrays["E"], arrays["pdgID"]
    mask = []
    for i in range(len(pid)):
        # count charged leptons
        n_lep = sum(1 for q in pid[i] if abs(q) in (11, 13))
        if (exact_n_leptons is not None) and (n_lep != exact_n_leptons):
            mask.append(False)
            continue

        # trigger OR
        fired = False
        for j in range(len(pid[i])):
            apid = abs(pid[i][j])
            if apid not in (11, 13):
                continue
            pt  = _pt(px[i][j], py[i][j])
            eta = abs(_eta(px[i][j], py[i][j], pz[i][j], E[i][j]))
            if apid == 11 and (pt >= pt_e) and (eta <= eta_max):
                fired = True
                break
            if apid == 13 and (pt >= pt_mu) and (eta <= eta_max):
                fired = True
                break

        mask.append(fired)
    return mask

def _apply_event_mask(arrays, mask):
    """Filter variable-length per-event arrays by a boolean mask."""
    keep_idx = [i for i, m in enumerate(mask) if m]
    out = {}
    for key in arrays.keys():
        out[key] = [arrays[key][i] for i in keep_idx]
    return out
# -----------------------
# Leading kinematics
# -----------------------
def compute_leading_pt_electron(tree):
    arrays = build_vectors(tree)
    electrons = get_particles(arrays, ELECTRON_ID)
    pts = [max([v.pt for v in e]) if e else np.nan for e in electrons]
    return np.array([p for p in pts if not np.isnan(p)], dtype=float)


def compute_leading_pt_muon(tree):
    arrays = build_vectors(tree)
    muons = get_particles(arrays, MUON_ID)
    pts = [max([v.pt for v in m]) if m else np.nan for m in muons]
    return np.array([p for p in pts if not np.isnan(p)], dtype=float)


def compute_leading_eta_electron(tree):
    arrays = build_vectors(tree)
    lead = _leading_vector_by_pt(arrays, ELECTRON_ID)
    return np.array([v.eta for v in lead if v is not None], dtype=float)


def compute_leading_eta_muon(tree):
    arrays = build_vectors(tree)
    lead = _leading_vector_by_pt(arrays, MUON_ID)
    return np.array([v.eta for v in lead if v is not None], dtype=float)


def compute_leading_phi_electron(tree):
    arrays = build_vectors(tree)
    lead = _leading_vector_by_pt(arrays, ELECTRON_ID)
    return np.array([v.phi for v in lead if v is not None], dtype=float)


def compute_leading_phi_muon(tree):
    arrays = build_vectors(tree)
    lead = _leading_vector_by_pt(arrays, MUON_ID)
    return np.array([v.phi for v in lead if v is not None], dtype=float)


# -----------------------
# Invariant masses (existing style)
# -----------------------
#def compute_mass_ee(tree):
#    arrays = build_vectors(tree)
#    px = arrays["px"]; py = arrays["py"]; pz = arrays["pz"]; E = arrays["E"]; pdg = arrays["pdgID"]

#    masses = []
#    for i in range(len(pdg)):
#        e_plus = []
#        e_minus = []
#        for j in range(len(pdg[i])):
#            if pdg[i][j] == -11:  # e+
#                e_plus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
#            elif pdg[i][j] == 11:  # e−
#                e_minus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
#        if e_plus and e_minus:
#            masses.append((e_plus[0] + e_minus[0]).mass)
#    return np.array(masses, dtype=float)

def compute_mass_ee_closestZ(tree, mZ=91.1876):
    arrays = build_vectors(tree)
    e_plus, e_minus, _, _ = _split_leptons_with_charge(arrays)  # uses the helper that sorts by pT
    # build all e+e- masses per event
    ee = _all_masses(e_plus, e_minus)
    # pick the one closest to mZ per event
    pick = _pick_closest_to_target_per_event(ee, mZ)
    return np.array([x for x in pick if x is not None], dtype=float)



def compute_mass_mumu(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]; py = arrays["py"]; pz = arrays["pz"]; E = arrays["E"]; pdg = arrays["pdgID"]

    masses = []
    for i in range(len(pdg)):
        mu_plus = []
        mu_minus = []
        for j in range(len(pdg[i])):
            if pdg[i][j] == -13:  # μ+
                mu_plus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
            elif pdg[i][j] == 13:  # μ−
                mu_minus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
        if mu_plus and mu_minus:
            masses.append((mu_plus[0] + mu_minus[0]).mass)
    return np.array(masses, dtype=float)


def compute_mass_epmup(tree):
    """Mass of e+ with μ+ (kept for backward compatibility with your naming)."""
    arrays = build_vectors(tree)
    px = arrays["px"]; py = arrays["py"]; pz = arrays["pz"]; E = arrays["E"]; pdg = arrays["pdgID"]

    masses = []
    for i in range(len(pdg)):
        e_plus = []
        mu_plus = []
        for j in range(len(pdg[i])):
            if pdg[i][j] == -11:  # e+
                e_plus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
            elif pdg[i][j] == -13:  # μ+
                mu_plus.append(vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j]))
        if e_plus and mu_plus:
            masses.append((e_plus[0] + mu_plus[0]).mass)
    return np.array(masses, dtype=float)


# -----------------------
# One-fill-per-event representative mass
# -----------------------
def compute_mass_event_representative(tree, mZ=91.1876):
    """
    Pick exactly one mass per event:
      1) If an SFOS (ee or μμ) pair exists, choose the pair mass closest to mZ.
      2) Else if a DFOS pair exists (e+μ− or μ+e−), take the first such mass.
      3) Else drop the event.
    Returns a 1D numpy array with one value per kept event.
    """
    arrays = build_vectors(tree)
    e_plus, e_minus, mu_plus, mu_minus = _split_leptons_with_charge(arrays)

    # SFOS candidates → pick closest to mZ per event
    ee = _all_masses(e_plus, e_minus)
    mm = _all_masses(mu_plus, mu_minus)
    pick_ee = _pick_closest_to_target_per_event(ee, mZ)
    pick_mm = _pick_closest_to_target_per_event(mm, mZ)
    sfos = [a if a is not None else b for a, b in zip(pick_ee, pick_mm)]

    # DFOS candidates (e+μ− and μ+e−) → use if SFOS missing
    em_os1 = _all_masses(e_plus, mu_minus)
    em_os2 = _all_masses(mu_plus, e_minus)
    pick_em1 = _pick_first_per_event(em_os1)
    pick_em2 = _pick_first_per_event(em_os2)

    rep = []
    for s, a, b in zip(sfos, pick_em1, pick_em2):
        if s is not None:
            rep.append(s)
        elif a is not None:
            rep.append(a)
        elif b is not None:
            rep.append(b)
        else:
            rep.append(None)

    return np.array([x for x in rep if x is not None], dtype=float)


# -----------------------
# ΔR and Δφ utilities
# -----------------------
def compute_delta_r(vecs1, vecs2):
    drs = []
    for v1, v2 in zip(vecs1, vecs2):
        if v1 and v2:
            drs.append(v1[0].deltaR(v2[0]))
    return np.array(drs, dtype=float)


def compute_delta_phi(vecs1, vecs2):
    dphis = []
    for v1, v2 in zip(vecs1, vecs2):
        if v1 and v2:
            delta_phi = v1[0].phi - v2[0].phi
            delta_phi = np.arctan2(np.sin(delta_phi), np.cos(delta_phi))  # wrap to [-π, π]
            dphis.append(np.degrees(np.abs(delta_phi)))
    return np.array(dphis, dtype=float)


def compute_dr_mumu(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]; py = arrays["py"]; pz = arrays["pz"]; E = arrays["E"]; pdg = arrays["pdgID"]

    dr_vals = []
    for i in range(len(pdg)):
        mu_plus = []
        mu_minus = []
        for j in range(len(pdg[i])):
            v = vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j])
            if pdg[i][j] == -13:  # μ+
                mu_plus.append(v)
            elif pdg[i][j] == 13:  # μ−
                mu_minus.append(v)

        if mu_plus and mu_minus:
            dr_vals.append(mu_plus[0].deltaR(mu_minus[0]))

    return np.array(dr_vals, dtype=float)


def compute_dr_epmup(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]; py = arrays["py"]; pz = arrays["pz"]; E = arrays["E"]; pdg = arrays["pdgID"]

    dr_vals = []
    for i in range(len(pdg)):
        e_plus = []
        mu_plus = []
        for j in range(len(pdg[i])):
            v = vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j])
            if pdg[i][j] == -11:  # e+
                e_plus.append(v)
            elif pdg[i][j] == -13:  # μ+
                mu_plus.append(v)

        if e_plus and mu_plus:
            dr_vals.append(e_plus[0].deltaR(mu_plus[0]))

    return np.array(dr_vals, dtype=float)


def compute_dphi_mumu(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]; py = arrays["py"]; pz = arrays["pz"]; E = arrays["E"]; pdg = arrays["pdgID"]

    dphi_list = []
    for i in range(len(pdg)):
        mu_plus = []
        mu_minus = []
        for j in range(len(pdg[i])):
            v = vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j])
            if pdg[i][j] == -13:  # μ+
                mu_plus.append(v)
            elif pdg[i][j] == 13:  # μ−
                mu_minus.append(v)

        if mu_plus and mu_minus:
            delta_phi = mu_plus[0].phi - mu_minus[0].phi
            delta_phi = np.arctan2(np.sin(delta_phi), np.cos(delta_phi))  # wrap to [-π, π]
            dphi_list.append(np.degrees(np.abs(delta_phi)))

    return np.array(dphi_list, dtype=float)


def compute_dphi_epmup(tree):
    arrays = build_vectors(tree)
    px = arrays["px"]; py = arrays["py"]; pz = arrays["pz"]; E = arrays["E"]; pdg = arrays["pdgID"]

    dphi_list = []
    for i in range(len(pdg)):
        ep = []
        mup = []
        for j in range(len(pdg[i])):
            v = vector.obj(px=px[i][j], py=py[i][j], pz=pz[i][j], E=E[i][j])
            if pdg[i][j] == -11:  # e+
                ep.append(v)
            elif pdg[i][j] == -13:  # μ+
                mup.append(v)

        if ep and mup:
            delta_phi = ep[0].phi - mup[0].phi
            delta_phi = np.arctan2(np.sin(delta_phi), np.cos(delta_phi))
            dphi_list.append(np.degrees(np.abs(delta_phi)))

    return np.array(dphi_list, dtype=float)


# -----------------------
# Rapidity and Δy
# -----------------------
def compute_rapidity(vecs):
    y_vals = []
    for v in vecs:
        if v:
            y_vals.append(v[0].rapidity)
    return np.array(y_vals, dtype=float)


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
    return np.array(dy_vals, dtype=float)


def compute_dypp(tree):
    arrays = build_vectors(tree)
    electrons = get_particles(arrays, ELECTRON_ID)
    muons = get_particles(arrays, MUON_ID)
    dy_vals = []
    for e, m in zip(electrons, muons):
        if e and m:
            dy_vals.append(np.abs(e[0].rapidity - m[0].rapidity))
    return np.array(dy_vals, dtype=float)


# -----------------------
# MET (from neutrinos)
# -----------------------
def compute_met(tree):
    arrays = build_vectors(tree)
    pdg_ids = arrays["pdgID"]; pxs = arrays["px"]; pys = arrays["py"]

    met = []
    for ids, px, py in zip(pdg_ids, pxs, pys):
        met_px = sum(px[i] for i, pid in enumerate(ids) if pid in NEUTRINO_IDS)
        met_py = sum(py[i] for i, pid in enumerate(ids) if pid in NEUTRINO_IDS)
        met_val = np.sqrt(met_px**2 + met_py**2)
        met.append(met_val)

    return np.array(met, dtype=float)


def compute_met_phi(tree):
    arrays = build_vectors(tree)
    pdg_ids = arrays["pdgID"]; pxs = arrays["px"]; pys = arrays["py"]

    met_phi = []
    for ids, px, py in zip(pdg_ids, pxs, pys):
        met_px = sum(px[i] for i, pid in enumerate(ids) if pid in NEUTRINO_IDS)
        met_py = sum(py[i] for i, pid in enumerate(ids) if pid in NEUTRINO_IDS)
        if met_px == 0 and met_py == 0:
            continue  # Skip zero-MET events
        phi = np.arctan2(met_py, met_px)
        met_phi.append(phi)

    return np.array(met_phi, dtype=float)




def _cos_theta_star_positive_lepton_from_Z(arrays, mZ=91.1876):
    """
    cos(theta*) for the Z decay positive lepton in the Z rest frame.
    Reference axis: Z direction in the WZ rest frame (helicity frame).
    Returns per-event tuple: (flavor, cos_theta_star) or (None, None).
    """
    zc = _pick_best_z_pair(arrays, mZ=mZ)
    nu = _sum_neutrinos(arrays)
    wlep = _pick_w_lepton(arrays, zc)

    out = []
    for zinfo, nuv, lv in zip(zc, nu, wlep):
        if zinfo is None or nuv is None or lv is None:
            out.append((None, None)); continue

        Z = zinfo["Z"]
        W = lv + nuv
        WZ = W + Z

        # Z direction in WZ rest frame
        p4Z_in_WZ = _to_rest_frame(
            (Z.E, Z.px, Z.py, Z.pz),
            (WZ.E, WZ.px, WZ.py, WZ.pz),
        )
        z_axis = _unit(p4Z_in_WZ[1:])
        if z_axis is None:
            out.append((None, None)); continue

        # positive lepton momentum in Z rest frame
        lp = zinfo["lp"]
        p4lp_in_Z = _to_rest_frame(
            (lp.E, lp.px, lp.py, lp.pz),
            (Z.E,  Z.px,  Z.py,  Z.pz),
        )
        p_lp = p4lp_in_Z[1:]
        u_lp = _unit(p_lp)
        if u_lp is None:
            out.append((None, None)); continue

        out.append((zinfo["flavor"], float(np.dot(u_lp, z_axis))))
    return out

def compute_cthep(tree):
    arrays = build_vectors(tree)
    vals = _cos_theta_star_positive_lepton_from_Z(arrays)
    # keep only Z->ee events
    out = [c for fl, c in vals if fl == "ee" and c is not None]
    return np.array(out, dtype=float)

def compute_cthmm(tree):
    arrays = build_vectors(tree)
    vals = _cos_theta_star_positive_lepton_from_Z(arrays)
    # keep only Z->mumu events
    out = [c for fl, c in vals if fl == "mumu" and c is not None]
    return np.array(out, dtype=float)

def compute_cscatt(tree):
    arrays = build_vectors(tree)
    zc = _pick_best_z_pair(arrays)
    nu = _sum_neutrinos(arrays)
    wlep = _pick_w_lepton(arrays, zc)

    out = []
    for zinfo, nuv, lv in zip(zc, nu, wlep):
        if zinfo is None or nuv is None or lv is None:
            continue
        Z = zinfo["Z"]
        W = lv + nuv
        WZ = W + Z

        p4W_in_WZ = _to_rest_frame(
            (W.E, W.px, W.py, W.pz),
            (WZ.E, WZ.px, WZ.py, WZ.pz),
        )
        pW = p4W_in_WZ[1:]
        n = np.linalg.norm(pW)
        if n == 0:
            continue
        cos_theta_W = pW[2] / n  # z-axis is beam
        out.append(abs(float(cos_theta_W)))

    return np.array(out, dtype=float)


def compute_cosuu(tree):
    """
    cos Δθ_{μ+μ-} : opening angle between μ+ and μ− (from the Z candidate if present).
    """
    arrays = build_vectors(tree)
    zc = _pick_best_z_pair(arrays)

    out = []
    for zinfo in zc:
        if zinfo is None or zinfo["flavor"] != "mumu":
            continue
        lp = zinfo["lp"]  # μ+
        lm = zinfo["lm"]  # μ-
        c = _cos_opening_angle_3d(lp, lm)
        if c is not None:
            out.append(c)
    return np.array(out, dtype=float)

def compute_cospp(tree):
    """
    cos Δθ_{e+μ+} : opening angle between leading e+ and leading μ+ in the event.
    (Matches your table label intent; adjust if you want a different pairing rule.)
    """
    arrays = build_vectors(tree)
    e_plus, _, mu_plus, _ = _split_leptons_with_charge(arrays)

    out = []
    for ep, mp in zip(e_plus, mu_plus):
        if not ep or not mp:
            continue
        c = _cos_opening_angle_3d(ep[0], mp[0])  # leading by pT (already sorted)
        if c is not None:
            out.append(c)
    return np.array(out, dtype=float)

##new functions per event

def _pick_z_pair_vectors(arrays, lep_id=11, mZ=91.1876):
    # lep_id=11 -> electrons, 13 -> muons
    e_plus, e_minus, mu_plus, mu_minus = _split_leptons_with_charge(arrays)
    plus, minus = (e_plus, e_minus) if lep_id==11 else (mu_plus, mu_minus)
    all_m = _all_masses(plus, minus)
    picked = _pick_closest_to_target_per_event(all_m, mZ)
    pair = []
    for P, M, mstar in zip(plus, minus, picked):
        if mstar is None:
            pair.append(None)
            continue
        # find the specific pair whose mass matches mstar (first match)
        chosen = None
        for a in P:
            for b in M:
                if abs((a+b).mass - mstar) < 1e-6:
                    chosen = (a, b); break
            if chosen: break
        pair.append(chosen)
    return pair  # list of (v_plus, v_minus) or None

def compute_dphi_ee_event(tree):
    arrays = build_vectors(tree)
    pairs = _pick_z_pair_vectors(arrays, lep_id=11)
    out = []
    for pr in pairs:
        if pr is None: continue
        a,b = pr
        dphi = np.arctan2(np.sin(a.phi-b.phi), np.cos(a.phi-b.phi))
        out.append(np.degrees(abs(dphi)))
    return np.array(out, dtype=float)

def compute_dr_ee_event(tree):
    arrays = build_vectors(tree)
    pairs = _pick_z_pair_vectors(arrays, lep_id=11)
    return np.array([a.deltaR(b) for pr in pairs if pr is not None for a,b in [pr]], dtype=float)

def compute_dphi_mumu_event(tree):
    arrays = build_vectors(tree)
    pairs = _pick_z_pair_vectors(arrays, lep_id=13)
    out = []
    for pr in pairs:
        if pr is None: continue
        a,b = pr
        dphi = np.arctan2(np.sin(a.phi-b.phi), np.cos(a.phi-b.phi))
        out.append(np.degrees(abs(dphi)))
    return np.array(out, dtype=float)

def compute_dr_mumu_event(tree):
    arrays = build_vectors(tree)
    pairs = _pick_z_pair_vectors(arrays, lep_id=13)
    return np.array([a.deltaR(b) for pr in pairs if pr is not None for a,b in [pr]], dtype=float)


def compute_dphi_emu_event(tree):
    arrays = build_vectors(tree)
    e_lead = _leading_vector_by_pt(arrays, 11)
    m_lead = _leading_vector_by_pt(arrays, 13)
    out = []
    for e,mu in zip(e_lead, m_lead):
        if e is None or mu is None: continue
        dphi = np.arctan2(np.sin(e.phi-mu.phi), np.cos(e.phi-mu.phi))
        out.append(np.degrees(abs(dphi)))
    return np.array(out, dtype=float)

def compute_dr_emu_event(tree):
    arrays = build_vectors(tree)
    e_lead = _leading_vector_by_pt(arrays, 11)
    m_lead = _leading_vector_by_pt(arrays, 13)
    return np.array([e.deltaR(mu) for e,mu in zip(e_lead, m_lead) if e is not None and mu is not None], dtype=float)



# -----------------------
# Event category counters (call from main if needed)
# -----------------------
def count_event_categories(tree, weight_branch=None):
    """
    Returns dict with unweighted and weighted counts:
      - SFOS: has ee OS or μμ OS
      - DF:   has at least one e and one μ (any charges)
      - DFOS: has an opposite-sign eμ pair
    If weight_branch is given (e.g. "weight"), multiply counts by weights.
    """
    arrays = build_vectors(tree)
    e_plus, e_minus, mu_plus, mu_minus = _split_leptons_with_charge(arrays)

    has_e  = [bool(ep or em) for ep, em in zip(e_plus, e_minus)]
    has_mu = [bool(mp or mm) for mp, mm in zip(mu_plus, mu_minus)]

    has_ee_os = [bool(ep) and bool(em) for ep, em in zip(e_plus, e_minus)]
    has_mm_os = [bool(mp) and bool(mm) for mp, mm in zip(mu_plus, mu_minus)]
    SFOS = np.array([a or b for a, b in zip(has_ee_os, has_mm_os)], dtype=bool)
    DF   = np.array([a and b for a, b in zip(has_e, has_mu)], dtype=bool)
    DFOS = np.array([ (bool(ep) and bool(mm)) or (bool(mp) and bool(em))
                      for ep, em, mp, mm in zip(e_plus, e_minus, mu_plus, mu_minus) ], dtype=bool)

    if weight_branch is not None:
        w = tree.arrays([weight_branch], library="np")[weight_branch]
        # If weights are stored per-event as variable-length arrays, take the first; else 1.0
        weights = np.array([ (w[i][0] if len(w[i])>0 else 1.0) for i in range(len(w)) ], dtype=float)
    else:
        weights = np.ones(len(SFOS), dtype=float)

    def wsum(mask):
        return float(np.sum(weights[mask]))

    return {
        "unweighted": {
            "N_SFOS": int(SFOS.sum()),
            "N_DF":   int(DF.sum()),
            "N_DFOS": int(DFOS.sum()),
        },
        "weighted": {
            "N_SFOS": wsum(SFOS),
            "N_DF":   wsum(DF),
            "N_DFOS": wsum(DFOS),
        },
    }


# -----------------------
# Dispatch & labels
# -----------------------
dispatch_map = {
    #"M_ee": compute_mass_ee,
    "M_mumu": compute_mass_mumu,
    "M_epmup": compute_mass_epmup,
    "M_event_rep": compute_mass_event_representative,  # NEW: one fill per event
    "M_ee_Z": compute_mass_ee_closestZ,
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
    "eta_leading_electron": compute_leading_eta_electron,  # fixed: true leading
    "eta_leading_muon": compute_leading_eta_muon,          # fixed: true leading
    "phi_leading_electron": compute_leading_phi_electron,  # fixed: true leading
    "phi_leading_muon": compute_leading_phi_muon,          # fixed: true leading

    "MET": compute_met,
    "MET_phi": compute_met_phi,
    
    "dphi_ee_event": compute_dphi_ee_event,
    "DR_ee_event":   compute_dr_ee_event,
    "dphi_mumu_event": compute_dphi_mumu_event,
    "DR_mumu_event":   compute_dr_mumu_event,
    "dphi_emu_event":  compute_dphi_emu_event,
    "DR_emu_event":    compute_dr_emu_event,

    "cthep":  compute_cthep,   # cos θ*_{e+}
    "cthmm":  compute_cthmm,   # cos θ*_{μ+}
    "cscatt": compute_cscatt,  # |cos θ_W|
    "cosuu":  compute_cosuu,   # cos Δθ_{μ+μ-}
    "cospp":  compute_cospp,   # cos Δθ_{e+μ+}
    
    
    
}

label_map = {
    "M_ee": r"$M_{e^+e^-}$ [GeV]",
    "M_mumu": r"$M_{\mu^+\mu^-}$ [GeV]",
    "M_epmup": r"$M_{e^+\mu^+}$ [GeV]",
    "M_event_rep": r"$M_{\mathrm{event}}$ [GeV]",
    "M_ee_Z" : r"$M_{e^+e^-}\ (\text{closest to } m_Z)$ [GeV]",
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
    
    "dphi_ee_event":   r"$\Delta\phi(e^+, e^-)$ [deg] (per event)",
    "DR_ee_event":     r"$\Delta R(e^+, e^-)$ (per event)",
    "dphi_mumu_event": r"$\Delta\phi(\mu^+, \mu^-)$ [deg] (per event)",
    "DR_mumu_event":   r"$\Delta R(\mu^+, \mu^-)$ (per event)",
    "dphi_emu_event":  r"$\Delta\phi(e, \mu)$ [deg] (per event)",
    "DR_emu_event":    r"$\Delta R(e, \mu)$ (per event)",
    "cthep":  r"$\cos\theta^{*}_{e^{+}}$",
    "cthmm":  r"$\cos\theta^{*}_{\mu^{+}}$",
    "cscatt": r"$|\cos\theta_{W}|$",
    "cosuu":  r"$\cos\Delta\theta_{\mu^{+},\mu^{-}}$",
    "cospp":  r"$\cos\Delta\theta_{e^{+},\mu^{+}}$",

}

