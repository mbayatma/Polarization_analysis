"""Reco-level (Pythia8 + Delphes) reader and observables for the 0j+1j WZ samples.

The LHE-level path in `physics_variables.py` works on truth particle records
(`px, py, pz, E, pdgID` in `mytree`). The showered samples are a different
format entirely: a flat `events` tree of reconstructed objects
(`el_pt/el_eta/el_phi/el_charge`, the same for muons, plus jets, MET and HT),
so they need their own reader and their own observable definitions.

Two things differ from the LHE path by design:

* **Normalization comes from the file, not a CSV.** Each event carries
  `xsec_weight` (pb), and summing it over a sample reproduces that sample's
  cross section times the shower/reco acceptance -- which is what the
  producer's VALIDATION.txt checks. Histograms are therefore filled with
  `xsec_weight` directly instead of rescaling by a hand-maintained
  `xsections_all.csv` entry. A small fraction of events carry negative
  weights, which is normal here and is handled by summing weights rather
  than counting entries.
* **Everything is vectorized.** The LHE helpers loop in Python over events,
  which is fine for one variable but not for 52 samples x ~30 variables.
  After the exactly-3-lepton requirement the event shape is regular, so the
  leptons collapse to fixed (N, 3) arrays and every observable is plain numpy.
"""

import numpy as np
import uproot

MZ = 91.1876
MW = 80.377
M_EL = 0.000511
M_MU = 0.10566

FLAV_EL = 11
FLAV_MU = 13

TREE_NAME = "events"

BRANCHES = [
    "weight", "xsec_weight",
    "el_pt", "el_eta", "el_phi", "el_charge",
    "mu_pt", "mu_eta", "mu_phi", "mu_charge",
    "jet_pt", "jet_eta", "jet_phi", "jet_m", "jet_btag",
    "met", "met_phi", "ht",
]

# Baseline object selection, ATLAS-like and deliberately loose.
LEP_PT_MIN = 7.0
LEP_ETA_MAX = 2.5
JET_PT_MIN = 25.0
JET_ETA_MAX = 4.5

# WZ 3-lepton event selection
N_LEPTONS = 3
LEAD_LEP_PT_MIN = 20.0   # single-lepton-trigger-like requirement
MZ_WINDOW = 15.0         # |m_SFOS - mZ| < 15 GeV


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def load_arrays(path, tree_name=TREE_NAME, branches=BRANCHES):
    """Read the flat reco tree into an awkward record array."""
    with uproot.open(path) as f:
        return f[tree_name].arrays(branches)


def _padded(jagged, n, fill=np.nan):
    """Convert a jagged array to a dense (len, n) float64 numpy array."""
    import awkward as ak
    padded = ak.pad_none(jagged, n, axis=1, clip=True)
    return np.asarray(ak.to_numpy(ak.fill_none(padded, fill)), dtype=np.float64)


def build_leptons(ev):
    """Merge electrons and muons into one pT-ordered lepton collection.

    Object-level quality cuts are applied here, before any event-level
    counting, so the lepton multiplicity reflects *selected* leptons rather
    than everything Delphes reconstructed.
    """
    import awkward as ak

    el_ok = (ev.el_pt > LEP_PT_MIN) & (abs(ev.el_eta) < LEP_ETA_MAX)
    mu_ok = (ev.mu_pt > LEP_PT_MIN) & (abs(ev.mu_eta) < LEP_ETA_MAX)

    def _coll(pt, eta, phi, charge, flav, mass):
        ptf = ak.values_astype(pt, np.float64)
        return ak.zip({
            "pt": ptf,
            "eta": ak.values_astype(eta, np.float64),
            "phi": ak.values_astype(phi, np.float64),
            "charge": ak.values_astype(charge, np.float64),
            "flav": ak.full_like(ptf, flav),
            "mass": ak.full_like(ptf, mass),
        })

    el = _coll(ev.el_pt[el_ok], ev.el_eta[el_ok], ev.el_phi[el_ok],
               ev.el_charge[el_ok], FLAV_EL, M_EL)
    mu = _coll(ev.mu_pt[mu_ok], ev.mu_eta[mu_ok], ev.mu_phi[mu_ok],
               ev.mu_charge[mu_ok], FLAV_MU, M_MU)

    lep = ak.concatenate([el, mu], axis=1)
    return lep[ak.argsort(lep.pt, axis=1, ascending=False)]


# ---------------------------------------------------------------------------
# Kinematic helpers (all numpy, all vectorized)
# ---------------------------------------------------------------------------
def _p4(pt, eta, phi, mass):
    """(pt, eta, phi, m) -> four-vector array with last axis (px, py, pz, E)."""
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    E = np.sqrt(np.maximum(px**2 + py**2 + pz**2 + mass**2, 0.0))
    return np.stack([px, py, pz, E], axis=-1)


def _mass(p4):
    m2 = p4[..., 3] ** 2 - np.sum(p4[..., :3] ** 2, axis=-1)
    return np.sqrt(np.maximum(m2, 0.0))


def _pt(p4):
    return np.hypot(p4[..., 0], p4[..., 1])


def _phi(p4):
    return np.arctan2(p4[..., 1], p4[..., 0])


def _eta(p4):
    p = np.linalg.norm(p4[..., :3], axis=-1)
    pz = p4[..., 2]
    return 0.5 * np.log(np.clip((p + pz) / np.clip(p - pz, 1e-12, None), 1e-12, None))


def _rapidity(p4):
    E, pz = p4[..., 3], p4[..., 2]
    return 0.5 * np.log(np.clip((E + pz) / np.clip(E - pz, 1e-12, None), 1e-12, None))


def delta_phi(phi1, phi2):
    """Signed phi difference wrapped into (-pi, pi]."""
    return (phi1 - phi2 + np.pi) % (2 * np.pi) - np.pi


def delta_r(eta1, phi1, eta2, phi2):
    return np.hypot(eta1 - eta2, delta_phi(phi1, phi2))


def _boost_to_rest_frame(p4_child, p4_parent):
    """Boost `p4_child` into the rest frame of `p4_parent` (both (N, 4))."""
    E_par = np.clip(p4_parent[:, 3], 1e-12, None)
    beta = p4_parent[:, :3] / E_par[:, None]
    b2 = np.clip(np.sum(beta**2, axis=1), 0.0, 1.0 - 1e-12)
    gamma = 1.0 / np.sqrt(1.0 - b2)

    bp = np.sum(beta * p4_child[:, :3], axis=1)
    gamma2 = np.where(b2 > 0, (gamma - 1.0) / np.where(b2 > 0, b2, 1.0), 0.0)

    p_new = (p4_child[:, :3]
             + (gamma2 * bp)[:, None] * beta
             - (gamma * p4_child[:, 3])[:, None] * beta)
    E_new = gamma * (p4_child[:, 3] - bp)
    return np.concatenate([p_new, E_new[:, None]], axis=1)


def cos_theta_star(p4_lepton, p4_boson):
    """cos of the angle between the lepton -- boosted into the boson rest
    frame -- and the boson flight direction in the lab.

    This is the standard polarization observable: a longitudinally polarized
    boson gives a 1 - cos^2 distribution, a transverse one gives
    (1 +/- cos)^2 / 2.
    """
    lep_rest = _boost_to_rest_frame(p4_lepton, p4_boson)
    v_lep, v_bos = lep_rest[:, :3], p4_boson[:, :3]
    num = np.sum(v_lep * v_bos, axis=1)
    den = np.linalg.norm(v_lep, axis=1) * np.linalg.norm(v_bos, axis=1)
    out = np.where(den > 0, num / np.clip(den, 1e-12, None), np.nan)
    return np.clip(out, -1.0, 1.0)


def _neutrino_pz(p4_lep, nu_px, nu_py, mw=MW):
    """Reconstruct the neutrino pz from the on-shell W mass constraint.

    Standard quadratic solution; where the discriminant is negative (which
    happens for a sizeable fraction of events once MET is smeared) the real
    part is used, and of the two real roots the one with smaller |pz| is
    taken. Both are the usual conventions in WZ analyses.
    """
    lpx, lpy, lpz, lE = p4_lep[:, 0], p4_lep[:, 1], p4_lep[:, 2], p4_lep[:, 3]
    met = np.hypot(nu_px, nu_py)
    mu = 0.5 * mw**2 + lpx * nu_px + lpy * nu_py
    lpt2 = lpx**2 + lpy**2
    denom = np.clip(lpt2, 1e-12, None)

    a = mu * lpz / denom
    disc = (mu**2 * lpz**2) / denom**2 - (lE**2 * met**2 - mu**2) / denom

    real = disc < 0
    sqrt_disc = np.sqrt(np.maximum(disc, 0.0))
    pz1, pz2 = a + sqrt_disc, a - sqrt_disc
    # smaller |pz| root, or the degenerate root when the discriminant is < 0
    pz = np.where(np.abs(pz1) < np.abs(pz2), pz1, pz2)
    return np.where(real, a, pz)


# ---------------------------------------------------------------------------
# WZ 3-lepton reconstruction
# ---------------------------------------------------------------------------
_PAIRS = ((0, 1), (0, 2), (1, 2))


def select_wz_3l(ev):
    """Apply the WZ 3-lepton selection and reconstruct the Z and W candidates.

    Z candidate = the same-flavour opposite-sign pair whose mass is closest to
    mZ; the W lepton is the remaining one. With exactly three leptons the Z
    pair indices (i, j) and the W index k satisfy i + j + k = 3, which is how
    the W lepton is identified without a second combinatorial pass.

    Returns a dict of per-event numpy arrays (one entry per selected event),
    or None if no event survives.
    """
    import awkward as ak

    lep = build_leptons(ev)
    n_lep = np.asarray(ak.num(lep, axis=1))
    lead_pt = np.asarray(ak.fill_none(ak.firsts(lep.pt), -1.0))

    presel = (n_lep == N_LEPTONS) & (lead_pt > LEAD_LEP_PT_MIN)
    if not presel.any():
        return None

    lep3, ev3 = lep[presel], ev[presel]

    pt = _padded(lep3.pt, 3)
    eta = _padded(lep3.eta, 3)
    phi = _padded(lep3.phi, 3)
    charge = _padded(lep3.charge, 3)
    flav = _padded(lep3.flav, 3)
    mass = _padded(lep3.mass, 3)

    p4 = _p4(pt, eta, phi, mass)            # (N, 3, 4)
    n = p4.shape[0]

    # All three lepton pairs; only SFOS pairs can be the Z
    pair_mass = np.full((n, 3), np.nan)
    is_sfos = np.zeros((n, 3), dtype=bool)
    for k, (i, j) in enumerate(_PAIRS):
        pair_mass[:, k] = _mass(p4[:, i, :] + p4[:, j, :])
        is_sfos[:, k] = (flav[:, i] == flav[:, j]) & (charge[:, i] * charge[:, j] < 0)

    dist = np.where(is_sfos, np.abs(pair_mass - MZ), np.inf)
    best = np.argmin(dist, axis=1)
    keep = dist[np.arange(n), best] < MZ_WINDOW
    if not keep.any():
        return None

    rows = np.arange(n)[keep]
    i_z1 = np.array([_PAIRS[b][0] for b in best])[keep]
    i_z2 = np.array([_PAIRS[b][1] for b in best])[keep]
    i_w = 3 - i_z1 - i_z2                   # leftover index of {0, 1, 2}

    p4_z1, p4_z2 = p4[rows, i_z1, :], p4[rows, i_z2, :]
    p4_wl = p4[rows, i_w, :]
    q_z1 = charge[rows, i_z1]
    q_wl = charge[rows, i_w]

    # Order Z leptons by charge so cos(theta*) uses the negative lepton,
    # which is the usual convention and keeps the observable sign-consistent.
    z1_is_minus = (q_z1 < 0)[:, None]
    p4_zm = np.where(z1_is_minus, p4_z1, p4_z2)
    p4_z = p4_z1 + p4_z2

    met = np.asarray(ev3.met, dtype=np.float64)[keep]
    met_phi = np.asarray(ev3.met_phi, dtype=np.float64)[keep]
    ht = np.asarray(ev3.ht, dtype=np.float64)[keep]
    nu_px, nu_py = met * np.cos(met_phi), met * np.sin(met_phi)

    # W candidate: transverse mass needs no pz hypothesis; the full
    # four-vector (for cos theta*) uses the mW on-shell constraint.
    wl_pt, wl_phi = _pt(p4_wl), _phi(p4_wl)
    mt_w = np.sqrt(np.maximum(
        2.0 * wl_pt * met * (1.0 - np.cos(delta_phi(wl_phi, met_phi))), 0.0))

    nu_pz = _neutrino_pz(p4_wl, nu_px, nu_py)
    nu_E = np.sqrt(nu_px**2 + nu_py**2 + nu_pz**2)
    p4_nu = np.stack([nu_px, nu_py, nu_pz, nu_E], axis=1)
    p4_w = p4_wl + p4_nu

    p4_3l = p4[rows, 0, :] + p4[rows, 1, :] + p4[rows, 2, :]
    p4_wz = p4_z + p4_w

    # Selected jets
    jet_pt_j = ev3.jet_pt[(ev3.jet_pt > JET_PT_MIN) & (abs(ev3.jet_eta) < JET_ETA_MAX)]
    n_jet = np.asarray(ak.num(jet_pt_j, axis=1))[keep].astype(np.float64)
    lead_jet_pt = np.asarray(ak.fill_none(ak.firsts(jet_pt_j), 0.0), dtype=np.float64)[keep]

    z_eta, z_phi = _eta(p4_z), _phi(p4_z)
    w_eta, w_phi = _eta(p4_w), _phi(p4_w)
    wl_eta = _eta(p4_wl)

    out = {
        # --- Z candidate ---
        "M_Z": _mass(p4_z),
        "pt_Z": _pt(p4_z),
        "y_Z": _rapidity(p4_z),
        "cos_theta_star_Z": cos_theta_star(p4_zm, p4_z),
        "dphi_Z_leptons": np.abs(delta_phi(_phi(p4_z1), _phi(p4_z2))),
        "dR_Z_leptons": delta_r(_eta(p4_z1), _phi(p4_z1), _eta(p4_z2), _phi(p4_z2)),
        "pt_Z_lep1": _pt(p4_z1),
        "pt_Z_lep2": _pt(p4_z2),

        # --- W candidate ---
        "mT_W": mt_w,
        "pt_W": _pt(p4_w),
        "y_W": _rapidity(p4_w),
        "cos_theta_star_W": cos_theta_star(p4_wl, p4_w),
        "pt_lep_W": wl_pt,
        "eta_lep_W": wl_eta,
        "dphi_lepW_MET": np.abs(delta_phi(wl_phi, met_phi)),

        # --- WZ system / event ---
        "M_3l": _mass(p4_3l),
        "pt_3l": _pt(p4_3l),
        "M_WZ": _mass(p4_wz),
        "pt_WZ": _pt(p4_wz),
        "dphi_Z_W": np.abs(delta_phi(z_phi, w_phi)),
        "dy_Z_W": np.abs(_rapidity(p4_z) - _rapidity(p4_w)),
        "dR_lepW_Z": delta_r(wl_eta, wl_phi, z_eta, z_phi),
        "MET": met,
        "HT": ht,
        "n_jet": n_jet,
        "pt_lead_jet": lead_jet_pt,
        "pt_lead_lep": pt[rows, 0],
        "eta_lead_lep": eta[rows, 0],

        # --- bookkeeping (not features) ---
        "weight": np.asarray(ev3.weight, dtype=np.float64)[keep],
        "xsec_weight": np.asarray(ev3.xsec_weight, dtype=np.float64)[keep],
        "W_charge": q_wl,
        "Z_is_mumu": (flav[rows, i_z1] == FLAV_MU).astype(np.float64),
    }
    return out


# Observables plotted / exported. `weight`, `xsec_weight`, `W_charge` and
# `Z_is_mumu` are bookkeeping, so they are deliberately not in this list.
RECO_VARIABLES = [
    "M_Z", "pt_Z", "y_Z", "cos_theta_star_Z",
    "dphi_Z_leptons", "dR_Z_leptons", "pt_Z_lep1", "pt_Z_lep2",
    "mT_W", "pt_W", "y_W", "cos_theta_star_W",
    "pt_lep_W", "eta_lep_W", "dphi_lepW_MET",
    "M_3l", "pt_3l", "M_WZ", "pt_WZ",
    "dphi_Z_W", "dy_Z_W", "dR_lepW_Z",
    "MET", "HT", "n_jet", "pt_lead_jet",
    "pt_lead_lep", "eta_lead_lep",
]

# Axis labels and binning for the comparison plots.
VAR_SPECS = {
    "M_Z":              (r"$m_{\ell^+\ell^-}$ (Z cand.) [GeV]", np.linspace(76, 106, 46)),
    "pt_Z":             (r"$p_T^{Z}$ [GeV]", np.linspace(0, 400, 41)),
    "y_Z":              (r"$y_{Z}$", np.linspace(-3, 3, 41)),
    "cos_theta_star_Z": (r"$\cos\theta^{*}_{\ell^-}$ (Z frame)", np.linspace(-1, 1, 41)),
    "dphi_Z_leptons":   (r"$\Delta\phi(\ell^+,\ell^-)_Z$ [rad]", np.linspace(0, np.pi, 37)),
    "dR_Z_leptons":     (r"$\Delta R(\ell^+,\ell^-)_Z$", np.linspace(0, 5, 41)),
    "pt_Z_lep1":        (r"$p_T^{\ell_1, Z}$ [GeV]", np.linspace(0, 300, 41)),
    "pt_Z_lep2":        (r"$p_T^{\ell_2, Z}$ [GeV]", np.linspace(0, 200, 41)),
    "mT_W":             (r"$m_T^{W}$ [GeV]", np.linspace(0, 200, 41)),
    "pt_W":             (r"$p_T^{W}$ [GeV]", np.linspace(0, 400, 41)),
    "y_W":              (r"$y_{W}$", np.linspace(-4, 4, 41)),
    "cos_theta_star_W": (r"$\cos\theta^{*}_{\ell_W}$ (W frame)", np.linspace(-1, 1, 41)),
    "pt_lep_W":         (r"$p_T^{\ell_W}$ [GeV]", np.linspace(0, 300, 41)),
    "eta_lep_W":        (r"$\eta_{\ell_W}$", np.linspace(-3, 3, 41)),
    "dphi_lepW_MET":    (r"$\Delta\phi(\ell_W, E_T^{miss})$ [rad]", np.linspace(0, np.pi, 37)),
    "M_3l":             (r"$m_{3\ell}$ [GeV]", np.linspace(0, 600, 41)),
    "pt_3l":            (r"$p_T^{3\ell}$ [GeV]", np.linspace(0, 400, 41)),
    "M_WZ":             (r"$m_{WZ}$ [GeV]", np.linspace(0, 1200, 41)),
    "pt_WZ":            (r"$p_T^{WZ}$ [GeV]", np.linspace(0, 300, 41)),
    "dphi_Z_W":         (r"$\Delta\phi(Z, W)$ [rad]", np.linspace(0, np.pi, 37)),
    "dy_Z_W":           (r"$|\Delta y(Z, W)|$", np.linspace(0, 5, 41)),
    "dR_lepW_Z":        (r"$\Delta R(\ell_W, Z)$", np.linspace(0, 6, 41)),
    "MET":              (r"$E_T^{miss}$ [GeV]", np.linspace(0, 300, 41)),
    "HT":               (r"$H_T$ [GeV]", np.linspace(0, 800, 41)),
    "n_jet":            (r"$N_{jets}$ ($p_T > 25$ GeV)", np.arange(-0.5, 7.5, 1.0)),
    "pt_lead_jet":      (r"$p_T^{jet_1}$ [GeV]", np.linspace(0, 400, 41)),
    "pt_lead_lep":      (r"$p_T^{\ell_1}$ [GeV]", np.linspace(0, 400, 41)),
    "eta_lead_lep":     (r"$\eta_{\ell_1}$", np.linspace(-3, 3, 41)),
}
