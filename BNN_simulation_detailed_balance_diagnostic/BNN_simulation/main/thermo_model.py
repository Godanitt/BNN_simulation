"""Thermodynamic trajectories for the BBN notebooks.

The assignment gives three hard snapshots: (50 s, 1.7 GK, 2e-4 g/cm3),
(200 s, 1.0 GK, 2e-5 g/cm3), and (1000 s, 0.4 GK, 2e-6 g/cm3).

To avoid a purely piecewise log-linear interpolation, we use a radiation-dominated
baseline and multiply it by a smooth PCHIP correction constrained to be exactly 1
at the assignment nodes after rescaling. This keeps the requested points exactly,
while using the expected T ~ t^{-1/2} and rho_b ~ eta n_gamma(T) behavior as the
reference shape between them.
"""
from __future__ import annotations

from types import SimpleNamespace
import numpy as np
from scipy.interpolate import PchipInterpolator

ETA10_PLANCK_APPROX = 6.10
M_U_CGS = 1.66053906660e-24  # g


def default_snapshots():
    """Return the three assignment snapshots as a pandas DataFrame if pandas exists."""
    import pandas as pd
    return pd.DataFrame({
        "label": ["50 s", "200 s", "1000 s"],
        "t_s": [50.0, 200.0, 1000.0],
        "T9": [1.7, 1.0, 0.4],
        "rho_g_cm3": [2.0e-4, 2.0e-5, 2.0e-6],
    })


def n_gamma_cm3(T_K):
    """Blackbody photon number density in cm^-3: n_gamma = 20.28 T^3."""
    return 20.28 * np.asarray(T_K, dtype=float)**3


def T9_radiation_base(t, t_ref=200.0, T9_ref=1.0):
    """Radiation-dominated reference, normalized at (t_ref, T9_ref)."""
    return T9_ref * (np.asarray(t, dtype=float) / t_ref)**(-0.5)


def rho_eta_reference_from_T9(T9, eta10=ETA10_PLANCK_APPROX):
    """Baryon mass density from eta and blackbody photon density."""
    eta = eta10 * 1.0e-10
    T_K = np.asarray(T9, dtype=float) * 1.0e9
    return M_U_CGS * eta * n_gamma_cm3(T_K)


def make_thermo_functions(snapshots=None, eta10_ref=ETA10_PLANCK_APPROX):
    """Build smooth forced thermodynamic functions.

    T9(t) = T9_rad_base(t) * C_T(t), with C_T constrained so that T9(t_i)
    matches the assignment snapshots exactly.

    rho(t) = rho_eta_reference[T9(t), eta_ref] * C_rho(t), with C_rho
    constrained so that rho(t_i) matches the assignment snapshots exactly.
    """
    if snapshots is None:
        snapshots = default_snapshots()

    t_nodes = np.asarray(snapshots["t_s"], dtype=float)
    T9_nodes = np.asarray(snapshots["T9"], dtype=float)
    rho_nodes = np.asarray(snapshots["rho_g_cm3"], dtype=float)

    order = np.argsort(t_nodes)
    t_nodes = t_nodes[order]
    T9_nodes = T9_nodes[order]
    rho_nodes = rho_nodes[order]

    logt_nodes = np.log(t_nodes)

    T9_base_nodes = T9_radiation_base(t_nodes)
    log_T_corr_nodes = np.log(T9_nodes / T9_base_nodes)
    T_corr = PchipInterpolator(logt_nodes, log_T_corr_nodes, extrapolate=True)

    def T9_of_t(t):
        t = np.asarray(t, dtype=float)
        return T9_radiation_base(t) * np.exp(T_corr(np.log(t)))

    def T_of_t(t):
        return 1.0e9 * T9_of_t(t)

    rho_base_nodes = rho_eta_reference_from_T9(T9_nodes, eta10=eta10_ref)
    log_rho_corr_nodes = np.log(rho_nodes / rho_base_nodes)
    rho_corr = PchipInterpolator(logt_nodes, log_rho_corr_nodes, extrapolate=True)

    def rho_eta_reference_of_t(t):
        return rho_eta_reference_from_T9(T9_of_t(t), eta10=eta10_ref)

    def rho_of_t(t):
        t = np.asarray(t, dtype=float)
        return rho_eta_reference_of_t(t) * np.exp(rho_corr(np.log(t)))

    def eta_of_t(t):
        return (rho_of_t(t) / M_U_CGS) / n_gamma_cm3(T_of_t(t))

    def rho_correction_of_t(t):
        t = np.asarray(t, dtype=float)
        return np.exp(rho_corr(np.log(t)))

    def T_correction_of_t(t):
        t = np.asarray(t, dtype=float)
        return np.exp(T_corr(np.log(t)))

    return SimpleNamespace(
        eta10_ref=eta10_ref,
        t_nodes=t_nodes,
        T9_nodes=T9_nodes,
        rho_nodes=rho_nodes,
        T9_of_t=T9_of_t,
        T_of_t=T_of_t,
        rho_of_t=rho_of_t,
        eta_of_t=eta_of_t,
        n_gamma_cm3=n_gamma_cm3,
        T9_radiation_base=T9_radiation_base,
        rho_eta_reference_from_T9=lambda T9: rho_eta_reference_from_T9(T9, eta10=eta10_ref),
        rho_eta_reference_of_t=rho_eta_reference_of_t,
        T_correction_of_t=T_correction_of_t,
        rho_correction_of_t=rho_correction_of_t,
    )
