"""Photophysical rate constants.

k_r   : simplified Strickler-Berg (Einstein A) from the S1 oscillator strength f and energy E.
k_ISC : Fermi golden rule with a Marcus-Levich-Jortner Franck-Condon weighted density of states
        (one classical low-frequency mode + solvent lambda_s, and one effective quantum mode
        hbar*omega with Huang-Rhys factor S). It includes the S-T gap, which the |<T|H|S>|^2
        proxy ignores.
Yields: Phi_ISC = k_ISC/(k_r + k_ISC + k_IC), Phi_F = k_r/(...).
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from math import factorial

KB_EV = 8.617333262e-5
HBAR_EVS = 6.582119569e-16
CM1_TO_EV = 1 / 8065.54393734921


@dataclass
class MLJParams:
    lambda_s: float = 0.10     # eV, classical reorganization energy (soft modes + solvent)
    hw: float = 0.17           # eV, effective quantum mode (~1370 cm^-1, C=C/C=N stretches)
    S: float = 0.8             # Huang-Rhys factor of the effective mode
    T: float = 298.15
    nmax: int = 30


def fcwd_mlj(dE: float, p: MLJParams) -> float:
    """Franck-Condon weighted density of states (1/eV). dE = E_final - E_initial (eV, negative = downhill)."""
    kT = KB_EV * p.T
    pref = 1 / np.sqrt(4 * np.pi * p.lambda_s * kT)
    tot = 0.0
    for n in range(p.nmax + 1):
        w = np.exp(-p.S) * p.S ** n / factorial(n)
        tot += w * np.exp(-(dE + p.lambda_s + n * p.hw) ** 2 / (4 * p.lambda_s * kT))
    return pref * tot


def k_isc(V_cm: float, dE: float, p: MLJParams = MLJParams()) -> float:
    """s^-1. V in cm^-1, dE = E_T - E_S in eV."""
    V = V_cm * CM1_TO_EV
    return 2 * np.pi / HBAR_EVS * V ** 2 * fcwd_mlj(dE, p)


def k_rad(f: float, E_eV: float, n_refr: float = 1.0) -> float:
    """Einstein: k_r = n^2 f nu~^2 / 1.499 with nu~ in cm^-1 -> s^-1."""
    nu = E_eV * 8065.54393734921
    return n_refr ** 2 * f * nu ** 2 / 1.499


def photophysics(E_S1, f_S1, E_T, V_S1T, k_ic=1.5e8, n_refr=1.0, p=MLJParams(), shift_S=0.0, shift_T=0.0):
    """Compute k_r, the per-channel S1->Tn k_ISC, and the yields.
    shift_S/shift_T allow applying energy corrections (e.g. from a correlated backend)."""
    ES = E_S1 + shift_S
    channels = []
    for n, (ET, V) in enumerate(zip(E_T, V_S1T)):
        dE = (ET + shift_T) - ES
        channels.append({"T": n + 1, "dE_eV": float(dE), "V_cm": float(V), "k": float(k_isc(V, dE, p))})
    kisc = sum(c["k"] for c in channels)
    kr = k_rad(f_S1, ES, n_refr)
    ktot = kr + kisc + k_ic
    return {"k_r": kr, "k_isc": kisc, "k_ic": k_ic, "phi_F": kr / ktot, "phi_ISC": kisc / ktot,
            "tau_S_ns": 1e9 / ktot, "channels": channels}
