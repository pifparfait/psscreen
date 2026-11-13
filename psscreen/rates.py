"""Constantes de velocidad fotofísicas.

k_r   : Strickler-Berg simplificado (Einstein A) a partir de f y E de S1.
k_ISC : regla de oro de Fermi con densidad de estados Franck-Condon de Marcus-Levich-Jortner
        (un modo clásico de baja frecuencia + disolvente lambda_s, y un modo cuántico efectivo
        hbar*omega con factor Huang-Rhys S). Incluye el gap S-T, que el proxy |<T|H|S>|^2 ignora.
Rendimientos: Phi_ISC = k_ISC/(k_r + k_ISC + k_IC), Phi_F = k_r/(...).
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
    lambda_s: float = 0.10     # eV, reorganización clásica (modos blandos + disolvente)
    hw: float = 0.17           # eV, modo cuántico efectivo (~1370 cm^-1, estiramientos C=C/C=N)
    S: float = 0.8             # Huang-Rhys del modo efectivo
    T: float = 298.15
    nmax: int = 30


def fcwd_mlj(dE: float, p: MLJParams) -> float:
    """Densidad FC ponderada (1/eV). dE = E_final - E_inicial (eV, negativo = cuesta abajo)."""
    kT = KB_EV * p.T
    pref = 1 / np.sqrt(4 * np.pi * p.lambda_s * kT)
    tot = 0.0
    for n in range(p.nmax + 1):
        w = np.exp(-p.S) * p.S ** n / factorial(n)
        tot += w * np.exp(-(dE + p.lambda_s + n * p.hw) ** 2 / (4 * p.lambda_s * kT))
    return pref * tot


def k_isc(V_cm: float, dE: float, p: MLJParams = MLJParams()) -> float:
    """s^-1. V en cm^-1, dE = E_T - E_S en eV."""
    V = V_cm * CM1_TO_EV
    return 2 * np.pi / HBAR_EVS * V ** 2 * fcwd_mlj(dE, p)


def k_rad(f: float, E_eV: float, n_refr: float = 1.0) -> float:
    """Einstein: k_r = 2 n^2 f nu~^2 / (1.499) con nu~ en cm^-1 -> s^-1."""
    nu = E_eV * 8065.54393734921
    return n_refr ** 2 * f * nu ** 2 / 1.499


def photophysics(E_S1, f_S1, E_T, V_S1T, k_ic=1.5e8, n_refr=1.0, p=MLJParams(), shift_S=0.0, shift_T=0.0):
    """Calcula k_r, k_ISC por canal S1->Tn, y rendimientos.
    shift_S/shift_T permiten aplicar correcciones (p.ej. de un backend correlacionado)."""
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
