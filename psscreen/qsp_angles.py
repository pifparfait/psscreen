"""Ángulos QSP/QSVT de grado alto por iteración de punto fijo (Dong, Lin, Ni, 2022, arXiv:2209.10162).

El resolvedor iterativo de PennyLane parte de coeficientes en base monomial (mal condicionados a
grado > ~60) y es lento a grado alto. Aquí se trabaja directamente en base de Chebyshev, con fases
simétricas y la convención W_x; los ángulos se pasan a la convención QSVT de PennyLane con
`qml.transform_angles(..., "QSP", "QSVT")`. Verificado numéricamente contra `qml.QSVT` (tests).
"""
from __future__ import annotations
import numpy as np


def _nodes(dt):
    j = np.arange(1, dt + 1)
    return np.cos((2 * j - 1) * np.pi / (4 * dt))


def _qsp_im00(full_phases, x):
    """Im <0|U(x)|0> para U = e^{i f0 Z} prod_k [W(x) e^{i fk Z}], vectorizado en x."""
    th = np.arccos(x)
    c, s = np.cos(th), np.sin(th)
    # matrices 2x2 por punto: empezamos con e^{i f0 Z}
    a = np.exp(1j * full_phases[0]) * np.ones_like(x, dtype=complex)  # U = [[a, b],[-b*, a*]] SU(2)
    b = np.zeros_like(x, dtype=complex)
    for f in full_phases[1:]:
        # multiplicar por W = [[c, i s],[i s, c]]
        a, b = a * c + b * 1j * s, a * 1j * s + b * c
        # por e^{i f Z} = diag(e^{if}, e^{-if})
        a, b = a * np.exp(1j * f), b * np.exp(-1j * f)
    return a.imag, a


def qsp_phases_fpi(cheb_parity_coeffs, degree, max_iter=500, tol=1e-12):
    """Fases simétricas W_x tales que Im<0|U|0> = sum_k c_k T_{2k+p}(x).

    cheb_parity_coeffs: c_k de los T de la paridad de `degree` (longitud ceil((degree+1)/2)).
    """
    p = degree % 2
    dt = int(np.ceil((degree + 1) / 2))
    c = np.zeros(dt); c[:len(cheb_parity_coeffs)] = cheb_parity_coeffs[:dt]
    x = _nodes(dt)
    T = np.cos(np.outer(np.arccos(x), 2 * np.arange(dt) + p))  # (nodes, dt)
    Tinv = np.linalg.inv(T)

    def full(red):
        if degree % 2 == 1:   # d+1 par de fases
            return np.concatenate([red, red[::-1]])
        return np.concatenate([red, red[-2::-1]])

    red = np.zeros(dt)
    for it in range(max_iter):
        g, _ = _qsp_im00(full(red), x)
        r = Tinv @ g - c
        err = np.max(np.abs(r))
        if err < tol:
            break
        # Jacobiano en 0: la fase reducida k controla T_{d-2k} con peso 2 (peso 1 para la fase
        # central si d es par): orden inverso al de c.
        delta = r[::-1] / 2.0
        if degree % 2 == 0:
            delta[-1] = r[0]
        red = red - delta
    return full(red), err, it


def to_pennylane_qsvt(phases):
    """Fases W_x (Im<0|U|0> = f) -> ángulos para qml.QSVT con PCPhase, cuyo bloque tiene Re = f.
    El desplazamiento -pi/4 en los extremos se verificó numéricamente contra qml.QSVT."""
    import pennylane as qml
    ph = np.array(phases, dtype=float).copy()
    ph[0] -= np.pi / 4; ph[-1] -= np.pi / 4
    return qml.transform_angles(ph, "QSP", "QSVT")
