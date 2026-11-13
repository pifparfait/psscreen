"""High-degree QSP/QSVT angles by fixed-point iteration (Dong, Lin, Ni, 2022, arXiv:2209.10162).

PennyLane's iterative solver starts from monomial-basis coefficients (ill-conditioned above
degree ~60) and is slow at high degree. Here we work directly in the Chebyshev basis with
symmetric phases and the W_x convention; the angles are mapped to PennyLane's QSVT convention
with `qml.transform_angles(..., "QSP", "QSVT")`. Verified numerically against `qml.QSVT`.
"""
from __future__ import annotations
import numpy as np


def _nodes(dt):
    j = np.arange(1, dt + 1)
    return np.cos((2 * j - 1) * np.pi / (4 * dt))


def _qsp_im00(full_phases, x):
    """Im <0|U(x)|0> for U = e^{i f0 Z} prod_k [W(x) e^{i fk Z}], vectorized in x."""
    th = np.arccos(x)
    c, s = np.cos(th), np.sin(th)
    # one 2x2 matrix per point: start from e^{i f0 Z}
    a = np.exp(1j * full_phases[0]) * np.ones_like(x, dtype=complex)  # U = [[a, b],[-b*, a*]] SU(2)
    b = np.zeros_like(x, dtype=complex)
    for f in full_phases[1:]:
        # multiply by W = [[c, i s],[i s, c]]
        a, b = a * c + b * 1j * s, a * 1j * s + b * c
        # then by e^{i f Z} = diag(e^{if}, e^{-if})
        a, b = a * np.exp(1j * f), b * np.exp(-1j * f)
    return a.imag, a


def qsp_phases_fpi(cheb_parity_coeffs, degree, max_iter=500, tol=1e-12):
    """Symmetric W_x phases such that Im<0|U|0> = sum_k c_k T_{2k+p}(x).

    cheb_parity_coeffs: the c_k of the Chebyshev terms with the parity of `degree` (length ceil((degree+1)/2)).
    """
    p = degree % 2
    dt = int(np.ceil((degree + 1) / 2))
    c = np.zeros(dt); c[:len(cheb_parity_coeffs)] = cheb_parity_coeffs[:dt]
    x = _nodes(dt)
    T = np.cos(np.outer(np.arccos(x), 2 * np.arange(dt) + p))  # (nodes, dt)
    Tinv = np.linalg.inv(T)

    def full(red):
        if degree % 2 == 1:   # d+1 phases (even count)
            return np.concatenate([red, red[::-1]])
        return np.concatenate([red, red[-2::-1]])

    red = np.zeros(dt)
    for it in range(max_iter):
        g, _ = _qsp_im00(full(red), x)
        r = Tinv @ g - c
        err = np.max(np.abs(r))
        if err < tol:
            break
        # Jacobian at 0: reduced phase k controls T_{d-2k} with weight 2 (weight 1 for the
        # central phase when d is even): reverse order relative to c.
        delta = r[::-1] / 2.0
        if degree % 2 == 0:
            delta[-1] = r[0]
        red = red - delta
    return full(red), err, it


def to_pennylane_qsvt(phases):
    """W_x phases (Im<0|U|0> = f) -> angles for qml.QSVT with PCPhase, whose block has Re = f.
    The -pi/4 shift at both ends was verified numerically against qml.QSVT."""
    import pennylane as qml
    ph = np.array(phases, dtype=float).copy()
    ph[0] -= np.pi / 4; ph[-1] -= np.pi / 4
    return qml.transform_angles(ph, "QSP", "QSVT")
