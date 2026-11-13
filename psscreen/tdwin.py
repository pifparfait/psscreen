"""Explicit TDDFT/TDA in an orbital window (global hybrids, LDA/GGA).

Motivation: PySCF's Davidson re-evaluates the XC kernel on the grid for every trial vector
(~1 s/vector even for 20-atom molecules on 4 cores). Here the full A and B matrices are built
in the window's (i,a) pair subspace, evaluating the orbital products on the grid only once.
This is the standard approximation of truncating the excitation space by energy (as in sTDA,
but without simplifying the integrals). The truncation error is controlled by `ewin` and is
checked against full PySCF in tests/.

Convenciones (capa cerrada, orbitales reales):
  A_S = de + 2(ia|jb) - c_x (ij|ab) + (f_aa + f_ab)
  A_T = de            - c_x (ij|ab) + (f_aa - f_ab)
  B_S =      2(ia|jb) - c_x (ib|ja) + (f_aa + f_ab)
"""
from __future__ import annotations
import numpy as np
import scipy.linalg as sla
from pyscf import lib, df


def select_window(mo_energy, nocc, ewin_ev=8.0, max_occ=60, max_vir=80):
    e = mo_energy * 27.211386
    homo, lumo = e[nocc - 1], e[nocc]
    occ = [i for i in range(nocc) if e[i] > homo - ewin_ev][-max_occ:]
    vir = [a for a in range(nocc, len(e)) if e[a] < lumo + ewin_ev][:max_vir]
    return np.array(occ), np.array(vir)


def _mo_df(mf, Co, Cv, auxbasis):
    mol = mf.mol
    if getattr(mf, "with_df", None) is not None and getattr(mf.with_df, "_cderi", None) is not None:
        cderi = mf.with_df._cderi
    else:
        cderi = df.incore.cholesky_eri(mol, auxbasis=auxbasis)
    no, nv = Co.shape[1], Cv.shape[1]
    naux = cderi.shape[0]
    Lov = np.empty((naux, no, nv)); Loo = np.empty((naux, no, no)); Lvv = np.empty((naux, nv, nv))
    blk = max(1, int(4e7 // (mol.nao ** 2)))
    for p0 in range(0, naux, blk):
        L = lib.unpack_tril(np.asarray(cderi[p0:p0 + blk]))
        Lov[p0:p0 + blk] = lib.einsum("Lpq,pi,qa->Lia", L, Co, Cv)
        Loo[p0:p0 + blk] = lib.einsum("Lpq,pi,qj->Lij", L, Co, Co)
        Lvv[p0:p0 + blk] = lib.einsum("Lpq,pa,qb->Lab", L, Cv, Cv)
    return Lov, Loo, Lvv


def _xc_kernels(mf, Co, Cv, grid_level=None):
    """Return (K_S, K_T) in the pair basis, each (npair, npair)."""
    mol, ni, grids = mf.mol, mf._numint, mf.grids
    if grid_level is not None:
        from pyscf import dft
        grids = dft.gen_grid.Grids(mol); grids.level = grid_level; grids.build(with_non0tab=True)
    xctype = ni._xc_type(mf.xc)
    no, nv = Co.shape[1], Cv.shape[1]; npair = no * nv
    KS = np.zeros((npair, npair)); KT = np.zeros((npair, npair))
    if xctype == "HF":
        return KS, KT
    deriv = 0 if xctype == "LDA" else 1
    nvar = 1 if xctype == "LDA" else 4
    dm_a = mf.make_rdm1() * 0.5
    for ao, mask, weight, coords in ni.block_loop(mol, grids, mol.nao, deriv, max_memory=400):
        rho_a = ni.eval_rho(mol, ao, dm_a, mask, xctype, hermi=1)
        rho = np.stack([rho_a, rho_a])
        if xctype == "LDA":
            rho = rho.reshape(2, 1, -1)
        fxc = ni.eval_xc_eff(mf.xc, rho if xctype != "LDA" else rho[:, 0], spin=1, deriv=2, xctype=xctype)[2]
        fxc = fxc.reshape(2, nvar, 2, nvar, -1)
        fS = fxc[0, :, 0, :] + fxc[0, :, 1, :]
        fT = fxc[0, :, 0, :] - fxc[0, :, 1, :]
        ao = ao.reshape(nvar, -1, mol.nao) if xctype != "LDA" else ao.reshape(1, -1, mol.nao)
        po = np.einsum("kgm,mi->kgi", ao, Co); pv = np.einsum("kgm,ma->kga", ao, Cv)
        ng = po.shape[1]
        M = np.empty((nvar, ng, no, nv))
        M[0] = po[0][:, :, None] * pv[0][:, None, :]
        for x in range(1, nvar):
            M[x] = po[x][:, :, None] * pv[0][:, None, :] + po[0][:, :, None] * pv[x][:, None, :]
        M = M.reshape(nvar, ng, npair)
        for f, K in ((fS, KS), (fT, KT)):
            W = np.einsum("xyg,xgp->ygp", f * weight, M)
            K += M.reshape(nvar * ng, npair).T @ W.reshape(nvar * ng, npair)
    return KS, KT


def window_td(mf, ewin_ev=8.0, nstates=6, rpa_singlets=True, max_occ=60, max_vir=80,
              auxbasis="def2-universal-jkfit", grid_level=None):
    mol = mf.mol
    ni = mf._numint
    omega, alpha, hyb = ni.rsh_and_hybrid_coeff(mf.xc, mol.spin)
    if omega != 0:
        raise NotImplementedError("Only global hybrids are supported (B3LYP, PBE0, BHandHLYP...)")
    nocc = mol.nelectron // 2
    occ, vir = select_window(mf.mo_energy, nocc, ewin_ev, max_occ, max_vir)
    Co, Cv = mf.mo_coeff[:, occ], mf.mo_coeff[:, vir]
    no, nv = len(occ), len(vir); npair = no * nv
    de = (mf.mo_energy[vir][None, :] - mf.mo_energy[occ][:, None]).ravel()
    Lov, Loo, Lvv = _mo_df(mf, Co, Cv, auxbasis)
    J = np.einsum("Lp,Lq->pq", Lov.reshape(-1, npair), Lov.reshape(-1, npair))           # (ia|jb)
    Kx = lib.einsum("Lij,Lab->iajb", Loo, Lvv).reshape(npair, npair)                        # (ij|ab)
    Kb = lib.einsum("Lib,Lja->iajb", Lov, Lov).reshape(npair, npair)                        # (ib|ja)
    KS, KT = _xc_kernels(mf, Co, Cv, grid_level)
    A_S = np.diag(de) + 2 * J - hyb * Kx + KS
    A_T = np.diag(de) - hyb * Kx + KT
    # transition dipoles
    with mol.with_common_orig((0, 0, 0)):
        dip = mol.intor_symmetric("int1e_r", comp=3)
    mu = np.einsum("xpq,pi,qa->xia", dip, Co, Cv).reshape(3, npair)
    out = {"occ": occ, "vir": vir, "npair": npair}
    if rpa_singlets:
        B_S = 2 * J - hyb * Kb + KS
        AmB, ApB = A_S - B_S, A_S + B_S
        w, U = np.linalg.eigh(AmB)
        T = (U * np.sqrt(w)) @ U.T
        Tinv = (U / np.sqrt(w)) @ U.T
        w2, Z = np.linalg.eigh(T @ ApB @ T)
        om = np.sqrt(w2[:nstates]); Z = Z[:, :nstates]
        XpY = T @ Z / np.sqrt(om)
        XmY = (Tinv @ Z) * np.sqrt(om)
        X = 0.5 * (XpY + XmY)
        f = 2 / 3 * om * (np.einsum("xp,pn->xn", mu, XpY) ** 2).sum(0) * 2
        xs = X.T
    else:
        om, Xv = sla.eigh(A_S, subset_by_index=[0, nstates - 1])
        f = 2 / 3 * om * (np.einsum("xp,pn->xn", mu, Xv) ** 2).sum(0) * 2
        xs = Xv.T
    eT, XT = sla.eigh(A_T, subset_by_index=[0, nstates - 1])
    def embed(v):
        full = np.zeros((nocc, mf.mo_coeff.shape[1] - nocc))
        full[np.ix_(occ, vir - nocc)] = v.reshape(no, nv)
        return full / np.linalg.norm(full)
    out.update(e_S=om, f=f, x_S=[embed(v) for v in xs], e_T=eT, x_T=[embed(v) for v in XT.T])
    return out


# ---------------------------------------------------------------------------
# Singlet RPA (TDDFT) in the full (or nearly full) space with our own Davidson solver.
# Exact exchange in the MO basis with DF tensors (cheap), XC kernel through PySCF's AO route
# on a coarse secondary grid. Stratmann-Scuseria-Frisch algorithm.
# ---------------------------------------------------------------------------
def _fxc_singlet_op(mf, Co, Cv, grid_level=1):
    from pyscf import dft
    mol, ni = mf.mol, mf._numint
    g = dft.gen_grid.Grids(mol); g.level = grid_level; g.build(with_non0tab=True)
    rho0, vxc, fxc = ni.cache_xc_kernel(mol, g, mf.xc, mf.mo_coeff, mf.mo_occ, 1)
    fxc = fxc * 0.5  # as in pyscf.scf._response_functions (singlet)

    def op(X):  # X: (nvec, no, nv) -> (nvec, no, nv), devuelve (f_aa+f_ab) X
        dms = np.einsum("pi,nia,qa->npq", Co, X, Cv)
        dms = dms + dms.transpose(0, 2, 1)
        v = ni.nr_rks_fxc_st(mol, g, mf.xc, None, dms, 0, True, rho0, vxc, fxc)
        return np.einsum("pi,npq,qa->nia", Co, v, Cv)
    return op


def rpa_singlets_davidson(mf, nstates=3, ewin_ev=40.0, nfrozen_core=None, grid_level=1,
                          tol=1e-5, max_iter=40, auxbasis="def2-universal-jkfit", verbose=False):
    mol = mf.mol
    omega, alpha, hyb = mf._numint.rsh_and_hybrid_coeff(mf.xc, mol.spin)
    assert omega == 0
    nocc = mol.nelectron // 2
    occ, vir = select_window(mf.mo_energy, nocc, ewin_ev, max_occ=10 ** 6, max_vir=10 ** 6)
    Co, Cv = mf.mo_coeff[:, occ], mf.mo_coeff[:, vir]
    no, nv = len(occ), len(vir); npair = no * nv
    de = (mf.mo_energy[vir][None, :] - mf.mo_energy[occ][:, None])
    Lov, Loo, Lvv = _mo_df(mf, Co, Cv, auxbasis)
    fxc = _fxc_singlet_op(mf, Co, Cv, grid_level)

    def apply(X):  # X (n, no, nv) -> ((A+B)X, (A-B)X)
        LX = Lov.reshape(len(Lov), -1) @ X.reshape(len(X), -1).T              # (L, n)
        J4 = 4 * (Lov.reshape(len(Lov), -1).T @ LX).T.reshape(X.shape)
        K1 = np.empty_like(X); K2 = np.empty_like(X)
        naux = Lov.shape[0]
        blk = max(1, int(1.5e8 // (nv * nv)))   # at most ~1.2 GB per intermediate
        for n in range(len(X)):  # one vector at a time and in blocks of L: the (L,v,v) intermediates are large
            k1 = np.zeros((no, nv)); k2 = np.zeros((no, nv))
            for p0 in range(0, naux, blk):
                sl = slice(p0, p0 + blk)
                k1 += lib.einsum("Lij,Lja->ia", Loo[sl], lib.einsum("Lab,jb->Lja", Lvv[sl], X[n]))   # (ij|ab)X
                k2 += lib.einsum("Lib,Lba->ia", Lov[sl], lib.einsum("jb,Lja->Lba", X[n], Lov[sl]))   # (ib|ja)X
            K1[n] = k1; K2[n] = k2
        F = fxc(X)
        ApB = de * X + J4 + 2 * F - hyb * (K1 + K2)
        AmB = de * X - hyb * (K1 - K2)
        return ApB, AmB

    # initial vectors: smallest orbital energy differences
    nguess = min(npair, max(2 * nstates, nstates + 3))
    idx = np.argsort(de.ravel())[:nguess]
    V = np.zeros((nguess, npair)); V[np.arange(nguess), idx] = 1.0
    P_list, M_list, B_list = [], [], []
    conv = False
    for it in range(max_iter):
        # orthonormalize new vectors against the basis
        Bmat = np.array(B_list) if B_list else np.zeros((0, npair))
        new = []
        for v in V:
            if len(Bmat): v = v - Bmat.T @ (Bmat @ v)
            for w in new: v = v - w * (w @ v)
            nrm = np.linalg.norm(v)
            if nrm > 1e-6: new.append(v / nrm)
        if not new:
            break
        ApB, AmB = apply(np.array(new).reshape(-1, no, nv))
        B_list += new; P_list += list(ApB.reshape(len(new), -1)); M_list += list(AmB.reshape(len(new), -1))
        Bm = np.array(B_list); P = Bm @ np.array(P_list).T; M = Bm @ np.array(M_list).T
        P = 0.5 * (P + P.T); M = 0.5 * (M + M.T)
        w, U = np.linalg.eigh(M)
        T = (U * np.sqrt(np.maximum(w, 1e-12))) @ U.T
        w2, Z = np.linalg.eigh(T @ P @ T)
        om = np.sqrt(np.maximum(w2[:nstates], 0)); Z = Z[:, :nstates]
        xpy = T @ Z / np.sqrt(om); xmy = (np.linalg.solve(T, Z)) * np.sqrt(om)
        XpY = Bm.T @ xpy; XmY = Bm.T @ xmy
        R1 = np.array(P_list).T @ xpy - XmY * om      # (A+B)(X+Y) - w (X-Y)
        R2 = np.array(M_list).T @ xmy - XpY * om      # (A-B)(X-Y) - w (X+Y)
        rn = np.sqrt((R1 ** 2).sum(0) + (R2 ** 2).sum(0))
        if verbose:
            print(f"  rpa it {it} nvec {len(B_list)} w(eV) {np.round(om * 27.2114, 4)} |r| {rn.round(6)}", flush=True)
        if rn.max() < tol:
            conv = True; break
        dd = de.ravel()
        V = []
        for k in range(nstates):
            if rn[k] >= tol:
                V.append(R1[:, k] / (dd - om[k] + 1e-4))
                V.append(R2[:, k] / (dd - om[k] + 1e-4))
    with mol.with_common_orig((0, 0, 0)):
        dip = mol.intor_symmetric("int1e_r", comp=3)
    mu = np.einsum("xpq,pi,qa->xia", dip, Co, Cv).reshape(3, npair)
    f = 2 / 3 * om * (mu @ XpY) ** 2
    f = f.sum(0) * 2
    X = 0.5 * (XpY + XmY)
    xs = []
    for k in range(nstates):
        full = np.zeros((nocc, mf.mo_coeff.shape[1] - nocc))
        full[np.ix_(occ, vir - nocc)] = X[:, k].reshape(no, nv)
        xs.append(full / np.linalg.norm(full))
    return {"e_S": om, "f": f, "x_S": xs, "converged": conv, "npair": npair, "nvec": len(B_list)}
