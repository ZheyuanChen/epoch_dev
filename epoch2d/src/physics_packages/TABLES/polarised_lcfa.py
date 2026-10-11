#!/usr/bin/env python3
"""Spin- and polarisation-resolved LCFA NLCS and NBW in EPOCH's bases.

This derives, with sympy, every final-state quantity the spin- and
polarisation-resolved QED module needs, in EPOCH's conventions:

    NLCS:  w, g, F_0..F_3 (Qian's thesis, Appendix A.1, in our bases)
    NBW:   w_0, j_e, j_p, M, G_0..G_3 (Appendix A.2, in our bases)

starting from Torgrimsson, New J. Phys. 23, 065001 (2021) ("T21"),
Eqs. (7), (43)-(64). Only section 1 below transcribes physics; section 2
states EPOCH's conventions as substitutions; everything after that is
mechanical. Run with --show-input to print the transcribed building
blocks on T21's own axes, for a line-by-line comparison with the paper.

Section 1: T21 building blocks (the only hand transcription)
-----------------------------------------------------------
P = <R> + n_g.R_g + n_1.R_1 + n_0.R_0 + n_g.R_g1.n_1 + n_g.R_g0.n_0
    + n_1.R_01.n_0 + R_g01,ijk n_gi n_1j n_0k                      (7)

NLCS: Eqs. (43)-(50) with s_0 = 1, s_1 = 1 - lam, q_1 = lam.
NBW:  Eqs. (57)-(64) with q_1 = 1, s_2 = 1 - lam_p (electron),
      s_3 = lam_p (positron).

Contraction order. The three-vector blocks are contracted in the order
of their own subscripts, R_g01,ijk n_gi n_0j n_1k and R_g23,ijk n_gi
n_2j n_3k (R_01 and R_23 are symmetric). Eq. (7) labels the tensors
P_10, P_g10 with n_1 before n_0; read literally with P_g10 = R_g01
this gives negative probabilities for pure states (P ~ -P_unpol near
lam = 0.9), so P_g10,ijk = R_g01,ikj. Only the full swap restores
P >= 0; flipping the i sigma_2 term or the eps_2 dyad alone does not
(--check reproduces this). Qian's Appendix A uses the literal order.
S = eps_1 sigma_1^(3) + eps_3 sigma_3^(3), (i sigma_2^(3))_ij =
k_l eps_ijl, 1_par = k k, 1_perp = 1 - k k; Stokes vectors as (1)-(2).

Section 2: EPOCH's conventions (SPIN_DOCUMENTATION.md section 2)
----------------------------------------------------------------
Lepton axes (eps, beta, v): eps = rest-frame E direction (normalised
E_perp + v x B), beta = eps x v (rest-frame B direction), v = velocity.
Photon: xi_obs on the local axes, n_g = xi1 eps_E + xi2 eps_2
+ xi3 eps_EB (xi1 linear par/perp to E_hat, xi2 circular, xi3 diagonal).

    NLCS, emitter charge sign qs:  E_hat = qs eps, B_hat = qs beta,
                                   k_hat = v, n_0, n_1 = physical spins.
    NBW:  E_hat = -eps, B_hat = -beta, k_hat = v,
          n_2 = s_e, n_3 = -s_p (T21's positron vector is the
          v-spinor label).

Spin components are written on (eps, beta, v): si_e = s_i.eps etc.

Bessel form
-----------
With x = (2/3) xi^(3/2) (u' = 2u/(3 chi) for NLCS, u~' = 2u~/(3 chi_g)
for NBW), Ai1 = IntK_1/3/(pi sqrt 3), Ai'/xi = -K_2/3/(pi sqrt 3) and
Ai/sqrt(xi) = K_1/3/(pi sqrt 3). All outputs are divided by the common
positive factor K_2/3/(pi sqrt 3), which leaves

    R1 = K_1/3/K_2/3,   R2 = IntK_1/3/K_2/3   (bessel_k.table),

so every quantity is (rational function of lam) x {1, R1, R2} x
(spin and Stokes components). Ratios such as g/w and F/F_0 are then
exactly what the run-time code needs.

Checks (--check)
----------------
Consistency of the transcription:
  * S from sigma matrices equals Eq. (55); B_hat = E_hat x k_hat.
  * eps_E, eps_EB mean "linear along E_hat" and "along (E_hat+B_hat)/
    sqrt 2" by Eqs. (1)-(2).
  * The probability is independent of the field angle Omega (52) when
    written on local axes.
  * Positivity: sampled pure states on a (lam, x) grid; this is a
    numerical check, not a proof over the continuous state space.
  * Unpolarised spectra equal the textbook (Baier-Katkov) forms.
  * NBW charge-conjugation symmetry (e- <-> e+, lam_p <-> 1 - lam_p,
    field reversed).
One physical result per convention:
  * Sokolov-Ternov: e- polarise along -beta, e+ along +beta, degree
    8/(5 sqrt 3) as chi -> 0 (pins B_hat <-> beta and the charge rule).
  * Synchrotron polarisation: energy-weighted <xi1> -> 3/4 for
    unpolarised emitters (pins E_hat and xi1).
  * Helicity transfer at the endpoints (pins xi2, k_hat and n_3 = -s_p).
Cross-check: Qian's Appendix A (A.17-A.26, A.44-A.60), re-expressed on
his axes (eps_Q = E_hat, beta_Q = B_hat) and labels (xi1_Q = diagonal,
xi3_Q = linear), first against T21 in his reading (Eq. (7) order),
which tests his algebra, then listing what the contraction order
changes. These are reported as NOTEs, not failures.
Independent audit (polarised_lcfa_audit.py, also run by --check):
  * Exact full-kernel equality to Chen et al., PRD 105, 116013 (2022),
    (13)-(15), (22)-(23); magnetic channels against Seipt-King (2020).
  * Choi positivity, exact rank-three factorisation, joint sampling,
    fixed negative witnesses, table/tail positivity and export algebra.
  * Qian A.48/A.54/A.58 agree when transcribed as printed: the energy
    denominator multiplies IntK only. A.50 has a real sign error.

Usage
-----
    python3 polarised_lcfa.py                    # print the expressions
    python3 polarised_lcfa.py --show-input       # T21 blocks on T21 axes
    python3 polarised_lcfa.py --check            # all checks
    python3 polarised_lcfa.py --markdown out.md --fortran out.f90

Requires sympy, numpy, scipy, and gen_bessel_k_table.py (same
directory) for the numerical checks.
"""

import argparse
import os
import re
import sys

import numpy as np
import sympy as sp

# ---------------------------------------------------------------------
# 0. Symbols
# ---------------------------------------------------------------------

# Airy-function values at T21's xi: Ai1(xi), Ai'(xi)/xi, Ai(xi)/sqrt(xi)
AI1, AIPZ, AIS = sp.symbols('Ai1 AipZ AiS')
R1, R2 = sp.symbols('r1 r2', positive=True)
LAM = sp.Symbol('lam', positive=True)       # NLCS photon fraction
LAMP = sp.Symbol('lam_p', positive=True)    # NBW positron fraction
# Display and Fortran variables (T21's longitudinal fractions):
# s1 = 1 - lam (NLCS final lepton); s2 = 1 - lam_p, s3 = lam_p (NBW).
S1, S2, S3 = sp.symbols('s1 s2 s3', positive=True)
QS = sp.Symbol('qs')                        # NLCS emitter charge sign
XI = sp.symbols('xi1 xi2 xi3')              # xi_obs


def spin_symbols(prefix):
    return sp.symbols(' '.join(prefix + c for c in ('_e', '_b', '_v')))


SI, SF = spin_symbols('si'), spin_symbols('sf')     # NLCS initial, final
SE, SP = spin_symbols('se'), spin_symbols('sp')     # NBW electron, positron

TO_BESSEL = {AI1: R2, AIPZ: -1, AIS: R1}

# ---------------------------------------------------------------------
# 1. T21 building blocks (transcription of the paper)
# ---------------------------------------------------------------------

SIG1 = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
SIG3 = sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]])
# S = eps_1 sigma_1^(3) + eps_3 sigma_3^(3): one spin matrix per Stokes
# index (fixed Stokes basis eps_1, eps_2, eps_3).
S_T = [SIG1, sp.zeros(3, 3), SIG3]


def t21_frame(om):
    """T21's local vectors on the fixed axes (e1, e2, e3 = k), Eqs. 51-54."""
    c, s = sp.cos(om), sp.sin(om)
    c2, s2 = sp.cos(2 * om), sp.sin(2 * om)
    return dict(
        E=sp.Matrix([c, s, 0]),             # (52)
        B=sp.Matrix([s, -c, 0]),            # (52), = E x k (51)
        k=sp.Matrix([0, 0, 1]),
        eps_E=sp.Matrix([s2, 0, c2]),       # (53)
        eps_EB=sp.Matrix([-c2, 0, s2]),     # (54)
        eps_2=sp.Matrix([0, 1, 0]))


def i_sigma2(k):
    """(i sigma_2^(3))_ij = k_l eps_ijl."""
    return sp.Matrix(3, 3, lambda i, j: sum(sp.LeviCivita(i, j, l) * k[l]
                                            for l in range(3)))


def sandwich(a, b):
    """a . S . b, a Stokes vector."""
    return sp.Matrix([(a.T * S_T[i] * b)[0] for i in range(3)])


def s_dot(b):
    """S . b, a Stokes-valued spin vector (list over Stokes index)."""
    return [S_T[i] * b for i in range(3)]


def nlcs_blocks(fr, s0, s1, q1):
    """T21 Eqs. (43)-(50)."""
    E, B, k, e2 = fr['E'], fr['B'], fr['k'], fr['eps_2']
    kap = s0 / s1 + s1 / s0
    kapt = s0 / s1 - s1 / s0
    one = sp.eye(3)
    par = k * k.T
    perp = one - par
    ese, bse, sb = sandwich(E, E), sandwich(B, E), s_dot(B)
    jmat = i_sigma2(k)
    return dict(
        R=-AI1 - kap * AIPZ,                                       # (43)
        R0=q1 / s0 * AIS * B,                                      # (44)
        R1=q1 / s1 * AIS * B,                                      # (45)
        Rg=-AIPZ * ese,                                            # (46)
        R01=-AI1 * (perp + (kap - 1) * par)
            - AIPZ * (2 * perp + kap * par),                       # (47)
        Rg0=[-q1 / s1 * AIS * sb[i]
             - q1 / s0 * (AI1 + (1 + s0 / s1) * AIPZ) * e2[i] * k
             for i in range(3)],                                   # (48)
        Rg1=[-q1 / s0 * AIS * sb[i]
             - q1 / s1 * (AI1 + (1 + s1 / s0) * AIPZ) * e2[i] * k
             for i in range(3)],                                   # (49)
        Rg01=[-q1**2 / (2 * s0 * s1) * AI1 * S_T[i]
              - AIPZ * ese[i] * (kap / 2 * perp + par)
              - kapt / 2 * AIPZ * bse[i] * jmat
              + q1 * AIS * e2[i] * (k * B.T / s1 - B * k.T / s0)
              for i in range(3)])                                  # (50)


def nbw_blocks(fr, s2, s3, q1):
    """T21 Eqs. (57)-(64)."""
    E, B, k, e2 = fr['E'], fr['B'], fr['k'], fr['eps_2']
    kap = s2 / s3 + s3 / s2
    kapt = s2 / s3 - s3 / s2
    one = sp.eye(3)
    par = k * k.T
    perp = one - par
    ese, bse, sb = sandwich(E, E), sandwich(B, E), s_dot(B)
    jmat = i_sigma2(k)
    return dict(
        R=AI1 - kap * AIPZ,                                        # (57)
        R2=q1 / s2 * AIS * B,                                      # (58)
        R3=q1 / s3 * AIS * B,                                      # (59)
        Rg=AIPZ * ese,                                             # (60)
        R23=-AI1 * (perp + (kap + 1) * par)
            - AIPZ * (2 * perp + kap * par),                       # (61)
        Rg2=[q1 / s3 * AIS * sb[i]
             + q1 / s2 * (AI1 + (1 - s2 / s3) * AIPZ) * e2[i] * k
             for i in range(3)],                                   # (62)
        Rg3=[q1 / s2 * AIS * sb[i]
             - q1 / s3 * (AI1 + (1 - s3 / s2) * AIPZ) * e2[i] * k
             for i in range(3)],                                   # (63)
        Rg23=[q1**2 / (2 * s2 * s3) * AI1 * S_T[i]
              + AIPZ * ese[i] * (kap / 2 * perp + par)
              + kapt / 2 * AIPZ * bse[i] * jmat
              + q1 * AIS * e2[i] * (k * B.T / s3 - B * k.T / s2)
              for i in range(3)])                                  # (64)


def dot(a, b):
    return (a.T * b)[0]


def contract_v(ng, vlist, n):
    """n_g . V . n for a Stokes-valued spin vector V."""
    return sum(ng[i] * dot(vlist[i], n) for i in range(3))


def contract_m(ng, mlist, na, nb):
    """V_ijk n_gi na_j nb_k for a Stokes-valued spin matrix V."""
    return sum(ng[i] * (na.T * mlist[i] * nb)[0] for i in range(3))


def nlcs_t21(b, n0, n1, ng, block_order=True):
    """Eq. (7) for NLCS; R_g01 contracted n_g n_0 n_1 (block order)."""
    na, nb = (n0, n1) if block_order else (n1, n0)
    return (b['R'] + dot(ng, b['Rg']) + dot(n1, b['R1']) + dot(n0, b['R0'])
            + contract_v(ng, b['Rg1'], n1) + contract_v(ng, b['Rg0'], n0)
            + (na.T * b['R01'] * nb)[0] + contract_m(ng, b['Rg01'], na, nb))


def nbw_t21(b, n2, n3, ng, block_order=True):
    """Eq. (7) for NBW; R_g23 contracted n_g n_2 n_3 (block order)."""
    na, nb = (n2, n3) if block_order else (n3, n2)
    return (b['R'] + dot(ng, b['Rg']) + dot(n2, b['R2']) + dot(n3, b['R3'])
            + contract_v(ng, b['Rg2'], n2) + contract_v(ng, b['Rg3'], n3)
            + (na.T * b['R23'] * nb)[0] + contract_m(ng, b['Rg23'], na, nb))


# ---------------------------------------------------------------------
# 2. EPOCH's conventions (the dictionary)
# ---------------------------------------------------------------------

NBW_POSITRON_SIGN = -1      # n_3 = -s_p
BLOCK_ORDER = True          # contract R_g01 (R_g23) as n_g n_0 n_1 (n_2 n_3)


def vec(comps, axes):
    return sum((c * a for c, a in zip(comps, axes)), sp.zeros(3, 1))


def lepton_axes(fr, e_sign):
    """(eps, beta, v) with E_hat = e_sign * eps."""
    eps = e_sign * fr['E']
    v = fr['k']
    return eps, eps.cross(v), v


def photon_stokes(xi, fr):
    return xi[0] * fr['eps_E'] + xi[1] * fr['eps_2'] + xi[2] * fr['eps_EB']


def reduce_sign(expr, sym=QS):
    """Use sym**2 = 1."""
    expr = sp.expand(expr)
    return expr.replace(lambda e: e.is_Pow and e.base == sym,
                        lambda e: sym**(e.exp % 2))


def nlcs_probability(qs=QS, om=0, lam=LAM, axes_sign=None,
                     block_order=BLOCK_ORDER, bessel=True):
    """P(si, sf, xi_obs) for NLCS by a lepton of charge sign qs.

    axes_sign overrides E_hat = axes_sign * eps (used for Qian's axes).
    """
    fr = t21_frame(om)
    b = nlcs_blocks(fr, sp.Integer(1), 1 - lam, lam)
    axes = lepton_axes(fr, qs if axes_sign is None else axes_sign)
    p = nlcs_t21(b, vec(SI, axes), vec(SF, axes), photon_stokes(XI, fr),
                 block_order)
    if bessel:
        p = p.subs(TO_BESSEL)
    return reduce_sign(p)


def nbw_probability(om=0, lam_p=LAMP, axes_sign=-1,
                    positron_sign=NBW_POSITRON_SIGN,
                    block_order=BLOCK_ORDER, bessel=True):
    """P(se, sp, xi_obs) for NBW; spins physical, components on axes."""
    fr = t21_frame(om)
    b = nbw_blocks(fr, 1 - lam_p, lam_p, sp.Integer(1))
    axes = lepton_axes(fr, axes_sign)
    n2 = vec(SE, axes)
    n3 = positron_sign * vec(SP, axes)
    p = nbw_t21(b, n2, n3, photon_stokes(XI, fr), block_order)
    if bessel:
        p = p.subs(TO_BESSEL)
    return sp.expand(p)


# ---------------------------------------------------------------------
# 3. Derived quantities
# ---------------------------------------------------------------------

def zero(syms):
    return {s: 0 for s in syms}


def nlcs_quantities(p):
    """Thesis A.1 quantities from P(si, sf, xi)."""
    f0 = p.subs(zero(XI))
    f = [sp.diff(p, x) for x in XI]
    w = f0.subs(zero(SF))
    q = dict(w=w)
    for c, s in zip('ebv', SF):
        q['g_' + c] = sp.diff(f0, s)
    q['F0'] = f0
    for j in range(3):
        q['F%d' % (j + 1)] = f[j]
    # Photon given the initial spin only (final spin summed).
    for j in range(3):
        q['F%d_bar' % (j + 1)] = f[j].subs(zero(SF))
    return q


def nbw_quantities(p):
    """Thesis A.2 quantities from P(se, sp, xi), plus the pair matrix M."""
    g0 = p.subs(zero(XI))
    q = dict(w0=p.subs(zero(SE + SP)))
    for c, s in zip('ebv', SE):
        q['je_' + c] = sp.diff(p, s).subs(zero(SP))
    for c, s in zip('ebv', SP):
        q['jp_' + c] = sp.diff(p, s).subs(zero(SE))
    for a, sa in zip('ebv', SP):
        for c, sc in zip('ebv', SE):
            q['m_' + a + c] = sp.diff(p, sa, sc)
    q['G0'] = g0
    for j in range(3):
        q['G%d' % (j + 1)] = sp.diff(p, XI[j])
    return q


def tidy_coeff(c, display=False):
    """Collect a coefficient as A + B r1 + C r2, each part factorised.

    With display=True, write lam = 1 - s1 and lam_p = s3, 1 - lam_p = s2.
    """
    c = sp.expand(c)
    if display:
        c = sp.expand(c.subs({LAM: 1 - S1, LAMP: S3}))
    parts = sp.collect(c, [R1, R2], evaluate=False)
    out = sp.Integer(0)
    for k, v in parts.items():
        v = sp.factor(sp.cancel(v))
        if display:
            v = v.subs({S3 - 1: -S2, S1 - 1: -LAM})
            v = v.subs({2 * S3**2 - 2 * S3 + 1: S2**2 + S3**2,
                        2 * S3 - 1: S3 - S2})
        out += v * k
    return out


def terms(expr, gens, display=False):
    """[(monomial, tidy coefficient)] of expr in the spin/Stokes gens."""
    poly = sp.Poly(sp.expand(expr), *gens)
    out = []
    for monom, coeff in poly.terms():
        if coeff == 0:
            continue
        mono = sp.Mul(*[g**e for g, e in zip(gens, monom)])
        out.append((mono, tidy_coeff(coeff, display)))
    return out


def tidy(expr, gens, display=False):
    return sum((m * c for m, c in terms(expr, gens, display)),
               sp.Integer(0))


NLCS_GENS = list(XI) + list(SI) + list(SF)
NBW_GENS = list(XI) + list(SE) + list(SP)


def derive():
    nlcs = nlcs_quantities(nlcs_probability())
    nbw = nbw_quantities(nbw_probability())
    return ({k: tidy(v, NLCS_GENS) for k, v in nlcs.items()},
            {k: tidy(v, NBW_GENS) for k, v in nbw.items()})


# ---------------------------------------------------------------------
# 4. Output
# ---------------------------------------------------------------------

NLCS_ORDER = ['w', 'g_e', 'g_b', 'g_v', 'F1_bar', 'F2_bar', 'F3_bar',
              'F0', 'F1', 'F2', 'F3']
NBW_ORDER = ['w0', 'je_e', 'je_b', 'je_v', 'jp_e', 'jp_b', 'jp_v'] \
    + ['m_' + a + c for a in 'ebv' for c in 'ebv'] \
    + ['G0', 'G1', 'G2', 'G3']


def print_text(nlcs, nbw):
    print(__doc__.split('\n')[0])
    print('All quantities divided by K_2/3/(pi sqrt 3); r1 = K_1/3/K_2/3,'
          ' r2 = IntK_1/3/K_2/3.')
    print('\nNLCS (s1 = 1 - lam, lam = photon fraction; qs = emitter '
          'charge sign, e- = -1):')
    for name in NLCS_ORDER:
        print('\n  %s =' % name)
        for mono, coeff in terms(nlcs[name], NLCS_GENS, True) or [(0, '')]:
            print('    %-14s %s' % (mono, coeff))
    print('\nNBW (s3 = lam_p = positron fraction, s2 = 1 - lam_p; '
          'M: P contains sp.M.se):')
    for name in NBW_ORDER:
        print('\n  %s =' % name)
        for mono, coeff in terms(nbw[name], NBW_GENS, True) or [(0, '')]:
            print('    %-14s %s' % (mono, coeff))


LATEX_NAMES = {LAM: r'\lambda', LAMP: r'\tilde\lambda', QS: 'q',
               R1: 'R_1', R2: 'R_2', S1: 's_1', S2: 's_2', S3: 's_3'}
for _x, _n in zip(XI, ('1', '2', '3')):
    LATEX_NAMES[_x] = r'\xi_%s' % _n
for _syms, _p in ((SI, 'i'), (SF, 'f'), (SE, 'e'), (SP, 'p')):
    for _s, _c in zip(_syms, (r'\varepsilon', r'\beta', 'v')):
        LATEX_NAMES[_s] = r's_{%s,%s}' % (_p, _c)


def write_markdown(path, nlcs, nbw):
    axis_names = {'e': r'\varepsilon', 'b': r'\beta', 'v': 'v'}

    def label(name):
        if name.startswith('m_'):
            return r'M_{%s,%s}' % tuple(axis_names[c] for c in name[2:])
        if name.startswith(('g_', 'je_', 'jp_')):
            prefix, axis = name.split('_')
            prefix = {'g': 'g', 'je': r'j_e', 'jp': r'j_p'}[prefix]
            return r'(%s)_{%s}' % (prefix, axis_names[axis])
        if name.startswith(('F', 'G')):
            return (r'\bar F_%s' if name.endswith('_bar') else
                    name[0] + '_%s') % name[1]
        return 'w_0' if name == 'w0' else name

    def block(fh, name, expr, gens):
        rows = terms(expr, gens, True)
        fh.write('$$\n\\begin{aligned}\n%s &= ' % label(name))
        fh.write(' \\\\\n  &{}+ '.join(
            r'%s\,\left[%s\right]' % (
                sp.latex(m, symbol_names=LATEX_NAMES),
                sp.latex(c, symbol_names=LATEX_NAMES)) for m, c in rows)
            or '0')
        fh.write('\n\\end{aligned}\n$$\n\n')

    with open(path, 'w') as fh:
        fh.write('# Polarised LCFA expressions in EPOCH bases\n\n')
        fh.write('Generated by `polarised_lcfa.py`; do not edit by hand.\n'
                 'All quantities are divided by K_{2/3}/(pi sqrt 3); '
                 'R_1 = K_{1/3}/K_{2/3}, R_2 = IntK_{1/3}/K_{2/3}.\n'
                 'Spin components are on (eps, beta, v); xi = xi_obs.\n'
                 'NLCS: s_1 = 1 - lambda (lambda = photon fraction). '
                 'NBW: s_3 = positron fraction, s_2 = 1 - s_3.\n\n')
        fh.write('## NLCS (q = emitter charge sign)\n\n')
        for name in NLCS_ORDER:
            block(fh, name, nlcs[name], NLCS_GENS)
        fh.write('## NBW (P contains s_p . M . s_e)\n\n')
        for name in NBW_ORDER:
            block(fh, name, nbw[name], NBW_GENS)
    print('wrote %s' % path)


def fortran_lines(name, expr, width=79, indent=2):
    """Free-form Fortran assignment in EPOCH style (num kind, <= width).

    The printer's own line wrapping is bypassed (it can split tokens);
    lines are broken only at the spaces the printer puts around + and -.
    """
    from sympy.printing.fortran import FCodePrinter
    printer = FCodePrinter({'source_format': 'free', 'standard': 95})
    code = name + ' = ' + printer._print(expr)
    code = re.sub(r'(\d+)\.0d0', r'\1.0_num', code)
    code = code.replace('1.0_num/2.0_num', '0.5_num')
    words = code.split()
    lines, cur = [], ' ' * indent
    for word in words:
        if len(cur) + len(word) + 3 > width and cur.strip():
            lines.append(cur.rstrip() + ' &')
            cur = ' ' * (indent + 4)
        cur += word + ' '
    lines.append(cur.rstrip())
    return lines


def write_fortran(path, nlcs, nbw):
    with open(path, 'w') as fh:
        fh.write('! Generated by polarised_lcfa.py: reference expressions.\n'
                 '! Not yet wired into EPOCH.\n'
                 '! Inputs: s1 = 1 - lam (NLCS, lam = photon fraction); '
                 's3 = lam_p,\n! s2 = 1 - lam_p (NBW, lam_p = positron '
                 'fraction); r1, r2 (bessel_k.table);\n! qs (NLCS emitter '
                 'charge sign); spin components on (eps, beta, v);\n'
                 '! xi1..3 = xi_obs.\n'
                 '! All quantities are divided by K_2/3/(pi sqrt 3).\n\n')
        fh.write('  ! NLCS\n')
        for name in NLCS_ORDER:
            fh.write('\n'.join(fortran_lines(
                name.lower(), tidy(nlcs[name], NLCS_GENS, True))) + '\n')
        fh.write('\n  ! NBW\n')
        for name in NBW_ORDER:
            fh.write('\n'.join(fortran_lines(
                name.lower(), tidy(nbw[name], NBW_GENS, True))) + '\n')
    print('wrote %s' % path)


NLCS_KINDS = dict(R='scalar', R0='spin', R1='spin', Rg='stokes',
                  R01='spin2', Rg0='stokes_spin', Rg1='stokes_spin',
                  Rg01='stokes_spin2')
NBW_KINDS = dict(R='scalar', R2='spin', R3='spin', Rg='stokes',
                 R23='spin2', Rg2='stokes_spin', Rg3='stokes_spin',
                 Rg23='stokes_spin2')


def describe(obj, kind, fr):
    """Coefficients of a block on (E, B, k) and (eps_E, eps_2, eps_EB)."""
    axes = (('E', fr['E']), ('B', fr['B']), ('k', fr['k']))
    stokes = (('eps_E', fr['eps_E']), ('eps_2', fr['eps_2']),
              ('eps_EB', fr['eps_EB']))
    out = []

    def add(label, c):
        c = sp.factor(sp.simplify(c))
        if c != 0:
            out.append('%-18s %s' % (label, c))

    if kind == 'scalar':
        add('1', obj)
    elif kind == 'spin':
        for an, a in axes:
            add(an, dot(a, obj))
    elif kind == 'spin2':
        for an, a in axes:
            for bn, bv in axes:
                add(an + bn, (a.T * obj * bv)[0])
    elif kind == 'stokes':
        for sn, sv in stokes:
            add(sn, dot(sv, obj))
    elif kind == 'stokes_spin':
        for sn, sv in stokes:
            for an, a in axes:
                add(sn + ' x ' + an,
                    sum(sv[i] * dot(a, obj[i]) for i in range(3)))
    elif kind == 'stokes_spin2':
        for sn, sv in stokes:
            for an, a in axes:
                for bn, bv in axes:
                    add(sn + ' x ' + an + bn,
                        sum(sv[i] * (a.T * obj[i] * bv)[0]
                            for i in range(3)))
    return out


def show_input():
    s0, s1, q1, s2, s3 = sp.symbols('s0 s1 q1 s2 s3', positive=True)
    fr = t21_frame(0)
    print('T21 building blocks on T21 axes (Omega = 0). Spin matrices: '
          'EB means E_hat.M.B_hat.\nStokes-valued terms: "eps_E x B" is '
          'the coefficient of eps_E (x) B_hat. Use (55): S.B = -eps_E B +'
          ' eps_EB E,\nE.S.E = eps_E, B.S.E = eps_EB, and 1_perp = EE + '
          'BB. Ai1, AipZ = Ai\'/xi, AiS = Ai/sqrt(xi).')
    for title, b, kinds in (
            ('NLCS, Eqs. (43)-(50)', nlcs_blocks(fr, s0, s1, q1),
             NLCS_KINDS),
            ('NBW, Eqs. (57)-(64)', nbw_blocks(fr, s2, s3, q1), NBW_KINDS)):
        print('\n' + title)
        for name, kind in kinds.items():
            print('  %s:' % name)
            for line in describe(b[name], kind, fr):
                print('    ' + line)


# ---------------------------------------------------------------------
# 5. Checks
# ---------------------------------------------------------------------

RESULTS = []


def report(name, ok, detail='', expect=True):
    """Record a check; expect=False marks a diagnostic that must fail."""
    passed = ok == expect
    RESULTS.append((name, passed))
    print('  [%s] %s%s' % ('PASS' if passed else 'FAIL', name,
                          ('  (' + detail + ')') if detail else ''))


def note(name, same, detail=''):
    """A cross-check against another source: never a failure."""
    print('  [%s] %s%s' % ('SAME' if same else 'DIFF', name,
                          ('  (' + detail + ')') if detail else ''))


def bessel_ratios(x):
    """R1, R2 at x > 0, with the tails documented in gen_bessel_k_table."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    sys.dont_write_bytecode = True
    import gen_bessel_k_table as gk
    from scipy.special import kve
    if x < 1e-6:
        return gk.tail_small(x)
    if x > 10.0**2.5:
        return gk.tail_large(x)
    k23 = kve(2.0 / 3.0, x)
    return kve(1.0 / 3.0, x) / k23, gk.intk_scaled(x) / k23


def k23_abs(x):
    from scipy.special import kve
    return kve(2.0 / 3.0, x) * np.exp(-x)


def random_unit(rng, n):
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1)[:, None]


def check_frame():
    print('\nTranscription and frame:')
    om = sp.Symbol('om', real=True)
    fr = t21_frame(om)
    e, bv = fr['E'], fr['B']
    s55 = [fr['eps_E'][i] * (e * e.T - bv * bv.T)
           + fr['eps_EB'][i] * (e * bv.T + bv * e.T) for i in range(3)]
    ok = all(sp.simplify(s55[i] - S_T[i]) == sp.zeros(3, 3)
             for i in range(3))
    report('S = eps1 sigma1 + eps3 sigma3 equals Eq. (55)', ok)
    report('B_hat = E_hat x k_hat (51)',
           sp.simplify(bv - e.cross(fr['k'])) == sp.zeros(3, 1))

    # Stokes vector of real linear polarisation at angle t, Eqs. (1)-(2):
    # rho/2 = t, lambda = 0 -> n = (sin 2t, 0, cos 2t).
    def stokes_linear(t):
        return np.array([np.sin(2 * t), 0.0, np.cos(2 * t)])
    worst = 0.0
    for omv in np.linspace(-3, 3, 13):
        f = {k: np.array(v.subs(om, omv), dtype=float).ravel()
             for k, v in fr.items()}
        d = (f['E'] + f['B']) / np.sqrt(2.0)
        worst = max(worst,
                    np.abs(stokes_linear(np.arctan2(f['E'][1], f['E'][0]))
                           - f['eps_E']).max(),
                    np.abs(stokes_linear(np.arctan2(f['B'][1], f['B'][0]))
                           + f['eps_E']).max(),
                    np.abs(stokes_linear(np.arctan2(d[1], d[0]))
                           - f['eps_EB']).max())
    report('xi1_obs = +1 along E_hat, -1 along B_hat; xi3_obs = +1 along '
           '(E_hat+B_hat)/sqrt2', worst < 1e-12, 'max error %.1e' % worst)

    rng = np.random.default_rng(1)
    worst = 0.0
    for label, fun, syms in (
            ('NLCS e-', lambda o: nlcs_probability(-1, o), NLCS_GENS),
            ('NLCS e+', lambda o: nlcs_probability(+1, o), NLCS_GENS),
            ('NBW', lambda o: nbw_probability(o), NBW_GENS)):
        lam = LAM if 'NLCS' in label else LAMP
        p0 = fun(sp.Integer(0))
        for _ in range(3):
            omv = sp.Float(rng.uniform(-3, 3))
            vals = {s: rng.uniform(-1, 1) for s in syms}
            vals.update({lam: rng.uniform(0.05, 0.95),
                         R1: rng.uniform(0, 1), R2: rng.uniform(0, 1)})
            worst = max(worst, abs(float(fun(omv).subs(vals)
                                         - p0.subs(vals))))
    report('P on local axes independent of the field angle Omega',
           worst < 1e-10, 'max difference %.1e' % worst)


def positivity(p, lam, gens_a, gens_b, label, rng, expect=True):
    """min over pure states of P = F0 - |F| (min over xi analytic)."""
    f0 = p.subs(zero(XI))
    f = [sp.diff(p, x) for x in XI]
    args = [lam, R1, R2] + list(gens_a) + list(gens_b)
    lf0 = sp.lambdify(args, f0, 'numpy')
    lf = [sp.lambdify(args, fi, 'numpy') for fi in f]
    xs = np.geomspace(1e-4, 1e3, 25)
    rr = np.array([bessel_ratios(x) for x in xs])
    lams = np.linspace(0.005, 0.995, 25)
    a = random_unit(rng, 400)
    b = random_unit(rng, 400)
    worst = np.inf
    for lv in lams:
        for r1v, r2v in rr:
            cols = [lv, r1v, r2v] + [a[:, i] for i in range(3)] \
                + [b[:, i] for i in range(3)]
            v0 = lf0(*cols) * np.ones(len(a))
            vf = np.sqrt(sum((lfi(*cols) * np.ones(len(a)))**2
                             for lfi in lf))
            # Scale by the unpolarised value so the test is relative.
            scale = abs(float(f0.subs(zero(list(gens_a) + list(gens_b)))
                              .subs({lam: lv, R1: r1v, R2: r2v})))
            worst = min(worst, np.min(v0 - vf) / scale)
    report('positivity, %s: min P / P_unpol over pure states' % label,
           worst > -1e-12, '%.3f' % worst, expect)
    return worst


def check_physics(nlcs_e, nlcs_p, nbw):
    print('\nPositivity and textbook limits:')
    rng = np.random.default_rng(2)
    positivity(nlcs_probability(-1), LAM, SI, SF, 'NLCS e-', rng)
    positivity(nlcs_probability(+1), LAM, SI, SF, 'NLCS e+', rng)
    positivity(nbw_probability(), LAMP, SE, SP, 'NBW', rng)
    positivity(nlcs_probability(-1, block_order=False), LAM, SI, SF,
               'NLCS e-, Eq. (7) order (must fail)', rng, expect=False)
    positivity(nbw_probability(block_order=False), LAMP, SE, SP,
               'NBW, Eq. (7) order (must fail)', rng, expect=False)

    kap = 1 / (1 - LAM) + (1 - LAM)
    rho = (1 - LAMP) / LAMP + LAMP / (1 - LAMP)
    w_unpol = nlcs_e['w'].subs(zero(SI))
    report('NLCS unpolarised spectrum = kappa K_2/3 - IntK_1/3 '
           '(Baier-Katkov)', sp.cancel(w_unpol - (kap - R2)) == 0)
    report('NLCS rate depends on the initial spin only through si_b',
           not (nlcs_e['w'].free_symbols & {SI[0], SI[2]}))
    w0_unpol = nbw['w0'].subs(zero(XI))
    report('NBW unpolarised spectrum = rho K_2/3 + IntK_1/3 '
           '(Baier-Katkov)', sp.cancel(w0_unpol - (rho + R2)) == 0)
    report('NBW rate depends on the photon only through xi1_obs',
           not (nbw['w0'].free_symbols & {XI[1], XI[2]}))

    print('\nNBW charge conjugation (e- <-> e+, lam_p <-> 1 - lam_p, '
          'field reversed):')
    a, b = SE, SP
    swap = {SE[0]: -b[0], SE[1]: -b[1], SE[2]: b[2],
            SP[0]: -a[0], SP[1]: -a[1], SP[2]: a[2], LAMP: 1 - LAMP}
    for sign in (-1, +1):
        for order in (True, False):
            p = nbw_probability(positron_sign=sign, block_order=order)
            ok = sp.cancel(p - p.subs(swap, simultaneous=True)) == 0
            adopted = sign == NBW_POSITRON_SIGN and order == BLOCK_ORDER
            report('n3 = %ss_p, %s order%s' % (
                '-' if sign < 0 else '+', 'block' if order else 'Eq. (7)',
                ' (adopted)' if adopted else
                (' (must fail)' if sign > 0 else '')), ok,
                expect=sign < 0)

    print('\nPhysical results that pin the conventions:')
    from scipy.integrate import quad

    def spectrum_integral(expr, lam, chi, weight=lambda l: 1.0):
        """int dlam K_2/3(x) expr(lam, R1, R2), x = 2u/(3 chi)."""
        f = sp.lambdify((lam, R1, R2), expr, 'numpy')

        def integrand(u):
            lv = u / (1.0 + u)
            x = 2.0 * u / (3.0 * chi)
            r1v, r2v = bessel_ratios(x)
            return (f(lv, r1v, r2v) * k23_abs(x) * weight(lv)
                    / (1.0 + u)**2)
        edges = [0.0] + [chi * 10.0**e for e in range(-3, 4)] + [np.inf]
        return sum(quad(integrand, lo, hi, limit=200, epsrel=1e-10)[0]
                   for lo, hi in zip(edges[:-1], edges[1:]))

    chi = 1e-3
    st = 8.0 / (5.0 * np.sqrt(3.0))
    for label, q, expect in (('e-', nlcs_e, -st), ('e+', nlcs_p, +st)):
        f0 = q['F0']

        def rate(si_b, sf_b):
            return spectrum_integral(
                f0.subs({SI[0]: 0, SI[1]: si_b, SI[2]: 0,
                         SF[0]: 0, SF[1]: sf_b, SF[2]: 0}), LAM, chi)
        up_dn, dn_up = rate(1, -1), rate(-1, 1)
        pol = (dn_up - up_dn) / (dn_up + up_dn)
        report('Sokolov-Ternov %s: equilibrium s.beta = %+.4f, expect '
               '%+.4f (chi = %g)' % (label, pol, expect, chi),
               abs(pol - expect) < 1e-2)

    num = spectrum_integral(nlcs_e['F1_bar'].subs(zero(SI)), LAM, chi,
                            weight=lambda l: l)
    den = spectrum_integral(nlcs_e['w'].subs(zero(SI)), LAM, chi,
                            weight=lambda l: l)
    report('synchrotron: energy-weighted <xi1_obs> = %.4f, expect +0.75'
           % (num / den), abs(num / den - 0.75) < 1e-2)

    def at(expr, lam, lv, chi_v, subs):
        x_of = {LAM: lambda l: 2 * l / (1 - l) / (3 * chi_v),
                LAMP: lambda l: 2 * (1 / l + 1 / (1 - l)) / (3 * chi_v)}
        r1v, r2v = bessel_ratios(x_of[lam](lv))
        return float(expr.subs(subs).subs({lam: lv, R1: r1v, R2: r2v}))

    for label, q in (('e-', nlcs_e), ('e+', nlcs_p)):
        hel = {SI[0]: 0, SI[1]: 0, SI[2]: 1}
        val = at(q['F2_bar'], LAM, 1 - 1e-6, 1.0, hel) \
            / at(q['w'], LAM, 1 - 1e-6, 1.0, hel)
        report('helicity transfer NLCS %s: s = +v, lam -> 1 gives '
               '<xi2> = %+.4f, expect +1' % (label, val),
               abs(val - 1) < 1e-3)
    circ = {XI[0]: 0, XI[1]: 1, XI[2]: 0}
    val_e = at(nbw['je_v'], LAMP, 1e-6, 1.0, circ) \
        / at(nbw['w0'], LAMP, 1e-6, 1.0, circ)
    val_p = at(nbw['jp_v'], LAMP, 1 - 1e-6, 1.0, circ) \
        / at(nbw['w0'], LAMP, 1 - 1e-6, 1.0, circ)
    report('helicity transfer NBW: xi2 = +1, hard electron <s_e.v> = %+.4f'
           ', hard positron <s_p.v> = %+.4f, expect +1' % (val_e, val_p),
           abs(val_e - 1) < 1e-3 and abs(val_p - 1) < 1e-3)


def qian_nlcs(fr):
    """Thesis A.17, A.18, A.21, A.23, A.25 on his axes (Bessel form)."""
    k13, k23, ik = R1, sp.Integer(1), R2
    lam = LAM
    u = lam / (1 - lam)
    kap = 1 / (1 - lam) + (1 - lam)
    kapt = 1 / (1 - lam) - (1 - lam)
    eq, bq, kq = fr['E'], fr['B'], fr['k']
    si, sf = vec(SI, (eq, bq, kq)), vec(SF, (eq, bq, kq))
    w = -ik + kap * k23 + lam * dot(si, bq) * k13                 # A.17
    g = (u * k13 * bq + (-ik + 2 * k23) * si
         + lam * u * (-ik + k23) * dot(si, kq) * kq)              # A.18
    f1 = (-k13 * (u * dot(si, eq) + lam * dot(sf, eq))
          - lam * u / 2 * ik * (dot(sf, eq) * dot(si, bq)
                                + dot(sf, bq) * dot(si, eq))
          + kapt / 2 * k23 * dot(sf.cross(si), kq))               # A.21
    f2 = (kapt / 2 * k13 * dot(sf.cross(si), eq)
          - (lam * ik - kapt * k23) * dot(si, kq)
          - (u * ik - kapt * k23) * dot(sf, kq)
          + lam * u / 2 * k13 * (dot(sf, kq) * dot(si, bq)
                                 + dot(sf, bq) * dot(si, kq)))    # A.23
    f3 = (k23 + (u * dot(si, bq) + lam * dot(sf, bq)) * k13
          - lam * u / 2 * ik * (dot(sf, eq) * dot(si, eq)
                                - dot(sf, bq) * dot(si, bq))
          + kap / 2 * k23 * dot(sf, si)
          - lam * u / 2 * k23 * dot(sf, kq) * dot(si, kq))        # A.25
    return dict(w=w, g_e=dot(g, eq), g_b=dot(g, bq), g_v=dot(g, kq),
                F_diag=f1, F_circ=f2, F_lin=f3)


def qian_nbw(fr):
    """Thesis A.44, A.46, A.48, A.50, A.53, A.54, A.58 on his axes."""
    k13, k23, ik = R1, sp.Integer(1), R2
    lp = LAMP
    ut = 1 / lp + 1 / (1 - lp)
    rho = (1 - lp) / lp + lp / (1 - lp)
    rhot = (1 - lp) / lp - lp / (1 - lp)
    eq, bq, kq = fr['E'], fr['B'], fr['k']
    se, spv = vec(SE, (eq, bq, kq)), vec(SP, (eq, bq, kq))
    xd, xc, xl = XI[2], XI[1], XI[0]    # Qian's xi1, xi2, xi3
    g0 = (ik + rho * k23 + dot(se, bq) * k13 / (1 - lp)
          - dot(spv, bq) * k13 / lp + (ik - 2 * k23) * dot(spv, se)
          + (rho * ik - (rho - 2) * k23) * dot(spv, kq) * dot(se, kq))
    g1 = (k13 * (dot(se, eq) / lp - dot(spv, eq) / (1 - lp))
          - ut / 2 * ik * (dot(spv, eq) * dot(se, bq)
                           + dot(spv, bq) * dot(se, eq))
          + rhot / 2 * k23 * dot(spv.cross(se), kq))
    # A.48: the denominators multiply IntK only, not the rhot K term.
    # Verified visually on thesis p. 194 (PDF page 216).
    g2 = (-ut / 2 * k13 * dot(spv.cross(se), eq)
          + (ik / (1 - lp) + rhot * k23) * dot(se, kq)
          + (ik / lp - rhot * k23) * dot(spv, kq)
          - rhot / 2 * k13 * (dot(spv, kq) * dot(se, bq)
                              + dot(spv, bq) * dot(se, kq)))
    g3 = (-k23 + (-dot(se, bq) / lp + dot(spv, bq) / (1 - lp)) * k13
          - ut / 2 * ik * (dot(spv, eq) * dot(se, eq)
                           - dot(spv, bq) * dot(se, bq))
          + rho / 2 * k23 * dot(spv, se)
          - (1 - rho / 2) * k23 * dot(spv, kq) * dot(se, kq))
    w0 = ik + rho * k23 - xl * k23
    je = ((1 / (1 - lp) - xl / lp) * k13 * bq + xd * k13 / lp * eq
          + xc * (ik / (1 - lp) + rhot * k23) * kq)
    jp = ((-1 / lp + xl / (1 - lp)) * k13 * bq - xd * k13 / (1 - lp) * eq
          + xc * (ik / lp - rhot * k23) * kq)
    out = dict(G0=g0, G_diag=g1, G_circ=g2, G_lin=g3, w0=w0)
    for c, a in zip('ebv', (eq, bq, kq)):
        out['je_' + c] = dot(je, a)
        out['jp_' + c] = dot(jp, a)
    return out


def check_qian():
    print('\nCross-check with Qian\'s Appendix A (his axes eps_Q = E_hat, '
          'beta_Q = B_hat;\nhis F1/G1 = diagonal = our xi3_obs, F3/G3 = '
          'linear = our xi1_obs).')
    fr = t21_frame(0)
    nlcs_eq = dict(w='17', g_e='18', g_b='18', g_v='18', F_lin='25',
                   F_circ='23', F_diag='21')
    nbw_eq = dict(G0='44', G_lin='50', G_circ='48', G_diag='46', w0='53',
                  je_e='54', je_b='54', je_v='54', jp_e='58', jp_b='58',
                  jp_v='58')

    def relabel(q, f):
        q.update({f + '_diag': q[f + '3'], f + '_circ': q[f + '2'],
                  f + '_lin': q[f + '1']})
        return q

    def nlcs_on_his_axes(order):
        return relabel(nlcs_quantities(
            nlcs_probability(axes_sign=+1, block_order=order)), 'F')

    def nbw_on_his_axes(order):
        return relabel(nbw_quantities(
            nbw_probability(axes_sign=+1, block_order=order)), 'G')

    def compare(title, a, b, eqno, gens):
        print(title)
        for name in eqno:
            diff = tidy(sp.cancel(sp.expand(a[name] - b[name])), gens)
            note('%-6s (A.%s)' % (name, eqno[name]), diff == 0,
                 '' if diff == 0 else 'difference = %s' % diff)

    compare('T21 in Qian\'s reading (Eq. (7) order) minus Qian, NLCS:',
            nlcs_on_his_axes(False), qian_nlcs(fr), nlcs_eq, NLCS_GENS)
    compare('T21 in Qian\'s reading (Eq. (7) order) minus Qian, NBW:',
            nbw_on_his_axes(False), qian_nbw(fr), nbw_eq, NBW_GENS)
    compare('Effect of the contraction order (block minus Eq. (7)), NLCS:',
            nlcs_on_his_axes(True), nlcs_on_his_axes(False), nlcs_eq,
            NLCS_GENS)
    compare('Effect of the contraction order (block minus Eq. (7)), NBW:',
            nbw_on_his_axes(True), nbw_on_his_axes(False), nbw_eq,
            NBW_GENS)


def run_checks(nlcs, nbw):
    RESULTS.clear()
    check_frame()
    nlcs_e = {k: reduce_sign(v.subs(QS, -1)) for k, v in nlcs.items()}
    nlcs_p = {k: reduce_sign(v.subs(QS, 1)) for k, v in nlcs.items()}
    check_physics(nlcs_e, nlcs_p, nbw)
    check_qian()
    from polarised_lcfa_audit import run_audit
    run_audit(report)
    failed = [n for n, ok in RESULTS if not ok]
    print('\n%d checks, %d failed' % (len(RESULTS), len(failed)))
    for n in failed:
        print('  FAILED: ' + n)
    return not failed


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--show-input', action='store_true',
                   help='print the T21 building blocks on T21 axes')
    p.add_argument('--check', action='store_true', help='run all checks')
    p.add_argument('--markdown', help='write the expressions as Markdown')
    p.add_argument('--fortran', help='write reference Fortran assignments')
    p.add_argument('--quiet', action='store_true',
                   help='do not print the expressions')
    args = p.parse_args()

    if args.show_input:
        show_input()
        return 0
    nlcs, nbw = derive()
    if not args.quiet:
        print_text(nlcs, nbw)
    if args.markdown:
        write_markdown(args.markdown, nlcs, nbw)
    if args.fortran:
        write_fortran(args.fortran, nlcs, nbw)
    if args.check:
        return 0 if run_checks(nlcs, nbw) else 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
