#!/usr/bin/env python3
"""Independent, reproducible audit of the polarised LCFA reference.

Run with system python3, directly or through polarised_lcfa.py --check.
Chen's scalar formulae below are transcribed independently of T21's
tensor blocks. See POLARISED_LCFA_AUDIT.md for sources, conventions,
normalisation, numerical coverage and limitations.
"""

from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import kve

import polarised_lcfa as lcfa


def dot(a, b):
    return (a.T * b)[0]


def physical_vectors(symbols):
    # Cartesian axes: eps = x, beta = -y, v = z. The physical spin
    # component triad is LEFT-handed, so beta's Pauli matrix is -sigma_y.
    return sp.Matrix([symbols[0], -symbols[1], symbols[2]])


def chen_kernels():
    """Chen (2022), (13)-(15), (22)-(23), common prefactors removed.

    Both processes: s = electron acceleration direction = -eps;
    b = v x s = beta; zeta_minus/plus are PHYSICAL pair spins.
    xi_Chen = (-xi3_obs, xi2_obs, xi1_obs). Chen (20) explicitly gives
    the positron's two-component spinor label as -zeta_plus.
    F2 = i(W12-W21): its outer i reverses the bracket in (13).
    This transcription uses vector products, not T21's dyad matrices.
    """
    a, b, v = (sp.Matrix([-1, 0, 0]), sp.Matrix([0, -1, 0]),
               sp.Matrix([0, 0, 1]))
    si, sf = map(physical_vectors, (lcfa.SI, lcfa.SF))
    se, positron = map(physical_vectors, (lcfa.SE, lcfa.SP))
    r, t = lcfa.R1, lcfa.R2
    l, lp = lcfa.LAM, lcfa.LAMP
    ep, em = 1 - l, 1 - lp
    kap, kt, h = (1 + ep**2) / ep, (1 - ep**2) / ep, l**2 / ep
    f0 = (kap - t + (2 - t) * dot(si, sf)
          - r * (l * dot(si, b) + l / ep * dot(sf, b))
          + h * (1 - t) * dot(si, v) * dot(sf, v))
    f1 = (kt / 2 * dot(v, sf.cross(si))
          + r * (l / ep * dot(si, a) + l * dot(sf, a))
          - h / 2 * t * (dot(si, a) * dot(sf, b)
                         + dot(si, b) * dot(sf, a)))
    f2 = -(kt / 2 * r * dot(a, sf.cross(si))
           + (-kt + l * t) * dot(si, v)
           + (-kt + l / ep * t) * dot(sf, v)
           + h / 2 * r * (dot(si, v) * dot(sf, b)
                          + dot(si, b) * dot(sf, v)))
    f3 = (1 + kap / 2 * dot(si, sf)
          - r * (l / ep * dot(si, b) + l * dot(sf, b))
          + h / 2 * (-dot(si, v) * dot(sf, v)
                     + t * (dot(si, b) * dot(sf, b)
                            - dot(si, a) * dot(sf, a))))

    rho = (lp**2 + em**2) / (lp * em)
    delta = (lp**2 - em**2) / (lp * em)
    u = 1 / (lp * em)
    g0 = (t + rho + (t - 2) * dot(se, positron)
          + r * (dot(positron, b) / lp - dot(se, b) / em)
          + (rho * t - (lp - em)**2 / (lp * em))
          * dot(se, v) * dot(positron, v))
    g1 = (-delta / 2 * dot(v, positron.cross(se))
          + r * (dot(positron, a) / em - dot(se, a) / lp)
          - u / 2 * t * (dot(se, b) * dot(positron, a)
                         + dot(se, a) * dot(positron, b)))
    g2 = (-u / 2 * r * dot(a, se.cross(positron))
          + (t / lp + delta) * dot(positron, v)
          + (t / em - delta) * dot(se, v)
          - delta / 2 * r * (dot(se, v) * dot(positron, b)
                             + dot(se, b) * dot(positron, v)))
    g3 = (-1 + rho / 2 * dot(se, positron)
          - r * (dot(positron, b) / em - dot(se, b) / lp)
          - (lp - em)**2 / (2 * lp * em)
          * dot(se, v) * dot(positron, v)
          + u / 2 * t * (dot(se, b) * dot(positron, b)
                         - dot(se, a) * dot(positron, a)))
    x1, x2, x3 = lcfa.XI
    return (sp.expand(f0 - x3 * f1 + x2 * f2 + x1 * f3),
            sp.expand(g0 - x3 * g1 + x2 * g2 + x1 * g3))


def zero_difference(a, b):
    return sp.cancel(sp.expand(a - b)) == 0


def check_papers(report):
    c, bw = chen_kernels()
    for label, a, b in [('NLCS', c, lcfa.nlcs_probability(-1)),
                         ('NBW', bw, lcfa.nbw_probability())]:
        # Equality of the full multivariate polynomial includes every
        # nonzero spin/Stokes coefficient, not just selected states.
        report('Chen full %s kernel: all coefficients exactly equal' % label,
               zero_difference(a, b))
    charge_flip = {lcfa.SI[0]: -lcfa.SI[0],
                   lcfa.SI[1]: -lcfa.SI[1],
                   lcfa.SF[0]: -lcfa.SF[0],
                   lcfa.SF[1]: -lcfa.SF[1]}
    report('NLCS positron: charge-transformed Chen kernel exactly equal',
           zero_difference(c.subs(charge_flip, simultaneous=True),
                           lcfa.nlcs_probability(+1)))

    # SK20 (38), (59): remove alpha/(4b), convert Airy functions.
    # h_dot positive: sigma_p/q refer to physical beta spin components.
    # h_dot negative: local beta reverses, so both components reverse.
    i, f, e, p, tau = sp.symbols('sk_i sk_f sk_e sk_p sk_tau')
    l, lp, r, t = lcfa.LAM, lcfa.LAMP, lcfa.R1, lcfa.R2
    g = 1 + l**2 / (2 * (1 - l))
    gt = 1 - 1 / (2 * lp * (1 - lp))
    for sign in (-1, 1):
        sk_c = (-t * (1 + i*f + tau*i*f*(1-g))
                - sign*r * (l*i + l/(1-l)*f
                            + tau*(l/(1-l)*i + l*f))
                + 2 * (g + i*f + tau*(1+g*i*f)/2))
        sk_b = (t * (1 + p*e + tau*p*e*(1-gt))
                + sign*r * (p/lp - e/(1-lp)
                            + tau*(e/lp - p/(1-lp)))
                - 2 * (gt + p*e + tau*(1+gt*p*e)/2))
        subs_c = dict(zip(lcfa.SI + lcfa.SF + lcfa.XI,
                          (0, sign*i, 0, 0, sign*f, 0, tau, 0, 0)))
        subs_b = dict(zip(lcfa.SE + lcfa.SP + lcfa.XI,
                          (0, sign*e, 0, 0, sign*p, 0, tau, 0, 0)))
        report('SK20 (38): all NLCS magnetic channels, field sign %d' % sign,
               zero_difference(c.subs(subs_c), sk_c))
        report('SK20 (59): all NBW magnetic channels, field sign %d' % sign,
               zero_difference(bw.subs(subs_b), sk_b))

    # These enforce the actual printed parentheses in Qian, independently
    # of the informational comparisons in the main reference.
    fr = lcfa.t21_frame(0)
    qb = lcfa.qian_nbw(fr)
    tb = lcfa.nbw_quantities(lcfa.nbw_probability(axes_sign=1,
                                                    block_order=False))
    for q, actual in [('G_circ', 'G2'), ('je_v', 'je_v'), ('jp_v', 'jp_v')]:
        report('Qian printed %s: no circular marginal slip' % q,
               zero_difference(qb[q], tb[actual]))
    rho = (1-lp)/lp + lp/(1-lp)
    expected = (2-rho) * lcfa.SE[2] * lcfa.SP[2]
    report('Qian A.50: exact remaining algebra difference',
           zero_difference(tb['G1'] - qb['G_lin'], expected))


def independent_ratios(x):
    """Scaled cosh integral, independent of the generator's primary method."""
    hi = 2 * np.arcsinh(np.sqrt(45 / (2*x)))
    integral = quad(lambda u: np.exp(-2*x*np.sinh(u/2)**2)
                    * np.cosh(u/3) / np.cosh(u), 0, hi,
                    epsabs=0, epsrel=3e-12)[0]
    k23 = kve(2/3, x)
    return kve(1/3, x) / k23, integral / k23


IDENTITY = sp.eye(2)
SX = sp.Matrix([[0, 1], [1, 0]])
SY = sp.Matrix([[0, -sp.I], [sp.I, 0]])
SZ = sp.diag(1, -1)
SPIN_PAULI = [IDENTITY, SX, -SY, SZ]
PHOTON_PAULI = [IDENTITY, SZ, SY, -SX]


def choi_matrix(expr, groups, bases):
    """8x8 Choi operator; first group is input and its matrices transpose.

    P is a bilinear analyser expectation on the output and a density
    expectation on the input. Transposing the INPUT Pauli matrices is
    essential; the untransposed coefficient operator need not be PSD.
    Overall positive factors and the conventional factor 1/4 are omitted.
    """
    gens = sum((list(g) for g in groups), [])
    out = sp.zeros(8)
    bases = [[m.T for m in bases[0]]] + list(bases[1:])
    for powers, coeff in sp.Poly(sp.expand(expr), *gens).terms():
        indices = []
        for start in range(0, 9, 3):
            block = powers[start:start+3]
            if sum(block) > 1:
                raise ValueError('kernel is not affine in each state')
            indices.append(next((j+1 for j, power in enumerate(block)
                                 if power), 0))
        out += coeff * sp.kronecker_product(
            *(basis[index] for basis, index in zip(bases, indices)))
    return out.applyfunc(sp.factor)


def check_quantum_maps(report):
    r, t, l, lp = lcfa.R1, lcfa.R2, lcfa.LAM, lcfa.LAMP
    discriminant = t*(3-t) - 2*r**2
    # Symbolic factorisation covers ALL fractions and state directions.
    # The numerical part then checks the physical scalar Bessel domain.
    for label, expr, lam, groups, bases, pivots, determinant in [
            ('NLCS e-', lcfa.nlcs_probability(-1), l,
             [lcfa.SI, lcfa.SF, lcfa.XI],
             [SPIN_PAULI, SPIN_PAULI, PHOTON_PAULI], [0, 1, 2],
             4*l**2/(1-l)*(1-t)*discriminant),
            ('NLCS e+', lcfa.nlcs_probability(1), l,
             [lcfa.SI, lcfa.SF, lcfa.XI],
             [SPIN_PAULI, SPIN_PAULI, PHOTON_PAULI], [0, 1, 2],
             4*l**2/(1-l)*(1-t)*discriminant),
            ('NBW', lcfa.nbw_probability(), lp,
             [lcfa.XI, lcfa.SE, lcfa.SP],
             [PHOTON_PAULI, SPIN_PAULI, SPIN_PAULI], [0, 1, 5],
             4/(lp*(1-lp))*(1-t)*discriminant)]:
        j = choi_matrix(expr, groups, bases)
        h = j.extract(pivots, pivots)
        seed = {r: 0, t: sp.Rational(1, 2)}
        transform = (j[:, pivots].subs(seed) * h.subs(seed).inv())
        transform = transform.applyfunc(sp.factor)
        residual = j - transform*h*transform.conjugate().T
        report('%s Choi: exact rank-three factorisation' % label,
               all(sp.cancel(v) == 0 for v in residual)
               and zero_difference(h.det(), determinant))
        f = sp.lambdify((lam, r, t), j, 'numpy')
        worst = np.inf
        for x in np.geomspace(1e-12, 1e6, 49):
            r1, r2 = independent_ratios(x)
            for fraction in (1e-6, .005, .1, .37, .5, .9, .995, 1-1e-6):
                m = f(fraction, r1, r2)
                if not np.allclose(m, m.conjugate().T):
                    raise ValueError('non-Hermitian Choi matrix')
                ev = np.linalg.eigvalsh(m)
                worst = min(worst, ev.min() / ev.max())
        report('%s Choi PSD on independent Bessel grid' % label,
               worst >= -2e-12, 'min eigenvalue / max = %.2e' % worst)


def check_table(report):
    path = Path(lcfa.__file__).with_name('bessel_k.table')
    with path.open() as fh:
        count, lo, hi = map(float, fh.readline().split())
        table = np.loadtxt(fh)
    shape_ok = (table.shape == (int(count), 6)
                and np.all(np.isfinite(table))
                and np.all(np.diff(table[:, 0]) > 0)
                and abs(table[0, 0]-lo) < 1e-12
                and abs(table[-1, 0]-hi) < 1e-12)
    report('on-disk Bessel table: header, shape and finite entries', shape_ok)
    if not shape_ok:
        return
    nodes_r, nodes_t = 10**table[:, 1], 10**table[:, 2]
    nodes_d = nodes_t*(3-nodes_t)-2*nodes_r**2
    # a=log R1, b=log R2: the allowed region is a <= f(b),
    # f(b)=(b+log(3-exp b)-log 2)/2, with
    # f''(b)=-3 exp(b)/(2(3-exp(b))**2) < 0. Its hypograph is convex.
    # Thus every log-interpolation segment is allowed if its endpoints
    # are allowed. This certifies the entire cells, not just a dense grid.
    report('all table nodes satisfy the convex log-interpolation domain',
           np.all((nodes_r > 0) & (nodes_r < 1)
                  & (nodes_t > 0) & (nodes_t < 1) & (nodes_d > 0)),
           'minimum node discriminant = %.2e' % np.min(nodes_d))
    import gen_bessel_k_table as gk
    e, d = gk.E1, gk.D1
    y_max, z_max = (10**lo)**(1/3), 1/10**hi
    # Exact discriminant of the small-x polynomial is y^2 times the
    # bracket below. A conservative coefficient bound covers 0<y<y_max.
    small_bound = (3*d-2*e*e - abs(-d*d+1.5*e)*y_max**2
                   - abs(3*d*e-4.5)*y_max**4 - 2.25*e*e*y_max**6)
    # Exact large-x discriminant is
    # z^2 (288 + 3624z - 2641z^2)/2592, z=1/x.
    large_bound = 288-2641*z_max*z_max
    report('analytic tail bounds preserve positivity for every outside x',
           small_bound > 0 and large_bound > 0 and y_max <= .01+1e-15
           and z_max < 2/3,
           'small bracket >= %.6f; large bracket >= %.6f'
           % (small_bound, large_bound))
    # Use the ACTUAL documented interpolation: linear in log10 R,
    # at 24 interior points per cell, plus every node and both tails.
    logx = np.linspace(lo, hi, 25*(len(table)-1)+1)
    r = 10**np.interp(logx, table[:, 0], table[:, 1])
    t = 10**np.interp(logx, table[:, 0], table[:, 2])
    xs = np.concatenate((np.geomspace(1e-12, 10**lo, 101),
                         np.geomspace(10**hi, 1e8, 101)))
    tails = np.array([lcfa.bessel_ratios(x) for x in xs])
    r, t = np.concatenate((r, tails[:, 0])), np.concatenate((t, tails[:, 1]))
    d = t*(3-t) - 2*r*r
    report('table interpolation and tails: quantum positivity domain',
           np.all((r > 0) & (r <= 1) & (t > 0) & (t <= 1))
           and np.min(d) >= -2e-12, 'min discriminant = %.2e' % np.min(d))
    error = 0.
    for logv in np.linspace(lo+.005, hi-.005, 151):
        interpolated = 10**np.array([
            np.interp(logv, table[:, 0], table[:, col]) for col in (1, 2)])
        exact = independent_ratios(10**logv)
        error = max(error, np.max(np.abs(interpolated / exact - 1)))
    report('on-disk interpolated ratios vs independent cosh integral',
           error < 1e-5, 'max relative error = %.2e' % error)


def check_witnesses(report):
    # Fixed axes, no random search: x=1 is a valid physical argument
    # (NLCS chi=6 at lambda=.9; NBW chi_gamma=25/6 at lp=.2).
    r1, r2 = independent_ratios(1.)
    for label, lam, lv, ga, gb, a, b, fun in [
            ('NLCS', lcfa.LAM, .9, lcfa.SI, lcfa.SF,
             (-1, 0, 0), (0, 1, 0), lambda order:
             lcfa.nlcs_probability(-1, block_order=order)),
            ('NBW', lcfa.LAMP, .2, lcfa.SE, lcfa.SP,
             (0, 0, -1), (0, -1, 0), lambda order:
             lcfa.nbw_probability(block_order=order))]:
        values = dict(zip(ga+gb, a+b))
        values.update({lam: lv, lcfa.R1: r1, lcfa.R2: r2})
        minima = []
        for order in (False, True):
            p = fun(order)
            f0 = float(p.subs(lcfa.zero(lcfa.XI)).subs(values))
            f = np.array([float(sp.diff(p, x).subs(values))
                          for x in lcfa.XI])
            scale = float(p.subs(lcfa.zero(lcfa.XI+ga+gb)).subs(values))
            minima.append((f0 - np.linalg.norm(f)) / scale)
        report('%s fixed witness: literal negative, block nonnegative' % label,
               minima[0] < -.8 and minima[1] >= -2e-12,
               'literal %.6f; block %.6f (P/P_unpol)' % tuple(minima))
        literal, correct = fun(False), fun(True)
        change = correct-literal
        # Two independent antisymmetric sectors: diagonal (i sigma2)
        # and circular (k B - B k). Correcting either alone is insufficient.
        axes = np.concatenate((np.eye(3), -np.eye(3)))
        worst = []
        for component in (lcfa.XI[2], lcfa.XI[1]):
            partial = literal + component*sp.diff(change, component)
            f0 = partial.subs(lcfa.zero(lcfa.XI))
            f = [sp.diff(partial, x) for x in lcfa.XI]
            evaluate = sp.lambdify(list(ga+gb),
                                   [v.subs({lam: lv, lcfa.R1: r1,
                                            lcfa.R2: r2}) for v in [f0]+f],
                                   'numpy')
            minimum = np.inf
            for aa, bb in product(axes, repeat=2):
                fv = np.array(evaluate(*aa, *bb))
                minimum = min(minimum, fv[0]-np.linalg.norm(fv[1:]))
            worst.append(minimum)
        report('%s: correcting either antisymmetric sector alone fails' % label,
               max(worst) < -1e-3)


def check_sampling(report):
    rng = np.random.default_rng(21)
    r1, r2 = independent_ratios(1.)
    for label, expr, lam, initial, ga, gb in [
            ('NLCS e-', lcfa.nlcs_probability(-1), lcfa.LAM,
             lcfa.SI, lcfa.SF, lcfa.XI),
            ('NLCS e+', lcfa.nlcs_probability(1), lcfa.LAM,
             lcfa.SI, lcfa.SF, lcfa.XI),
            ('NBW', lcfa.nbw_probability(), lcfa.LAMP,
             lcfa.XI, lcfa.SE, lcfa.SP)]:
        worst, error = np.inf, 0.
        for _ in range(16):
            state = rng.normal(size=3)
            state *= rng.uniform(0, 1) / np.linalg.norm(state)
            a, b = rng.normal(size=(2, 3))
            a, b = a/np.linalg.norm(a), b/np.linalg.norm(b)
            vals = dict(zip(initial, state))
            vals.update({lam: rng.uniform(.01, .99), lcfa.R1: r1,
                         lcfa.R2: r2})
            p = expr.subs(vals)
            w = float(p.subs(lcfa.zero(ga+gb)))
            bins = {}
            for sign_a, sign_b in product((-1, 1), repeat=2):
                subs = dict(zip(ga+gb, tuple(sign_a*a)+tuple(sign_b*b)))
                bins[sign_a, sign_b] = float(p.subs(subs)) / (4*w)
            worst = min(worst, min(bins.values()))
            error = max(error, abs(sum(bins.values()) - 1))
            for syms, direction, index in [(ga, a, 0), (gb, b, 1)]:
                marginal = sum(float(sp.diff(p, s).subs(lcfa.zero(ga+gb)))
                               * d for s, d in zip(syms, direction)) / w
                sampled = sum(signs[index]*prob for signs, prob in bins.items())
                error = max(error, abs(marginal-sampled))
            # Sequential collapse: first output analyser, then the
            # conditional second analyser. Their product MUST equal P/4w.
            for (sa, sb), prob in bins.items():
                q = p.subs(dict(zip(ga, sa*a)))
                normalisation = float(q.subs(lcfa.zero(gb)))
                sequential = normalisation/(2*w) * float(
                    q.subs(dict(zip(gb, sb*b))))/(2*normalisation)
                error = max(error, abs(sequential-prob))
        report('%s four-bin/conditional sampler: positivity and marginals' % label,
               worst >= -2e-12 and error < 2e-12,
               'min bin %.3e; max error %.2e' % (worst, error))

    # A mean output state is not a measurement outcome. Substituting g/w
    # as the final analyser before evaluating F/F0 changes the photon
    # marginal even with the correct kernel. Keep the two marginals in
    # mean mode, or use a genuine sampled analyser for conditioning.
    q, _ = lcfa.derive()
    vals = {lcfa.QS: -1, lcfa.LAM: .5, lcfa.R1: r1, lcfa.R2: r2,
            **lcfa.zero(lcfa.SI)}
    w = float(q['w'].subs(vals))
    mean_spin = [float(q['g_'+c].subs(vals))/w for c in 'ebv']
    selected = {**vals, **dict(zip(lcfa.SF, mean_spin))}
    wrong = float(q['F1'].subs(selected)/q['F0'].subs(selected))
    marginal = float(q['F1_bar'].subs(vals))/w
    report('mean-mode counterexample: conditional analyser is not the marginal',
           abs(wrong-marginal) > .005,
           'conditional %.8f; correct %.8f' % (wrong, marginal))


def check_exports(report):
    c, b = lcfa.derive()
    for label, quantities, gens in [('NLCS', c, lcfa.NLCS_GENS),
                                    ('NBW', b, lcfa.NBW_GENS)]:
        ok = True
        for expr in quantities.values():
            displayed = lcfa.tidy(expr, gens, display=True)
            restored = displayed.subs({lcfa.S1: 1-lcfa.LAM,
                                       lcfa.S2: 1-lcfa.LAMP,
                                       lcfa.S3: lcfa.LAMP})
            ok &= zero_difference(expr, restored)
        report('%s: every exported coefficient equals the kernel' % label, ok)


def check_ptarmigan(report):
    # Independent native-code fixtures, pinned to local Ptarmigan 1.6.1,
    # 9e7caf793645b3fc6bef34ddd8808e367c7f7439. Ptarmigan has no lepton
    # spin here: these check absolute spectra, NOT the three-way terms.
    def k_values(x):
        if x > 700:
            return 0., 0.
        _, r2 = independent_ratios(x)
        k23 = kve(2/3, x)*np.exp(-x)
        return k23, r2*k23

    from scipy.constants import alpha, hbar, m_e, c
    tau = hbar/(m_e*c*c)
    error = 0.
    for chi, target in [(0.026, 2.07935e14), (3.5, 1.58485e16),
                        (9.98, 3.45844e16), (12.4, 4.04647e16)]:
        def integrand(u):
            k, integral = k_values(2*u/(3*chi))
            return ((1+u+1/(1+u))*k-integral)/(1+u)**2
        value = alpha/(np.sqrt(3)*np.pi*tau*1000) * quad(
            integrand, 0, np.inf, epsabs=0, epsrel=2e-8, limit=200)[0]
        error = max(error, abs(value/target-1))
    report('absolute NLCS rates vs four Ptarmigan native fixtures',
           error < 1e-3, 'max relative error = %.2e' % error)
    error = 0.
    for chi, target in [(0.042, 6.07753899492993e-29),
                        (.105, 2.1082097875655205e-12),
                        (.42, .0003779613236658133),
                        (1.05, .015977478443872017),
                        (4.2, .08917816786414409),
                        (12., .10884579479913804),
                        (42., .09266735324318656)]:
        def integrand(s):
            k, integral = k_values(2/(3*chi*s*(1-s)))
            return ((s/(1-s)+(1-s)/s)*k+integral)
        value = 2/(np.sqrt(3)*np.pi*chi) * quad(
            integrand, 0, .5, epsabs=0, epsrel=2e-8)[0]
        error = max(error, abs(value/target-1))
    report('absolute NBW rates vs seven Ptarmigan native fixtures',
           error < 1e-7, 'max relative error = %.2e' % error)


def run_audit(report=lcfa.report):
    print('\nIndependent theory, quantum maps, sampling and exports:')
    check_papers(report)
    check_quantum_maps(report)
    check_table(report)
    check_witnesses(report)
    check_sampling(report)
    check_exports(report)
    check_ptarmigan(report)


if __name__ == '__main__':
    lcfa.RESULTS.clear()
    run_audit()
    failures = [name for name, passed in lcfa.RESULTS if not passed]
    print('\n%d independent checks, %d failed' % (
        len(lcfa.RESULTS), len(failures)))
    raise SystemExit(bool(failures))
