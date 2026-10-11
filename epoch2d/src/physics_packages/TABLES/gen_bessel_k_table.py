#!/usr/bin/env python3
"""Generate bessel_k.table: K_1/3, K_2/3, IntK_1/3 and their ratios.

Every spin- and polarisation-resolved LCFA final-state quantity in Qian's
thesis (Appendix A: w and g, F_{0-3}, w_0, j_e, j_p, G_{0-3}) needs the
following building blocks

    K_1/3(x),   K_2/3(x),   IntK_1/3(x) = int_x^inf K_1/3(t) dt,

all evaluated at one argument:

    NLCS:  u'  = 2 u / (3 chi),         u = lambda / (1 - lambda),
    NBW:   u~' = 2 u~ / (3 chi_gamma),  u~ = 1/lambda~ + 1/(1 - lambda~),

with chi_gamma the conventional photon parameter (EPOCH's photon chi is
chi_gamma / 2). In Torgrimsson's (T21) Airy notation x = (2/3) z^(3/2).

The code only needs ratios of such sums (<s_f> = g/w,
<xi> = F/F_0, <s_e> = j_e/w_0, the collapse probabilities, the pair-spin
bin probabilities G/sum(G)). Here sum(G) = 4 w_0 for four spin outcomes;
the photon-independent joint G_0 is not that normalisation. Dividing
through by K_2/3 leaves two functions,

    R1(x) = K_1/3(x) / K_2/3(x)    = -sqrt(z) Ai(z)  / Ai'(z),
    R2(x) = IntK_1/3(x) / K_2/3(x) = -z       Ai1(z) / Ai'(z),

both smooth and in (0, 1), so there is no exp(-x) underflow at large x
and no divergence at small x. For example w/K_2/3 = kappa - R2
+ lambda s_beta R1 (thesis A.17). The absolute functions are tabulated
for offline use and diagnostics.

Numerics
--------
K_nu uses scipy's exponentially scaled kve, so log10 K is accurate
at any x. IntK_1/3 is computed as

    x <  1 :  pi/sqrt(3) - int_0^x K_1/3(t) dt   (t = w^3 removes the
                                                  t^(-1/3) singularity)
    x >= 1 :  exp(-x) int_0^inf kve_1/3(x + s) exp(-s) ds

and checked (--check) against the independent representation

    IntK_1/3(x) = int_0^inf exp(-x cosh t) cosh(t/3) / cosh t dt,

K_nu against scipy's (scaled) Airy functions, and IntK_1/3 against
the derivative identity d IntK_1/3 / dx = -K_1/3.

Range, interpolation and tails
------------------------------
The table spans log10 x in [-6, 2.5] at 100 points per decade. Columns
are interpolated linearly in log10 x; storing log10 R (rather than R)
makes the small-x power laws exact and keeps the relative error of R1
and R2 below 1e-5 (the worst case is near x = 0.65). The absolute
functions are stored as log10(exp(x) K), which interpolates as well;
multiply by exp(-x) to recover K.

Outside the table use the two-term expansions

    x < x_min :  R1 = e1 x^(1/3) - (3/2) x,
                 R2 = d1 x^(2/3) - (3/2) e1 x^(4/3),
    x > x_max :  R1 = 1 - 1/(6x) + 7/(72 x^2),
                 R2 = 1 - 2/(3x) + 1/x^2,

    e1 = Gamma(1/3) / (2^(1/3) Gamma(2/3)),
    d1 = 2 pi / (sqrt(3) 2^(2/3) Gamma(2/3)).

The small-x errors fall as x^(4/3) (1e-8 relative at x = 1e-6) and the
large-x errors as x^(-3) (1e-7 at x = 10^2.5), so the tails join the
table without rescaling. This script prints e1 and d1 for the Fortran
constants.

Format (header as hsokolov.table; EPOCH reads it with list-directed I/O)
------
    <n> <log10 x_min> <log10 x_max>
    <log10 x> <log10 R1> <log10 R2>
        <log10 e^x K_1/3> <log10 e^x K_2/3> <log10 e^x IntK_1/3>
    ...

(one row per point; the second line above continues the first).

Usage
-----
    python gen_bessel_k_table.py            # writes the table
    python gen_bessel_k_table.py --check    # verifies it

Requires numpy and scipy.
"""

import argparse
import os
import warnings

import numpy as np
from scipy.integrate import quad, IntegrationWarning
from scipy.special import airye, gamma, kv, kve

LN10 = np.log(10.0)
SQRT3 = np.sqrt(3.0)
INTK0 = np.pi / SQRT3
E1 = gamma(1.0 / 3.0) / (2.0**(1.0 / 3.0) * gamma(2.0 / 3.0))
D1 = 2.0 * np.pi / (SQRT3 * 2.0**(2.0 / 3.0) * gamma(2.0 / 3.0))
QUAD = dict(epsabs=0.0, epsrel=1e-13, limit=400)


def intk_scaled(x):
    """exp(x) IntK_1/3(x), the primary method."""
    if x < 1.0:
        part = quad(lambda w: 3.0 * w * w * kv(1.0 / 3.0, w**3),
                    0.0, x**(1.0 / 3.0), **QUAD)[0]
        return np.exp(x) * (INTK0 - part)
    return quad(lambda s: kve(1.0 / 3.0, x + s) * np.exp(-s),
                0.0, np.inf, **QUAD)[0]


def intk_scaled_cosh(x):
    """exp(x) IntK_1/3(x) from the cosh representation (independent)."""
    # The scaled integrand is below exp(-40) beyond t_max, and below
    # exp(-2t/3) < 1e-17 beyond t = 60 whatever x is.
    t_max = min(60.0, np.arccosh(1.0 + 40.0 / x))
    return quad(lambda t: np.exp(-x * (np.cosh(t) - 1.0))
                * np.cosh(t / 3.0) / np.cosh(t), 0.0, t_max, **QUAD)[0]


def columns(x):
    """log10 of R1, R2 and exp(x) times K_1/3, K_2/3, IntK_1/3 at x."""
    k13 = kve(1.0 / 3.0, x)
    k23 = kve(2.0 / 3.0, x)
    ik = intk_scaled(x)
    return np.log10([k13 / k23, ik / k23, k13, k23, ik])


def tail_small(x):
    return (E1 * x**(1.0 / 3.0) - 1.5 * x,
            D1 * x**(2.0 / 3.0) - 1.5 * E1 * x**(4.0 / 3.0))


def tail_large(x):
    return (1.0 - 1.0 / (6.0 * x) + 7.0 / (72.0 * x * x),
            1.0 - 2.0 / (3.0 * x) + 1.0 / (x * x))


def ratios_exact(x):
    """R1, R2 from the independent methods (Airy for K, cosh for IntK)."""
    z = (1.5 * x)**(2.0 / 3.0)
    eai, eaip = airye(z)[:2]
    return -np.sqrt(z) * eai / eaip, intk_scaled_cosh(x) / kve(2.0 / 3.0, x)


def check(log_x, table):
    failures = []

    def require(error, tolerance, label):
        if not np.isfinite(error) or error > tolerance:
            failures.append('%s: %.3g > %.3g' % (label, error, tolerance))

    names = ['R1', 'R2', 'K_1/3', 'K_2/3', 'IntK_1/3']
    x = 10.0**log_x

    # 1. Independent methods at every grid point.
    z = (1.5 * x)**(2.0 / 3.0)
    eai, eaip = airye(z)[:2]
    k13_airy = np.pi * np.sqrt(3.0 / z) * eai
    k23_airy = -np.pi * SQRT3 * eaip / z
    ik_cosh = np.array([intk_scaled_cosh(v) for v in x])
    k13 = kve(1.0 / 3.0, x)
    k23 = kve(2.0 / 3.0, x)
    ik = np.array([intk_scaled(v) for v in x])
    print('independent methods, max relative difference on the grid:')
    print('  K_1/3 vs Airy     %.1e' % np.max(np.abs(k13 / k13_airy - 1)))
    print('  K_2/3 vs Airy     %.1e' % np.max(np.abs(k23 / k23_airy - 1)))
    print('  IntK_1/3 vs cosh  %.1e' % np.max(np.abs(ik / ik_cosh - 1)))
    for label, values, reference in (
            ('K_1/3 vs Airy', k13, k13_airy),
            ('K_2/3 vs Airy', k23, k23_airy),
            ('IntK_1/3 vs cosh', ik, ik_cosh)):
        require(np.max(np.abs(values / reference - 1)), 1e-10, label)
    if not np.all(np.isfinite(table)):
        failures.append('table contains non-finite entries')
    if not np.all((table[:, :2] < 0.0)):
        failures.append('ratios must lie strictly between zero and one')

    # 2. d IntK_1/3 / dx = -K_1/3 (central difference).
    worst = 0.0
    for v in (1e-3, 0.1, 1.0, 10.0, 100.0):
        h = 1e-5 * min(v, 1.0)
        deriv = (np.exp(-h) * intk_scaled(v + h)
                 - np.exp(h) * intk_scaled(v - h)) / (2.0 * h)
        worst = max(worst, abs(-deriv / kve(1.0 / 3.0, v) - 1.0))
    print('  dIntK/dx = -K_1/3 %.1e (finite difference)' % worst)
    require(worst, 1e-6, 'derivative identity')
    # IntK_1/3 -> pi/sqrt(3) - (3/2) Gamma(1/3) 2^(-2/3) x^(2/3).
    v = 1e-9
    lead = INTK0 - 1.5 * gamma(1.0 / 3.0) * 2.0**(-2.0 / 3.0) * v**(2.0 / 3.0)
    print('  IntK_1/3 small-x limit %.1e'
          % abs(intk_scaled(v) * np.exp(-v) / lead - 1))

    # 3. Interpolation error at the midpoints.
    mid = 0.5 * (log_x[:-1] + log_x[1:])
    exact = np.array([columns(10.0**v) for v in mid])
    print('linear interpolation in log10 x, max error at midpoints:')
    for k, name in enumerate(names):
        interp = np.interp(mid, log_x, table[:, k])
        err = np.max(np.abs(10.0**(interp - exact[:, k]) - 1.0))
        print('  %-15s relative %.1e' % (name, err))
        require(err, 1e-5, name + ' interpolation')

    # 4. Tails at and beyond the table edges.
    print('tails, relative error against the independent methods:')
    for v in (x[0], 1e-8, 1e-10, x[-1], 1e3, 1e4):
        approx = tail_small(v) if v <= x[0] else tail_large(v)
        ref = ratios_exact(v)
        print('  x = %-9.3g R1 %.1e   R2 %.1e'
              % (v, abs(approx[0] / ref[0] - 1),
                 abs(approx[1] / ref[1] - 1)))
        require(max(abs(approx[0] / ref[0] - 1),
                    abs(approx[1] / ref[1] - 1)), 1e-6,
                'tail at x=%g' % v)
    if failures:
        raise RuntimeError('validation failed:\n' + '\n'.join(failures))
    print('all checks passed')


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--log-min', type=float, default=-6.0)
    p.add_argument('--log-max', type=float, default=2.5)
    p.add_argument('--per-decade', type=int, default=100)
    p.add_argument('--output', default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'bessel_k.table'))
    p.add_argument('--check', action='store_true',
                   help='verify methods, interpolation error and tails')
    args = p.parse_args()
    if (not np.isfinite(args.log_min) or not np.isfinite(args.log_max)
            or args.log_min >= args.log_max or args.per_decade < 1):
        p.error('require finite log-min < log-max and per-decade >= 1')
    # Reject powers outside the floating-point range before quadrature.
    if args.log_min < -300.0 or args.log_max > 300.0:
        p.error('logarithmic bounds must lie within [-300, 300]')

    # quad occasionally reports round-off at these tight tolerances; the
    # independent checks show it does not affect the result.
    warnings.simplefilter('ignore', IntegrationWarning)

    n = int(round((args.log_max - args.log_min) * args.per_decade)) + 1
    if n < 2:
        p.error('the requested grid must contain at least two points')
    log_x = np.linspace(args.log_min, args.log_max, n)
    table = np.array([columns(10.0**v) for v in log_x])
    if not np.all(np.isfinite(table)):
        raise RuntimeError('non-finite table entries; output was not written')
    if args.check:
        check(log_x, table)

    with open(args.output, 'w') as fh:
        fh.write('%d\t%.16g\t%.16g\n' % (n, args.log_min, args.log_max))
        for v, row in zip(log_x, table):
            fh.write('%.16f\t' % v
                     + '\t'.join('%.16f' % c for c in row) + '\n')
    print('wrote %s (%d points, log10 x %g..%g)'
          % (args.output, n, args.log_min, args.log_max))
    print('tail constants: e1 = %.16e  d1 = %.16e' % (E1, D1))
    print('R1(1) = %.10f   R2(1) = %.10f'
          % tuple(10.0**columns(1.0)[:2]))

if __name__ == '__main__':
    main()
