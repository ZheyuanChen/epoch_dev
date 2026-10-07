#!/usr/bin/env python3
"""Generate anomalous_moment.table: the chi-dependent anomalous moment.

The one-loop anomalous magnetic moment of an electron in a locally
constant crossed field depends on the quantum parameter chi (Ritus 1970;
Seipt & Thomas, arXiv:2307.02114, eq. 37):

    a(chi) = (alpha/chi) int_0^1 dl l Gi(z) / sqrt(z),
    z = (l / (chi (1 - l)))^(2/3),

with Gi the Scorer function. The table holds the ratio

    f(chi) = a(chi) / (alpha / 2 pi),

which tends to 1 as chi -> 0, so EPOCH can use a_eff = a_species * f(chi)
with the species' own (deck or identify) moment.

Numerics
--------
Substituting u = l / (1 - l) gives a form that is well behaved for all chi:

    f(chi) = 2 pi chi^(-2/3)
             * int_0^inf u^(2/3) (1 + u)^(-3) Gi((u/chi)^(2/3)) du.

Gi(z), z >= 0, uses the exponentially damped representation

    Gi(z) = -(1/pi) int_0^inf exp(-t^3/3 - z t/2)
                              * cos(sqrt(3) z t/2 + 2 pi/3) dt,

which has no cancellation at large z (Gi = Bi - Hi loses ~1e-6 at z = 10);
the two agree to 12 digits at moderate z.

Range and tail
--------------
The table spans log10(chi) in [-3, 4] at 20 points per decade (0.05 dex,
the spacing of hsokolov.table). Linear interpolation of log10 f in
log10 chi is then accurate to 3e-4. Outside the table EPOCH uses:

    chi < 1e-3 :  f = 1                    (1 - f = 6e-5 at chi = 1e-3)
    chi > 1e4  :  f = s y (c0 + c1 y),   y = chi^(-2/3)

The tail is the two-term large-chi expansion of the integral,

    c0 = 2 pi Gi(0)  B(5/3, 4/3),   Gi(0)  = 1 / (3^(7/6) Gamma(2/3)),
    c1 = 2 pi Gi'(0) B(7/3, 2/3),   Gi'(0) = 1 / (3^(5/6) Gamma(1/3)),

(B the Beta function), with s chosen at run time so that the tail meets the
last table entry. The unscaled expansion is off by 8e-5 at chi = 1e4 and its
error falls as chi^(-4/3) ln chi; after scaling, the tail error tends to
|s - 1| = 8e-5 and stays below it for all chi > 1e4 (checked to 1e6).
This script prints c0 and c1 for the Fortran constants.

Note that the one-loop result itself is not expected to hold beyond
alpha chi^(2/3) ~ 1 (chi ~ 1600, the Ritus-Narozhny conjecture), so the 
upper truncation is the least thing to worry about.

Format (as hsokolov.table; EPOCH reads the header with list-directed I/O)
------
    <n> <log10 chi_min> <log10 chi_max>
    <log10 chi_1> <log10 f_1>
    ...

Usage
-----
    python gen_anomalous_moment_table.py            # writes the table 
    python gen_anomalous_moment_table.py --check    # verifies it 

The quadrature has a noise floor of ~1e-6 in f at chi <~ 1e-4, below the
1 - f it would resolve there; another reason the table starts at 1e-3.

Requires numpy and scipy.
"""

import argparse
import os
import warnings

import numpy as np
from scipy.integrate import quad, IntegrationWarning
from scipy.special import beta, gamma

SQRT3 = np.sqrt(3.0)
GI0 = 1.0 / (3.0**(7.0 / 6.0) * gamma(2.0 / 3.0))
GIP0 = 1.0 / (3.0**(5.0 / 6.0) * gamma(1.0 / 3.0))
C0 = 2.0 * np.pi * GI0 * beta(5.0 / 3.0, 4.0 / 3.0)
C1 = 2.0 * np.pi * GIP0 * beta(7.0 / 3.0, 2.0 / 3.0)


def scorer_gi(z):
    """Scorer function Gi(z) for real z >= 0 (damped integral)."""
    val = quad(lambda t: np.exp(-t**3 / 3.0 - 0.5 * z * t)
               * np.cos(0.5 * SQRT3 * z * t + 2.0 * np.pi / 3.0),
               0.0, np.inf, limit=400, epsabs=0.0, epsrel=1e-13)[0]
    return -val / np.pi


def ratio(chi):
    """f(chi) = a(chi) / (alpha / 2 pi)."""
    def integrand(u):
        return (u**(2.0 / 3.0) / (1.0 + u)**3
                * scorer_gi((u / chi)**(2.0 / 3.0)))
    # Split where the integrand changes character: Gi's argument passes
    # through 1 at u = chi, and (1 + u)^-3 turns over at u = 1.
    edges = sorted({0.0, 1e-2 * chi, chi, 1e2 * chi, 1.0, np.inf})
    total = sum(quad(integrand, a, b, limit=800, epsabs=0.0,
                     epsrel=1e-12)[0] for a, b in zip(edges[:-1], edges[1:]))
    return 2.0 * np.pi * chi**(-2.0 / 3.0) * total


def tail(chi):
    """Two-term large-chi expansion (unscaled)."""
    y = chi**(-2.0 / 3.0)
    return y * (C0 + C1 * y)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--log-min', type=float, default=-3.0)
    p.add_argument('--log-max', type=float, default=4.0)
    p.add_argument('--per-decade', type=int, default=20)
    p.add_argument('--output', default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'anomalous_moment.table'))
    p.add_argument('--check', action='store_true',
                   help='verify limits, tail and interpolation error')
    args = p.parse_args()

    # The tolerances above are far tighter than the table needs; quad
    # occasionally reports round-off on the tiny outer segments, which
    # does not affect the result (checked against an independent Gi).
    warnings.simplefilter('ignore', IntegrationWarning)

    n = int(round((args.log_max - args.log_min) * args.per_decade)) + 1
    log_chi = np.linspace(args.log_min, args.log_max, n)
    log_f = np.log10([ratio(10.0**x) for x in log_chi])

    with open(args.output, 'w') as fh:
        fh.write('%d\t%.1f\t%.1f\n' % (n, args.log_min, args.log_max))
        for x, y in zip(log_chi, log_f):
            fh.write('%.16f\t%.16f\n' % (x, y))
    print('wrote %s (%d points, log10 chi %g..%g)'
          % (args.output, n, args.log_min, args.log_max))
    print('tail constants: c0 = %.16e  c1 = %.16e' % (C0, C1))
    print('f(chi_min) = %.10f   f(1) = %.6f   tail/table at chi_max = %.8f'
          % (10**log_f[0], ratio(1.0), tail(10**args.log_max) / 10**log_f[-1]))

    if args.check:
        mid = 0.5 * (log_chi[:-1] + log_chi[1:])
        exact = np.array([ratio(10.0**x) for x in mid])
        interp = 10.0**np.interp(mid, log_chi, log_f)
        print('max |interp/exact - 1| at midpoints: %.2e'
              % np.max(np.abs(interp / exact - 1.0)))
        s = 10**log_f[-1] / tail(10**args.log_max)
        for chi in (3e4, 1e5, 1e6):
            print('tail check chi=%g: |s*tail/exact - 1| = %.2e'
                  % (chi, abs(s * tail(chi) / ratio(chi) - 1.0)))
        for chi in (1e-3, 3e-4, 1e-4):
            print('small-chi check chi=%g: 1 - f = %.2e'
                  % (chi, 1 - ratio(chi)))


if __name__ == '__main__':
    main()
