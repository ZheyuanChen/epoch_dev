#!/usr/bin/env python3
"""Generate audited NLCS/NBW rates and inverse CDFs for EPOCH.

Run with system python3. Requires NumPy and SciPy. The new, versioned
interface is documented in POLARISED_LCFA_TABLES.md; production EPOCH
does not yet read it. Legacy tables are never overwritten.

NLCS components are h = qs*s_beta = +/-1, with final states summed.
NBW components are Q = xi1_obs = +/-1, with pair spins summed.
Totals include 1/(sqrt(3)*pi), but exclude alpha/(tau_c*energy).
The NBW grid uses EPOCH chi = conventional chi_gamma/2. Exponential
scaling retains tiny pair rates without imposing the legacy rate floor.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import alpha, c, hbar, m_e
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator
from scipy.special import gamma, kve

from gen_bessel_k_table import columns, intk_scaled, intk_scaled_cosh

HERE = Path(__file__).resolve().parent
NORM = np.pi * np.sqrt(3.0)
LN10 = np.log(10.0)
TAIL = 120.0
SIGNS = np.array([1.0, -1.0])


class SpecialFunctions:
    """Direct scaled K functions, smooth log interpolation of scaled I.

    Unlike a final-state table, the integration backend is constructed
    afresh from quadrature. Its interpolation is checked independently.
    """

    def __init__(self, log_min=-14.0, log_max=6.0):
        self.log_x = np.linspace(log_min, log_max,
                                 int(np.ceil((log_max-log_min)*200))+1)
        values = [intk_scaled(10.0**z) for z in self.log_x]
        self.interp = PchipInterpolator(self.log_x, np.log(values),
                                        extrapolate=False)

    def __call__(self, x):
        x = np.asarray(x)
        lx = np.log10(x)
        if np.any(x <= 0.0) or np.any(lx > self.log_x[-1]):
            raise ValueError('special-function argument outside backend grid')
        ik = np.exp(self.interp(np.maximum(lx, self.log_x[0])))
        # Adaptive quadrature can probe arbitrarily close to the regular
        # NLCS endpoint. The first omitted term here is O(x^(4/3)).
        xs = np.minimum(x, 10.0**self.log_x[0])
        small = np.exp(xs)*(np.pi/np.sqrt(3.0)
                            - 1.5*gamma(1/3)*2.0**(-2/3)*xs**(2/3))
        ik = np.where(lx < self.log_x[0], small, ik)
        return kve(1/3, x), kve(2/3, x), ik


def coordinate(process, chi, v):
    """Fraction, Jacobian, argument, and removable exponent offset.

    Both transformed densities are regular and resolved over the entire
    default legacy chi grids. NBW integrates the positive half and uses
    its exact sharing symmetry to reconstruct the full distribution.
    """
    sh, ch = np.sinh(v), np.cosh(v)
    if process == 'nlcs':
        a = 1.5 * chi / (1.0 + chi)
        u = a * sh**3
        fraction = u / (1.0 + u)
        jac = 3.0 * a * sh**2 * ch / (1.0 + u)**2
        x = sh**3 / (1.0 + chi)
        return fraction, jac, x, 0.0
    x0 = 4.0 / (3.0 * chi)  # chi_gamma = 2 chi_EPOCH
    f = 1.0 / (1.0 + x0)
    z = np.sqrt(f) * sh
    fraction = 0.5 * (1.0 + z / np.sqrt(1.0 + z*z))
    jac = np.sqrt(f) * ch / (2.0 * (1.0 + z*z)**1.5)
    x = x0 * (1.0 + z*z)
    return fraction, jac, x, x0


def upper_limit(process, chi):
    if process == 'nlcs':
        return np.arcsinh((TAIL * (1.0+chi))**(1/3))
    x0 = 4.0 / (3.0 * chi)
    return np.arcsinh(np.sqrt(TAIL * (1.0+x0) / x0))


def density(process, chi, v, special):
    """Two positive component densities, with exp(-x0) removed for NBW."""
    l, jac, x, x0 = coordinate(process, chi, np.asarray(v))
    k13, k23, ik = special(x)
    if process == 'nlcs':
        # Compute s1 without subtracting nearly equal numbers.
        s1 = 1.0 / (1.0 + 1.5*chi/(1.0+chi)*np.sinh(v)**3)
        base, correction = (s1+1.0/s1)*k23-ik, l*k13
    else:
        # rho = 2 + 4*z^2 avoids loss of precision near equal sharing.
        z2 = np.sinh(v)**2 / (1.0+x0)
        base, correction = (2.0+4.0*z2)*k23+ik, -k23
    # For NBW form x-x0 directly: subtracting x0 loses small exponents.
    exponent = x if process == 'nlcs' else x0/(1.0+x0)*np.sinh(v)**2
    return ((base[..., None] + correction[..., None]*SIGNS)
            * (np.exp(-exponent)*jac/NORM)[..., None])


def make_row(process, chi, probabilities, special, cells=512):
    """Integrate positive cell masses and invert their cumulative sums."""
    edges = np.linspace(0.0, upper_limit(process, chi), cells+1)
    nodes, weights = np.polynomial.legendre.leggauss(8)
    mid, half = (edges[1:]+edges[:-1])/2, np.diff(edges)/2
    v = mid[:, None] + half[:, None]*nodes
    values = density(process, chi, v, special)
    masses = np.einsum('cni,n,c->ci', values, weights, half)
    if not np.all(np.isfinite(masses)) or np.any(masses < 0.0):
        raise RuntimeError('non-finite or negative integration cell mass')
    cumulative = np.vstack((np.zeros(2), np.cumsum(masses, axis=0)))
    total = cumulative[-1].copy()
    if np.any(total <= 0.0):
        raise RuntimeError('non-positive component total')
    cdf = cumulative / total
    quantiles = np.empty((len(probabilities), 2))
    for component in range(2):
        keep = np.r_[True, np.diff(cdf[:, component]) > 0.0]
        inverse = PchipInterpolator(cdf[keep, component], edges[keep],
                                     extrapolate=False)
        if process == 'nlcs':
            vq = inverse(probabilities[1:-1])
            quantiles[1:-1, component] = coordinate(process, chi, vq)[0]**(1/3)
        else:
            half_p = np.abs(2.0*probabilities[1:-1]-1.0)
            vq = inverse(half_p)
            high = coordinate(process, chi, vq)[0]
            quantiles[1:-1, component] = np.where(
                probabilities[1:-1] >= 0.5, high, 1.0-high)
    quantiles[0], quantiles[-1] = 0.0, 1.0
    if process == 'nbw':
        # Symmetrise round-off in mirrored probabilities. The table
        # probability grid is itself exactly symmetric to round-off.
        quantiles = 0.5*(quantiles+1.0-quantiles[::-1])
    if (not np.all(np.isfinite(quantiles))
            or np.any(np.diff(quantiles, axis=0) < -2e-15)):
        raise RuntimeError('invalid inverse CDF')
    if process == 'nbw':
        total *= 2.0
    # Natural logs avoid exponent underflow; files store decimal logs.
    log_scaled = np.log(total)
    x0 = 0.0 if process == 'nlcs' else 4.0/(3.0*chi)
    return log_scaled-x0, log_scaled, quantiles


def generate(process, log_chi, probabilities, special, cells=512):
    rows = [make_row(process, 10.0**z, probabilities, special, cells)
            for z in log_chi]
    logs, scaled, quantiles = map(np.asarray, zip(*rows))
    return dict(process=process, log_chi=np.asarray(log_chi),
                probabilities=probabilities, log_rates=logs,
                log_scaled_rates=scaled, quantiles=quantiles)


def bracket(grid, value):
    if not np.isfinite(value) or value < grid[0] or value > grid[-1]:
        raise ValueError('lookup outside table domain; no extrapolation')
    i = min(np.searchsorted(grid, value, side='right')-1, len(grid)-2)
    i = max(i, 0)
    return i, (value-grid[i])/(grid[i+1]-grid[i])


def lookup(table, chi, state, choose_u, energy_u):
    """Executable reader reference: (log dimensionless rate, fraction).

    NLCS state = qs*s_beta; NBW state = xi1_obs. Component choice and
    inverse-CDF sampling require TWO independent uniform variates.
    """
    if not np.isfinite(chi) or chi <= 0.0 or not -1.0 <= state <= 1.0:
        raise ValueError('require chi > 0 and physical state in [-1,1]')
    if not 0.0 <= choose_u < 1.0 or not 0.0 <= energy_u < 1.0:
        raise ValueError('random variates must be in [0,1)')
    i, f = bracket(table['log_chi'], np.log10(chi))
    rates = ((1.0-f)*table['log_scaled_rates'][i]
             + f*table['log_scaled_rates'][i+1])
    x0 = 0.0 if table['process'] == 'nlcs' else 4.0/(3.0*chi)
    weights = np.array([(1.0+state)/2.0, (1.0-state)/2.0])
    log_weights = np.full(2, -np.inf)
    np.log(weights, out=log_weights, where=weights > 0.0)
    log_components = rates+log_weights
    log_total_scaled = np.logaddexp(*log_components)
    plus_probability = np.exp(log_components[0]-log_total_scaled)
    component = 0 if choose_u < plus_probability else 1
    p, g = bracket(table['probabilities'], energy_u)
    q = table['quantiles']
    value = ((1-f)*((1-g)*q[i, p, component]+g*q[i, p+1, component])
             + f*((1-g)*q[i+1, p, component]+g*q[i+1, p+1, component]))
    fraction = value**3 if table['process'] == 'nlcs' else value
    return log_total_scaled-x0, fraction


def reference(process, chi, fraction=None):
    """Independent adaptive quadrature, using the cosh I representation."""
    def exact(x):
        return kve(1/3, x), kve(2/3, x), intk_scaled_cosh(float(x))

    end = upper_limit(process, chi)
    totals = np.array([quad(lambda v: density(process, chi, v, exact)[c],
                            0.0, end, epsabs=0.0, epsrel=2e-10,
                            limit=250)[0] for c in range(2)])
    if process == 'nbw':
        totals *= 2.0
    if fraction is None:
        return totals
    if fraction <= 0.0:
        return np.zeros(2)
    if fraction >= 1.0:
        return np.ones(2)
    if process == 'nlcs':
        a = 1.5*chi/(1.0+chi)
        end_cdf = np.arcsinh((fraction/(a*(1.0-fraction)))**(1/3))
        sign, offset = 1.0, 0.0
    else:
        x0 = 4.0/(3.0*chi)
        t = abs(2.0*fraction-1.0)
        end_cdf = np.arcsinh(t*np.sqrt((1.0+x0)/(1.0-t*t)))
        sign, offset = (1.0 if fraction >= 0.5 else -1.0), 0.5
    part = np.array([quad(lambda v: density(process, chi, v, exact)[c],
                          0.0, min(end_cdf, end), epsabs=0.0,
                          epsrel=2e-9, limit=250)[0] for c in range(2)])
    return offset+sign*part/totals


def validate(tables, special):
    """Off-grid rates and quantiles, extremes, symmetries and zero weights."""
    worst_rate = worst_cdf = 0.0
    x = np.logspace(special.log_x[0]+0.013,
                    special.log_x[-1]-0.013, 101)
    error = max(abs(float(special(v)[2])/intk_scaled_cosh(v)-1.0) for v in x)
    if error > 2e-7:
        raise RuntimeError('special-function interpolation error %.3g' % error)
    for table in tables:
        grid = table['log_chi']
        indices = np.unique(np.linspace(0, len(grid)-2, 8, dtype=int))
        probes = np.r_[grid[0], grid[-1], (grid[indices]+grid[indices+1])/2]
        for z in probes:
            chi = 10.0**z
            exact = reference(table['process'], chi)
            x0 = 0.0 if table['process'] == 'nlcs' else 4.0/(3.0*chi)
            for state in (-1.0, 0.0, 0.37, 1.0):
                log_ref = np.log(np.dot([(1+state)/2, (1-state)/2], exact))-x0
                log_rate, _ = lookup(table, chi, state, 0.31, 0.43)
                worst_rate = max(worst_rate, abs(np.expm1(log_rate-log_ref)))
            # Pure components permit a direct inverse-CDF check. This is
            # stronger than merely checking monotonicity of written data.
            for component, state in enumerate((1.0, -1.0)):
                for p in (1e-5, 0.01, 0.13, 0.5, 0.87, 0.99, 1-1e-5):
                    _, fraction = lookup(table, chi, state, 0.31, p)
                    actual = reference(table['process'], chi, fraction)[component]
                    worst_cdf = max(worst_cdf, abs(actual-p))
        if not np.all(np.isfinite(table['log_rates'])):
            raise RuntimeError('non-finite logarithmic rates')
        q = table['quantiles']
        if table['process'] == 'nbw' and np.max(abs(q+q[:, ::-1]-1.0)) > 2e-13:
            raise RuntimeError('NBW sharing symmetry failed')
    # Native Ptarmigan fixtures, independently pinned during the audit.
    # These prevent a shared prefactor/chi error in the reference and
    # generation paths from passing a self-consistent quadrature check.
    tau_c = hbar/(m_e*c*c)
    fixtures = [('nlcs', chi, target*tau_c*1000/alpha) for chi, target in
                [(0.026, 2.07935e14), (3.5, 1.58485e16),
                 (9.98, 3.45844e16), (12.4, 4.04647e16)]]
    fixtures += [('nbw', chi/2, target*chi) for chi, target in
                 [(0.042, 6.07753899492993e-29),
                  (0.105, 2.1082097875655205e-12),
                  (0.42, 0.0003779613236658133),
                  (1.05, 0.015977478443872017),
                  (4.2, 0.08917816786414409),
                  (12.0, 0.10884579479913804),
                  (42.0, 0.09266735324318656)]]
    fixture_error = 0.0
    fixture_count = 0
    for process, chi, target in fixtures:
        table = next(t for t in tables if t['process'] == process)
        if table['log_chi'][0] <= np.log10(chi) <= table['log_chi'][-1]:
            log_rate, _ = lookup(table, chi, 0.0, 0.31, 0.43)
            fixture_error = max(fixture_error,
                                abs(np.expm1(log_rate-np.log(target))))
            fixture_count += 1
    print('independent checks: special functions %.3g, rate %.3g, CDF %.3g'
          % (error, worst_rate, worst_cdf), flush=True)
    print('Ptarmigan absolute fixtures: %d checked, max relative error %.3g'
          % (fixture_count, fixture_error), flush=True)
    if worst_rate > 1e-3 or worst_cdf > 3e-4:
        raise RuntimeError('lookup accuracy failed; refine chi/probability grids')
    if fixture_error > 1e-3:
        raise RuntimeError('independent absolute normalisation check failed')
    return dict(special_relative_error=error, rate_relative_error=worst_rate,
                cdf_absolute_error=worst_cdf,
                native_fixture_count=fixture_count,
                native_fixture_relative_error=fixture_error)


def legacy_comparison(tables):
    """Compare absolute unpolarised totals; never fit to legacy artefacts."""
    result = {}
    for table, filename in zip(tables, ('hsokolov.table', 'pairprod.table')):
        old = np.loadtxt(HERE/filename, skiprows=1)
        samples = []
        for chi in ((0.01, 0.1, 1.0, 10.0, 100.0) if table['process'] == 'nlcs'
                    else (0.05, 0.1, 0.5, 1.0, 10.0)):
            if not table['log_chi'][0] <= np.log10(chi) <= table['log_chi'][-1]:
                continue
            log_rate, _ = lookup(table, chi, 0.0, 0.3, 0.5)
            if table['process'] == 'nlcs':
                legacy_log = np.interp(np.log10(chi), old[:, 0], old[:, 1])
                legacy_log = LN10*legacy_log+np.log(chi*np.sqrt(3)/(2*np.pi))
            else:
                legacy_log = LN10*np.interp(np.log10(chi), old[:, 0], old[:, 2])
                legacy_log += np.log(chi)
            samples.append(dict(chi_epoch=chi,
                                new_over_legacy=float(np.exp(log_rate-legacy_log))))
        result[table['process']] = samples
        print('%s unpolarised / legacy: %s' % (table['process'],
              ', '.join('%g: %.6f' % (v['chi_epoch'], v['new_over_legacy'])
                        for v in samples)), flush=True)
    return result


def write_tables(output, tables, special, checks, legacy):
    output.mkdir(parents=True, exist_ok=True)
    files = []
    for table in tables:
        name = table['process']
        grid, p = table['log_chi'], table['probabilities']
        path = output/(name+'_rates.table')
        with path.open('w') as f:
            f.write('%d %.17e %.17e\n' % (len(grid), grid[0], grid[-1]))
            np.savetxt(f, np.column_stack((grid, table['log_rates']/LN10,
                                           table['log_scaled_rates']/LN10)),
                       fmt='%.17e')
        files.append(path)
        path = output/(name+'_quantiles.table')
        with path.open('w') as f:
            f.write('%d %d %.17e %.17e\n' % (len(grid), len(p), grid[0], grid[-1]))
            # One probability per record; followed by one record per
            # (chi, probability), with plus then minus components.
            np.savetxt(f, p[:, None], fmt='%.17e')
            np.savetxt(f, table['quantiles'].reshape(-1, 2), fmt='%.17e')
        files.append(path)
    path = output/'bessel_k.table'
    grid = np.linspace(-6.0, 2.5, 851)
    bessel = np.array([columns(10.0**z) for z in grid])
    with path.open('w') as f:
        f.write('851 -6 2.5\n')
        np.savetxt(f, np.column_stack((grid, bessel)), fmt='%.17e')
    files.append(path)
    manifest = dict(format_version=1, kernels='POLARISED_LCFA_AUDIT.md',
                    rate_factor='alpha/(tau_c*gamma_or_photon_energy_mc2)',
                    nlcs_state='qs*s_beta', nbw_state='xi1_obs',
                    nbw_chi='chi_EPOCH=chi_gamma/2',
                    nlcs_quantile='lambda^(1/3)', nbw_quantile='positron_fraction',
                    tail_exponent=TAIL, checks=checks,
                    legacy_comparison=legacy,
                    special_log_bounds=special.log_x[[0, -1]].tolist(),
                    sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in files},
                    generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


def read_table(output, process):
    """List-directed-compatible round-trip reader used by the checks."""
    output = Path(output)
    if process not in ('nlcs', 'nbw'):
        raise ValueError('unknown process')
    with (output/(process+'_rates.table')).open() as f:
        rn, rlo, rhi = map(float, f.readline().split())
        rates = np.loadtxt(f)
    with (output/(process+'_quantiles.table')).open() as f:
        n, nq, lo, hi = map(float, f.readline().split())
        n, nq = int(n), int(nq)
        probabilities = np.array([float(f.readline()) for _ in range(nq)])
        q = np.loadtxt(f).reshape(n, nq, 2)
    if (n < 2 or nq < 2 or rn != n or rlo != lo or rhi != hi
            or rates.shape != (n, 5) or rates[0, 0] != lo
            or rates[-1, 0] != hi or not np.all(np.isfinite(rates))
            or not np.all(np.diff(rates[:, 0]) > 0.0)
            or not np.all(np.isfinite(probabilities))
            or probabilities[0] != 0.0 or probabilities[-1] != 1.0
            or not np.all(np.diff(probabilities) > 0.0)
            or not np.all(np.isfinite(q)) or np.any(q < 0) or np.any(q > 1)
            or np.any(q[:, 0] != 0) or np.any(q[:, -1] != 1)
            or np.any(np.diff(q, axis=1) < 0)):
        raise ValueError('inconsistent table header or probability grid')
    return dict(process=process, log_chi=rates[:, 0], probabilities=probabilities,
                log_rates=rates[:, 1:3]*LN10,
                log_scaled_rates=rates[:, 3:5]*LN10, quantiles=q)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--output-dir', type=Path, default=HERE/'polarised')
    parser.add_argument('--chi-per-decade', type=int,
                        help='replace legacy chi grids with uniform grids')
    parser.add_argument('--nlcs-log-min', type=float, default=-5.0)
    parser.add_argument('--nlcs-log-max', type=float, default=5.0)
    parser.add_argument('--nbw-log-min', type=float, default=-5.0)
    parser.add_argument('--nbw-log-max', type=float, default=8.0)
    parser.add_argument('--quantiles', type=int, default=1025)
    parser.add_argument('--cells', type=int, default=512)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.output_dir.resolve() == HERE or (args.output_dir/'hsokolov.table').exists():
        parser.error('use a separate output directory to preserve legacy tables')
    if args.quantiles < 65 or args.cells < 64 or (args.chi_per_decade is not None
                                               and args.chi_per_decade < 1):
        parser.error('require quantiles >= 65, cells >= 64, chi-per-decade >= 1')
    grids = []
    for process, legacy in (('nlcs', 'hsokolov.table'), ('nbw', 'pairprod.table')):
        lo, hi = getattr(args, process+'_log_min'), getattr(args, process+'_log_max')
        if not np.isfinite(lo+hi) or not -8.0 <= lo < hi <= 10.0:
            parser.error('require -8 <= log-min < log-max <= 10')
        old = np.loadtxt(HERE/legacy, skiprows=1)[:, 0]
        if args.chi_per_decade is None and lo == old[0] and hi == old[-1]:
            grid = old
        else:
            ppd = args.chi_per_decade or 100
            grid = np.linspace(lo, hi, int(np.ceil((hi-lo)*ppd))+1)
        grids.append(grid)
    # Bounds include all GL nodes and the highest NBW central argument.
    special = SpecialFunctions(log_max=max(6.0, -grids[1][0]+1.0))
    probabilities = np.sin(np.linspace(0.0, np.pi/2, args.quantiles))**2
    probabilities[0], probabilities[-1] = 0.0, 1.0
    tables = []
    for process, grid in zip(('nlcs', 'nbw'), grids):
        print('generating %s: %d chi rows, %d quantiles' %
              (process, len(grid), args.quantiles), flush=True)
        tables.append(generate(process, grid, probabilities, special, args.cells))
    checks = validate(tables, special) if args.check else None
    legacy = legacy_comparison(tables)
    write_tables(args.output_dir, tables, special, checks, legacy)
    for table in tables:
        restored = read_table(args.output_dir, table['process'])
        if not np.array_equal(restored['quantiles'], table['quantiles']):
            raise RuntimeError('written inverse-CDF table failed round trip')
    print('wrote five tables and manifest to %s' % args.output_dir)


if __name__ == '__main__':
    main()
