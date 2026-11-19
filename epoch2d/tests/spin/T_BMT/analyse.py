#!/usr/bin/env python3
"""Compare EPOCH's T-BMT spin precession with independent references.

For every test*/run_<label>/ directory the script reads the dumps, takes
the initial momentum and spin of each species from dump 0, the uniform
fields from the grid, and the anomalous moment from the run's input.deck.
It then integrates the Lorentz force and the T-BMT equation together
(DOP853, rtol 1e-12) as the reference,

    du/dt = q/(m c) (E + c beta x B),            u = p/(m c)
    dS/dt = (q/m) S x W,
    W = (a + 1/g) B - a g/(g+1) (beta.B) beta - (a + 1/(g+1)) beta x E / c,

and overlays the closed-form predictions of each test.

When the deck has use_qed = T, species identified as electrons or
positrons use the chi-dependent moment a f(chi) in EPOCH
(TABLES/anomalous_moment.table). The reference then uses a f(chi0), with
f from the exact quadrature in TABLES/gen_anomalous_moment_table.py (not
the table) and chi0 from the initial state. That is exact only while chi
is constant, as in test 5 (v perp uniform B, E = 0).

Usage:
    python analyse.py              # figures for the hc/boris runs
    python analyse.py --scan       # also the dt-convergence figure
Figures are written to figures/.
"""

import glob
import os
import re
import sys

import numpy as np
import sdf
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# EPOCH's own values (src/constants.F90)
Q0 = 1.602176565e-19
M0 = 9.10938291e-31
C = 2.99792458e8
QE, ME = -Q0, M0          # every test species is an electron

# Schwinger field e_s (V/m), EPOCH's value
E_S = 1.323285417001326061279735961512150e18

HERE = os.path.dirname(os.path.abspath(__file__))
TABLES = os.path.join(HERE, '..', '..', '..', 'src', 'physics_packages',
                      'TABLES')
FIGDIR = os.path.join(HERE, 'figures')
PUSHERS = {'hc': ('Higuera-Cary', 'C0'), 'boris': ('Boris', 'C3')}


# ------------------------------------------------------------------ input

def deck_moments(deck):
    """Anomalous moment of each species, evaluating the constant block."""
    ns = {'pi': np.pi, 'sqrt': np.sqrt, 'micron': 1e-6, 'femto': 1e-15,
          'pico': 1e-12, 'me': M0, 'qe': Q0, 'c': C}
    text = open(deck).read()
    for blk in re.findall(r'begin:constant(.*?)end:constant', text, re.S):
        for line in blk.splitlines():
            line = line.split('#')[0].strip()
            if '=' in line:
                k, v = (s.strip() for s in line.split('=', 1))
                # EPOCH writes powers as ^
                ns[k] = eval(v.replace('^', '**'), {}, ns)
    out = {}
    for blk in re.findall(r'begin:species(.*?)end:species', text, re.S):
        name = re.search(r'^\s*name\s*=\s*(\S+)', blk, re.M).group(1)
        a = re.search(r'^\s*anomalous_magnetic_moment\s*=\s*(.+)$', blk,
                      re.M).group(1)
        out[name] = float(eval(a, {}, ns))
    return out


def deck_chi_dependent(deck):
    """Species that EPOCH gives a(chi): use_qed on and identified as an
    electron or positron (species_type electron/positron)."""
    text = open(deck).read()
    qed = re.search(r'^\s*use_qed\s*=\s*T', text, re.M | re.I) is not None
    out = {}
    for blk in re.findall(r'begin:species(.*?)end:species', text, re.S):
        name = re.search(r'^\s*name\s*=\s*(\S+)', blk, re.M).group(1)
        lepton = re.search(r'^\s*identify\s*:\s*\S*(electron|positron)\s*$',
                           blk, re.M | re.I) is not None
        out[name] = qed and lepton
    return out


def chi_of(u, E, B):
    """Quantum parameter chi for normalised momentum u (electron mass)."""
    g = np.sqrt(1.0 + u @ u)
    b = u / g
    er = E + C * np.cross(b, B)
    return g / E_S * np.sqrt(max(er @ er - (b @ E)**2, 0.0))


def f_exact(chi):
    """a(chi)/a(0) by quadrature, from the table generator."""
    sys.dont_write_bytecode = True
    sys.path.insert(0, TABLES)
    import gen_anomalous_moment_table as gen
    return gen.ratio(chi)


def f_table(chi):
    """a(chi)/a(0) as EPOCH looks it up (log-log linear in the table)."""
    tab = np.loadtxt(os.path.join(TABLES, 'anomalous_moment.table'),
                     skiprows=1)
    return 10.0**np.interp(np.log10(chi), tab[:, 0], tab[:, 1])


def load_run(rdir):
    """Time series of u and S (species means) plus the uniform fields."""
    files = sorted(glob.glob(os.path.join(rdir, '[0-9]*.sdf')))
    moments = deck_moments(os.path.join(rdir, 'input.deck'))
    steps, times = [], []
    u = {s: [] for s in moments}
    S = {s: [] for s in moments}
    spread = 0.0
    for f in files:
        d = sdf.read(f, dict=True)
        steps.append(d['Header']['step'])
        times.append(d['Header']['time'])
        for s in moments:
            uu = np.array([d['Particles/P%s/%s' % (c, s)].data
                           for c in 'xyz']) / (M0 * C)
            ss = np.array([d['Particles/Spin_%s/%s' % (c, s)].data
                           for c in 'xyz'])
            spread = max(spread, np.ptp(ss, axis=1).max())
            u[s].append(uu.mean(axis=1))
            S[s].append(ss.mean(axis=1))
        if f == files[0]:
            E = np.array([d['Electric Field/E%s' % c].data.mean()
                          for c in 'xyz'])
            B = np.array([d['Magnetic Field/B%s' % c].data.mean()
                          for c in 'xyz'])
    steps = np.array(steps)
    # Dump 0 holds the t = 0 state but its header says dt/2; later dumps
    # have time = step * dt, so rebuild every time from the step count
    k = np.argmax(steps)
    dt = times[k] / steps[k]
    # Effective moment for the reference: a f(chi0) where EPOCH uses a(chi)
    chi_dep = deck_chi_dependent(os.path.join(rdir, 'input.deck'))
    chi0 = {s: chi_of(np.array(u[s][0]), E, B) for s in moments}
    a_eff = {s: moments[s] * (f_exact(chi0[s]) if chi_dep[s] else 1.0)
             for s in moments}
    return dict(t=steps * dt, dt=dt, E=E, B=B, a=moments, a_eff=a_eff,
                chi_dep=chi_dep, chi0=chi0, spread=spread,
                u={s: np.array(v) for s, v in u.items()},
                S={s: np.array(v) for s, v in S.items()})


# -------------------------------------------------------------- reference

def reference(run, s):
    """Exact (to ODE tolerance) u(t), S(t) for species s of a run."""
    E, B, a = run['E'], run['B'], run['a_eff'][s]

    def rhs(_, y):
        u, S = y[:3], y[3:]
        g = np.sqrt(1.0 + u @ u)
        b = u / g
        du = QE / (ME * C) * (E + C * np.cross(b, B))
        W = ((a + 1.0 / g) * B - a * g / (g + 1.0) * (b @ B) * b
             - (a + 1.0 / (g + 1.0)) * np.cross(b, E) / C)
        return np.concatenate([du, QE / ME * np.cross(S, W)])

    y0 = np.concatenate([run['u'][s][0], run['S'][s][0]])
    sol = solve_ivp(rhs, (0.0, run['t'][-1]), y0, method='DOP853',
                    t_eval=run['t'], rtol=1e-12, atol=1e-14)
    return sol.y[:3].T, sol.y[3:].T


def runs_for(test, prefix=''):
    out = {}
    for label in PUSHERS:
        rdir = os.path.join(HERE, test, 'run_' + prefix + label)
        if glob.glob(os.path.join(rdir, '0000.sdf')):
            out[label] = load_run(rdir)
    return out


def gamma(u):
    return np.sqrt(1.0 + np.sum(u**2, axis=-1))


def plane_angle(S, axis, e1):
    """Unwrapped angle of S about `axis`, measured from e1."""
    e2 = np.cross(axis, e1)
    return np.unwrap(np.arctan2(S @ e2, S @ e1))


# ------------------------------------------------------------- plot helpers

def plot_components(ax, run, s, title):
    _, Sr = reference(run, s)
    t = run['t'] * 1e12
    for i, c in enumerate('xyz'):
        ax.plot(t, Sr[:, i], color='C%d' % i, lw=1.2, label='$S_%s$ ref' % c)
        ax.plot(t[::4], run['S'][s][::4, i], 'o', ms=3, mfc='none',
                color='C%d' % i)
    ax.set_xlabel('t (ps)')
    ax.set_ylabel('S')
    ax.set_ylim(-1.1, 1.1)
    ax.set_title(title, fontsize=10)


def plot_errors(ax, runs, species):
    styles = ['-', '--', ':', '-.']
    for label, run in runs.items():
        name, col = PUSHERS[label]
        for j, s in enumerate(species):
            _, Sr = reference(run, s)
            err = np.linalg.norm(run['S'][s] - Sr, axis=1)
            norm = np.abs(np.linalg.norm(run['S'][s], axis=1) - 1.0)
            # HC drawn wider underneath, so coincident curves stay visible
            ax.semilogy(run['t'][1:] * 1e12, np.maximum(err[1:], 1e-17),
                        styles[j % 4], color=col,
                        lw=3.0 if label == 'hc' else 1.2,
                        alpha=0.5 if label == 'hc' else 1.0,
                        label='%s, %s' % (name, s))
            if j == 0:
                ax.semilogy(run['t'][1:] * 1e12,
                            np.maximum(norm[1:], 1e-17), styles[j % 4],
                            color=col, lw=0.6, alpha=0.4)
    ax.set_xlabel('t (ps)')
    ax.set_ylabel(r'$|S_{EPOCH} - S_{ref}|$')
    ax.set_ylim(1e-17, 1)
    ax.set_title(r'error vs reference (faint: $||S|-1|$)', fontsize=10)
    ax.legend(fontsize=7)


def summary(test, runs, species):
    for label, run in runs.items():
        for s in species:
            _, Sr = reference(run, s)
            err = np.linalg.norm(run['S'][s] - Sr, axis=1).max()
            drift = np.linalg.norm(run['S'][s] - run['S'][s][0],
                                   axis=1).max()
            print('  %-6s %-13s a=%-10.4g max|S-S_ref|=%.2e  '
                  'max|S-S0|=%.3f  spread=%.1e'
                  % (label, s, run['a_eff'][s], err, drift, run['spread']))


# ------------------------------------------------------------------ tests

def test1():
    runs = runs_for('test1_B_perp_v')
    run = runs['hc']
    t = run['t']
    s = 'transverse'
    B = run['B'][2]
    g0 = gamma(run['u'][s][0])
    a = run['a'][s]
    omega = QE / ME * (a + 1.0 / g0) * B        # spin, about +z
    omega_c = QE * B / (g0 * ME)                # momentum, about +z
    z = np.array([0.0, 0.0, 1.0])
    phi_s = plane_angle(run['S'][s], z, np.array([1.0, 0, 0]))
    phi_u = plane_angle(run['u'][s], z, np.array([1.0, 0, 0]))

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    plot_components(ax[0, 0], run, s,
                    r'transverse spin, $S_0=\hat x$ (HC; lines ref, '
                    r'circles EPOCH)')
    ax[0, 0].plot(t * 1e12, run['S']['antiparallel'][:, 2], 'k-', lw=2,
                  label=r'$S_z$ antiparallel ($S_0=-\hat z$)')
    ax[0, 0].axvline(np.pi / abs(omega) * 1e12, color='grey', ls=':')
    ax[0, 0].legend(fontsize=7, loc='lower left')

    ax[0, 1].plot(t * 1e12, -omega * t, 'k-',
                  label=r'$-\Omega t$, $\Omega=(q/m)(a+1/\gamma)B$')
    ax[0, 1].plot(t[::4] * 1e12, phi_s[::4], 'o', ms=3, mfc='none',
                  label='EPOCH spin angle')
    ax[0, 1].axhline(np.pi, color='grey', ls=':')
    ax[0, 1].text(0.02, np.pi + 0.15, r'$\pi$: in-plane spin reversed',
                  fontsize=8)
    ax[0, 1].set_xlabel('t (ps)')
    ax[0, 1].set_ylabel('spin angle about z (rad)')
    ax[0, 1].set_title('precession angle is linear in time', fontsize=10)
    ax[0, 1].legend(fontsize=8)

    ax[1, 0].plot(t * 1e12, -(omega - omega_c) * t * 1e3, 'k-',
                  label=r'$-(q/m)\,aB\,t$')
    ax[1, 0].plot(t[::4] * 1e12, (phi_s - phi_u)[::4] * 1e3, 'o', ms=3,
                  mfc='none', label='EPOCH: spin angle - momentum angle')
    ax[1, 0].set_xlabel('t (ps)')
    ax[1, 0].set_ylabel('mrad')
    ax[1, 0].set_title('g-2: spin advances on the momentum at (q/m) a B',
                       fontsize=10)
    ax[1, 0].legend(fontsize=8)

    plot_errors(ax[1, 1], runs, ['transverse', 'antiparallel'])
    fig.suptitle(r'Test 1: $v \perp B$ (B = %g T $\hat z$, $\gamma$ = %g)'
                 % (B, g0))
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'test1_B_perp_v.png'), dpi=130)

    slope = np.polyfit(t, phi_s, 1)[0]
    print('Test 1: v perp B')
    print('  spin rate  EPOCH %.10e  analytic %.10e  rel diff %.1e'
          % (slope, -omega, slope / -omega - 1))
    slope = np.polyfit(t, phi_s - phi_u, 1)[0]
    print('  g-2 rate   EPOCH %.10e  analytic %.10e  rel diff %.1e'
          % (slope, -(omega - omega_c), slope / -(omega - omega_c) - 1))
    summary('test1', runs, ['transverse', 'antiparallel'])


def test2():
    runs = runs_for('test2_B_par_v')
    run = runs['hc']
    t = run['t']
    species = ['a_minus1', 'a_electron', 'a_one']
    Bv = run['B']
    B0 = np.linalg.norm(Bv)
    d = Bv / B0
    g0 = gamma(run['u']['a_one'][0])

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    plot_components(ax[0, 0], run, 'a_one',
                    r'$a=1$ (HC; lines ref, circles EPOCH)')
    ax[0, 0].legend(fontsize=7, loc='lower left')
    plot_components(ax[0, 1], run, 'a_minus1',
                    r'$a=-1$: first two terms cancel, S must not move')

    print('Test 2: v parallel B')
    for j, s in enumerate(species):
        a = run['a'][s]
        e1 = run['S'][s][0] / np.linalg.norm(run['S'][s][0])
        phi = plane_angle(run['S'][s], d, e1)
        omega = QE / ME * (1.0 + a) * B0 / g0
        ax[1, 0].plot(t * 1e12, -omega * t, '-', color='C%d' % j, lw=1)
        ax[1, 0].plot(t[::4] * 1e12, phi[::4], 'o', ms=3, mfc='none',
                      color='C%d' % j,
                      label=r'%s ($a$=%.4g)' % (s, a))
        if omega != 0.0:
            slope = np.polyfit(t, phi, 1)[0]
            print('  %-10s rate EPOCH %.10e  (q/m)(1+a)B/gamma %.10e  '
                  'rel diff %.1e' % (s, slope, -omega, slope / -omega - 1))
        else:
            print('  %-10s rate EPOCH %.3e  (expected 0)'
                  % (s, np.polyfit(t, phi, 1)[0]))
    ax[1, 0].set_xlabel('t (ps)')
    ax[1, 0].set_ylabel('spin angle about B (rad)')
    ax[1, 0].set_title(r'lines: $-(q/m)(1+a)B\,t/\gamma$; circles: EPOCH',
                       fontsize=10)
    ax[1, 0].legend(fontsize=8)
    plot_errors(ax[1, 1], runs, species)
    fig.suptitle(r'Test 2: $v \parallel B \parallel (1,2,2)/3$ '
                 r'(|B| = %g T, $\gamma$ = %g)' % (B0, g0))
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'test2_B_par_v.png'), dpi=130)
    summary('test2', runs, species)


def test3():
    runs = runs_for('test3_E_par_v')
    run = runs['hc']
    t = run['t']
    species = ['a_electron', 'a_one']
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    plot_components(ax[0], run, 'a_one',
                    r'$a=1$ (HC; lines ref, circles EPOCH)')
    ax[0].legend(fontsize=7)
    for s in species:
        ax[1].plot(t * 1e12, gamma(run['u'][s]), label=s)
    ax[1].set_xlabel('t (ps)')
    ax[1].set_ylabel(r'$\gamma$')
    ax[1].set_title('the electron is accelerated along E', fontsize=10)
    plot_errors(ax[2], runs, species)
    fig.suptitle(r'Test 3: $B=0$, $v \parallel E \parallel (1,2,2)/3$')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'test3_E_par_v.png'), dpi=130)
    print('Test 3: v parallel E')
    summary('test3', runs, species)


def test4_quadrature(run, s, n=200001):
    """Spin angle about z from the closed-form trajectory, by quadrature."""
    Ey = run['E'][1]
    a = run['a'][s]
    ux, uy0 = run['u'][s][0][:2]
    tt = np.linspace(0.0, run['t'][-1], n)
    uy = uy0 + QE * Ey * tt / (ME * C)
    g = np.sqrt(1.0 + ux**2 + uy**2)
    # Omega_z = -(q/m)(a + 1/(g+1)) (beta x E)_z / c; dphi/dt = -Omega_z
    rate = QE / ME * (a + 1.0 / (g + 1.0)) * (ux / g) * Ey / C
    phi = np.concatenate([[0.0], np.cumsum(0.5 * (rate[1:] + rate[:-1])
                                           * np.diff(tt))])
    return np.interp(run['t'], tt, phi)


def test4():
    runs = runs_for('test4_E_perp_v')
    run = runs['hc']
    t = run['t']
    species = ['a_electron', 'a_one', 'axial']
    z = np.array([0.0, 0.0, 1.0])
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    plot_components(ax[0, 0], run, 'a_one',
                    r'$a=1$, $S_0=\hat x$ (HC; lines ref, circles EPOCH)')
    ax[0, 0].plot(t * 1e12, run['S']['axial'][:, 2], 'k-', lw=2,
                  label=r'$S_z$ axial ($S_0=\hat z$)')
    ax[0, 0].legend(fontsize=7, loc='lower left')

    print('Test 4: v perp E')
    for j, s in enumerate(['a_electron', 'a_one']):
        phi = plane_angle(run['S'][s], z, np.array([1.0, 0, 0]))
        quad = test4_quadrature(run, s)
        ax[0, 1].plot(t * 1e12, quad, '-', color='C%d' % j, lw=1)
        ax[0, 1].plot(t[::4] * 1e12, phi[::4], 'o', ms=3, mfc='none',
                      color='C%d' % j, label=r'%s ($a$=%.4g)'
                      % (s, run['a'][s]))
        print('  %-10s final angle EPOCH %.10f  quadrature %.10f  '
              'diff %.1e' % (s, phi[-1], quad[-1], phi[-1] - quad[-1]))
    ax[0, 1].set_xlabel('t (ps)')
    ax[0, 1].set_ylabel('spin angle about z (rad)')
    ax[0, 1].set_title(r'lines: $\int (q/m)(a+\frac{1}{\gamma+1})'
                       r'\beta_x E_y/c\,dt$; circles: EPOCH', fontsize=10)
    ax[0, 1].legend(fontsize=8)

    ax[1, 0].plot(t * 1e12, gamma(run['u']['a_one']), 'k-',
                  label=r'$\gamma$')
    ax[1, 0].plot(t * 1e12, run['u']['a_one'][:, 0] /
                  gamma(run['u']['a_one']), 'k--', label=r'$\beta_x$')
    ax[1, 0].set_xlabel('t (ps)')
    ax[1, 0].set_title(r'E accelerates the electron: $\gamma$, $\beta$ '
                       r'vary, so the rate is not constant', fontsize=10)
    ax[1, 0].legend(fontsize=8)
    plot_errors(ax[1, 1], runs, species)
    fig.suptitle(r'Test 4: $B=0$, $v \perp E$ (E = %g V/m $\hat y$, '
                 r'$\gamma_0$ = %g)' % (run['E'][1],
                                        gamma(run['u']['a_one'][0])))
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'test4_E_perp_v.png'), dpi=130)
    summary('test4', runs, species)


def test5():
    runs = runs_for('test5_chi_dependent_moment')
    run = runs['hc']
    t = run['t']
    species = ['lepton', 'generic']
    z = np.array([0.0, 0.0, 1.0])
    x = np.array([1.0, 0.0, 0.0])
    B = run['B'][2]
    a = run['a']['lepton']
    chi0 = run['chi0']['lepton']
    fe, ft = f_exact(chi0), f_table(chi0)

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    plot_components(ax[0, 0], run, 'lepton',
                    r'lepton, $a_e f(\chi)$ (HC; lines ref, circles EPOCH)')
    ax[0, 0].legend(fontsize=7, loc='lower left')

    print('Test 5: chi-dependent anomalous moment')
    print('  chi0 = %.6f   f exact = %.6f   f table = %.6f'
          % (chi0, fe, ft))
    rates = {}
    for j, s in enumerate(species):
        g2 = (plane_angle(run['S'][s], z, x)
              - plane_angle(run['u'][s], z, x))
        rates[s] = np.polyfit(t, g2, 1)[0]
        expect = -QE / ME * run['a_eff'][s] * B
        ax[0, 1].plot(t * 1e15, -QE / ME * run['a_eff'][s] * B * t,
                      '-', color='C%d' % j, lw=1)
        ax[0, 1].plot(t[::4] * 1e15, g2[::4], 'o', ms=3, mfc='none',
                      color='C%d' % j, label='%s (%s)' % (
                          s, r'$a_e f(\chi)$' if run['chi_dep'][s]
                          else r'$a_e$'))
        print('  %-8s g-2 rate EPOCH %.6e  -(q/m) a_eff B %.6e  '
              'rel diff %.1e' % (s, rates[s], expect, rates[s] / expect - 1))
    ratio = rates['lepton'] / rates['generic']
    print('  rate ratio lepton/generic %.6f  vs f(chi0) %.6f  '
          'rel diff %.1e' % (ratio, fe, ratio / fe - 1))
    ax[0, 1].set_xlabel('t (fs)')
    ax[0, 1].set_ylabel('spin angle - momentum angle (rad)')
    ax[0, 1].set_title(r'g-2 phase; lines: $-(q/m)\,a_{eff}B\,t$',
                       fontsize=10)
    ax[0, 1].legend(fontsize=8)

    tab = np.loadtxt(os.path.join(TABLES, 'anomalous_moment.table'),
                     skiprows=1)
    ax[1, 0].semilogx(10**tab[:, 0], 10**tab[:, 1], 'k-', lw=1,
                      label='anomalous_moment.table')
    ax[1, 0].semilogx([chi0], [ratio], 'o', color='C3', ms=7,
                      label='EPOCH: g-2 rate ratio (%.5f)' % ratio)
    ax[1, 0].set_xlabel(r'$\chi$')
    ax[1, 0].set_ylabel(r'$a(\chi)/a(0)$')
    ax[1, 0].set_title('rate ratio lepton/generic vs the table', fontsize=10)
    ax[1, 0].legend(fontsize=8)
    plot_errors(ax[1, 1], runs, species)
    fig.suptitle(r'Test 5: $a(\chi)$ with QED on ($\chi$ = %.3f, '
                 r'$\gamma$ = %g, B = %.3g T)'
                 % (chi0, gamma(run['u']['lepton'][0]), B))
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'test5_chi_dependent_moment.png'),
                dpi=130)
    summary('test5', runs, species)


# ------------------------------------------------------------ dt scan

def scan():
    """Final-time spin error against dt for the non-trivial cases."""
    cases = [('test1_B_perp_v', 'transverse', 'test 1'),
             ('test2_B_par_v', 'a_one', r'test 2, $a=1$'),
             ('test4_E_perp_v', 'a_electron', r'test 4, $a=a_e$'),
             ('test4_E_perp_v', 'a_one', r'test 4, $a=1$'),
             ('test5_chi_dependent_moment', 'lepton',
              r'test 5, $a_e f(\chi)$')]
    fig, axs = plt.subplots(1, len(cases), figsize=(4 * len(cases), 4),
                            sharey=True)
    print('dt scan: final-time |S - S_ref|')
    for ax, (test, s, title) in zip(axs, cases):
        for label, (name, col) in PUSHERS.items():
            dts, errs = [], []
            for rdir in glob.glob(os.path.join(HERE, test,
                                               'run_%s_dx*' % label)):
                run = load_run(rdir)
                _, Sr = reference(run, s)
                dts.append(run['dt'])
                errs.append(np.linalg.norm(run['S'][s][-1] - Sr[-1]))
            if not dts:
                continue
            o = np.argsort(dts)
            dts, errs = np.array(dts)[o], np.array(errs)[o]
            p = np.polyfit(np.log(dts), np.log(errs), 1)[0]
            ax.loglog(dts * 1e15, errs, 'o-', color=col,
                      lw=3.0 if label == 'hc' else 1.2,
                      ms=9 if label == 'hc' else 5,
                      alpha=0.5 if label == 'hc' else 1.0,
                      label='%s (slope %.2f)' % (name, p))
            print('  %-15s %-11s %-6s slope %.2f  errors %s'
                  % (test, s, label, p,
                     ' '.join('%.1e' % e for e in errs)))
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('dt (fs)')
        ax.legend(fontsize=8)
    axs[0].set_ylabel(r'$|S - S_{ref}|$ at $t_{end}$')
    fig.suptitle('Convergence with time step (second order expected)')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'convergence.png'), dpi=130)


if __name__ == '__main__':
    os.makedirs(FIGDIR, exist_ok=True)
    test1()
    test2()
    test3()
    test4()
    test5()
    if '--scan' in sys.argv:
        scan()
