#!/usr/bin/env python3
"""Check the spin restart test (see README.md).

For each restarted run (run_<label>/restart2, restart1), every dump it
wrote (0006-0010) is compared with the same dump of the uninterrupted run
(run_<label>/full), particle by particle. Particles carry no IDs, so the
two dumps are matched by sorting each species on position; the run is
deterministic, so the positions agree to round-off and are unique.

Checks, for every species and every dump after the restart:
  R1  matched positions agree           max |dx| / cell     < 1e-9
  R2  momenta agree                     max |du|            < 1e-9 (*)
  R3  spins/polarisations agree         max |dS|            < 1e-12
and once:
  R0  the restart dump holds Spin_x/y/z for every species
  R4  the lepton spins at the restart dump are spread out (ptp > 0.5 in
      each component), so a restart that reset them to any common value
      would fail R3 by more than 0.25
(*) u = p / (m c) for leptons and p c / (10 MeV) for photons.

Usage:
    python analyse.py LABEL
Writes run_<label>/summary.txt and figures/restart_<label>.png; exits 1 if
any check fails.
"""

import glob
import os
import sys

import numpy as np
import sdf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# EPOCH's own values (src/constants.F90)
Q0 = 1.602176565e-19
M0 = 9.10938291e-31
C = 2.99792458e8

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, 'figures')
SPECIES = ('electron', 'positron', 'pol', 'unpol')
LEPTONS = ('electron', 'positron')
P_SCALE = {'electron': M0 * C, 'positron': M0 * C,
           'pol': 1.0e7 * Q0 / C, 'unpol': 1.0e7 * Q0 / C}
CELL = 0.1e-6
RESTART_DUMP = '0005.sdf'


def own_dumps(rdir):
    """Dumps a run wrote itself (not the linked restart dump)."""
    return [f for f in sorted(glob.glob(os.path.join(rdir, '[0-9]*.sdf')))
            if not os.path.islink(f)]


def species_state(d, sp):
    """Positions (2, n), momenta (3, n) and spins (3, n), sorted by
    position so that two runs can be compared particle by particle."""
    x, y = d['Grid/Particles/%s' % sp].data
    p = np.array([d['Particles/P%s/%s' % (c, sp)].data for c in 'xyz'])
    s = np.array([d['Particles/Spin_%s/%s' % (c, sp)].data for c in 'xyz'])
    order = np.lexsort((y, x))
    return (np.array([x, y])[:, order], p[:, order] / P_SCALE[sp],
            s[:, order])


class Checks:
    def __init__(self):
        self.rows = []

    def add(self, cid, text, value, ok):
        if not isinstance(value, str):
            value = '%.2e' % value
        self.rows.append((cid, text, value, bool(ok)))

    def report(self):
        lines = ['%-3s %-62s %-10s %s' % ('', 'check', 'value', 'result')]
        for cid, text, value, ok in self.rows:
            lines.append('%-3s %-62s %-10s %s'
                         % (cid, text, value, 'PASS' if ok else 'FAIL'))
        npass = sum(r[3] for r in self.rows)
        lines.append('%d/%d checks pass' % (npass, len(self.rows)))
        return '\n'.join(lines)


def check_restart_dump(full, ck):
    d = sdf.read(os.path.join(full, RESTART_DUMP), dict=True)
    missing = [sp for sp in SPECIES for c in 'xyz'
               if 'Particles/Spin_%s/%s' % (c, sp) not in d]
    ck.add('R0', '%s holds Spin_x/y/z for all species' % RESTART_DUMP,
           'missing %d' % len(missing) if missing else 'all', not missing)
    spread = min(float(np.ptp(species_state(d, sp)[2], axis=1).min())
                 for sp in LEPTONS)
    ck.add('R4', 'lepton spins at the restart dump: min component ptp',
           spread, spread > 0.5)
    return d


def compare(full, rdir, ck, tag):
    """R1-R3 for one restarted run; returns max |dS| per dump and species."""
    files = own_dumps(rdir)
    if not files:
        ck.add('R1', '%s: no dumps written' % tag, 'missing', False)
        return {}
    worst = {k: 0.0 for k in ('dx', 'du', 'ds')}
    history = {sp: [] for sp in SPECIES}
    counts_ok = True
    for f in files:
        dr = sdf.read(f, dict=True)
        df = sdf.read(os.path.join(full, os.path.basename(f)), dict=True)
        for sp in SPECIES:
            xr, ur, sr = species_state(dr, sp)
            xf, uf, sf = species_state(df, sp)
            if xr.shape != xf.shape:
                counts_ok = False
                history[sp].append(np.nan)
                continue
            worst['dx'] = max(worst['dx'], np.abs(xr - xf).max() / CELL)
            worst['du'] = max(worst['du'], np.abs(ur - uf).max())
            ds = float(np.abs(sr - sf).max())
            worst['ds'] = max(worst['ds'], ds)
            history[sp].append(ds)
    span = '%s-%s' % (os.path.basename(files[0])[:4],
                      os.path.basename(files[-1])[:4])
    ck.add('R1', '%s, dumps %s: matched positions, max |dx|/cell'
           % (tag, span), worst['dx'], counts_ok and worst['dx'] < 1e-9)
    ck.add('R2', '%s: momenta, max |du|' % tag, worst['du'],
           counts_ok and worst['du'] < 1e-9)
    ck.add('R3', '%s: spins/polarisations, max |dS|' % tag, worst['ds'],
           counts_ok and worst['ds'] < 1e-12)
    return history


def figure(histories, d5, full, label):
    os.makedirs(FIGDIR, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    floor = 1e-18
    for (tag, hist), ls in zip(histories.items(), ('-', '--')):
        for sp, col in zip(SPECIES, ('C0', 'C1', 'C2', 'C3')):
            v = np.maximum(np.array(hist[sp], dtype=float), floor)
            ax1.semilogy(range(6, 6 + len(v)), v, ls, marker='o',
                         color=col, label='%s %s' % (tag, sp))
    ax1.axhline(1e-12, color='k', lw=0.5, ls=':')
    ax1.set_xticks(range(6, 11))
    ax1.set_xlabel('dump')
    ax1.set_ylabel('max |S_restart - S_full|  (floor 1e-18 = exact)')
    ax1.set_title('restarted vs uninterrupted run')
    ax1.legend(fontsize=7, ncol=2)
    d10 = sdf.read(own_dumps(full)[-1], dict=True)
    bins = np.linspace(-1, 1, 41)
    for d, name, ls in ((d5, 'dump 5 (restart)', 'step'),
                        (d10, 'dump 10', 'stepfilled')):
        s = species_state(d, 'positron')[2]
        ax2.hist(s[1], bins=bins, histtype=ls, alpha=0.6, label=name)
    ax2.set_xlabel('positron S_y (all start at +1)')
    ax2.set_ylabel('count')
    ax2.set_title('spins are spread out by the restart time')
    ax2.legend(fontsize=8)
    fig.suptitle('Spin restart test [%s]' % label)
    fig.tight_layout()
    path = os.path.join(FIGDIR, 'restart_%s.png' % label)
    fig.savefig(path, dpi=130)
    return path


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    label = sys.argv[1]
    top = os.path.join(HERE, 'run_%s' % label)
    full = os.path.join(top, 'full')
    ck = Checks()
    d5 = check_restart_dump(full, ck)
    histories = {}
    for n in (2, 1):
        tag = 'restart%d' % n
        histories[tag] = compare(full, os.path.join(top, tag), ck, tag)
    text = ck.report()
    print(text)
    with open(os.path.join(top, 'summary.txt'), 'w') as fh:
        fh.write(text + '\n')
    print('figure: %s' % figure(histories, d5, full, label))
    sys.exit(0 if all(r[3] for r in ck.rows) else 1)


if __name__ == '__main__':
    main()
