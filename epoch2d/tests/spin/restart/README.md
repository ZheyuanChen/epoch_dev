# Spin restart test

Status (7 October 2026): **PASS** with the restart read-back plus the
three `it_spin_*` iterators. The negative control fails, as designed.

## Purpose

A run restarted from a restart dump must continue exactly as the
uninterrupted run would have. With `-DSPIN`, that requires two things:

- the restart dump holds each particle's `spin(3)`;
- `restart_data` (`housekeeping/setup.F90`) reads it back for every
  species: lepton spin and photon polarisation alike.

The test is convention-agnostic. It compares a restarted run with the
uninterrupted one, so it holds whatever the deck check does to photon
polarisation vectors.

## Design

- **Domain:** 2D, 32 × 32 cells of 0.1 µm, periodic, 0–20 fs.
- **Field:** uniform B = (1000, 0, 2000) T. The precession axis is tilted,
  so every spin component changes.
- **QED:** on, so leptons take the a(χ) path, but with no emission and no
  recoil. The run is therefore deterministic.
- **Species:**

| species | particles | initial state |
|---|---|---|
| `electron` | 2000 | thermal, 100 keV; isotropic random spins (`spin = uniform`) |
| `positron` | 2000 | thermal, 100 keV; every spin along +y |
| `pol` | 1000 | photons at 10 MeV along (1, 1, 0), with a polarisation vector |
| `unpol` | 1000 | the same photons, with no polarisation settings |

Each lepton precesses at its own rate, so by the restart dump (dump 5,
t = 10 fs) the spins are spread over the sphere. A restart that lost them,
or reset them to any common value, cannot reproduce the later dumps.

`run_test.sh` makes three runs:

| run | ranks | what it is |
|---|---|---|
| `full` | 2 | dumps every 2 fs; 0000, 0005 and 0010 are restart dumps |
| `restart2` | 2 | `restart_snapshot = 5`, same decomposition |
| `restart1` | 1 | `restart_snapshot = 5`, so every particle is read onto another rank |

Particles carry no IDs. `analyse.py` therefore matches them by sorting
each species on position: the run is deterministic, so positions agree to
round-off and are unique.

## Checks (`analyse.py`)

| | check | pass |
|---|---|---|
| R0 | dump 0005 holds `Spin_x/y/z` for every species | — |
| R1 | matched positions, max \|Δx\|/cell, for every species and every dump 0006–0010 | < 1e-9 |
| R2 | momenta, max \|Δu\| | < 1e-9 |
| R3 | spins and polarisations, max \|ΔS\| | < 1e-12 |
| R4 | the lepton spins at dump 5 are spread out: the smallest per-component range is > 0.5, so a reset to a common value would fail R3 by > 0.25 | > 0.5 |

R1–R3 are judged separately for `restart2` and `restart1`.

## Running

```bash
cd epoch2d/tests/spin/restart
./run_test.sh hc                     # default binary ../../../bin/epoch2d
~/Desktop/Project_EPOCH/.venv_EPOCH/bin/python analyse.py hc
```

This takes a few seconds and writes:

- `run_hc/summary.txt`;
- `figures/restart_hc.png`, which shows max \|ΔS\| per dump and species,
  and the positron S_y histograms at dumps 5 and 10.

## Results (7 October 2026, scratch builds of `spin` + working tree)

| binary | R0 | R4 | R1–R3 `restart2` | R1–R3 `restart1` |
|---|---|---|---|---|
| read-back + `it_spin_*` | PASS | 0.62 | 0, 0, 0 (bit-identical) | 7e-14, 1.6e-15, 1.0e-15 |
| negative control: `setup.F90` at HEAD (no read-back) | PASS | 0.62 | 0, 0, **1.0** | 7e-14, 1.6e-15, **1.0** |

The negative control's R3 = 1.0 is the restored spins sitting at the
`init_particle` default (0, 0, 0). Positions and momenta still match,
which confirms that only spin is lost.

## Known issue: restart at 3–4 ranks (not spin)

On this deck, restart dumps written by 3 or 4 ranks hold corrupted field
blocks (Ex … Bz, Jx). Most cells hold zeros or denormal garbage, and
occasionally a NaN.

- Normal dumps at the same time are correct, and so is dump 0000.
- The restarted run then crashes (a NaN position at the cell lookup) or
  diverges.

This happens with and without `-DSPIN`, and also in an unmodified
`upstream/4.21-devel` build (default Makefile). Restarting at 2 ranks is
bit-identical there. The test therefore uses 1–2 ranks; the cause has not
been investigated yet.
