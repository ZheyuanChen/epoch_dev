# Pair-production macro-particle boosting — integration report

**Branch:** `laser-injection-and-pair-boost`
**Date:** 7 July 2026
**Base:** `upstream-pr-custom-laser-injection` (custom laser profile
injection) plus the CLAUDE-related project files.

This branch adds the `boost_pairs` / `boost_muons` macro-particle
splitting feature from the upstream `stu/muons_and_trident` branch of
[epochpic/epoch](https://github.com/epochpic/epoch), together with the
nuclear-trident and muon pair-production physics modules it depends on.
The feature plays the same role for the optical-depth pair-production
channels (nonlinear Breit-Wheeler, Bethe-Heitler, field trident,
nuclear trident) that `amplify_LBW_factor` already plays for the linear
Breit-Wheeler channel: it trades a few high-weight macro-pairs for many
low-weight ones without changing the expected physical yield.

---

## 1. What was integrated

Seven commits were cherry-picked from `upstream/stu/muons_and_trident`
(original authorship preserved):

| Commit (new) | Upstream | Description |
|---|---|---|
| `a4c6038a` | `5512ea4b` | Nuclear trident and muon pair production (epoch2d) |
| `ad9e1302` | `d2dd407d` | `ignore_dt_corrections` control-block key |
| `4e14c7b1` | `941ab6c9` | Port of the above to epoch1d and epoch3d |
| `d55cb864` | `b8183075` | **`boost_pairs` / `boost_muons` macro-particle boosting** |
| `3e9f3121` | `5eaf4846` | Fix: photon deallocation in boosted Breit-Wheeler |
| `fbb9a945` | `a610508d` | Fix: muon gamma sampling |
| `92552de6` | `8485e692` | Fix: optical-depth accumulation at high boost (epoch3d) |

Two follow-up commits of our own:

| Commit | Description |
|---|---|
| `a5716caf` | Port the optical-depth accumulation fix (`92552de6`, upstream PR #818) to epoch1d and epoch2d, which upstream had only applied to epoch3d |
| `3c373673` | Guard the `identify:muon` / `identify:antimuon` handlers in `deck_species_block.F90` with `#ifdef BREMSSTRAHLUNG`: they referenced species variables that only exist under that define, breaking any build without it. Upstream never noticed because their Makefile force-enables the define. |

### Conflict resolutions and deviations from upstream

- **Constants renumbered** (`constants.F90`, all three codes): upstream
  assigned `c_def_brem_trident = 2**28` and dump IDs 74–75, which
  collide with the `4.20-devel` additions already on this branch
  (`c_def_transition_rates = 2**28`, dump IDs 74–78). The new entries
  now use `2**29`/`2**30` and dump IDs 79–80. These are internal
  identifiers only; no user-visible behaviour changes.
- **Makefile defines left commented out**: upstream enabled
  `BREMSSTRAHLUNG`, `BREM_TRIDENT` and `BREM_MUON` by default. This
  branch keeps them commented out, consistent with `PHOTONS` /
  `TRIDENT_PHOTONS` and with our convention that optional physics is
  enabled per build (each define adds per-particle storage, so
  default-off avoids a silent memory-footprint increase for laser runs
  that do not use these packages).
- **SDF submodule pointer unchanged**: the upstream port commit moved
  the SDF submodule; that change was excluded.
- The `boost_pairs` deck keys were merged after the
  `use_*_recombination` keys added by `4.20-devel` in
  `deck_control_block.F90`.

### Build verification

All three codes were built successfully with `gfortran` (exit 0,
binaries produced):

```
make COMPILER=gfortran DEFINE="-DPHOTONS -DTRIDENT_PHOTONS \
    -DBREMSSTRAHLUNG -DBREM_TRIDENT -DBREM_MUON"
```

and epoch2d was rebuilt cleanly with no optional defines to confirm the
default configuration is unaffected.

---

## 2. User documentation: pair-yield boosting

### 2.1 The sampling problem

The pair-production cross sections are small. With physical rates, a
typical run produces very few pair events, and each event converts a
whole (often heavy) macro-photon or macro-electron into a single pair
of equal weight. The result is a handful of high-weight macro-pairs and
extremely poor statistics on positron spectra and angular
distributions.

Two independent variance-reduction controls now exist:

| Deck key | Block | Processes affected |
|---|---|---|
| `amplify_LBW_factor` | `qed` | Linear (two-photon) Breit-Wheeler only |
| `boost_pairs` | `control` | Nonlinear Breit-Wheeler, field trident, Bethe-Heitler, nuclear trident |
| `boost_muons` | `control` | Bethe-Heitler muon pair production |

They act on disjoint code paths and **do not compound**: using
`amplify_LBW_factor` together with `boost_pairs` is safe — each factor
tunes only the sampling of its own channel, and no event is ever
boosted by both.

### 2.2 How `boost_pairs` works

For every affected process, with `boost_pairs = N`:

1. the optical-depth decay rate of the parent particle is multiplied by
   `N`, so pair events trigger ~`N` times more often;
2. each event creates an electron and a positron with
   `weight = parent_weight / N` instead of the full parent weight;
3. the parent is destroyed (photon channels) or recoiled (nuclear
   trident) only with probability `1/N`; otherwise it survives with its
   remaining optical depth and can convert again.

The expected pair yield is therefore unchanged, but it is carried by up
to `N` times more macro-particles. Total weight is conserved on
average, not exactly per event — over-splitting a small number of
parent particles adds noise of its own, so increase `boost_pairs`
until the positron statistics converge rather than starting huge.
`boost_pairs = 1` (the default) reproduces the unboosted algorithm
exactly.

### 2.3 Building

Enable the packages you need at compile time (or uncomment the
corresponding `DEFINES` lines in the Makefile):

```bash
cd epoch2d          # or epoch1d / epoch3d
make COMPILER=gfortran DEFINE="-DPHOTONS -DTRIDENT_PHOTONS \
    -DBREMSSTRAHLUNG -DBREM_TRIDENT -DBREM_MUON"
```

| Define | Enables |
|---|---|
| `PHOTONS` | QED synchrotron emission + nonlinear Breit-Wheeler pairs |
| `TRIDENT_PHOTONS` | Field (virtual-photon) trident process |
| `BREMSSTRAHLUNG` | Bremsstrahlung photons + Bethe-Heitler pairs |
| `BREM_TRIDENT` | Nuclear trident pairs (needs `BREMSSTRAHLUNG`) |
| `BREM_MUON` | Bethe-Heitler muon pairs (needs `BREMSSTRAHLUNG`) |

`BREM_TRIDENT`/`BREM_MUON` each add one optical-depth REAL per
particle, so leave them out of builds that do not use them.

### 2.4 Deck configuration

**Control block** — the boost factors live here, *not* in the `qed`
block:

```
begin:control
  ...
  boost_pairs = 100        # e+e- channels (default 1 = off)
  boost_muons = 10         # muon channel  (default 1 = off)
  ignore_dt_corrections = T  # optional: skip the automatic dt
                             # reduction from the plasma frequency /
                             # laser period (useful for very large
                             # domains, e.g. converter-target runs)
end:control
```

**Nonlinear Breit-Wheeler** (photon + strong field → e⁺e⁻) — as
before, via the `qed` block:

```
begin:qed
  use_qed = T
  produce_photons = T
  produce_pairs = T          # NBW pair creation
  photon_dynamics = T        # required when produce_pairs = T
  qed_table_location = src/physics_packages/TABLES
end:qed
```

with `identify:breit_wheeler_electron` / `identify:breit_wheeler_positron`
(or `bw_electron`/`bw_positron`) species; plain `electron`/`positron`
species are used as a fallback.

**Field trident** (electron + strong field → e⁺e⁻, one-step) —
compile with `TRIDENT_PHOTONS`; pairs go to
`identify:trident_electron` / `identify:trident_positron` species
(fallback: first electron/positron species).

**Bethe-Heitler and nuclear trident** (nuclear-field processes in
high-Z targets) — via the `bremsstrahlung` block:

```
begin:bremsstrahlung
  use_bremsstrahlung = T
  produce_photons = T
  photon_energy_min = 100 * kev
  photon_dynamics = T        # photons must move to reach the nuclei
  use_bethe_heitler = T      # photon + nucleus -> e+e-
  use_brem_trident = T       # electron + nucleus -> e- + e+e-
  use_brem_muon = F          # photon + nucleus -> mu+mu-
end:bremsstrahlung
```

The target species must have `atomic_number` set (see
`epoch2d/example_decks/bethe_heitler_trident.deck` for a complete
100 MeV beam-into-gold example). Pair species:
`identify:bh_electron` / `identify:bh_positron` for Bethe-Heitler,
`identify:brem_trident_electron` / `identify:brem_trident_positron`
for nuclear trident (both fall back to the first electron/positron
species). Muons use `identify:muon` / `identify:antimuon` (or
`bh_muon`/`bh_antimuon`), see `epoch2d/example_decks/muons.deck`.

**Output block** — two new per-particle dump keys are available when
the corresponding define is enabled:

```
begin:output
  brem_trident_optical_depth = always
  bh_muon_optical_depth = always
end:output
```

### 2.5 Choosing a boost factor

- Start at `boost_pairs = 1` to establish the unboosted baseline, then
  raise it (10, 100, 1000, ...) until the positron yield/spectrum stops
  changing within noise. The example decks use 1–10; our LBW work used
  `amplify_LBW_factor = 1000`.
- The boost multiplies the per-step optical-depth decrement. Keep
  `N * delta_opdep` per timestep well below 1, otherwise a parent can
  be "due" several conversions in one step and the discrete algorithm
  saturates — the accumulation fix (`92552de6`/`a5716caf`) makes this
  regime more robust for Bethe-Heitler but cannot remove the bias
  entirely.
- Watch particle counts: every event adds two macro-particles and
  parents survive with probability `1 - 1/N`, so memory and load
  imbalance grow with the boost.

### 2.6 Physics references

- Nuclear trident: Bhabha, Proc. R. Soc. A **152**, 559 (1935);
  implementation follows Martinez *et al.*, Phys. Plasmas **26**,
  103109 (2019), Section II.D.
- Bethe-Heitler sampling uses Geant4 algorithms (see
  `bethe_heitler.F90`).

---

## 3. Notes and caveats

- **Muon modules**: `BREM_MUON`/`boost_muons` came bundled with the
  same commits and are fully integrated, but are opt-in at compile time
  and default-off at runtime; they add nothing to a build without the
  define.
- **Survival-path optical depth**: the accumulation fix is applied to
  Bethe-Heitler in all three codes. The equivalent survival paths in
  nonlinear Breit-Wheeler (`photons.F90`) and the muon module still
  redraw a fresh optical depth, as upstream wrote them; at very high
  `boost_pairs` this slightly under-counts multiple conversions per
  step for those channels.
- **`.vscode/`** remains untracked (editor-local configuration).
- The custom laser profile injection feature carried over unchanged
  from `upstream-pr-custom-laser-injection`; see
  `DOCUMENTATION_LASER_INJECTION.tex`.
