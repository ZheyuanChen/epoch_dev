# Interior laser antenna (Huygens/TFSF injection plane) -- Stage 2A

**Status**: implemented and validated in `epoch2d` on branch
`interior-laser-antenna` (diverged from `laser-injection-and-pair-boost`).
Builds clean (`gfortran`), all Stage 2A validation targets met. epoch3d,
oblique incidence, a paraxial Gaussian beam and a full SCiPIC comparison
are explicit follow-on stages, not attempted here.

## 1. Motivation

EPOCH's existing laser injection (`simple_laser` / the custom-file
spatiotemporal pipeline in `laser.f90` + `custom_laser.f90`) works at a
domain **boundary**: it prescribes tangential E and constructs B itself
via a fixed one-way relation. This is a structural ceiling for
near-unity-NA 3D work: the boundary source term is effectively
first-order in incidence angle, and it cannot accept a genuinely
independent B/H field or reproduce a real longitudinal field component
(e.g. the TM01 longitudinal field a companion project, SCiPIC, already
computes).

Investigation before implementing found:

- Upstream EPOCH's `5.0-devel` branch has a `begin:antenna` block, but it
  is a **plain volumetric current source** (`generate_antennae_currents`
  adds parser-expression `jx/jy/jz` every step, called once per step
  alongside particle-current deposition). It radiates symmetrically in
  both directions and has **no backward-wave cancellation** -- it is not
  a Huygens surface and does not solve this problem on its own, though
  its deck-block scaffolding pattern was a useful reference.
- No TFSF/Huygens-surface code exists anywhere in this repo's history or
  any branch.
- This project's own prior high-NA-injection discussion (2026-07-02,
  recorded in memory and `TIGHT_FOCUSING_INJECTION_PROMPT.md`)
  independently reached the same conclusion: a current-sheet antenna is
  the structural fix, not a higher-order rewrite of the boundary BC.

The implemented mechanism is a genuine total-field/scattered-field (TFSF)
correction applied directly to the interior Yee-grid curl updates at one
plane `x = x0`, not a volumetric current source and not a copy of the
boundary characteristic-BC mechanism.

## 2. Design summary

- New, distinct deck block `begin:laser_antenna` (not an extension of
  `begin:laser`): a plane position, a direction (`x_min`/`x_max`), and
  **four independent tangential incident-field files**
  (`ey_inc_file`/`ez_inc_file`/`by_inc_file`/`bz_inc_file`), each a raw
  binary spatiotemporal `(n_transverse_points x n_t_points)` array --
  same convention as the existing boundary-laser profile/phase files
  (`access='stream'`, no embedded header, shape/bounds deck-declared).
- New `antenna_block` type (`shared_data.F90`), a linked list `antennas`
  (mirrors `lasers`).
- New module `laser_antenna.f90`: setup/loading, MPI ownership
  (rank-local index arithmetic, generalising correctly to a plane that
  straddles two ranks), and the two correction subroutines.
- `custom_laser.f90` was refactored (behaviour-preserving) to extract two
  generic routines -- `load_spatiotemporal_matrix` (allocate + load) and
  `sample_spatiotemporal_matrix` (bilinear (pos, t) interpolation) -- so
  the antenna's four channels reuse the exact same, already-validated
  file-loading/interpolation code the boundary laser uses, rather than
  duplicating it. Verified behaviour-preserving before adding any antenna
  code: clean rebuild, then an existing-feature smoke test (two-file
  custom-profile boundary laser) loaded and ran identically.
- Corrections are hooked into `fields.f90` immediately after every
  `CALL update_e_field` / `CALL update_b_field`, in **both**
  `update_eb_fields_half` and `update_eb_fields_final` (each is called
  twice per full timestep in EPOCH's symmetric-split scheme) -- placed
  *before* the following `efield_bcs`/`bfield_bcs`/`bfield_final_bcs`
  call, so corrected values reach neighbour ranks' ghost cells through
  the ordinary halo exchange with **no new communication**.
- The correction was derived from exactly one branch of
  `update_e_field`/`update_b_field`: `field_order = 2`,
  `maxwell_solver = yee`, `cpml_boundaries = F` (the plain nearest-
  neighbour Yee curl). `field_order = 4/6` reach wider neighbours; the
  Lehe/Lehe_x/Lehe_y/Pukhov/Cowan solvers and the `cpml_boundaries`
  branch all use different stencil coefficients (`alphax`/`betaxy`/
  `deltax` cross terms) the correction knows nothing about. None of these
  combinations would crash on their own -- they would silently
  under-cancel or misapply the correction. `laser_antenna_block_check`
  aborts explicitly (via a direct `abort_code` call, not the usual IOR'd
  errcode return -- see below) if any `laser_antenna` is present with
  `field_order /= 2`, `maxwell_solver /= yee`, or `cpml_boundaries = T`.
  Both defaults already match what the correction assumes, so this only
  bites a deck that deliberately opts into a different solver (e.g. Lehe,
  commonly chosen specifically to suppress numerical Cherenkov radiation
  from a relativistic particle beam -- worth remembering if a future
  campaign needs both that and this antenna at once; extending the
  correction to cover Lehe is a bounded, derivable follow-on, not
  attempted here). Implemented as a direct `abort_code` call rather than
  the usual IOR'd errcode return, because tracing `handle_deck_element`
  showed that `check_compulsory_blocks`' aggregate errcode is never
  upgraded to the hard-terminate bit the way per-element parse errors
  are -- an informational-only return here would print the error and
  silently continue running with an uncancelled seam. Verified: test
  decks with `field_order = 4` and separately with `maxwell_solver =
  lehe` both abort cleanly with the expected message; the normal
  (`yee`, order 2, no CPML) case is unaffected (re-ran test 6, identical
  0.86% residual).

### Files created/modified (epoch2d only)

- **New**: `src/laser_antenna.f90`, `src/deck/deck_laser_antenna_block.f90`.
- **Modified**: `src/user_interaction/custom_laser.f90` (generic loader/
  sampler extraction), `src/shared_data.F90` (`antenna_block` type,
  `antennas`/`n_antennas`), `src/housekeeping/setup.F90` (`NULLIFY
  (antennas)`), `src/fields.f90` (correction hooks, `antenna_t_mid`),
  `src/deck/deck.F90` (six mirrored dispatch sites), `src/housekeeping/
  finish.f90` (`CALL deallocate_antennas`), `Makefile` (new source files
  + dependency lines).

## 3. The correction equations -- independently re-derived and verified

Yee staggering (confirmed by direct inspection of `fields.f90`'s
`update_e_field`/`update_b_field`, not assumed): `ey`/`ez` sit at integer
x-grid points; `by`/`bz` sit at half-integer ("B-grid") x-points, with
`by(k,:)`/`bz(k,:)` at physical position `x_min + k*dx` (derived directly
from `setup_grid`'s `x_global`/`dx` construction -- **not** from the
`xb`/`xb_global` arrays, whose own index convention turned out to be
offset from this by a full cell and would have given a silently-wrong
snap if trusted blindly).

`x0` is snapped to the nearest such B-grid point (`i0_global-1`); the
first E-grid point strictly east of the seam is `i0_global`.

For `direction = x_max` (total field for `x > x0`, scattered-only for
`x <= x0`, with the seam's own B-value designated scattered by a fixed
convention):

```fortran
! After update_b_field, at the seam (hdtx = 0.5*dt/dx):
by(i0-1, iy) -= hdtx * Ez_inc(x0, y(iy), t_mid)
bz(i0-1, iy) += hdtx * Ey_inc(x0, y(iy), t_mid)

! After update_e_field, at i0 (cnx = hdtx*c**2):
ey(i0, iy) += cnx * Bz_inc(x0, y(iy), t_e)
ez(i0, iy) -= cnx * By_inc(x0, y(iy), t_e)
```

**`direction = x_min` does not just flip every sign at the same
indices** -- this was checked carefully rather than assumed by symmetry,
and the naive "just flip signs" version is wrong. Re-deriving from the
same base equations for `direction = x_min` (total field for `x < x0`)
shows the B-correction stays at the same seam index with flipped signs,
but the **E-correction moves to the seam's west-adjacent E-grid point**
(`i0-1`, not `i0`), because that point's total/scattered designation
flips with direction while the seam B-value's own designation (always
"scattered" by the fixed convention) does not. Implemented correctly in
`apply_antenna_e_correction`; not independently validated by a test this
pass (only `x_max` was exercised).

**No separate Ex correction or incident file is needed**: `ex(i0-1,iy)`
differences `bz(i0-1,iy)` at adjacent `iy`, both already carrying the
B-correction above -- the transverse derivative of the incident field
falls out automatically, which is exactly the mechanism by which a
genuine longitudinal field develops self-consistently rather than being
injected directly.

**Time sampling**: `update_e_field`/`update_b_field` are each called
twice per full step, and the module `time` variable is `t_n` at the
first call but `t_n + dt` (not `t_n + dt/2`) at the second (confirmed
directly against `epoch2d.F90`'s main loop: two separate `time = time +
dt/2` statements straddle `output_routines` in between). Both
B-corrections in a step must use the *same* `t_n + dt/2` -- reading
`time` directly at the second call site would silently sample the wrong
half-step's incident field. Fixed by stashing `antenna_t_mid = time +
hdt` at the top of `update_eb_fields_half` and reusing it (not `time`)
at both B-correction call sites.

## 4. The half-cell file-generation convention (important, verified by test)

Strictly, per the Yee stagger, the E-correction should sample the
incident fields at `x0_actual + dx/2` (the E-grid plane) while the
B-correction samples at `x0_actual` (the B-grid plane) -- **not both at
the same nominal x0**. Stage 2A's Fortran code makes no attempt to
enforce or correct for this: it samples all four files at whatever
position the caller passes in, and treats "x0" as a single label. This
means the physical correctness of the half-cell distinction is entirely
a **file-generation convention**, not a Fortran-enforced one:

> **`ey_inc_file`/`ez_inc_file` must be generated at `x0_actual + dx/2`;
> `by_inc_file`/`bz_inc_file` must be generated at `x0_actual` exactly.**
> Getting this wrong (generating all four at the same nominal plane) does
> not crash or warn -- it reintroduces a resolution-dependent backward-
> wave residual, quantified in §5 below.

This was not obvious going in and was only found by running the
backward-cancellation validation test itself (exactly the point of that
test, per the plan). For a normal-incidence plane wave the fix is a
simple time-shift (`E(x0+dx/2, t) = E(x0, t - dx/(2c))` for a
forward-propagating wave); for the general case (an arbitrary incident
field, e.g. a focused beam), the generating tool (e.g. the Thiele
propagator) should simply be run twice, once outputting at each of the
two staggered planes, since it is a continuous model and not restricted
to a single plane the way a TFSF surface's *consumed* data is.

## 5. Validation results

All runs: `epoch2d`, single rank, vacuum unless noted, normal incidence,
`Ey`/`Bz` polarisation (`Ez`/`By` = 0), antenna at `x0 = -5 um`,
`direction = x_max`. "Residual" = peak `|Ey|` several cells west of the
antenna (scattered-only "quiet zone") relative to the peak transmitted
`|Ey|` east of the antenna.

| Test | dx | lambda | Half-cell file convention | East peak vs analytic E0 | West residual |
|---|---|---|---|---|---|
| 1 | 0.1 um | 0.8 um | No (naive, same-plane) | +5.3% | 19.7% |
| 2 | 0.05 um | 0.8 um | No | +0.6% | 9.4% |
| 3 | 0.1 um | 0.8 um | No (smooth Gaussian envelope) | +5.5% | 19.1% |
| 4 | 0.1 um | 4.0 um | No | +0.17% | 3.6% |
| 5 | 0.05 um | 4.0 um | No | -0.03% | 1.8% |
| 6 | 0.1 um | 0.8 um | **Yes** | -0.35% | **0.86%** |

Tests 1-2 and 4-5 (resolution doubled at fixed wavelength, then again at
a 5x longer wavelength) show the residual scaling essentially linearly
with `dx/lambda` (equivalently `k*dx`) in every case -- the signature of
a genuine, single, well-understood discretisation effect (the half-cell
offset in §4), not a fixed sign/logic bug (which would not shrink
with resolution, or would shrink discontinuously/inconsistently). Test 6
confirms the fix: applying the correct half-cell file-generation
convention collapses the residual from 19.7% to **0.86%** at the same,
otherwise-unremarkable resolution (dx = lambda/8) -- meeting the plan's
<1% target without needing any Fortran change, only the correct file
convention.

**Genuine-backscatter pass-through** (test 9): a second, independent,
ordinary boundary laser (`begin:laser`, `bc_x_max = simple_laser`)
injected a known -x-travelling vacuum wave from `x_max`, timed to start
well after the antenna's own pulse had finished. It arrived at the west
probe at almost exactly the predicted light-travel time (46.3 fs
observed vs 46.7 fs predicted) with amplitude 8.5-8.7e10 V/m against a
predicted 8.67e10 V/m (within a few percent) -- confirming the antenna's
correction does not filter, attenuate or otherwise interfere with a real,
independent backward-going field; it only cancels its own self-generated
backward copy. (An initial attempt to demonstrate this with a genuine
overdense-plasma reflector, test 7/8, hit a numerical instability from
the plasma being drastically under-resolved at 50x critical density --
skin depth ~18 nm against a 100 nm grid spacing -- an orthogonal, well-
understood PIC-resolution problem unrelated to the antenna; replaced
with the cleaner vacuum-wave test above. A properly-resolved plasma-
reflector cross-check is a reasonable follow-up.)

Not attempted this pass (explicit follow-on, per the approved plan):
oblique incidence, a paraxial Gaussian beam, SCiPIC's 2D f/2 field, and
comparison against the Thiele-propagator source-free prediction.

## 6. Known limitations / follow-on work

- **`direction = x_min` is implemented but not independently validated**
  (see §3) -- only `x_max` was exercised by the test suite above.
- **Half-cell file-generation convention (§4) is documentation-only,
  not enforced.** A deck using an antenna generated without respecting
  it will run to completion with a silently elevated (but bounded,
  resolution-dependent) backward residual, not an error.
- **epoch2d only; every rank loads the full incident-field files** (no
  per-rank slab windowing) -- matches epoch2d's existing boundary-laser
  behaviour, but will not scale to epoch3d's large 2D transverse grids.
  The epoch3d port should adopt epoch3d's existing `local_slab_window`-
  style per-rank slab loading and the `balance.F90` re-slab hook already
  built for the boundary laser, applied to the four antenna files.
- **MPI multi-rank plane-split case is implemented (rank-local index
  ownership, no new communication) but not exercised by a multi-rank
  test this pass** -- all validation above ran on a single rank.
- Oblique incidence, Gaussian beam, SCiPIC comparison, and the epoch3d
  port are Stage 2B+ as originally scoped.

## 7. Reproducing the validation

Test artefacts (Python generators + decks) live in the session scratchpad,
not in the repo. Regenerating test 6 (the headline result):

```python
# Ey_inc/Ez_inc sampled at x0_actual + dx/2 (time-shifted by dx/(2c) for
# a normal-incidence wave); By_inc/Bz_inc sampled at x0_actual exactly.
# See section 4 above for why the two channels are NOT sampled at the
# same nominal plane.
```

Deck: `begin:laser_antenna` with `direction = x_max`, `x0 = -5 * micron`,
the four `*_inc_file` binaries, `n_t = 300`, `n_y = 3`,
`y_min`/`y_max = -2/+2 * micron`, `t_start = 0`, `t_end = 15e-15`; domain
`nx = 200`, `x_min/x_max = -10/+10 um`, `bc_x_min = bc_x_max = open`.
Probe `Electric_Field_Ey` a few cells either side of the antenna via the
standard SDF output (`Project_EPOCH/.venv_EPOCH`'s `sdf` reader).
