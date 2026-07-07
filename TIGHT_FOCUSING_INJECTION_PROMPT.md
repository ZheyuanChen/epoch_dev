# epoch_dev session prompt — vectorial injection + binary file format

## Status (updated after the 2026-07-01/02 epoch_dev session)

### Released

| Tag | Commit | Content | Compatibility |
|-----|--------|---------|---------------|
| `v1.1.0` | `945fd060` | Per-laser profile/phase storage fix | Backward-compatible |
| `v2.0.0` | `9b14ed4b` | Binary spatiotemporal profile/phase format | **Breaking** |
| `v2.1.0` | `537b0445` | epoch3d port of the binary injection pipeline | **Breaking** (epoch3d spatial text decks) |

Both tags are GitHub Releases with full release notes at
`github.com/ZheyuanChen/epoch_dev/releases`. `my-epoch-mods` has been
merged into `main`; both branches point to `9b14ed4b`.

### epoch2d — complete

**Post-v2.1.0 addition (commit `cb57b7c1`)**: the static spatial path
(`use_spatiotemporal_profile = F`) was converted from the legacy 1D text
format to the same raw binary convention (**breaking** for old static
text decks — belongs in the next release's notes), and gained static
phase-from-file support, matching epoch3d. Both epoch2d paths are now
binary. Verified: spatiotemporal regression, static field-shape check
(encoded off-centre Gaussian reproduced within one cell), y_min-boundary
load, file-size abort.

Issue 1 is resolved for epoch2d: the `pol_angle`/custom-file interaction is
confirmed correct, and the second-independent-channel question (open
question 3, below) was decided **yes** and implemented — see "Resolution"
under Issue 1. Issue 2 (binary file format) is **also resolved and
implemented** for epoch2d — see the updated Issue 2 section below; the
final spec diverged from the original header-based recommendation after
reading EPOCH's own documented binary-file convention (shape/bounds belong
in the deck, not in the file).

### epoch3d — ported (2026-07-02 session)

The full binary injection pipeline is now implemented in epoch3d
(committed as `537b0445`, released as `v2.1.0`, both branches pushed):

| Feature | epoch3d status |
|---------|---------------|
| 3D spatiotemporal profile E(tr1, tr2, t) | **Implemented** — raw binary, trilinear O(1) sampling |
| 3D spatiotemporal phase from file | **Implemented** — same format/grid as amplitude |
| Static spatial profile (no time) | **Converted to raw binary** — the old text format with embedded coordinates is no longer read (**breaking**) |
| Static spatial phase from file | **Implemented** (new — never existed in the text era) |
| Per-laser storage | Built per-laser from the start (no module-level singletons) |

See "epoch3d port — implementation report" at the end for the deck
elements, file format, design decisions and validation results.

### AELP tutorial — broken

The binary format change **breaks** the existing text-format writer at
`AELP/tutorial/temporal_spatial_gaussian_beam/numerical_input/
generate_spatial_temporal_profile.py` — needs updating to emit raw binary
before that tutorial will run against current epoch2d again.

**Independently verified (separate session, not the build agent's own
check):** built old (committed, pre-fix) and new (working-tree, fixed)
binaries locally via `git stash`/`stash pop`, ran an A/B comparison.
Regression (single custom-file laser, amp+phase): old vs new agree to
~1.8e-15 relative — clean no-op. New targeted test (two `x_min` lasers,
`pol=0`/`pol=90`, same-sized but differently-shaped profile files, chosen
specifically to avoid tripping the grid-mismatch abort and instead expose
the silent-aliasing path): old binary's `Ez` output is a bit-for-bit reuse
of laser A's `Ey` shape (same peak location, same FWHM) instead of laser
B's; fixed binary gives independently correct `Ey`/`Ez` matching each
laser's own file (width ratio and amplitude ratio both match the inputs to
within a few %). Bug confirmed real, fix confirmed effective. Full
methodology and numbers in `Project_EPOCH/daily_log/2026/06/30.md` under
"Validated the epoch_dev per-laser-storage fix".

## Context

This is a design/decision session, not yet a full implementation pass. It
follows up on work done in `Project_EPOCH` (local, sister repo) building
towards the supervisor's next-goal: study an **f/1 tightly-focused laser**
(NA = 0.5) colliding with an electron bunch, and its effect on the
**nonlinear Breit-Wheeler positron spectrum**. The existing custom laser
injector (amplitude + phase from a `.dat` file, validated for 2D at moderate
NA) isn't enough on its own at f/1: LASY (the Python tool generating the
`.dat` files) is fundamentally **scalar** — one polarisation envelope, no
`Ez`, no `B`-field — and its `GaussianProfile` seeds the field directly at
focus using the standard *paraxial* formula, which is exactly the
approximation that breaks down at this NA.

Two things need deciding before the injection pipeline can be extended for
this regime: (1) whether/how a genuinely vectorial field source can feed
into EPOCH's existing boundary condition, and (2) what file format the
`.dat`-replacement should use, given the mod is intended to eventually merge
into public EPOCH.

## Issue 1 — vectorial injection: what EPOCH's boundary condition actually does

**Already established by reading the source (epoch2d, `my-epoch-mods`) —
please verify this independently before relying on it:**

`outflow_bcs_x_min` in `src/laser.f90` (lines ~457-560) builds two source
terms per grid line:

```fortran
t_env = laser_time_profile(current) * current%amp
DO i = 0,ny
  base = t_env * current%profile(i) &
    * SIN(current%current_integral_phase + current%phase(i))
  source1(i) = source1(i) + base * COS(current%pol_angle)
  source2(i) = source2(i) + base * SIN(current%pol_angle)
END DO
```

`source1` drives `bz` (→ the `Ey`-like transverse channel), `source2` drives
`by` (→ the `Ez`-like transverse channel, i.e. the second polarisation
component, **not** a longitudinal-along-propagation field — `y` and `z` are
both transverse to the `x`-propagating beam in EPOCH2D's x-min boundary).
The split is controlled entirely by the single deck-level `pol_angle`
constant.

**Key implication:** EPOCH's boundary framework already supports two
*transverse* polarisation channels (`source1`/`source2`), but both are
currently driven from **one** scalar `profile`/`phase` pair via
`cos(pol_angle)`/`sin(pol_angle)` — not two independently-specified
channels. `custom_laser.f90` does not reference `pol_angle` at all (grepped,
zero hits), so the custom-file path presumably inherits whatever `pol_angle`
the deck sets, same as the analytic path — **worth confirming this isn't
silently broken or ignored when `use_custom_profile = T`.**

**Also important:** there is no mechanism here to inject a genuine
**longitudinal** field (along the propagation/boundary-normal direction,
i.e. `Ex` in this geometry) at the boundary at all — and that's expected,
not a gap to fix. Boundary conditions specify the tangential (transverse)
field; the longitudinal component is not an independent degree of freedom
to inject — it develops self-consistently from `∇·E = 0` as the wave
propagates into the bulk via EPOCH's own Maxwell solver. This reframes what
a vector-field-accurate generator (e.g. a properly-seeded Thiele propagator
— see `Project_EPOCH/src/epoch_tools/field_propagator/thiele_injector.py`,
recently audited and found mathematically sound: round-trip self-consistent
to ~1e-7 relative error, energy-conserving to ~1e-13) is actually useful
for at the injection stage: **getting the *transverse* `Ey`/`Ez` profile
itself accurate at high NA**, not feeding EPOCH a pre-computed longitudinal
field or `B`. The longitudinal structure EPOCH itself produces in the bulk
becomes a natural cross-check: if the transverse injection is accurate,
the `Ez`/`B` that develop a few cells into the box should match what Thiele
predicts at that plane.

### Open questions for this session — resolved 1-3, 4 still open

1. **Resolved — confirmed.** Re-read `outflow_bcs_x_min`/`outflow_bcs_x_max`
   in epoch2d and all six boundary subroutines in epoch3d
   (`x_min/max`, `y_min/max`, `z_min/max`): the `source1`/`source2`
   `cos`/`sin(pol_angle)` split pattern is identical everywhere.
2. **Resolved — `pol_angle` is correctly respected on the custom-file path,
   not silently broken.** `custom_laser.f90` has zero references to
   `pol_angle` because it doesn't need any: it only populates
   `laser%profile`/`laser%phase`. The `cos`/`sin(pol_angle)` split happens
   downstream in `laser.f90`'s boundary routines, applied uniformly
   regardless of whether `profile`/`phase` came from the analytic deck
   functions or the file loader. `pol_angle` itself is set unconditionally
   from the deck (`deck_laser_block.f90`), independent of
   `use_custom_profile`.
3. **Resolved — decided yes, implemented and verified this session.**
   While checking whether two custom-file laser blocks (`pol_angle` = 0 and
   π/2) could already give independent channels "for free" via the existing
   `source1`/`source2` summation over all lasers on a boundary, found that
   they could **not**: `custom_laser.f90`'s profile/phase storage
   (`file_field_matrix`, `file_phase_matrix`, the coordinate arrays, and the
   `profile_loaded`/`phase_loaded` guards) was module-level `SAVE` — a
   global singleton shared across every laser block, not per-laser. A
   second laser block requesting a different file would hit the
   already-set `loaded` guard and silently keep reusing the first laser's
   data (or abort on a grid-size mismatch), with no warning. Given this was
   a real silent-bug risk regardless of whether f/1 LP needs two channels,
   fixed it: moved the grid/data storage onto the `laser_block` type itself
   (`shared_data.F90`) and refactored `ensure_shared_grid`,
   `load_temporal_spatial_profile`, `load_phase_profile`,
   `custom_laser_profile`, `custom_laser_phase` (`custom_laser.f90`) to
   operate on the passed-in `laser` instance. `deallocate_laser`
   (`laser.f90`) frees the new per-laser pointer fields. Also renamed
   `file_y_coords`/`n_y_points` to `file_transverse_coords`/
   `n_transverse_points`, since the same field serves `y_min`/`y_max`
   boundary lasers too, where the transverse axis is physically `x`, not
   `y` — the old name was misleading once formalised into the type.
   Rebuilt epoch2d clean and verified end-to-end: two `x_min` laser blocks
   (`pol_angle` 0 and π/2) with deliberately different-sized synthetic
   profile/phase files (6×4 and 9×7) both loaded independently with the
   correct distinct grid sizes, no aliasing, no abort, ran to completion
   under 2 MPI ranks. Two custom-file lasers with distinct `pol_angle` and
   distinct files now genuinely give independent amplitude+phase per
   transverse polarisation channel via existing deck syntax — no new file
   format or deck element needed for the elliptical/orthogonal-channel
   case.
   ~~**Not yet done: port this same fix to epoch3d**~~ **Done in the
   2026-07-02 session** — the epoch3d port was built per-laser from the
   start (see the implementation report at the end).
4. **Still open.** If a future workflow uses Thiele (or another vector
   method) to generate the transverse field, what's the cleanest place to
   validate it inside EPOCH — e.g. dumping `ex`/`ey`/`ez`/`bz`/`by` a few
   cells downstream of the boundary and comparing against the
   externally-computed prediction at that plane?

## Issue 2 — file format for the profile/phase data — RESOLVED (epoch2d)

**Independently verified (separate session):** read the current source
directly (not the handoff doc's description) and ran three targeted tests.
(1) Array ordering — the highest-risk item: a deliberately non-separable,
asymmetric profile (a peak moving diagonally across y as a function of t,
`n_t != n_y` so a swap would be structurally loud) written as a numpy
`(n_t, n_y)` array via `.tofile()` with no transpose, exactly as specced.
Measured at the boundary cell, the injected field's peak trajectory matched
the encoded diagonal to within grid resolution (0.01-0.13 um) across 5
snapshots — confirms the no-transpose claim is actually correct, not just
internally consistent reasoning. (2) File-size-mismatch abort — retriggered
independently with a wrong `n_t`, got a clean `MPI_ABORT` with the exact
expected/actual byte counts. (3) Two-channel independence — re-ran the
Issue-1 polarisation-independence test under the new binary format (not the
old text format): both lasers loaded independently, `Ey`/`Ez` shapes and
amplitude ratios matched their respective generator inputs, same numbers as
the original text-format test. All three passed. Full methodology and
numbers in `Project_EPOCH/daily_log/2026/06/30.md` under "Validated the new
raw-binary profile/phase format".

**Final decision: binary, not OpenPMD/HDF5** — same reasoning as the
original recommendation below, now implemented.

Why binary-not-OpenPMD, given the mod is meant to merge upstream:
- EPOCH's own native I/O (SDF) is a deliberately lightweight, dependency-free
  binary format — no HDF5 link anywhere in the build. Requiring HDF5 (even
  optionally) for a laser-injection feature is a much harder sell for a
  public-EPOCH merge than something that fits the existing minimal-dependency
  philosophy and portability story across HPC clusters with inconsistent
  HDF5 module availability.
- This isn't introducing a new I/O pattern to EPOCH — `custom_laser.f90`
  already loads external profile data from a file; this is extending an
  established convention (text → binary), not inventing a new one. Smaller,
  more reviewable diff for an eventual upstream PR.
- OpenPMD's real value (self-describing metadata, cross-code/tool
  interoperability) matters when data crosses between different tools or
  groups. Here the format is a narrow, private interface between a Python
  preprocessing step (under our control) and this Fortran reader (also under
  our control) — the interoperability case is weak.
- Keeps epoch_dev decoupled from LASY specifically at the file-format level,
  per the constraint that the public-mergeable epoch-mod should not be tied
  to any one Python generation tool.

**The original "small in-file header" layout (below) was superseded** after
re-reading EPOCH's own documented "Binary files" convention (covers
`particles_from_file` and the other file-reading deck blocks): EPOCH's
binary files carry **no shape metadata at all** — "do not contain any
information about the shape of the arrays... must be supplied using the
input deck" — and the documented example confirms Fortran `access='stream'`
column-major output (first index fastest-varying), one array per file.
Matching that exactly, rather than inventing an in-file header, was judged
the better fit for an eventual upstream merge.

**Implemented spec** (`custom_laser.f90`, `deck_laser_block.f90`,
`shared_data.F90`):
- `profile_data_file` / `phase_data_file` are raw binary, `access='stream'`,
  `form='unformatted'`: just `n_t_points × n_transverse_points` values of
  EPOCH's `REAL(num)` (always 8 bytes — `num = KIND(1.d0)`, confirmed in
  `constants.F90`), no header, no coordinate arrays.
- **Array ordering, pinned down explicitly** (the exact bug class flagged
  below): Fortran array declared `(n_transverse_points, n_t_points)`
  (transverse fastest-varying, matching the documented convention) written
  via a single un-looped `WRITE(unit) array`. On the Python side, this is
  written/read as a numpy array of shape `(n_t, n_transverse)` in default
  C (row-major) order — `E.tofile(...)` / `np.fromfile(...).reshape(n_t,
  n_transverse)` round-trips with **no transpose needed**, because
  numpy's row axis (t) lines up with Fortran's slower-varying axis. Verified
  by actually generating files this way and loading them in epoch2d (see
  below) — not just reasoned about.
- **Shape and bounds are deck-declared, not file-embedded**, per the
  official convention: new laser-block elements `n_t_points`/`n_t`,
  `n_transverse_points`/`n_y`, `profile_transverse_min`/`y_min`,
  `profile_transverse_max`/`y_max` (each pair are aliases of each other).
  The temporal grid extent reuses the laser's existing `t_start`/`t_end`
  rather than adding a separate pair of elements, since the laser is only
  ever active in that window anyway.
- **Uniform grid only** — matches the existing requirement already baked
  into `custom_laser_profile`/`custom_laser_phase`'s O(1) bilinear lookup
  (the non-uniform search path was already dead code). This is what makes
  dropping the coordinate arrays safe: the grid is fully reconstructible
  from `n_*_points` + bounds + `t_start`/`t_end`.
- **Defensive file-size check**: aborts with a clear error if
  `file size ≠ n_transverse_points * n_t_points * 8 bytes`, since this is
  the only available shape sanity check once there's no embedded header.
  Verified this actually fires (deliberately wrong `n_t_points` in a test
  deck → clean abort with the expected/actual byte counts printed, not a
  garbled read).
- `ensure_shared_grid` (the cross-check from the Issue 1 fix, validating
  that a laser's amplitude and phase files agreed on grid size) is now
  **unnecessary and removed** — both files share one deck-declared grid by
  construction, not by runtime validation.

Verified end-to-end: two `x_min` laser blocks (`pol_angle` 0 and π/2) with
deliberately different-shaped raw binary files (6×4 and 9×7, generated via
numpy) both loaded correctly; a deliberately-wrong `n_t_points` triggered
the file-size abort as expected.

**Known breaking consequence**: the old text `.dat` format is no longer
read for the spatiotemporal path. `AELP/tutorial/temporal_spatial_gaussian_
beam/numerical_input/generate_spatial_temporal_profile.py` (separate repo)
still writes the old format and needs updating before that tutorial works
against current epoch2d again — not done this session.

~~**Not done this session**: porting any of this (binary format, per-laser
storage) to epoch3d.~~ **Done in the 2026-07-02 session** — see the
implementation report at the end.

<details>
<summary>Original recommendation text (superseded, kept for context)</summary>

**Two implementation details to nail down explicitly, not leave implicit**
(both are exactly the class of bug this project has repeatedly hit this
month — silent convention mismatches):

1. **Use Fortran `access='stream'`**, not default sequential unformatted
   I/O. Stream files are a flat byte sequence directly readable by
   `np.fromfile`/`numpy.frombuffer` with no surprises. Default unformatted
   I/O inserts compiler-dependent record-length markers that silently
   corrupt cross-language reads (numpy has no idea those markers are there).
2. **Pin down array ordering (row-major vs column-major) explicitly in the
   format spec, in writing.** Fortran is column-major, numpy is row-major by
   default. This exact bug class (implicit transposition / axis mismatch)
   has bitten this project multiple times already (see the carrier-offset
   and phase-unwrap bugs documented in `Project_EPOCH`'s daily logs from
   late June 2026) — don't let the binary format introduce a new instance of
   it by leaving the ordering convention unstated.

A reasonable minimal layout (open to revision): small fixed-size header
(format version tag, `n_t`, `n_y`, and `n_z` if extending to a true 3D
spatial+temporal profile), then the coordinate arrays (float64), then the
data block(s), all written/read via stream access with the ordering
convention documented at the top of both the Fortran reader and whatever
Python writer eventually replaces the current `.dat` generators.

</details>

### What to produce this session — status

- ~~A confirmed (or corrected) account of how `pol_angle` interacts with the
  custom-file injection path.~~ **Done** — see open question 2 above.
- ~~A decision on whether the second-channel question (3, above) is in scope
  now or deferred.~~ **Done — in scope, implemented and verified for
  epoch2d.** See open question 3 above.
- ~~A concrete binary format spec (header layout, byte order, array ordering,
  versioning) that both a Fortran reader and a future Python writer can be
  built against.~~ **Done — implemented and verified for epoch2d.** See the
  "RESOLVED" Issue 2 section above. Note the final spec has **no header**
  (shape/bounds are deck-declared instead), which is a deliberate departure
  from the original recommendation, not an oversight.

## epoch3d port — implementation report (2026-07-02 session)

**Status: implemented, built clean (gfortran), functionally and
physically validated. Committed as `537b0445` and released as `v2.1.0`
(GitHub release with notes; `my-epoch-mods` and `main` both point at
it).** Steps 1-3
of the previous roadmap are all done (including the "optional" step 3,
plus a static spatial phase-from-file path that never existed in the text
era). Step 4 (AELP tutorial writer) remains outstanding.

### What was implemented

1. **3D spatiotemporal amplitude + phase** — exactly the roadmap spec:
   Fortran array `REAL(num), DIMENSION(n_tr1, n_tr2, n_t)`, tr1
   fastest-varying, time slowest, `access='stream'`,
   `form='unformatted'`, no header, single `READ(unit)`. Python side
   writes numpy C-order `(n_t, n_tr2, n_tr1)` via `.tofile()` with **no
   transpose** — the direct extension of the verified 2D rule. Sampled
   per boundary cell per step by an O(1) trilinear interpolator
   (`sample_file_matrix`) on the uniform grid reconstructed from the
   deck-declared bounds/counts and `t_start`/`t_end`.
2. **Static spatial path converted to raw binary** — `(n_tr1, n_tr2)`
   values, same convention minus the time axis, bilinearly interpolated
   onto `laser%profile` once at setup. **Breaking**: the old text
   `spatial_profile.dat` (integer-count header + coordinate rows) is no
   longer read; old decks abort loudly with the missing-grid-declaration
   error.
3. **Static spatial phase-from-file** (new) — same binary plane format,
   interpolated onto `laser%phase` once at setup; `laser_update_phase`
   early-returns on this path so a deck `phase = ...` expression can
   never overwrite it.
4. **Per-laser storage from the start** — grid declaration and both data
   matrices live on `laser_block` (`shared_data.F90`), freed in
   `deallocate_laser`; the epoch2d singleton bug was never reproduced.
5. **Deadlock-safe loading** — all file loads happen at deck-parse time
   (`attach_laser` → `custom_laser_spatial_setup`), so the rank-0 read +
   `MPI_BCAST` involves all ranks. Reading is centralised in one
   `load_binary_file` helper (missing-file abort; file-size abort
   printing expected/actual byte counts — the only shape sanity check
   available with headerless files).

### Deck elements (laser block)

Canonical boundary-agnostic names, a short form, and axis-named aliases:

| Canonical | Short | Axis-named alias |
|-----------|-------|-----------------|
| `n_transverse1_points` | `n_tr1` | `n_x_points`/`n_x`, `n_y_points`/`n_y`, `n_z_points`/`n_z` |
| `n_transverse2_points` | `n_tr2` | (same six — resolved per boundary) |
| `n_t_points` | `n_t` | — |
| `profile_transverse1_min`/`_max` | `tr1_min`/`tr1_max` | `x_min`/`x_max`, `y_min`/`y_max`, `z_min`/`z_max` |
| `profile_transverse2_min`/`_max` | `tr2_min`/`tr2_max` | (same six — resolved per boundary) |

Plus `use_phase_from_file` / `phase_data_file`, identical to epoch2d.
Temporal extent reuses `t_start`/`t_end` as in 2D. Defaults filenames:
`temporal_spatial_profile.dat`, `spatial_profile.dat`,
`phase_profile.dat`.

**Aliasing decision** (the open naming question): the axis-named aliases
are **boundary-aware**, not blind synonyms — `n_y_points` maps to tr1 on
an `x_min`/`x_max` laser but tr2 on a `z_min`/`z_max` laser, per the
roadmap's axis table. Naming the propagation axis (e.g. `n_x_points` on
an `x_min` laser) is a deck **error** with a message pointing at the
transverse axes and the boundary-agnostic names. This makes decks read
naturally per boundary while keeping `n_tr1`/`n_tr2` available for
boundary-independent tooling.

Validation at load time: each point count must be >= 2 and each bounds
pair must satisfy max > min (and `t_start < t_end` on the spatiotemporal
path), else a clear `MPI_ABORT` listing the required elements.

**Deliberate 2D/3D difference kept**: epoch3d's `use_spatiotemporal`
default remains `F` (epoch2d defaults to `T`), preserving epoch3d's
historical spatial-only default for `use_custom_profile = T` decks.

### Files changed (epoch3d only)

- `src/shared_data.F90` — `laser_block` gains `use_phase_from_file`,
  `phase_data_file`, `profile_loaded`/`phase_loaded`, `n_t_points`,
  `n_tr1_points`/`n_tr2_points`, `profile_tr1_min`/`max`,
  `profile_tr2_min`/`max`, 3D `file_field_matrix`/`file_phase_matrix`.
- `src/user_interaction/custom_laser.f90` — rewritten: dispatcher
  (`custom_laser_spatial_setup`), `load_binary_file`,
  `check_file_grid_declared`, spatiotemporal loaders + trilinear
  sampler, binary spatial loader (`load_spatial_fields` +
  `interp_plane_to_boundary`); the `custom_laser_profile_3d` stub is
  gone.
- `src/deck/deck_laser_block.f90` — new elements above;
  `transverse_axis_index` / `set_transverse_count` /
  `set_transverse_bound` helpers implementing the boundary-aware
  aliases.
- `src/laser.f90` — phase-from-file branch in `laser_update_phase`
  (+ static-path early return); stub call replaced with
  `custom_laser_profile` in `laser_update_profile`; six boundary
  routines call `laser_update_phase` when `use_phase_from_file`;
  per-laser matrices freed in `deallocate_laser`; unconditional
  setup-time load call in `attach_laser`.

### Validation (all with the gfortran build, 2 MPI ranks)

1. **Spatiotemporal amp + phase load** — 6(y)×5(z)×4(t) files on an
   `x_min` laser declared entirely via axis aliases (`n_y_points`,
   `n_z_points`, `y_min`...`z_max`): both loaded with the correct
   distinct sizes, run to completion.
2. **Static spatial amp + phase load** — 7(x)×9(y) files on a `z_min`
   laser declared via mixed canonical/short names (`n_tr1`,
   `n_transverse2_points`, `tr1_min`, `profile_transverse2_min`, ...):
   both loaded, run to completion.
3. **File-size abort** — deliberately wrong `n_t` → clean `MPI_ABORT`:
   "amp3d.dat is 960 bytes; expected 1200".
4. **Propagation-axis alias error** — `n_x_points` on an `x_min` laser →
   clean deck error naming the element and line number.
5. **Array-ordering test** (the transpose bug class, measured physically,
   not just reasoned): 48(y)×40(z)×16(t) profile — an asymmetric
   Gaussian (w_y = 2.5 µm, w_z = 5 µm), centre fixed at z = -2 µm but
   drifting in y from -3 to +3 µm across the time window — written from
   numpy with no transpose. Measured `|Ey|` on the plane one cell inside
   `x_min` across 6 snapshots: peak y tracked the encoded drift to
   within one grid cell (0.33 µm) at every snapshot, peak z stayed at
   -2 µm, and the width asymmetry came out in the correct orientation
   (~2.6 µm vs ~4.4 µm; the wider axis reads slightly low because the
   second-moment estimate is domain-truncated). With n_y ≠ n_z and a
   y-only drift, a tr1/tr2 swap or time-axis mixup would have been
   structurally loud. Time, tr1 and tr2 axes all confirmed.

Test artefacts live in the session scratchpad (synthetic numpy
generators + minimal decks), not in the repo.

### Known non-bug behavioural divergence — out-of-file-grid edge handling
(found during independent post-port verification, `Project_EPOCH` session,
2 July 2026)

If a laser's declared file grid (`profile_tr1_min/max`, `profile_tr2_min/max`
in the deck) is **smaller** than the simulation's transverse box, the two
injection paths in `custom_laser.f90` fall back differently for points
outside the file grid — both intentional per their own local logic, but
divergent from each other:
- **Static spatial path** (`interp2d`, lines ~359-395): clamps the
  interpolation fraction to `[0,1]`, so it keeps returning the file's edge
  value indefinitely outward — a flat plateau with no decay.
- **Spatiotemporal path** (`sample_file_matrix`, lines ~405-459): returns a
  hard `0.0` immediately outside `[tr1_min,tr1_max] x [tr2_min,tr2_max]` —
  an abrupt cutoff, no interpolation at all.

Measured with a Gaussian file grid declared only to ±2 µm inside a ±6 µm
box: static path stays at ~87-91% of its edge value all the way to 5 µm;
spatiotemporal path drops three orders of magnitude within ~1 µm of the
edge. Two qualitatively different silent artefacts from the same mistake,
not a subtle numerical difference.

**Relevance to the f/1 tight-focusing goal**: whichever path is used for
the LASY-generated file, its declared transverse extent must comfortably
cover the simulation's transverse box with the Gaussian tails already
negligible at the file edge — not just match the nominal spot size.
Otherwise: spatiotemporal path → artificial hard-edged aperture (could be
mistaken for real diffraction structure); static path → a spurious
low-level illumination skirt across the rest of the boundary. No error or
warning fires either way. Worth an explicit check (deck assertion or a
Python-side pre-flight check comparing file extent to box extent) before
trusting an f/1 run's near-field pattern at the domain edges.

### Memory scaling — per-rank slab storage (added post-v2.1.0)

Implemented in a follow-up commit on `my-epoch-mods`: each rank stores
only the window of the spatiotemporal file grid covering its own
transverse boundary patch (+1 file cell margin); ranks not owning the
laser's boundary face store an empty slab; loading streams the file one
time-slice at a time on rank 0 (so no rank ever holds the full 3D
array); the slab keeps its global index bounds via pointer allocation so
the sampler is unchanged apart from two guards (empty-slab early return;
loud `MPI_ABORT` on out-of-slab access). Full-plane fallback on boundary
ranks when the local patch can move: dynamic load balancing, or moving
window with a `y`/`z`-boundary laser (whose tr1 axis is x). Verified
bit-for-bit against the pre-slab v2.1.0 fields at 2 ranks, and to ~2e-15
at 4 ranks with `nprocx = 2` (two ranks off the boundary, exercising the
empty-slab path); z-boundary laser and abort/spatial tests re-passed.
**Future option if pulses get long: temporal chunking** (time is the
slowest file axis, so a rolling window is one contiguous positioned
stream read; timestep synchronisation makes the collective refresh
deadlock-safe). Not implemented.

**Follow-up (`aa605cca`, see below): the slab window is only used when
`use_pre_balance = F` and `use_balance = F` are both set explicitly.**
Both default to values that trigger the full-plane fallback
(`use_pre_balance` defaults `.TRUE.`), so out of the box this commit's
memory saving does not apply to a typical deck — only to one that has
deliberately frozen its domain decomposition. Worth revisiting: a
structural fix (deferring the slab load until after `pre_load_balance`
completes) would recover the saving under default settings, at the cost
of restructuring the startup sequence; not attempted here.

### MPI slab-vs-load-balancing bug — root-caused and fixed (`aa605cca`)

Follow-up to the "Known MPI limitation" reported 2 July 2026
(`Project_EPOCH` session, against `9d50a6a6`): the per-rank-slab
optimisation aborted with `ERROR: custom laser profile sampled outside
the stored per-rank slab` under certain multi-rank decompositions.

**Root cause, confirmed by reading the epoch3d startup sequence
(`epoch3d.F90`) and reproducing directly:** `local_slab_window` computes
each rank's transverse window at laser setup (`attach_laser`, during the
second deck-parse pass). But `use_pre_balance` **defaults to `.TRUE.`**
(`deck_control_block.F90`), and the one-off startup load balance it
triggers (`pre_load_balance`, `setup.F90`) runs *after* that second
`read_deck` call and redistributes cell ownership between ranks based on
particle load — invalidating the already-sized-and-filled slab whenever
the initial particle distribution is uneven. Reproduced directly: a
2-rank `x_min` spatiotemporal laser with a y-graded electron density
(`number_density = 1e24 * (1.5 + y/4e-6)`) printed `Redistributing.
Balance: 0.881, after: 0.993 (pre-load balance)` immediately after laser
setup, then aborted on the very next boundary update — matching the
originally reported symptom exactly. (The 320x320x3 / `1 1 4` repro in
the original report likely triggered the same path via
`use_optimal_layout`, also default `.TRUE.`, reshaping the processor
grid itself post-setup — same root cause, different trigger.) The static
spatial path was correctly identified as unaffected: it has no slab
logic, full plane on every rank always.

**Fix**: extended `local_slab_window`'s existing "can the local patch
change after this?" fallback (already covering `use_balance` and a
moving window on y/z-boundary lasers) to also cover `use_pre_balance`.
Since that defaults to `.TRUE.`, **the full per-rank-plane fallback is
now the default for typical decks** — correctness restored, no aborts —
and the slab-window memory saving from `9d50a6a6` now applies only once
a deck explicitly sets `use_pre_balance = F` (and `use_balance = F` for
the dynamic case). This is a real reduction in the memory feature's
day-to-day benefit versus what was originally claimed; see the follow-up
note under "Memory scaling" above.

Verified: the particle-imbalance repro now completes with no abort; the
2-rank ordering test and all four epoch3d custom-laser regression decks
(spatiotemporal, static spatial, z-boundary spatiotemporal, file-size
abort) re-pass; re-running with `use_pre_balance = F` explicitly still
exercises the real windowed slab path and reproduces the reference
fields bit-for-bit, confirming the slab logic itself was correct once
the domain is actually static — the bug was purely about *when* the
window was computed relative to EPOCH's own startup rebalancing, not the
window arithmetic.

**Action for the postponed two-channel-polarisation test**
(`Viking_results/3D/test2_two_channel_polarisation/`, per the handover
doc): safe to re-attempt now at the target rank count. If it uses
default settings it will silently take the full-plane fallback (correct
but no memory saving); set `use_pre_balance = F` explicitly if the
memory saving specifically needs to be exercised.

### Still outstanding

- **Step 4 below (AELP tutorial writer)** — now needed for both 2D and
  3D.
- ~~**Commit + release tag**~~ **Done** — committed as `537b0445`,
  released as `v2.1.0` with the breaking change (epoch3d spatial
  text→binary) called out in the release notes.
- **3D field-level benchmark** — quantitative validation (focusing
  Gaussian: waist position/size/energy vs theory; two-channel test in
  3D) to be run separately.
- Open question 4 (Thiele validation plane) — unchanged.

### Step 4 — update AELP tutorial

Update `AELP/tutorial/temporal_spatial_gaussian_beam/numerical_input/
generate_spatial_temporal_profile.py` to emit the raw binary format used by
epoch2d (and eventually epoch3d). Currently writes old text `.dat` which
epoch2d can no longer read.

### Open question 4 (still unaddressed)

If a future workflow uses Thiele (or another vector method) to generate the
transverse field, what's the cleanest place to validate it inside EPOCH —
e.g. dumping `ex`/`ey`/`ez`/`bz`/`by` a few cells downstream of the
boundary and comparing against the externally-computed prediction at that
plane?

## High-NA injection fidelity — options (2026-07-02 discussion)

Parked for a later session; recorded so the reasoning isn't lost. The
concern: EPOCH's characteristic-based laser boundary source term is
**first-order in incidence angle** — it injects as if each field component
propagates normal to the boundary. An f/1 beam carries plane-wave
components out to ~30°, which the BC will inject with angle-dependent
amplitude/phase errors *no matter how accurate the profile/phase files
are*. Moving the boundary further from focus does **not** help: the
angular spectrum is invariant under propagation.

**Agreed plan, in order of preference:**

1. **Measure before engineering.** Run the downstream-plane comparison
   (open question 4) with the Thiele reference at f/1: dump
   `ex`/`ey`/`ez`/`by`/`bz` a few cells inside the boundary, compare with
   the Thiele prediction at that plane. A Gaussian f/1 beam's energy at
   its extreme angles is limited — the error may be a few per cent and
   acceptable for the Breit-Wheeler spectrum study. Don't build anything
   until this number exists.
2. **Current-sheet "laser antenna" injection** — the structural fix if
   the measured error is unacceptable, and the community-proven route
   (Smilei's antenna, WarpX's laser antenna). Prescribe a surface current
   `K(y, z, t) = -2 E_target / Z0` on one plane of cells *inside* the
   domain; the Maxwell solver itself radiates every angular component
   exactly (no one-way approximation), and the backward-propagating twin
   simply exits through the now purely-absorbing `x_min` boundary.
   Crucially this is **additive, not an overhaul**: a new laser-block
   mode (e.g. `inject_plane = <x0>`) that reuses the entire v2.1.0
   file-loading/interpolation machinery and deposits into `jy`/`jz` at
   that plane each step instead of driving the boundary source terms.
   Existing decks and boundary routines untouched. (Note: an interior
   injection plane means interior ranks need profile data — interacts
   with the per-rank slab storage; the antenna loader should compute its
   own window at the injection plane.)
3. **Python-side angular pre-compensation** — cheap experiment, zero
   Fortran changes: measure the BC's angular transfer function
   `T(k_y, k_z)` numerically (inject single-k components, fit amplitude
   and phase downstream), then pre-divide the generated field's angular
   spectrum by `T` before writing the files. Risk: ill-conditioning
   where `T` becomes small at the largest angles; cannot recover
   genuinely missing physics.
4. **Volume initialisation** (set `E` and `B` over the whole grid at
   t = 0 from Thiele): near-zero code but the box must contain the whole
   pulse — long box, expensive in 3D. Useful as a one-off cross-check of
   options 2-3, not for production.

**Rejected:** rewriting the boundary source term to higher order in angle
— a genuine overhaul of `laser.f90` with regression risk across all six
boundaries, when option 2 achieves the same accuracy additively.