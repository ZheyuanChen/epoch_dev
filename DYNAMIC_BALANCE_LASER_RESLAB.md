# Continuous re-slabbing of the spatiotemporal custom-laser file
# cache under `use_balance = T`

**Status**: implemented and locally verified in `epoch_dev` (branch
`laser-injection-and-pair-boost`, uncommitted). Touches
`epoch3d/src/housekeeping/balance.F90`,
`epoch3d/src/user_interaction/custom_laser.f90`, and a comment-only
change in `epoch3d/src/laser.f90`. Builds directly on `6ac93129`
("Defer spatiotemporal laser slab load past startup balance settle"),
which fixed the equivalent problem for the one-off `use_pre_balance`
startup rebalance but explicitly left `use_balance` (continuous dynamic
load balancing) on the full-plane fallback path.

---

## 1. The problem this addresses

Motivating case: a laser-solid deck with a pre-plasma target, where
hole-boring, target expansion, and relativistic transparency are all
expected to shift the particle-load distribution substantially over
the course of the run. That calls for `use_balance = T` (continuous
dynamic rebalancing), not just the one-off `use_pre_balance` startup
settle — the domain needs to keep moving to track the target, not just
settle once at t = 0.

`local_slab_window` (`custom_laser.f90`) stores only the transverse
window of a spatiotemporal laser's injected file that a rank's own
patch of the boundary actually needs, rather than the full file on
every rank — this is what makes large injected files (e.g. a high-NA
focusing profile at fine transverse resolution) affordable at scale.
But its fallback condition currently reads:

```fortran
IF (use_balance .OR. (use_pre_balance .AND. .NOT. startup_balance_done) &
    .OR. (move_window .AND. x_is_transverse)) THEN
  ! full plane on every rank
```

`use_balance` unconditionally forces the full-plane fallback, for the
entire run, with no path back to the windowed slab — unlike
`use_pre_balance`, which only forces it until `startup_balance_done`
is set. The result: a deck that needs both real dynamic load balancing
*and* a large injected file cannot have either without compromise.
Per the Viking report that prompted this: 64 ranks × 2 files × 5.76 GB
≈ 737 GB if the file is coarsened enough to survive that, resolution
that was explicitly presented for validation is lost; if it isn't
coarsened, the run does not fit on a node at all. Neither outcome is
acceptable for a target-physics campaign where the injected profile's
resolution is part of what is being validated.

## 2. Why `use_balance` was never given the `use_pre_balance` treatment

`6ac93129`'s fix works by deferring the *first* load of a spatiotemporal
laser's file past deck-parse time, until the one-off startup rebalance
has settled (`startup_balance_done`), then loading it once against the
final, static domain. That works for `use_pre_balance` because the
domain genuinely stops moving after that point — there is a single
"has it settled" instant to wait for.

`use_balance = T` has no such instant: the domain can move again at any
later step (governed by `dlb_threshold`, `dlb_maximum_interval`,
`dlb_force_interval` in `balance_workload`). Deferring the first load
until "the domain is done moving" is not an option — it never is.
`custom_laser_spatial_setup`'s existing comment makes exactly this
argument for why it did not extend the deferred-load treatment to
`use_balance`:

```fortran
! use_balance keeps the domain moving for the entire run, so
! deferring buys nothing there; only the one-off use_pre_balance
! startup settle is worth waiting for.
```

That reasoning is correct as far as it goes — deferring the first load
alone buys nothing — but it treats "wait once" and "never track" as the
only two options. The gap is a third option: keep the window correct
by re-deriving it every time the domain actually changes, rather than
either waiting for it to stop changing (impossible) or storing enough
to never need updating (the current, expensive, fallback).

## 3. Strategy: re-slab on every real redistribution, not just at startup

### 3.1 Where "the domain actually changed" is already decided, collectively

`balance_workload` (`balance.F90`) already computes, identically on
every rank (via `MPI_ALLREDUCE`), whether this step's rebalance check
results in an actual redistribution:

```fortran
IF (balance_improvement > 0.05_num) THEN
  use_redistribute_domain = .TRUE.
  ...
END IF
...
IF (use_redistribute_domain) THEN
  old_comm = comm
  old_coordinates(:) = coordinates(:)
  CALL redistribute_domain
END IF
```

`use_redistribute_domain` is a plain `LOGICAL`, not communicated
separately — it is already guaranteed identical across ranks because
`balance_frac`/`balance_frac_final` (which it is derived from) come
from `MPI_ALLREDUCE(..., MPI_MAX, ...)` and `MPI_ALLREDUCE(...,
MPI_SUM, ...)` calls a few lines above. Any code placed inside that
`IF (use_redistribute_domain) THEN` block after `CALL
redistribute_domain` is, by construction, entered by every rank
together or by none — the same invariant `redistribute_domain` itself,
and the existing `laser%profile`/`laser%phase` slice remap inside
`redistribute_fields`, already depend on. This is the natural, already
paid-for hook: no new consensus mechanism is needed.

`balance.F90` already `USE`s `boundary`, which `USE`s `laser`, which
`USE`s `custom_laser` — none of the three modules declare a module-level
`PRIVATE`, so everything in `custom_laser` is already transitively
visible from `balance.F90`. A new public subroutine in `custom_laser`
is callable from `balance_workload` with no new `USE` statement and no
circular-dependency risk (`balance → boundary → laser → custom_laser`
is one-directional).

### 3.2 The re-slab routine: reuse, not new I/O machinery

Add `reslab_custom_laser_files` to `custom_laser.f90`, called from that
hook. For every laser with `use_custom_profile .AND.
use_spatiotemporal .AND. profile_loaded` (i.e. one that has already
been loaded at least once — covers both the startup load and any prior
reslab):

1. Deallocate the current `file_field_matrix` (and `file_phase_matrix`
   if `phase_loaded`), and clear `profile_loaded`/`phase_loaded`.
2. Call the existing `custom_laser_spatial_setup(laser)` (no
   `allow_defer` argument, so it runs immediately, not deferred).

That routine already resolves the profile/phase filenames and calls
`load_temporal_spatial_profile` → `load_spatiotemporal_file` →
`local_slab_window`, exactly the path used for the first load. No new
file-reading code is needed — `reslab_custom_laser_files` only has to
clear the stale slab and re-enter the existing load path, which will
naturally size the new slab against the now-current `x`/`y`/`z` grid
(see §3.3 on ordering) and re-stream just that window from disk on rank
0, broadcasting to all ranks exactly as before.

Every rank participates in this unconditionally, whether or not *its
own* window actually changed this time — deliberately, since the
collective `MPI_BCAST` calls inside `load_spatiotemporal_file` require
symmetric participation. Skipping the call on ranks whose local window
happens not to have moved (which would need its own inter-rank
consensus check to do safely) is not worth the complexity: a real
redistribution event already means *some* rank's ownership changed, so
in practice this "wasted" reload on an unaffected rank is rare, and the
whole re-slab is already throttled by the same `dlb_threshold` /
`dlb_maximum_interval` / `dlb_force_interval` machinery that limits how
often `use_redistribute_domain` goes true in the first place.

### 3.3 `local_slab_window`'s fallback condition

Change:

```fortran
IF (use_balance .OR. (use_pre_balance .AND. .NOT. startup_balance_done) &
    .OR. (move_window .AND. x_is_transverse)) THEN
```

to:

```fortran
IF (((use_balance .OR. use_pre_balance) .AND. .NOT. startup_balance_done) &
    .OR. (move_window .AND. x_is_transverse)) THEN
```

`use_balance` no longer forces the full-plane fallback outright — it
joins `use_pre_balance` in only forcing it before the *first* settle
(`startup_balance_done`, set once, pre-loop, exactly as today). After
that, the window computed is trusted, and `reslab_custom_laser_files`
keeps it correct across every subsequent real redistribution. The
`move_window` clause (a moving-window x-boundary laser under a
y/z-boundary laser) is untouched — it is a different mechanism
(injection-window motion, not domain-decomposition motion) and is not
addressed by this change; see §5.

### 3.4 `custom_laser_spatial_setup`'s defer condition

Change:

```fortran
IF (defer_ok .AND. laser%use_spatiotemporal .AND. use_pre_balance &
    .AND. .NOT. use_balance) RETURN
```

to:

```fortran
IF (defer_ok .AND. laser%use_spatiotemporal &
    .AND. (use_pre_balance .OR. use_balance)) RETURN
```

At deck-parse time, a spatiotemporal laser now defers its first load
whenever *either* flag applies, not just `use_pre_balance` alone. Both
converge on the same first real load, performed once by
`finalize_custom_laser_domain` after the startup rebalance settles —
which already applies to `use_balance` decks too, since its loop
condition only checks `.NOT. laser%profile_loaded`, not which flag is
set. This removes the current wasted full-plane read at deck-parse
time for `use_balance` decks as a side effect, before any reslab logic
even runs.

## 4. Ordering check (why this cannot race the sampler)

Concern to verify before implementing: could the reslab run too late,
after a step's field solve has already sampled the (now stale, wrong
domain) per-rank slab?

Traced the actual per-step call sequence in `epoch3d.F90`'s main loop:

```
epoch3d.F90:220   CALL update_eb_fields_half   ! samples the current slab
epoch3d.F90:222   CALL run_injectors
epoch3d.F90:224   IF (use_balance) CALL balance_workload(.FALSE.)
                    -> may CALL redistribute_domain, x/y/z updated
epoch3d.F90:225   CALL push_particles
...
epoch3d.F90:283   CALL update_eb_fields_final  ! samples the slab again
```

and inside `redistribute_domain` (`balance.F90`):

```
CALL redistribute_fields(domain)   ! uses the OLD x/y/z + old/new cell
                                    ! bounds together to remap
                                    ! laser%profile/laser%phase
DEALLOCATE(x, y, z)
ALLOCATE(x(...), y(...), z(...))
CALL setup_grid_x / setup_grid_y / setup_grid_z   ! x/y/z now reflect
                                                    ! the NEW domain
```

`local_slab_window` reads `y(0), y(ny)` etc. directly, so it must run
*after* `setup_grid_*` inside `redistribute_domain`, not before. The
planned hook — inside `balance_workload`, immediately after `CALL
redistribute_domain` returns (not from inside `redistribute_domain`
itself) — satisfies this by construction: by the time
`reslab_custom_laser_files` runs, `x`/`y`/`z` already reflect the new
domain, and it completes at line ≈224 of the main loop, strictly before
the next sample at `update_eb_fields_final` (line 283). No reordering
fix is needed; the natural hook point is already correctly placed.

Separately, `redistribute_fields`' existing remap of `laser%profile`/
`laser%phase` (the small, *derived*, per-timestep 2D slice the field
solver actually consumes) is independent data from the raw per-rank
file cache this change re-slabs — it already gets carried across a
redistribution correctly today, by interpolating old-decomposition
values onto the new one. This change only fixes what backs the *next*
time that derived slice gets re-sampled from the raw file, not the
value already in flight for the step during which redistribution
happens.

## 5. Scope and residual limitations

- **Resolution and per-rank memory are fully preserved.** Same file,
  same `n_tr1_points`/`n_tr2_points`, same window-sizing arithmetic —
  the only change is *when* the window gets recomputed (on every real
  redistribution, not never). No coarsening, no extra hardware.
- **Cost is bounded by EPOCH's existing rebalance throttling** —
  `dlb_threshold`, `dlb_maximum_interval`, `dlb_force_interval`, and
  the doubling backoff in `balance_workload` already limit how often
  `use_redistribute_domain` goes true. Each real redistribution now
  additionally re-streams one file window per spatiotemporal custom
  laser; this is the same per-event cost the startup load already
  pays, just potentially repeated a bounded number of times through
  the run.
- **`move_window` is not addressed.** A moving injection window on a
  transverse-boundary laser under a y/z-boundary configuration still
  forces the full-plane fallback unconditionally; that is a distinct
  mechanism (the injection window itself moving, not the domain
  decomposition) and out of scope here.
- **epoch3d only** — the per-rank-slab optimisation this builds on does
  not exist in epoch2d.

## 6. Validation plan

- Clean `gfortran` build of `epoch3d`.
- A ≥2-rank repro deck with `use_balance = T`, a spatiotemporal
  boundary laser, and a graded/moving particle density, engineered to
  force at least one genuine post-startup redistribution
  (`Redistributing. Balance: ... after: ...` printed mid-run, not just
  at "(initial setup)").
- Confirm no `ERROR: custom laser profile sampled outside the stored
  per-rank slab` abort after the mid-run redistribution — this is the
  exact failure `aa605cca` originally fixed for the startup case, and
  is what a wrong or stale window after this change would reproduce.
- Confirm the per-rank slab is narrower than the full plane
  post-redistribution (print or debug the `LBOUND`/`UBOUND` of
  `file_field_matrix`), not just that no abort occurred.
- Existing epoch3d custom-laser regression decks (spatiotemporal,
  static spatial, z-boundary spatiotemporal, file-size abort,
  `use_pre_balance = F`) re-passed to confirm no regression to the
  paths this change does not touch.

## 7. Implementation

Matches the strategy above exactly, in `custom_laser.f90`:

- `custom_laser_spatial_setup`'s defer condition now reads
  `defer_ok .AND. laser%use_spatiotemporal .AND. (use_pre_balance .OR.
  use_balance)`.
- `local_slab_window`'s fallback condition now reads
  `((use_balance .OR. use_pre_balance) .AND. .NOT. startup_balance_done)
  .OR. (move_window .AND. x_is_transverse)` — `use_balance` no longer
  forces the full plane unconditionally.
- New `reslab_custom_laser_files`: for every laser with
  `use_custom_profile .AND. use_spatiotemporal .AND. profile_loaded`,
  deallocates `file_field_matrix` (and `file_phase_matrix` if
  `phase_loaded`), clears the `*_loaded` flags, and re-enters
  `custom_laser_spatial_setup(laser)` (no `allow_defer`, so it runs
  immediately) — reusing the existing load path rather than duplicating
  filename-resolution or I/O code.
- `balance.F90`'s `balance_workload`, right after `CALL
  redistribute_domain` inside the existing `IF (use_redistribute_domain)
  THEN` block, now also calls `CALL reslab_custom_laser_files`. No new
  `USE` statement needed — `balance` already reaches `custom_laser`
  transitively via `USE boundary` → `USE laser` → `USE custom_laser`,
  and none of the three declares a module-level `PRIVATE`.

## 8. Validation performed

- **Ordering check (§4) confirmed by direct code inspection**, not just
  reasoning: traced `redistribute_domain` in `balance.F90` and confirmed
  `redistribute_fields` (which remaps the small derived
  `laser%profile`/`laser%phase`) runs first, against the *old* grid,
  followed by `setup_grid_x`/`_y`/`_z`, which is what actually updates
  `x`/`y`/`z` to the new domain — so the hook placed after
  `redistribute_domain` returns sees the correct, already-updated grid.
  Confirmed against the real `epoch3d.F90` main-loop ordering that
  `balance_workload(.FALSE.)` (line 224) completes before the next
  sample of the raw file cache at `update_eb_fields_final` (line 283) —
  same iteration, hook strictly precedes next use.
- **Clean `gfortran` build** of `epoch3d`, no warnings.
- **2-rank repro** (`nprocx=1, nprocy=2, nprocz=1`, `x_min`
  spatiotemporal laser, `n_y_points=6, n_z_points=4, n_t_points=4`, one
  electron "bunch" species concentrated near the rank boundary with
  `drift_py = 0.99 * me * c` so it physically crosses ranks during the
  run, `dlb_threshold = 0.9` with a short `dlb_maximum_interval`/
  `dlb_force_interval` to force frequent rebalance checks). Confirmed
  via a temporary debug print (removed before finalising) that:
  - The startup pre-load balance narrows the window asymmetrically
    per rank (rank 0: `[1,5]`, rank 1: `[2,6]`, of `[1,6]` full) —
    the pre-existing `6ac93129` behaviour, unaffected.
  - At least ten further `Redistributing.` events fired through the
    run as the bunch drifted (steps 3, 6, 7, 10, 11, 14, 17, 18, 21,
    24, 27, 30, ...).
  - At step 17 the window actually changed
    (rank 0: `[1,5]→[1,6]`, rank 1: `[2,6]→[3,6]`) — direct proof
    `reslab_custom_laser_files` re-derived and reloaded the slab against
    the new domain, not just skipped or silently kept stale data.
  - Windows stayed narrower than the full `[1,6]` plane for the whole
    run except where genuinely appropriate (rank 0 legitimately owns
    the full transverse range from step 17 onward given how the domain
    split at that point) — no permanent full-plane fallback.
  - Zero `ERROR: custom laser profile sampled outside the stored
    per-rank slab` aborts, zero MPI hangs, clean completion to
    `t_end`.
- **Untouched-path spot check**: same repro with `use_pre_balance = F`
  and `dlb_threshold = -1` (`use_balance = F`) still loads eagerly at
  deck-parse (`Loaded Successfully` before the timestep loop, no
  `Redistributing.` lines, since neither flag is set), completes
  cleanly, no errors — confirms the pre-existing eager-load path is
  untouched by the defer-condition change.
- **Full existing epoch3d custom-laser regression suite** re-run against
  the rebuilt binary (`Project_EPOCH/dev_test/laser_profile_injection/
  epoch3d_verification/test{1..7}`, all default `use_pre_balance = T,
  use_balance = F`, i.e. exercising `finalize_custom_laser_domain`'s
  path with the updated-but-equivalent defer condition): test1
  (spatiotemporal y_min boundary), test2 (alias-error aborts, all 4
  sub-cases), test3 (two-laser independence), test4 (interpolation
  accuracy), test5 (static + dynamic edge behaviour), test6 (static
  phase override, all 3 sub-cases), test7 (file-size abort, all 4
  sub-cases) — all passed with no change in behaviour.

## 9. Still outstanding

- ~~Debug instrumentation was temporary and has been removed~~ —
  **done**: `load_spatiotemporal_file` now carries a permanent,
  load-time-only print of the per-rank slab bounds (fires once per
  load/reslab event, never per timestep), and a permanent regression
  test, `epoch3d_verification/test8_dynamic_balance_reslab`
  (Project_EPOCH), parses it to assert the window stays narrower than
  the full plane and actually changes across a reslab event. Checked
  to discriminate: fails against a reverted `local_slab_window`
  condition (every window reports full-plane), passes against the
  current one. This closes the exact gap this section flagged — a
  future silent reversion to full-plane is now a caught memory
  regression, not a blind spot.
- `move_window` (moving injection window on a transverse-boundary laser
  under a y/z-boundary configuration) still forces the full-plane
  fallback unconditionally — out of scope here, see §5.
- Not yet run on Viking / at the scale that originally motivated this
  (large `N_TR`, many ranks, real hole-boring/pre-plasma target deck).
  The 2-rank repro proves the mechanism is correct; it does not
  establish the actual memory saving at production scale, though that
  follows directly from the per-rank window arithmetic being identical
  to the already-proven `6ac93129` case.
- epoch3d only; epoch2d has no per-rank-slab optimisation to extend.
