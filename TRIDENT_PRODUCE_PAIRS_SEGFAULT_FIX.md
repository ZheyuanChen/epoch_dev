# Fix: TRIDENT_PHOTONS segfault when `produce_pairs = F`

**Status**: fixed in `photons.F90`, all three dimensions (`epoch1d`,
`epoch2d`, `epoch3d`). Originally found and root-caused in
`Project_EPOCH` (2026-06-26 daily log) while running the GPR intensity
scan tutorial deck; carried as an open `epoch_dev` issue in
`Project_EPOCH/TODO.md` §3 until reproduced and fixed here.

## 1. The bug

Compiling with `-DTRIDENT_PHOTONS` and running a deck with
`produce_pairs = F` (the default) segfaults as soon as any electron's
trident optical depth first drops below zero:

```
#3  __partlist_MOD_add_particle_to_partlist
      at src/housekeeping/partlist.F90:403
#4  __photons_MOD_generate_pair_tri
      at src/physics_packages/photons.F90:1109
#5  __photons_MOD_qed_update_optical_depth
      at src/physics_packages/photons.F90:588
```

## 2. Root cause

`check_qed_variables` (photons.F90) returns early when
`.NOT. produce_pairs`, *before* it reaches the code that defaults
`trident_electron_species`/`trident_positron_species` to
`first_electron`/`first_positron` when the deck did not declare
`identify:trident_electron`/`identify:trident_positron` explicitly:

```fortran
! If you're not producing pairs then you don't have to designate special
! electron or positron species so just return
IF (.NOT.produce_pairs) RETURN
...
#ifdef TRIDENT_PHOTONS
    IF (trident_positron_species < 0) THEN
      ...
      trident_positron_species = first_positron
    END IF
#endif
```

So with `produce_pairs = F` and no explicit trident species, both
stay at their unset default of `-1` (`shared_data.F90`).

`qed_update_optical_depth`'s electron/positron branch, however,
decrements `optical_depth_tri` and calls `generate_pair_tri`
**unconditionally** — it was never gated on `produce_pairs`, unlike
the photon branch a few lines below it (`... .AND. produce_pairs`),
which correctly is:

```fortran
#ifdef TRIDENT_PHOTONS
          IF (current%optical_depth_tri < 0.0_num) THEN
            CALL generate_pair_tri(current, trident_electron_species, &
                trident_positron_species)
```

`generate_pair_tri` calls `add_particle_to_partlist` on
`species_list(ielectron)`/`species_list(ipositron)` with `ielectron =
ipositron = -1` — an out-of-bounds array access, hence the segfault.
This is a pure omission bug: `produce_pairs` is EPOCH's single global
pair-production toggle (used identically for the Breit-Wheeler photon
branch), and `check_qed_variables` already treats trident as covered
by it; `qed_update_optical_depth`'s trident block was just never
updated to match.

## 3. Fix

Gate both the optical-depth decrement and the generate-and-reset call
on `produce_pairs`, mirroring the photon branch's existing pattern:

```fortran
#ifdef TRIDENT_PHOTONS
          IF (produce_pairs) THEN
            current%optical_depth_tri = current%optical_depth_tri &
                - delta_optical_depth_tri(eta, gamma_rel)
          END IF
#endif
          ...
#ifdef TRIDENT_PHOTONS
          IF (produce_pairs .AND. current%optical_depth_tri < 0.0_num) THEN
            CALL generate_pair_tri(current, trident_electron_species, &
                trident_positron_species)
```

Applied identically to `epoch1d`, `epoch2d`, and `epoch3d` (the
surrounding `qed_update_optical_depth` code is otherwise identical
across all three).

## 4. Verification

- **Reproduced** on unmodified `epoch2d` (default Makefile already
  builds `-DPHOTONS -DTRIDENT_PHOTONS`): a single-rank deck with an
  ultra-relativistic electron (`drift_px = 3000 * me * c`) driven by
  a `1e23 W/cm^2` boundary laser, `identify:electron` +
  `identify:photon` species only, `produce_pairs = F` (no trident
  species declared) — segfaults (signal 11) within 40 fs, exact stack
  above, confirmed via a temporary debug print showing
  `trident_electron_species = -1, trident_positron_species = -1` at
  the crash site.
- **Fix confirmed**: same deck, same binary rebuilt with the fix —
  runs to completion, exit code 0, no crash.
- **Positive path unaffected**: same deck with `produce_pairs = T`
  and an added `identify:positron` species — a temporary debug print
  confirmed `generate_pair_tri` still fires, with valid species
  indices (`trident_electron_species = 1`,
  `trident_positron_species = 3`, matching the deck's declared
  species order), and the run completes cleanly. Confirms the fix
  only removes the invalid-index path, not trident production itself.
- **Clean `gfortran` build** of all three dimensions with the fix.
- Debug prints used only for verification, removed before finalising;
  not part of the committed diff.

## 5. Scope note

Unrelated to the `boost_pairs`/`boost_muons` pair-boost integration
(`PAIR_BOOST_INTEGRATION.md`) — this is a pre-existing upstream-style
omission in the trident optical-depth update, not part of that
upscaling work. `produce_pairs = F` trident decks previously needed
the workaround documented in `Project_EPOCH`'s GPR tutorial (a dummy
`identify:positron` species purely to keep `check_qed_variables`
happy); that workaround is no longer necessary after this fix, though
it remains harmless if left in place.
