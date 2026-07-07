# CLAUDE.md — epoch_dev

## Project overview

This is a fork of the [EPOCH](https://github.com/epochpic/epoch) particle-in-cell plasma physics code (Fortran 90 + MPI). The repo contains `epoch1d`, `epoch2d`, and `epoch3d` sub-projects, each with its own Makefile and source tree. The `SDF` submodule provides the I/O library.

The `my-epoch-mods` branch carries custom modifications on top of upstream EPOCH, focused on custom laser profile injection and boundary condition changes.

## Key custom modifications

### epoch2d (mature — more extensively modified)

- **Custom spatio-temporal laser profile injection**: load laser profiles from external files instead of built-in analytic profiles. Deck-level configuration via `use_custom_profile`, `use_spatiotemporal_profile`, `profile_data_file`, `use_phase_from_file`, and `phase_data_file`.
- **Binary profile/phase file format**: `profile_data_file`/`phase_data_file` are raw binary (`access='stream'`, no embedded header), matching EPOCH's own documented binary-file convention — shape and bounds are declared in the deck, not in the file. Required deck elements (only when `use_spatiotemporal_profile = T`): `n_t_points`/`n_t`, `n_transverse_points`/`n_y`, `profile_transverse_min`/`y_min`, `profile_transverse_max`/`y_max` (aliases either way); the temporal grid extent reuses the laser's existing `t_start`/`t_end`. The grid is assumed uniform. A file-size check aborts with a clear error if the file doesn't match the declared shape. This replaced an earlier text `.dat` format (header line + coordinate rows + data rows) that is no longer supported for the spatiotemporal path.
- **Per-laser file-profile storage**: the loaded amplitude/phase data matrices live on each `laser_block` instance (not module-level singletons), so two custom-file laser blocks on the same boundary — e.g. one per transverse polarisation channel via distinct `pol_angle` values (0 and π/2 driving `source1`/`source2` independently) — load genuinely independent files instead of the second silently aliasing the first's data.
- **Half-cell coordinate offset fix**: corrected a half-cell offset in the custom laser profile sampling grid.
- **Parser and photon package changes**: modifications to `shunt.F90` and `photons.F90`.
- **Makefile customisation**.
- Modified files: `src/laser.f90`, `src/shared_data.F90`, `src/deck/deck_laser_block.f90`, `src/user_interaction/custom_laser.f90`, `src/parser/shunt.F90`, `src/physics_packages/photons.F90`, `Makefile`.

### epoch3d

- **Custom spatio-temporal laser profile injection**: same feature as epoch2d, ported across.
- **Outflow boundary conditions**: modified outflow BCs on all four boundaries (x_min, x_max, y_min, y_max).
- **Makefile customisation**.
- Modified files: `src/laser.f90`, `src/shared_data.F90`, `src/deck/deck_laser_block.f90`, `src/user_interaction/custom_laser.f90`, `Makefile`.

### epoch1d

No custom modifications.

## Build

```bash
cd epoch3d   # or epoch2d
make COMPILER=gfortran
# or: make COMPILER=intel
```

The binary lands in `<subproject>/bin/epoch<N>d`. The `obj/` and `bin/` directories are git-ignored per sub-project.

## Language and style

- Fortran 90 free-form (`.f90` / `.F90`).
- Use British English in all comments, commit messages, and documentation.
- Follow the existing EPOCH coding style (see `CODING_STYLE` in repo root) and
  `CONTRIBUTING.md`. When writing or editing Fortran source (and Makefiles),
  apply these rules as you go rather than fixing them up after the fact:
  - All lines, including comments, ≤80 columns.
  - No trailing whitespace, no tab characters.
  - Indentation: 2 spaces; continuation lines: 4 spaces deeper than the
    statement they continue.
  - `END IF` / `END DO`, not `ENDIF` / `ENDDO`; F90 operators (`==`, `<`,
    etc.), not F77 (`.EQ.`, `.LT.`, etc.).
  - Fortran keywords/intrinsics uppercase (`CALL`, `IF`, `SUBROUTINE`, ...);
    user-defined names lowercase. MPI constants/routines are the exception
    (uppercase).
  - Don't put a long explanation in a trailing `! comment` on a code line if
    it pushes the line over 80 columns — put it on its own line(s) above the
    statement instead.
  - Commit messages: subject line ≤50 chars, blank line, body wrapped at
    ≤72 chars per line.

## Daily logging

At the end of each session (or roughly every 5 messages), append a concise bullet-point summary of work done to:

```
~/Desktop/Project_EPOCH/daily_log/YYYY/MM/DD.md
```

under an `## epoch_dev` section.

- Create the file and parent directories if they don't exist yet.
- New files start with `# Daily Log — DD Month YYYY`.
- If other sections already exist, append — never alter them.
- Use British English, bullet points, not paragraphs.
- Only write to the daily log path — don't modify anything else in `Project_EPOCH`.
- Catch-up fallback: if today's file lacks an `## epoch_dev` section but work has already been done, backfill it.
