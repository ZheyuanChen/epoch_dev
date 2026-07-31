! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE custom_laser

  USE shared_data
  IMPLICIT NONE

  INTEGER, PARAMETER :: custom_laser_lu = 150

  ! EPOCH's "num" kind is hardcoded to KIND(1.d0) (constants.F90), so the
  ! binary profile/phase files are always 8 bytes per value. Not derived via
  ! STORAGE_SIZE (F2008) since the rest of this codebase targets F2003.
  INTEGER, PARAMETER :: real_bytes = 8

CONTAINS

  ! For now, we return a constant value for the time profile, but this could be
  ! extended to read from a file or compute a more complex function of time.
  FUNCTION custom_laser_time_profile(laser)
    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num) :: custom_laser_time_profile
    custom_laser_time_profile = 1.0_num
  END FUNCTION custom_laser_time_profile



  ! Entry point called from attach_laser for every laser at deck-parse
  ! time. Dispatches to the spatiotemporal or static spatial loader; both
  ! read raw binary files (see load_binary_file). Loading happens here so
  ! that the MPI_BCAST calls run during setup when ALL ranks participate,
  ! avoiding the deadlock that occurs if loading is deferred to the
  ! per-boundary-cell timestepping loop.
  SUBROUTINE custom_laser_spatial_setup(laser)

    TYPE(laser_block), INTENT(INOUT) :: laser
    CHARACTER(LEN=c_max_path_length) :: filename

    IF (.NOT. laser%use_custom_profile) RETURN

    IF (laser%use_spatiotemporal) THEN
      IF (LEN_TRIM(laser%profile_data_file) > 0) THEN
        filename = laser%profile_data_file
      ELSE
        filename = 'temporal_spatial_profile.dat'
      END IF
      CALL load_temporal_spatial_profile(laser, filename)

      IF (laser%use_phase_from_file) THEN
        IF (LEN_TRIM(laser%phase_data_file) > 0) THEN
          filename = laser%phase_data_file
        ELSE
          filename = 'phase_profile.dat'
        END IF
        CALL load_phase_profile(laser, filename)
      END IF

      RETURN
    END IF

    CALL load_spatial_fields(laser)

  END SUBROUTINE custom_laser_spatial_setup



  ! Abort with a clear error unless the deck declared a valid uniform grid
  ! for this laser's binary profile/phase files. Required because the files
  ! carry no embedded shape header (per EPOCH's documented binary-file
  ! convention) -- n_transverse_points, profile_transverse_min/max and (on
  ! the spatiotemporal path) n_t_points with t_start/t_end together fully
  ! determine the uniform grid. need_time selects the spatiotemporal
  ! variant.
  SUBROUTINE check_file_grid_declared(laser, need_time)

    TYPE(laser_block), INTENT(IN) :: laser
    LOGICAL, INTENT(IN) :: need_time
    LOGICAL :: ok
    INTEGER :: mpi_err

    ok = laser%n_transverse_points >= 2 &
        .AND. laser%profile_transverse_max > laser%profile_transverse_min
    IF (need_time) THEN
      ok = ok .AND. laser%n_t_points >= 2 .AND. laser%t_end > laser%t_start
    END IF
    IF (ok) RETURN

    IF (rank == 0) THEN
      IF (need_time) THEN
        PRINT *, 'ERROR: use_spatiotemporal_profile = T requires ' &
            // 'n_t_points and n_transverse_points (each >= 2), ' &
            // 'profile_transverse_min < profile_transverse_max and ' &
            // 't_start < t_end to be set in the laser block.'
      ELSE
        PRINT *, 'ERROR: use_custom_profile = T with ' &
            // 'use_spatiotemporal_profile = F requires ' &
            // 'n_transverse_points (>= 2) and profile_transverse_min ' &
            // '< profile_transverse_max to be set in the laser block.'
      END IF
    END IF
    CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)

  END SUBROUTINE check_file_grid_declared



  ! Read a headerless raw binary array (access='stream', n_values of
  ! EPOCH's 8-byte REAL(num)) on rank 0 and broadcast it to all ranks.
  ! Aborts with a clear error if the file is missing or its size does not
  ! match the deck-declared shape -- the only available sanity check, since
  ! the file embeds no shape metadata. Absolute paths are used directly;
  ! relative paths are resolved from data_dir.
  SUBROUTINE load_binary_file(filename, array, n_values)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: n_values
    REAL(num), DIMENSION(n_values), INTENT(OUT) :: array
    INTEGER :: io_err, mpi_err
    INTEGER(KIND=8) :: expected_bytes, actual_bytes
    LOGICAL :: file_exists
    CHARACTER(LEN=c_max_path_length) :: full_filename

    IF (filename(1:1) == '/') THEN
      full_filename = TRIM(filename)
    ELSE
      full_filename = TRIM(data_dir) // '/' // TRIM(filename)
    END IF

    IF (rank == 0) THEN
      expected_bytes = INT(n_values, 8) * INT(real_bytes, 8)

      INQUIRE(FILE=TRIM(full_filename), EXIST=file_exists, SIZE=actual_bytes)
      IF (.NOT. file_exists) THEN
        PRINT *, 'ERROR: Could not find ', TRIM(full_filename)
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF
      IF (actual_bytes /= expected_bytes) THEN
        PRINT *, 'ERROR: ', TRIM(full_filename), ' is ', actual_bytes, &
            ' bytes; expected ', expected_bytes, &
            ' (product of the deck-declared point counts * 8 bytes)'
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF

      OPEN(UNIT=custom_laser_lu, FILE=TRIM(full_filename), STATUS='OLD', &
          ACCESS='STREAM', FORM='UNFORMATTED', ACTION='READ', IOSTAT=io_err)
      IF (io_err /= 0) THEN
        PRINT *, 'ERROR: Could not open ', TRIM(full_filename)
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF

      READ(custom_laser_lu) array
      CLOSE(custom_laser_lu)
    END IF

    CALL MPI_BCAST(array, n_values, mpireal, 0, mpi_comm_world, mpi_err)

  END SUBROUTINE load_binary_file



  ! Generic allocate-and-load for a 2D spatiotemporal matrix (n1 x n2,
  ! first axis fastest-varying) from a raw binary file, via load_binary_file.
  ! Shared by the boundary-laser profile/phase loaders below and by the
  ! interior laser-antenna module (laser_antenna.f90), so the file-format
  ! convention only needs to be right in one place.
  SUBROUTINE load_spatiotemporal_matrix(filename, n1, n2, matrix)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: n1, n2
    REAL(num), DIMENSION(:,:), POINTER :: matrix

    ALLOCATE(matrix(n1, n2))
    CALL load_binary_file(filename, matrix, n1 * n2)

  END SUBROUTINE load_spatiotemporal_matrix



  ! Generic bilinear interpolation of a 2D spatiotemporal matrix (as loaded
  ! by load_spatiotemporal_matrix) at transverse position 'pos' and sample
  ! time 't_sample', given the deck-declared uniform grid bounds/counts.
  ! Returns 0 outside the declared grid in either axis. Shared by
  ! custom_laser_profile/custom_laser_phase below and by laser_antenna.f90
  ! -- 't_sample' is passed explicitly (not read from the module 'time'
  ! variable) so callers that need a specific half-step time (rather than
  ! the ambient one) can supply it directly.
  REAL(num) FUNCTION sample_spatiotemporal_matrix(matrix, pos, t_sample, &
      pos_min, pos_max, n_pos, t_start, t_end, n_t)

    REAL(num), DIMENSION(:,:), INTENT(IN) :: matrix
    REAL(num), INTENT(IN) :: pos, t_sample, pos_min, pos_max
    REAL(num), INTENT(IN) :: t_start, t_end
    INTEGER, INTENT(IN) :: n_pos, n_t
    INTEGER :: idx_pos, idx_t
    REAL(num) :: dpos, dt, pos0, t0, u, v, q11, q12, q21, q22

    sample_spatiotemporal_matrix = 0.0_num

    IF (pos < pos_min .OR. pos > pos_max) RETURN
    IF (t_sample < t_start .OR. t_sample > t_end) RETURN

    dpos = (pos_max - pos_min) / REAL(n_pos - 1, num)
    dt = (t_end - t_start) / REAL(n_t - 1, num)

    idx_pos = INT((pos - pos_min) / dpos) + 1
    idx_t = INT((t_sample - t_start) / dt) + 1

    ! Clamp to valid interpolation range [1, n-1]
    idx_pos = MAX(1, MIN(idx_pos, n_pos - 1))
    idx_t = MAX(1, MIN(idx_t, n_t - 1))

    pos0 = pos_min + REAL(idx_pos - 1, num) * dpos
    t0 = t_start + REAL(idx_t - 1, num) * dt
    u = (pos - pos0) / dpos
    v = (t_sample - t0) / dt

    q11 = matrix(idx_pos,   idx_t)      ! Bottom-Left
    q21 = matrix(idx_pos+1, idx_t)      ! Top-Left
    q12 = matrix(idx_pos,   idx_t+1)    ! Bottom-Right
    q22 = matrix(idx_pos+1, idx_t+1)    ! Top-Right

    sample_spatiotemporal_matrix = (1.0_num - u) * (1.0_num - v) * q11 &
                                 + u * (1.0_num - v) * q21             &
                                 + (1.0_num - u) * v * q12             &
                                 + u * v * q22

  END FUNCTION sample_spatiotemporal_matrix



  ! Load a 2D spatiotemporal amplitude profile from a raw binary file into
  ! the given laser block: access='stream', no embedded header, column-major
  ! (transverse axis fastest-varying), n_transverse_points * n_t_points
  ! values of EPOCH's REAL(num) (always 8-byte, see real_bytes above).
  ! Only loads once per laser (guarded by laser%profile_loaded).
  SUBROUTINE load_temporal_spatial_profile(laser, profile_filename)

    TYPE(laser_block), INTENT(INOUT) :: laser
    CHARACTER(LEN=*), INTENT(IN) :: profile_filename

    IF (laser%profile_loaded) RETURN

    CALL check_file_grid_declared(laser, .TRUE.)

    CALL load_spatiotemporal_matrix(profile_filename, &
        laser%n_transverse_points, laser%n_t_points, laser%file_field_matrix)

    laser%profile_loaded = .TRUE.

    IF (rank == 0) THEN
      PRINT *, '>>> Custom 2D Spatiotemporal Profile Loaded Successfully! <<<'
      PRINT *, '    Grid Size: ', laser%n_transverse_points, &
          ' (Spatial) x ', laser%n_t_points, ' (Temporal)'
    END IF

  END SUBROUTINE load_temporal_spatial_profile



  ! Load a 2D spatiotemporal phase profile from a raw binary file into the
  ! given laser block. Identical file convention to the amplitude profile
  ! (see load_temporal_spatial_profile); shares the same deck-declared grid
  ! via laser%n_t_points/n_transverse_points/profile_transverse_min/max and
  ! t_start/t_end. Phase values are read as-is -- the Python-side (LASY)
  ! converter writes them already in EPOCH's sign/offset convention
  ! (phase = -phi + pi/2).
  SUBROUTINE load_phase_profile(laser, phase_filename)

    TYPE(laser_block), INTENT(INOUT) :: laser
    CHARACTER(LEN=*), INTENT(IN) :: phase_filename

    IF (laser%phase_loaded) RETURN

    CALL check_file_grid_declared(laser, .TRUE.)

    CALL load_spatiotemporal_matrix(phase_filename, &
        laser%n_transverse_points, laser%n_t_points, laser%file_phase_matrix)

    laser%phase_loaded = .TRUE.

    IF (rank == 0) THEN
      PRINT *, '>>> Custom 2D Spatiotemporal Phase Profile Loaded ' &
          // 'Successfully! <<<'
      PRINT *, '    Grid Size: ', laser%n_transverse_points, &
          ' (Spatial) x ', laser%n_t_points, ' (Temporal)'
    END IF

  END SUBROUTINE load_phase_profile



  ! Static spatial path: load a 1D amplitude line (and optionally a phase
  ! line) from raw binary files and linearly interpolate them onto
  ! laser%profile / laser%phase once at setup. The file grid is declared
  ! in the deck exactly as on the spatiotemporal path, minus the temporal
  ! elements; the files themselves are headerless (n_transverse_points
  ! values each). This replaced the earlier text format (point count
  ! header plus coordinate-value pairs), which is no longer supported.
  SUBROUTINE load_spatial_fields(laser)

    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num), ALLOCATABLE, DIMENSION(:) :: line, coords
    CHARACTER(LEN=c_max_path_length) :: filename
    INTEGER :: i, n
    REAL(num) :: d

    CALL check_file_grid_declared(laser, .FALSE.)

    n = laser%n_transverse_points
    ALLOCATE(line(n), coords(n))

    ! Reconstruct the uniform file-grid axis from the deck declaration
    d = (laser%profile_transverse_max - laser%profile_transverse_min) &
        / REAL(n - 1, num)
    DO i = 1, n
      coords(i) = laser%profile_transverse_min + REAL(i - 1, num) * d
    END DO

    IF (LEN_TRIM(laser%profile_data_file) > 0) THEN
      filename = laser%profile_data_file
    ELSE
      filename = 'spatial_profile.dat'
    END IF
    CALL load_binary_file(filename, line, n)
    CALL interp_line_to_boundary(laser, coords, line, laser%profile)

    IF (rank == 0) THEN
      PRINT *, '>>> Custom 1D Spatial Profile Loaded Successfully! <<<'
      PRINT *, '    Grid Size: ', n
    END IF

    IF (laser%use_phase_from_file) THEN
      IF (LEN_TRIM(laser%phase_data_file) > 0) THEN
        filename = laser%phase_data_file
      ELSE
        filename = 'phase_profile.dat'
      END IF
      CALL load_binary_file(filename, line, n)
      CALL interp_line_to_boundary(laser, coords, line, laser%phase)

      IF (rank == 0) THEN
        PRINT *, '>>> Custom 1D Spatial Phase Profile Loaded ' &
            // 'Successfully! <<<'
        PRINT *, '    Grid Size: ', n
      END IF
    END IF

    DEALLOCATE(line, coords)

  END SUBROUTINE load_spatial_fields



  ! Linearly interpolate a file-grid line onto the local section of the
  ! boundary this laser is attached to, writing into 'dest' (laser%profile
  ! or laser%phase, passed as a pointer so its original bounds are kept).
  ! The transverse axis is y for an x_min/x_max laser and x for a
  ! y_min/y_max laser. Cell-centre coordinates are used, consistent with
  ! the analytical evaluator which resolves deck variables at y(i) / x(i).
  SUBROUTINE interp_line_to_boundary(laser, coords, vals, dest)

    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num), DIMENSION(:), INTENT(IN) :: coords, vals
    REAL(num), DIMENSION(:), POINTER :: dest
    INTEGER :: i, n

    n = SIZE(coords)

    SELECT CASE(laser%boundary)

      CASE(c_bd_x_min, c_bd_x_max)
        DO i = 0, ny
          dest(i) = interp1d(y(i), coords, vals, n)
        END DO

      CASE(c_bd_y_min, c_bd_y_max)
        DO i = 0, nx
          dest(i) = interp1d(x(i), coords, vals, n)
        END DO

    END SELECT

  END SUBROUTINE interp_line_to_boundary



  ! Linear interpolation on a 1D regular grid. Returns the interpolated
  ! value at p. Clamps to boundary values for points outside the data
  ! range.
  REAL(num) FUNCTION interp1d(p, c, vals, n)

    REAL(num), INTENT(IN) :: p
    INTEGER, INTENT(IN) :: n
    REAL(num), DIMENSION(n), INTENT(IN) :: c, vals
    INTEGER :: i1
    REAL(num) :: u

    ! Direct index calculation (O(1)) -- valid for uniform grids
    i1 = INT((p - c(1)) / (c(2) - c(1))) + 1

    ! Clamp to valid interpolation range [1, n-1]
    i1 = MAX(1, MIN(i1, n - 1))

    u = (p - c(i1)) / (c(i1+1) - c(i1))

    ! Clamp fractional position for points outside the data range
    u = MAX(0.0_num, MIN(u, 1.0_num))

    interp1d = (1.0_num - u) * vals(i1) + u * vals(i1+1)

  END FUNCTION interp1d



  REAL(num) FUNCTION custom_laser_profile(laser, pos)
    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num), INTENT(IN) :: pos
    CHARACTER(LEN=c_max_path_length) :: fname

    ! Ensure this laser's 2D profile data is loaded into memory on first
    ! call. Use the deck-specified filename if given, otherwise the legacy
    ! default.
    IF (.NOT. laser%profile_loaded) THEN
      IF (LEN_TRIM(laser%profile_data_file) > 0) THEN
        fname = laser%profile_data_file
      ELSE
        fname = 'temporal_spatial_profile.dat'
      END IF
      CALL load_temporal_spatial_profile(laser, fname)
    END IF

    custom_laser_profile = sample_spatiotemporal_matrix( &
        laser%file_field_matrix, pos, time, &
        laser%profile_transverse_min, laser%profile_transverse_max, &
        laser%n_transverse_points, laser%t_start, laser%t_end, &
        laser%n_t_points)

  END FUNCTION custom_laser_profile

  ! Bilinear interpolation of this laser's spatiotemporal phase profile at
  ! spatial position 'pos' and the current simulation 'time'. Mirrors
  ! custom_laser_profile exactly, but reads from laser%file_phase_matrix
  ! (loaded by load_phase_profile).
  REAL(num) FUNCTION custom_laser_phase(laser, pos)
    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num), INTENT(IN) :: pos
    CHARACTER(LEN=c_max_path_length) :: fname

    ! Ensure this laser's 2D phase data is loaded into memory on first call.
    ! Use the deck-specified filename if given, otherwise the default.
    IF (.NOT. laser%phase_loaded) THEN
      IF (LEN_TRIM(laser%phase_data_file) > 0) THEN
        fname = laser%phase_data_file
      ELSE
        fname = 'phase_profile.dat'
      END IF
      CALL load_phase_profile(laser, fname)
    END IF

    custom_laser_phase = sample_spatiotemporal_matrix( &
        laser%file_phase_matrix, pos, time, &
        laser%profile_transverse_min, laser%profile_transverse_max, &
        laser%n_transverse_points, laser%t_start, laser%t_end, &
        laser%n_t_points)

  END FUNCTION custom_laser_phase

END MODULE custom_laser
