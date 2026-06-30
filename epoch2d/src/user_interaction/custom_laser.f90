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

  ! This subroutine reads a spatial profile from a file and interpolates it onto
  ! the local grid of each MPI rank. The file is expected to have a simple
  ! format:
    ! the first line contains the number of points, followed by lines of
    ! coordinate-value pairs.
    ! The profile file must be named 'spatial_profile.dat' and located in the
    ! directory specified by 'USE_DATA_DIRECTORY'.
  SUBROUTINE custom_laser_spatial_setup(laser)
    TYPE(laser_block), INTENT(INOUT) :: laser

    CHARACTER(LEN=c_max_path_length) :: filename
    INTEGER :: file_unit, i, n_points, err
    REAL(num) :: pos, weight
    INTEGER :: idx_low, idx_high

    ! Arrays to hold the raw file data
    REAL(num), ALLOCATABLE, DIMENSION(:) :: file_coords, file_values

    IF (.NOT. laser%use_custom_profile) RETURN

    ! 2D spatiotemporal path: pre-load amplitude (and optionally phase) profiles
    ! here so that the MPI_BCAST calls happen during setup when ALL ranks
    ! participate, avoiding the deadlock that occurs if loading is deferred to
    ! the per-boundary-cell timestepping loop.
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

    ! Resolve the profile data filename:
    !   - If the user specified profile_data_file in the deck, use it
    !     (absolute paths used as-is, relative paths resolved from data_dir)
    !   - Otherwise fall back to the legacy default 'spatial_profile.dat'
    IF (LEN_TRIM(laser%profile_data_file) > 0) THEN
      IF (laser%profile_data_file(1:1) == '/') THEN
        filename = TRIM(laser%profile_data_file)
      ELSE
        filename = TRIM(data_dir) // '/' // TRIM(laser%profile_data_file)
      END IF
    ELSE
      filename = TRIM(data_dir) // '/' // 'spatial_profile.dat'
    END IF
    file_unit = custom_laser_lu

    ! --- 1. RANK 0 READS THE FILE ---
    IF (rank == 0) THEN
      OPEN(UNIT=file_unit, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', &
          IOSTAT=err)
      IF (err /= 0) THEN
        PRINT*, 'ERROR: Could not open laser profile file: ', TRIM(filename)
        CALL MPI_ABORT(mpi_comm_world, 1, err)
      END IF

      ! Read the number of points specified at the top of your data file
      READ(file_unit, *) n_points
      ALLOCATE(file_coords(n_points), file_values(n_points))

      ! Read coordinate-value pairs
        ! By the end of this loop, file_coords will hold the coordinates (e.g.,
        ! Y or X positions) and file_values will hold the corresponding laser
        ! profile values at those coordinates.
      DO i = 1, n_points
        READ(file_unit, *) file_coords(i), file_values(i)
      END DO
      CLOSE(file_unit)
    END IF

    ! --- 2. BROADCAST DATA TO ALL MPI RANKS ---
    CALL MPI_BCAST(n_points, 1, MPI_INTEGER, 0, mpi_comm_world, err)

    IF (rank /= 0) ALLOCATE(file_coords(n_points), file_values(n_points))

    CALL MPI_BCAST(file_coords, n_points, MPI_DOUBLE_PRECISION, 0, &
        mpi_comm_world, err)
    CALL MPI_BCAST(file_values, n_points, MPI_DOUBLE_PRECISION, 0, &
        mpi_comm_world, err)

    ! --- 3. INTERPOLATE ONTO THE LOCAL PROCESSOR GRID ---
    SELECT CASE(laser%boundary)

      CASE(c_bd_x_min, c_bd_x_max)
        ! Laser is on a vertical boundary, varying along the Y-axis
        DO i = 0, ny
          ! Use cell-centre coordinate y(i), consistent with the analytical
          ! evaluator which resolves 'y' at y(pack_iy).
          pos = y(i)

          ! Perform a simple linear search & interpolation from file data
          IF (pos <= file_coords(1)) THEN
            laser%profile(i) = file_values(1)
          ELSE IF (pos >= file_coords(n_points)) THEN
            laser%profile(i) = file_values(n_points)
          ELSE
            DO idx_low = 1, n_points - 1
              IF (pos >= file_coords(idx_low) &
                  .AND. pos <= file_coords(idx_low+1)) THEN
                idx_high = idx_low + 1
                EXIT
              END IF
            END DO
            weight = (pos - file_coords(idx_low)) &
                / (file_coords(idx_high) - file_coords(idx_low))
            laser%profile(i) = file_values(idx_low) &
                + weight * (file_values(idx_high) - file_values(idx_low))
          END IF
        END DO

      CASE(c_bd_y_min, c_bd_y_max)
        ! Laser is on a horizontal boundary, varying along the X-axis
        DO i = 0, nx
          pos = x(i)

          ! Perform a simple linear search & interpolation from file data
          IF (pos <= file_coords(1)) THEN
            laser%profile(i) = file_values(1)
          ELSE IF (pos >= file_coords(n_points)) THEN
            laser%profile(i) = file_values(n_points)
          ELSE
            DO idx_low = 1, n_points - 1
              IF (pos >= file_coords(idx_low) &
                  .AND. pos <= file_coords(idx_low+1)) THEN
                idx_high = idx_low + 1
                EXIT
              END IF
            END DO
            weight = (pos - file_coords(idx_low)) &
                / (file_coords(idx_high) - file_coords(idx_low))
            laser%profile(i) = file_values(idx_low) &
                + weight * (file_values(idx_high) - file_values(idx_low))
          END IF
        END DO

    END SELECT

    DEALLOCATE(file_coords, file_values)

  END SUBROUTINE custom_laser_spatial_setup

  ! Abort with a clear error if the deck didn't declare a valid spatiotemporal
  ! grid for this laser. Required because the binary profile/phase files
  ! carry no embedded shape header (per EPOCH's documented binary-file
  ! convention) -- n_t_points, n_transverse_points, profile_transverse_min/max
  ! and t_start/t_end together fully determine the uniform grid.
  SUBROUTINE check_spatiotemporal_grid_declared(laser)
    TYPE(laser_block), INTENT(IN) :: laser
    INTEGER :: mpi_err

    IF (laser%n_t_points > 0 .AND. laser%n_transverse_points > 0 &
        .AND. laser%profile_transverse_max > laser%profile_transverse_min &
        .AND. laser%t_end > laser%t_start) RETURN

    IF (rank == 0) THEN
      PRINT *, "ERROR: use_spatiotemporal_profile = T requires " // &
          "n_t_points, n_transverse_points, profile_transverse_min, " // &
          "profile_transverse_max (and t_start < t_end) to be set in " // &
          "the laser block."
    END IF
    CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)

  END SUBROUTINE check_spatiotemporal_grid_declared


  ! Load a 2D spatiotemporal amplitude profile from a raw binary file into
  ! the given laser block: access='stream', no embedded header, column-major
  ! (transverse axis fastest-varying), n_transverse_points * n_t_points
  ! values of EPOCH's REAL(num) (always 8-byte, see real_bytes above).
  ! Absolute paths are used directly; relative paths are resolved from
  ! data_dir. Only loads once per laser (guarded by laser%profile_loaded).
  SUBROUTINE load_temporal_spatial_profile(laser, profile_filename)
    TYPE(laser_block), INTENT(INOUT) :: laser
    CHARACTER(LEN=*), INTENT(IN) :: profile_filename
    INTEGER :: io_err, mpi_err
    INTEGER(KIND=8) :: expected_bytes, actual_bytes
    LOGICAL :: file_exists
    CHARACTER(LEN=c_max_path_length) :: full_filename

    IF (laser%profile_loaded) RETURN

    CALL check_spatiotemporal_grid_declared(laser)

    ! Resolve absolute vs relative path
    IF (profile_filename(1:1) == '/') THEN
      full_filename = TRIM(profile_filename)
    ELSE
      full_filename = TRIM(data_dir) // '/' // TRIM(profile_filename)
    END IF

    ALLOCATE(laser%file_field_matrix(laser%n_transverse_points, &
        laser%n_t_points))

    IF (rank == 0) THEN
      expected_bytes = INT(laser%n_transverse_points, 8) &
          * INT(laser%n_t_points, 8) * INT(real_bytes, 8)

      INQUIRE(FILE=TRIM(full_filename), EXIST=file_exists, SIZE=actual_bytes)
      IF (.NOT. file_exists) THEN
        PRINT *, "ERROR: Could not find ", TRIM(full_filename)
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF
      IF (actual_bytes /= expected_bytes) THEN
        PRINT *, "ERROR: ", TRIM(full_filename), " is ", actual_bytes, &
            " bytes; expected ", expected_bytes, &
            " (n_transverse_points * n_t_points * 8 bytes, from the deck)"
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF

      OPEN(UNIT=custom_laser_lu, FILE=TRIM(full_filename), STATUS='OLD', &
          ACCESS='STREAM', FORM='UNFORMATTED', ACTION='READ', IOSTAT=io_err)
      IF (io_err /= 0) THEN
        PRINT *, "ERROR: Could not open ", TRIM(full_filename)
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF

      READ(custom_laser_lu) laser%file_field_matrix
      CLOSE(custom_laser_lu)
    END IF

    CALL MPI_BCAST(laser%file_field_matrix, &
        laser%n_transverse_points * laser%n_t_points, &
        mpireal, 0, mpi_comm_world, mpi_err)

    laser%profile_loaded = .TRUE.

    IF (rank == 0) THEN
        PRINT *, ">>> Custom 2D Spatiotemporal Profile Loaded Successfully! <<<"
        PRINT *, "    Grid Size: ", laser%n_transverse_points, &
            " (Spatial) x ", laser%n_t_points, " (Temporal)"
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
    INTEGER :: io_err, mpi_err
    INTEGER(KIND=8) :: expected_bytes, actual_bytes
    LOGICAL :: file_exists
    CHARACTER(LEN=c_max_path_length) :: full_filename

    IF (laser%phase_loaded) RETURN

    CALL check_spatiotemporal_grid_declared(laser)

    ! Resolve absolute vs relative path
    IF (phase_filename(1:1) == '/') THEN
      full_filename = TRIM(phase_filename)
    ELSE
      full_filename = TRIM(data_dir) // '/' // TRIM(phase_filename)
    END IF

    ALLOCATE(laser%file_phase_matrix(laser%n_transverse_points, &
        laser%n_t_points))

    IF (rank == 0) THEN
      expected_bytes = INT(laser%n_transverse_points, 8) &
          * INT(laser%n_t_points, 8) * INT(real_bytes, 8)

      INQUIRE(FILE=TRIM(full_filename), EXIST=file_exists, SIZE=actual_bytes)
      IF (.NOT. file_exists) THEN
        PRINT *, "ERROR: Could not find ", TRIM(full_filename)
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF
      IF (actual_bytes /= expected_bytes) THEN
        PRINT *, "ERROR: ", TRIM(full_filename), " is ", actual_bytes, &
            " bytes; expected ", expected_bytes, &
            " (n_transverse_points * n_t_points * 8 bytes, from the deck)"
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF

      OPEN(UNIT=custom_laser_lu, FILE=TRIM(full_filename), STATUS='OLD', &
          ACCESS='STREAM', FORM='UNFORMATTED', ACTION='READ', IOSTAT=io_err)
      IF (io_err /= 0) THEN
        PRINT *, "ERROR: Could not open ", TRIM(full_filename)
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF

      READ(custom_laser_lu) laser%file_phase_matrix
      CLOSE(custom_laser_lu)
    END IF

    CALL MPI_BCAST(laser%file_phase_matrix, &
        laser%n_transverse_points * laser%n_t_points, &
        mpireal, 0, mpi_comm_world, mpi_err)

    laser%phase_loaded = .TRUE.

    IF (rank == 0) THEN
        PRINT *, ">>> Custom 2D Spatiotemporal Phase Profile Loaded " // &
            "Successfully! <<<"
        PRINT *, "    Grid Size: ", laser%n_transverse_points, &
            " (Spatial) x ", laser%n_t_points, " (Temporal)"
    END IF

  END SUBROUTINE load_phase_profile

  REAL(num) FUNCTION custom_laser_profile(laser, pos)
    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num), INTENT(IN) :: pos
    INTEGER :: idx_pos, idx_t
    REAL(num) :: dy, dt, pos0, t0, u, v, q11, q12, q21, q22
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

    ! Default return value if coordinates fall completely outside the
    ! deck-declared grid.
    custom_laser_profile = 0.0_num

    ! --- 1. Boundary & Guard Checks ---
    IF (pos < laser%profile_transverse_min &
        .OR. pos > laser%profile_transverse_max) RETURN
    IF (time < laser%t_start .OR. time > laser%t_end) RETURN

    ! --- 2. Locate the Bounding Cell Box ---
    ! The grid is uniform by construction (deck-declared bounds/counts), so
    ! the cell spacing and bounding indices are computed directly -- no
    ! stored coordinate array to search.
    dy = (laser%profile_transverse_max - laser%profile_transverse_min) &
        / REAL(laser%n_transverse_points - 1, num)
    dt = (laser%t_end - laser%t_start) / REAL(laser%n_t_points - 1, num)

    idx_pos = INT((pos - laser%profile_transverse_min) / dy) + 1
    idx_t = INT((time - laser%t_start) / dt) + 1

    ! Clamp to valid interpolation range [1, n-1]
    idx_pos = MAX(1, MIN(idx_pos, laser%n_transverse_points - 1))
    idx_t = MAX(1, MIN(idx_t, laser%n_t_points - 1))

    ! --- 3. Bilinear Interpolation Math ---
    pos0 = laser%profile_transverse_min + REAL(idx_pos - 1, num) * dy
    t0 = laser%t_start + REAL(idx_t - 1, num) * dt
    u = (pos - pos0) / dy
    v = (time - t0) / dt

    ! Grab the 4 surrounding pixel values from the data matrix
    q11 = laser%file_field_matrix(idx_pos,   idx_t)      ! Bottom-Left
    q21 = laser%file_field_matrix(idx_pos+1, idx_t)      ! Top-Left
    q12 = laser%file_field_matrix(idx_pos,   idx_t+1)    ! Bottom-Right
    q22 = laser%file_field_matrix(idx_pos+1, idx_t+1)    ! Top-Right

    ! Execute bilinear interpolation formula
    custom_laser_profile = (1.0_num - u) * (1.0_num - v) * q11 &
                         + u * (1.0_num - v) * q21             &
                         + (1.0_num - u) * v * q12             &
                         + u * v * q22

  END FUNCTION custom_laser_profile

  ! Bilinear interpolation of this laser's spatiotemporal phase profile at
  ! spatial position 'pos' and the current simulation 'time'. Mirrors
  ! custom_laser_profile exactly, but reads from laser%file_phase_matrix
  ! (loaded by load_phase_profile).
  REAL(num) FUNCTION custom_laser_phase(laser, pos)
    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num), INTENT(IN) :: pos
    INTEGER :: idx_pos, idx_t
    REAL(num) :: dy, dt, pos0, t0, u, v, q11, q12, q21, q22
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

    ! Default return value if coordinates fall completely outside the
    ! deck-declared grid. The amplitude envelope is likewise zero there, so
    ! the phase value is immaterial.
    custom_laser_phase = 0.0_num

    ! --- 1. Boundary & Guard Checks ---
    IF (pos < laser%profile_transverse_min &
        .OR. pos > laser%profile_transverse_max) RETURN
    IF (time < laser%t_start .OR. time > laser%t_end) RETURN

    ! --- 2. Locate the Bounding Cell Box ---
    dy = (laser%profile_transverse_max - laser%profile_transverse_min) &
        / REAL(laser%n_transverse_points - 1, num)
    dt = (laser%t_end - laser%t_start) / REAL(laser%n_t_points - 1, num)

    idx_pos = INT((pos - laser%profile_transverse_min) / dy) + 1
    idx_t = INT((time - laser%t_start) / dt) + 1

    ! Clamp to valid interpolation range [1, n-1]
    idx_pos = MAX(1, MIN(idx_pos, laser%n_transverse_points - 1))
    idx_t = MAX(1, MIN(idx_t, laser%n_t_points - 1))

    ! --- 3. Bilinear Interpolation Math ---
    pos0 = laser%profile_transverse_min + REAL(idx_pos - 1, num) * dy
    t0 = laser%t_start + REAL(idx_t - 1, num) * dt
    u = (pos - pos0) / dy
    v = (time - t0) / dt

    ! Grab the 4 surrounding pixel values from the phase matrix
    q11 = laser%file_phase_matrix(idx_pos,   idx_t)      ! Bottom-Left
    q21 = laser%file_phase_matrix(idx_pos+1, idx_t)      ! Top-Left
    q12 = laser%file_phase_matrix(idx_pos,   idx_t+1)    ! Bottom-Right
    q22 = laser%file_phase_matrix(idx_pos+1, idx_t+1)    ! Top-Right

    ! Execute bilinear interpolation formula
    custom_laser_phase = (1.0_num - u) * (1.0_num - v) * q11 &
                       + u * (1.0_num - v) * q21             &
                       + (1.0_num - u) * v * q12             &
                       + u * v * q22

  END FUNCTION custom_laser_phase

END MODULE custom_laser
