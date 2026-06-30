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

  ! Establish or validate the spatial/temporal grid shared by a single laser
  ! block's amplitude and phase profiles. The first loader to run for this
  ! laser (whichever of amplitude or phase is requested first) reads the
  ! coordinate arrays from its file into laser%file_transverse_coords /
  ! laser%file_t_coords. Any subsequent loader for the SAME laser instead
  ! validates that its file declares the same grid and consumes (discards)
  ! the coordinate lines so the reader advances to the data matrix.
  !
  ! Storage lives on the laser_block instance (not module-level), so two
  ! different laser blocks (e.g. one per transverse polarisation channel,
  ! distinguished by pol_angle) each load and keep their own independent
  ! grid and data, instead of colliding on shared state.
  !
  ! file_unit must be an open unit positioned just after the dimension line;
  ! only rank 0 performs file I/O.
  SUBROUTINE ensure_shared_grid(laser, file_unit, n_transverse_file, n_t_file)
    TYPE(laser_block), INTENT(INOUT) :: laser
    INTEGER, INTENT(IN) :: file_unit, n_transverse_file, n_t_file
    INTEGER :: mpi_err
    REAL(num), ALLOCATABLE, DIMENSION(:) :: scratch_pos, scratch_t

    IF (.NOT. ASSOCIATED(laser%file_transverse_coords)) THEN
      ! First profile to load for this laser: adopt this file's grid.
      laser%n_transverse_points = n_transverse_file
      laser%file_n_t_points = n_t_file
      ALLOCATE(laser%file_transverse_coords(laser%n_transverse_points))
      ALLOCATE(laser%file_t_coords(laser%file_n_t_points))

      IF (rank == 0) THEN
        READ(file_unit, *) laser%file_transverse_coords
        READ(file_unit, *) laser%file_t_coords
      END IF

      CALL MPI_BCAST(laser%file_transverse_coords, &
          laser%n_transverse_points, MPI_DOUBLE_PRECISION, 0, &
          mpi_comm_world, mpi_err)
      CALL MPI_BCAST(laser%file_t_coords, laser%file_n_t_points, &
          MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_err)
    ELSE
      ! Grid already established by this laser's other profile: the
      ! dimensions must match, since amplitude and phase are sampled on the
      ! same LASY grid for a given laser.
      IF (n_transverse_file /= laser%n_transverse_points &
          .OR. n_t_file /= laser%file_n_t_points) THEN
        IF (rank == 0) THEN
          PRINT *, "ERROR: phase and amplitude profile grids differ: ", &
                   n_transverse_file, " x ", n_t_file, " vs ", &
                   laser%n_transverse_points, " x ", laser%file_n_t_points
          CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
        END IF
      END IF

      ! Consume this file's coordinate lines on rank 0 to reach the matrix.
      IF (rank == 0) THEN
        ALLOCATE(scratch_pos(n_transverse_file), scratch_t(n_t_file))
        READ(file_unit, *) scratch_pos
        READ(file_unit, *) scratch_t
        DEALLOCATE(scratch_pos, scratch_t)
      END IF
    END IF

  END SUBROUTINE ensure_shared_grid


  ! Load a 2D spatiotemporal amplitude profile from the given filename into
  ! the given laser block. Absolute paths are used directly; relative paths
  ! are resolved from data_dir. Only loads once per laser (guarded by
  ! laser%profile_loaded).
  SUBROUTINE load_temporal_spatial_profile(laser, profile_filename)
    TYPE(laser_block), INTENT(INOUT) :: laser
    CHARACTER(LEN=*), INTENT(IN) :: profile_filename
    INTEGER :: io_err, j, mpi_err, n_t_file, n_transverse_file
    CHARACTER(LEN=c_max_path_length) :: full_filename

    IF (laser%profile_loaded) RETURN

    ! Resolve absolute vs relative path
    IF (profile_filename(1:1) == '/') THEN
      full_filename = TRIM(profile_filename)
    ELSE
      full_filename = TRIM(data_dir) // '/' // TRIM(profile_filename)
    END IF

    ! --- 1. RANK 0 READS DIMENSIONS ---
    IF (rank == 0) THEN
        OPEN(UNIT=custom_laser_lu, FILE=TRIM(full_filename), STATUS='OLD', &
             ACTION='READ', IOSTAT=io_err)
        IF (io_err /= 0) THEN
            PRINT *, "ERROR: Could not open ", TRIM(full_filename)
            CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
        END IF

        READ(custom_laser_lu, *) n_t_file, n_transverse_file
    END IF

    ! --- 2. BROADCAST DIMENSIONS SO ALL RANKS CAN ALLOCATE ---
    CALL MPI_BCAST(n_t_file, 1, MPI_INTEGER, 0, mpi_comm_world, mpi_err)
    CALL MPI_BCAST(n_transverse_file, 1, MPI_INTEGER, 0, mpi_comm_world, &
        mpi_err)

    ! --- 3. ESTABLISH/VALIDATE THE SHARED GRID, THEN ALLOCATE THE MATRIX ---
    CALL ensure_shared_grid(laser, custom_laser_lu, n_transverse_file, &
        n_t_file)
    ALLOCATE(laser%file_field_matrix(laser%n_transverse_points, &
        laser%file_n_t_points))

    ! --- 4. RANK 0 READS THE DATA MATRIX ---
    IF (rank == 0) THEN
        DO j = 1, laser%file_n_t_points
            READ(custom_laser_lu, *) laser%file_field_matrix(:, j)
        END DO

        CLOSE(custom_laser_lu)
    END IF

    ! --- 5. BROADCAST DATA TO ALL RANKS ---
    CALL MPI_BCAST(laser%file_field_matrix, &
        laser%n_transverse_points * laser%file_n_t_points, &
        MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_err)

    laser%profile_loaded = .TRUE.

    IF (rank == 0) THEN
        PRINT *, ">>> Custom 2D Spatiotemporal Profile Loaded Successfully! <<<"
        PRINT *, "    Grid Size: ", laser%n_transverse_points, &
            " (Spatial) x ", laser%file_n_t_points, " (Temporal)"
    END IF

  END SUBROUTINE load_temporal_spatial_profile


  ! Load a 2D spatiotemporal phase profile from the given filename into the
  ! given laser block. Identical file format to the amplitude profile
  ! (dimensions, transverse coords, t coords, then the matrix). The
  ! spatial/temporal grid is shared with this laser's amplitude profile via
  ! ensure_shared_grid, so only laser%file_phase_matrix is stored here. Phase
  ! values are read as-is — the Python-side (LASY) converter writes them
  ! already in EPOCH's sign/offset convention (phase = -phi + pi/2).
  SUBROUTINE load_phase_profile(laser, phase_filename)
    TYPE(laser_block), INTENT(INOUT) :: laser
    CHARACTER(LEN=*), INTENT(IN) :: phase_filename
    INTEGER :: io_err, j, mpi_err, n_t_file, n_transverse_file
    CHARACTER(LEN=c_max_path_length) :: full_filename

    IF (laser%phase_loaded) RETURN

    ! Resolve absolute vs relative path
    IF (phase_filename(1:1) == '/') THEN
      full_filename = TRIM(phase_filename)
    ELSE
      full_filename = TRIM(data_dir) // '/' // TRIM(phase_filename)
    END IF

    ! --- 1. RANK 0 READS DIMENSIONS ---
    IF (rank == 0) THEN
        OPEN(UNIT=custom_laser_lu, FILE=TRIM(full_filename), STATUS='OLD', &
             ACTION='READ', IOSTAT=io_err)
        IF (io_err /= 0) THEN
            PRINT *, "ERROR: Could not open ", TRIM(full_filename)
            CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
        END IF

        READ(custom_laser_lu, *) n_t_file, n_transverse_file
    END IF

    ! --- 2. BROADCAST DIMENSIONS SO ALL RANKS CAN ALLOCATE ---
    CALL MPI_BCAST(n_t_file, 1, MPI_INTEGER, 0, mpi_comm_world, mpi_err)
    CALL MPI_BCAST(n_transverse_file, 1, MPI_INTEGER, 0, mpi_comm_world, &
        mpi_err)

    ! --- 3. ESTABLISH/VALIDATE THE SHARED GRID, THEN ALLOCATE THE MATRIX ---
    CALL ensure_shared_grid(laser, custom_laser_lu, n_transverse_file, &
        n_t_file)
    ALLOCATE(laser%file_phase_matrix(laser%n_transverse_points, &
        laser%file_n_t_points))

    ! --- 4. RANK 0 READS THE DATA MATRIX ---
    IF (rank == 0) THEN
        DO j = 1, laser%file_n_t_points
            READ(custom_laser_lu, *) laser%file_phase_matrix(:, j)
        END DO

        CLOSE(custom_laser_lu)
    END IF

    ! --- 5. BROADCAST DATA TO ALL RANKS ---
    CALL MPI_BCAST(laser%file_phase_matrix, &
        laser%n_transverse_points * laser%file_n_t_points, &
        MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_err)

    laser%phase_loaded = .TRUE.

    IF (rank == 0) THEN
        PRINT *, ">>> Custom 2D Spatiotemporal Phase Profile Loaded " // &
            "Successfully! <<<"
        PRINT *, "    Grid Size: ", laser%n_transverse_points, &
            " (Spatial) x ", laser%file_n_t_points, " (Temporal)"
    END IF

  END SUBROUTINE load_phase_profile

  REAL(num) FUNCTION custom_laser_profile(laser, pos)
    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num), INTENT(IN) :: pos
    INTEGER :: idx_pos, idx_t
    REAL(num) :: u, v, q11, q12, q21, q22
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

    ! Default return value if coordinates fall completely outside our file
    ! scope.
    custom_laser_profile = 0.0_num

    ! --- 1. Boundary & Guard Checks ---
    ! Spatial check (pos = current transverse coordinate on grid)
    IF (pos < laser%file_transverse_coords(1) .OR. &
        pos > laser%file_transverse_coords(laser%n_transverse_points)) RETURN

    ! Temporal check (time = current global simulation time from shared_data)
    IF (time < laser%file_t_coords(1) &
        .OR. time > laser%file_t_coords(laser%file_n_t_points)) RETURN


    ! --- 2. Locate the Bounding Cell Box ---
    ! Find lower index spatial bounding point

    ! The commented code should apply for an arbitraryly spaced laser profile
    ! grid. However, this is an O(n) search and is not efficient for large
    ! grids.
    ! Instead, we can use a direct index calculation for uniform grids, which
    ! is O(1). The Python script guarantees uniform spacing using np.linspace.

    !!!!!! Must ensure that the laser profile file is generated with uniform
    !!!!!! spacing for this optimization to be valid. !!!!!!

    !idx_pos = 1
    !DO WHILE (laser%file_transverse_coords(idx_pos+1) < pos &
    !    .AND. idx_pos < laser%n_transverse_points - 1)
    !   idx_pos = idx_pos + 1
    !END DO

    ! Find lower index temporal bounding point
    !idx_t = 1
    !DO WHILE (laser%file_t_coords(idx_t+1) < time &
    !    .AND. idx_t < laser%file_n_t_points - 1)
    !   idx_t = idx_t + 1
    !END DO

    ! Direct index calculation (O(1)) — valid only for uniform grids. This is
    ! likely to be visited again if we decide to support non-uniform grids in
    ! the future.
    idx_pos = INT((pos - laser%file_transverse_coords(1)) &
        / (laser%file_transverse_coords(2) &
        - laser%file_transverse_coords(1))) + 1
    idx_t = INT((time - laser%file_t_coords(1)) &
        / (laser%file_t_coords(2) - laser%file_t_coords(1))) + 1

    ! Clamp to valid interpolation range [1, n-1]
    idx_pos = MAX(1, MIN(idx_pos, laser%n_transverse_points - 1))
    idx_t = MAX(1, MIN(idx_t, laser%file_n_t_points - 1))

    ! --- 3. Bilinear Interpolation Math ---
    ! Compute normalised fractional positions within the grid cell
    u = (pos - laser%file_transverse_coords(idx_pos)) &
        / (laser%file_transverse_coords(idx_pos+1) &
        - laser%file_transverse_coords(idx_pos))
    v = (time - laser%file_t_coords(idx_t)) &
        / (laser%file_t_coords(idx_t+1) - laser%file_t_coords(idx_t))

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
  ! (loaded by load_phase_profile) on the shared laser%file_transverse_coords
  ! / laser%file_t_coords grid.
  REAL(num) FUNCTION custom_laser_phase(laser, pos)
    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num), INTENT(IN) :: pos
    INTEGER :: idx_pos, idx_t
    REAL(num) :: u, v, q11, q12, q21, q22
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

    ! Default return value if coordinates fall completely outside our file
    ! scope. The amplitude envelope is likewise zero there, so the phase value
    ! is immaterial.
    custom_laser_phase = 0.0_num

    ! --- 1. Boundary & Guard Checks ---
    ! Spatial check (pos = current transverse coordinate on grid)
    IF (pos < laser%file_transverse_coords(1) .OR. &
        pos > laser%file_transverse_coords(laser%n_transverse_points)) RETURN

    ! Temporal check (time = current global simulation time from shared_data)
    IF (time < laser%file_t_coords(1) &
        .OR. time > laser%file_t_coords(laser%file_n_t_points)) RETURN

    ! --- 2. Locate the Bounding Cell Box ---
    ! Direct index calculation (O(1)) — valid only for uniform grids, which the
    ! Python generator guarantees via np.linspace (same grid as the amplitude).
    idx_pos = INT((pos - laser%file_transverse_coords(1)) &
        / (laser%file_transverse_coords(2) &
        - laser%file_transverse_coords(1))) + 1
    idx_t = INT((time - laser%file_t_coords(1)) &
        / (laser%file_t_coords(2) - laser%file_t_coords(1))) + 1

    ! Clamp to valid interpolation range [1, n-1]
    idx_pos = MAX(1, MIN(idx_pos, laser%n_transverse_points - 1))
    idx_t = MAX(1, MIN(idx_t, laser%file_n_t_points - 1))

    ! --- 3. Bilinear Interpolation Math ---
    ! Compute normalised fractional positions within the grid cell
    u = (pos - laser%file_transverse_coords(idx_pos)) &
        / (laser%file_transverse_coords(idx_pos+1) &
        - laser%file_transverse_coords(idx_pos))
    v = (time - laser%file_t_coords(idx_t)) &
        / (laser%file_t_coords(idx_t+1) - laser%file_t_coords(idx_t))

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
