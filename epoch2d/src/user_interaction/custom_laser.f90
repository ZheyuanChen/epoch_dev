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

  ! A global switch to toggle between the 1D spatial loader and  2D spatiotemporal matrix without rewriting code
  !LOGICAL, SAVE :: use_2d_spatiotemporal = .TRUE.

! Persistent storage arrays for the 2D spatiotemporal profile
  INTEGER, SAVE :: n_y_points = 0
  INTEGER, SAVE :: n_t_points = 0
  REAL(num), ALLOCATABLE, DIMENSION(:), SAVE :: file_y_coords
  REAL(num), ALLOCATABLE, DIMENSION(:), SAVE :: file_t_coords
  REAL(num), ALLOCATABLE, DIMENSION(:,:), SAVE :: file_field_matrix
  LOGICAL, SAVE :: profile_loaded = .FALSE.


CONTAINS

  ! For now, we return a constant value for the time profile, but this could be extended to read from a file or compute a more complex function of time.
  FUNCTION custom_laser_time_profile(laser)
    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num) :: custom_laser_time_profile
    custom_laser_time_profile = 1.0_num
  END FUNCTION custom_laser_time_profile

  ! This subroutine reads a spatial profile from a file and interpolates it onto the local grid of each MPI rank. The file is expected to have a simple format:
    ! the first line contains the number of points, followed by lines of coordinate-value pairs.
    ! The profile file must be named 'spatial_profile.dat' and located in the directory specified by 'USE_DATA_DIRECTORY'.
  SUBROUTINE custom_laser_spatial_setup(laser)
    TYPE(laser_block), INTENT(INOUT) :: laser
    
    CHARACTER(LEN=c_max_path_length) :: filename
    INTEGER :: file_unit, i, n_points, err
    REAL(num) :: pos, weight
    INTEGER :: idx_low, idx_high

    ! Arrays to hold the raw file data
    REAL(num), ALLOCATABLE, DIMENSION(:) :: file_coords, file_values

    ! Only proceed for the 1D spatial custom path (not 2D spatiotemporal)
    IF (.NOT. laser%use_custom_profile .OR. laser%use_spatiotemporal) RETURN

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
    file_unit = 99 ! Explicit assignment to satisfy strict -std=f2003 standard
    
    ! --- 1. RANK 0 READS THE FILE ---
    IF (rank == 0) THEN
      OPEN(UNIT=file_unit, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=err)
      IF (err /= 0) THEN
        PRINT*, 'ERROR: Could not open laser profile file: ', TRIM(filename)
        CALL MPI_ABORT(mpi_comm_world, 1, err)
      END IF
      
      ! Read the number of points specified at the top of your data file
      READ(file_unit, *) n_points
      ALLOCATE(file_coords(n_points), file_values(n_points))
      
      ! Read coordinate-value pairs
        ! By the end of this loop, file_coords will hold the coordinates (e.g., Y or X positions) and file_values will hold the corresponding laser profile values at those coordinates.
      DO i = 1, n_points
        READ(file_unit, *) file_coords(i), file_values(i)
      END DO
      CLOSE(file_unit)
    END IF
    
    ! --- 2. BROADCAST DATA TO ALL MPI RANKS ---
    CALL MPI_BCAST(n_points, 1, MPI_INTEGER, 0, mpi_comm_world, err)
    
    IF (rank /= 0) ALLOCATE(file_coords(n_points), file_values(n_points))
    
    CALL MPI_BCAST(file_coords, n_points, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, err)
    CALL MPI_BCAST(file_values, n_points, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, err)

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
              IF (pos >= file_coords(idx_low) .AND. pos <= file_coords(idx_low+1)) THEN
                idx_high = idx_low + 1
                EXIT
              END IF
            END DO
            weight = (pos - file_coords(idx_low)) / (file_coords(idx_high) - file_coords(idx_low))
            laser%profile(i) = file_values(idx_low) + weight * (file_values(idx_high) - file_values(idx_low))
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
              IF (pos >= file_coords(idx_low) .AND. pos <= file_coords(idx_low+1)) THEN
                idx_high = idx_low + 1
                EXIT
              END IF
            END DO
            weight = (pos - file_coords(idx_low)) / (file_coords(idx_high) - file_coords(idx_low))
            laser%profile(i) = file_values(idx_low) + weight * (file_values(idx_high) - file_values(idx_low))
          END IF
        END DO

    END SELECT

    DEALLOCATE(file_coords, file_values)

  END SUBROUTINE custom_laser_spatial_setup

  ! Load a 2D spatiotemporal profile from the given filename.
  ! Absolute paths are used directly; relative paths are resolved from data_dir.
  ! Only loads once (guarded by profile_loaded) — the first laser to trigger
  ! this call determines which file is read.
  SUBROUTINE load_temporal_spatial_profile(profile_filename)
    CHARACTER(LEN=*), INTENT(IN) :: profile_filename
    INTEGER :: io_err, i, j, mpi_err
    CHARACTER(LEN=c_max_path_length) :: full_filename

    IF (profile_loaded) RETURN

    ! Resolve absolute vs relative path
    IF (profile_filename(1:1) == '/') THEN
      full_filename = TRIM(profile_filename)
    ELSE
      full_filename = TRIM(data_dir) // '/' // TRIM(profile_filename)
    END IF

    ! --- 1. RANK 0 READS DIMENSIONS AND COORDINATES ---
    IF (rank == 0) THEN
        OPEN(UNIT=100, FILE=TRIM(full_filename), STATUS='OLD', &
             ACTION='READ', IOSTAT=io_err)
        IF (io_err /= 0) THEN
            PRINT *, "ERROR: Could not open ", TRIM(full_filename)
            CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
        END IF

        READ(100, *) n_t_points, n_y_points
    END IF

    ! --- 2. BROADCAST DIMENSIONS SO ALL RANKS CAN ALLOCATE ---
    CALL MPI_BCAST(n_t_points, 1, MPI_INTEGER, 0, mpi_comm_world, mpi_err)
    CALL MPI_BCAST(n_y_points, 1, MPI_INTEGER, 0, mpi_comm_world, mpi_err)

    ALLOCATE(file_y_coords(n_y_points))
    ALLOCATE(file_t_coords(n_t_points))
    ALLOCATE(file_field_matrix(n_y_points, n_t_points))

    ! --- 3. RANK 0 READS DATA ---
    IF (rank == 0) THEN
        READ(100, *) file_y_coords
        READ(100, *) file_t_coords

        DO j = 1, n_t_points
            READ(100, *) file_field_matrix(:, j)
        END DO

        CLOSE(100)
    END IF

    ! --- 4. BROADCAST DATA TO ALL RANKS ---
    CALL MPI_BCAST(file_y_coords,    n_y_points,              &
                   MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_err)
    CALL MPI_BCAST(file_t_coords,    n_t_points,              &
                   MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_err)
    CALL MPI_BCAST(file_field_matrix, n_y_points * n_t_points, &
                   MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_err)

    profile_loaded = .TRUE.

    IF (rank == 0) THEN
        PRINT *, ">>> Custom 2D Spatiotemporal Profile Loaded Successfully! <<<"
        PRINT *, "    Grid Size: ", n_y_points, " (Spatial) x ", n_t_points, " (Temporal)"
    END IF

  END SUBROUTINE load_temporal_spatial_profile

  REAL(num) FUNCTION custom_laser_profile(laser, pos)
    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num), INTENT(IN) :: pos
    INTEGER :: idx_y, idx_t
    REAL(num) :: u, v, q11, q12, q21, q22
    CHARACTER(LEN=c_max_path_length) :: fname

    ! Ensure the 2D profile data is loaded into memory on first call.
    ! Use the deck-specified filename if given, otherwise the legacy default.
    IF (.NOT. profile_loaded) THEN
      IF (LEN_TRIM(laser%profile_data_file) > 0) THEN
        fname = laser%profile_data_file
      ELSE
        fname = 'temporal_spatial_profile.dat'
      END IF
      CALL load_temporal_spatial_profile(fname)
    END IF

    ! Default return value if coordinates fall completely outside our file scope
    custom_laser_profile = 0.0_num

    ! --- 1. Boundary & Guard Checks ---
    ! Spatial check (pos = current y coordinate on grid)
    IF (pos < file_y_coords(1) .OR. pos > file_y_coords(n_y_points)) RETURN
    
    ! Temporal check (time = current global simulation time from shared_data)
    IF (time < file_t_coords(1) .OR. time > file_t_coords(n_t_points)) RETURN


    ! --- 2. Locate the Bounding Cell Box ---
    ! Find lower index spatial bounding point
    
    ! The commented code should apply for an arbitraryly spaced laser profile grid. However, this is an O(n) search and is not efficient for large grids. 
    ! Instead, we can use a direct index calculation for uniform grids, which is O(1). The Python script guarantees uniform spacing using np.linspace. 
    
    !!!!!! Must ensure that the laser profile file is generated with uniform spacing for this optimization to be valid. !!!!!!

    !idx_y = 1
    !DO WHILE (file_y_coords(idx_y+1) < pos .AND. idx_y < n_y_points - 1)
    !   idx_y = idx_y + 1
    !END DO

    ! Find lower index temporal bounding point
    !idx_t = 1
    !DO WHILE (file_t_coords(idx_t+1) < time .AND. idx_t < n_t_points - 1)
    !   idx_t = idx_t + 1
    !END DO

    ! Direct index calculation (O(1)) — valid only for uniform grids. This is likely to be visited again if we decide to support non-uniform grids in the future.
    idx_y = INT((pos  - file_y_coords(1)) / (file_y_coords(2) - file_y_coords(1))) + 1
    idx_t = INT((time - file_t_coords(1)) / (file_t_coords(2) - file_t_coords(1))) + 1

    ! Clamp to valid interpolation range [1, n-1]
    idx_y = MAX(1, MIN(idx_y, n_y_points - 1))
    idx_t = MAX(1, MIN(idx_t, n_t_points - 1))

    ! --- 3. Bilinear Interpolation Math ---
    ! Compute normalised fractional positions within the grid cell
    u = (pos - file_y_coords(idx_y)) / (file_y_coords(idx_y+1) - file_y_coords(idx_y))
    v = (time - file_t_coords(idx_t)) / (file_t_coords(idx_t+1) - file_t_coords(idx_t))

    ! Grab the 4 surrounding pixel values from the data matrix
    q11 = file_field_matrix(idx_y,     idx_t)      ! Bottom-Left
    q21 = file_field_matrix(idx_y+1,   idx_t)      ! Top-Left
    q12 = file_field_matrix(idx_y,     idx_t+1)    ! Bottom-Right
    q22 = file_field_matrix(idx_y+1,   idx_t+1)    ! Top-Right

    ! Execute bilinear interpolation formula
    custom_laser_profile = (1.0_num - u) * (1.0_num - v) * q11 &
                         + u * (1.0_num - v) * q21             &
                         + (1.0_num - u) * v * q12             &
                         + u * v * q22

  END FUNCTION custom_laser_profile

END MODULE custom_laser


