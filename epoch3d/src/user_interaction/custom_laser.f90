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

  FUNCTION custom_laser_time_profile(laser)

    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num) :: custom_laser_time_profile

    custom_laser_time_profile = 1.0_num

  END FUNCTION custom_laser_time_profile

  ! Read a 2D spatial profile from a data file and bilinearly interpolate
  ! it onto the local grid of each MPI rank.
  !
  ! The boundary determines which two coordinate axes the profile spans:
  !   x_min/x_max -> profile in (y, z)
  !   y_min/y_max -> profile in (x, z)
  !   z_min/z_max -> profile in (x, y)
  !
  ! Expected file format:
  !   Line 1:    n1  n2             (integer counts for coord1, coord2)
  !   Line 2:    coord1 values      (n1 values, e.g. y coordinates)
  !   Line 3:    coord2 values      (n2 values, e.g. z coordinates)
  !   Lines 4+:  n2 rows of n1 values each  (profile values)
  !              Row j contains profile(:, j) for coord2(j)
  SUBROUTINE custom_laser_spatial_setup(laser)
    TYPE(laser_block), INTENT(INOUT) :: laser

    CHARACTER(LEN=c_max_path_length) :: filename
    INTEGER :: file_unit, i, j, n1, n2, err, mpi_err
    REAL(num) :: pos1, pos2, u, v
    INTEGER :: i1, i2

    REAL(num), ALLOCATABLE, DIMENSION(:) :: file_c1, file_c2
    REAL(num), ALLOCATABLE, DIMENSION(:,:) :: file_vals

    ! Only proceed for the spatial custom path (not spatiotemporal)
    IF (.NOT. laser%use_custom_profile .OR. laser%use_spatiotemporal) RETURN

    ! Resolve the profile data filename
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
      OPEN(UNIT=file_unit, FILE=TRIM(filename), STATUS='OLD', &
           ACTION='READ', IOSTAT=err)
      IF (err /= 0) THEN
        PRINT *, 'ERROR: Could not open laser profile file: ', TRIM(filename)
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF

      READ(file_unit, *) n1, n2
    END IF

    ! --- 2. BROADCAST DIMENSIONS SO ALL RANKS CAN ALLOCATE ---
    CALL MPI_BCAST(n1, 1, MPI_INTEGER, 0, mpi_comm_world, mpi_err)
    CALL MPI_BCAST(n2, 1, MPI_INTEGER, 0, mpi_comm_world, mpi_err)

    ALLOCATE(file_c1(n1), file_c2(n2), file_vals(n1, n2))

    ! --- 3. RANK 0 READS COORDINATES AND DATA ---
    IF (rank == 0) THEN
      READ(file_unit, *) file_c1
      READ(file_unit, *) file_c2
      DO j = 1, n2
        READ(file_unit, *) file_vals(:, j)
      END DO
      CLOSE(file_unit)
    END IF

    ! --- 4. BROADCAST DATA TO ALL RANKS ---
    CALL MPI_BCAST(file_c1, n1, MPI_DOUBLE_PRECISION, 0, &
                   mpi_comm_world, mpi_err)
    CALL MPI_BCAST(file_c2, n2, MPI_DOUBLE_PRECISION, 0, &
                   mpi_comm_world, mpi_err)
    CALL MPI_BCAST(file_vals, n1 * n2, MPI_DOUBLE_PRECISION, 0, &
                   mpi_comm_world, mpi_err)

    IF (rank == 0) THEN
      PRINT *, '>>> Custom 2D Spatial Profile Loaded Successfully! <<<'
      PRINT *, '    Grid Size: ', n1, ' x ', n2
    END IF

    ! --- 5. BILINEARLY INTERPOLATE ONTO THE LOCAL PROCESSOR GRID ---
    ! The profile array index convention matches allocate_with_boundary:
    !   x_min/x_max -> profile(0:ny, 0:nz), coord1=y, coord2=z
    !   y_min/y_max -> profile(0:nx, 0:nz), coord1=x, coord2=z
    !   z_min/z_max -> profile(0:nx, 0:ny), coord1=x, coord2=y
    SELECT CASE(laser%boundary)

      CASE(c_bd_x_min, c_bd_x_max)
        DO j = 0, nz
          DO i = 0, ny
            ! Use cell-centre coordinates, consistent with the analytical
            ! evaluator which resolves deck variables at y(i), z(j).
            pos1 = y(i)
            pos2 = z(j)
            laser%profile(i, j) = interp2d(pos1, pos2, &
                file_c1, file_c2, file_vals, n1, n2)
          END DO
        END DO

      CASE(c_bd_y_min, c_bd_y_max)
        DO j = 0, nz
          DO i = 0, nx
            pos1 = x(i)
            pos2 = z(j)
            laser%profile(i, j) = interp2d(pos1, pos2, &
                file_c1, file_c2, file_vals, n1, n2)
          END DO
        END DO

      CASE(c_bd_z_min, c_bd_z_max)
        DO j = 0, ny
          DO i = 0, nx
            pos1 = x(i)
            pos2 = y(j)
            laser%profile(i, j) = interp2d(pos1, pos2, &
                file_c1, file_c2, file_vals, n1, n2)
          END DO
        END DO

    END SELECT

    DEALLOCATE(file_c1, file_c2, file_vals)

  END SUBROUTINE custom_laser_spatial_setup

  ! Bilinear interpolation on a 2D regular grid.
  ! Returns the interpolated value at (p1, p2). Clamps to boundary values
  ! for points outside the data range.
  REAL(num) FUNCTION interp2d(p1, p2, c1, c2, vals, n1, n2)
    REAL(num), INTENT(IN) :: p1, p2
    INTEGER, INTENT(IN) :: n1, n2
    REAL(num), DIMENSION(n1), INTENT(IN) :: c1
    REAL(num), DIMENSION(n2), INTENT(IN) :: c2
    REAL(num), DIMENSION(n1, n2), INTENT(IN) :: vals

    INTEGER :: i1, i2
    REAL(num) :: u, v, q11, q21, q12, q22

    ! Direct index calculation (O(1)) — valid for uniform grids
    i1 = INT((p1 - c1(1)) / (c1(2) - c1(1))) + 1
    i2 = INT((p2 - c2(1)) / (c2(2) - c2(1))) + 1

    ! Clamp to valid interpolation range [1, n-1]
    i1 = MAX(1, MIN(i1, n1 - 1))
    i2 = MAX(1, MIN(i2, n2 - 1))

    u = (p1 - c1(i1)) / (c1(i1+1) - c1(i1))
    v = (p2 - c2(i2)) / (c2(i2+1) - c2(i2))

    ! Clamp fractional positions for points outside the data range
    u = MAX(0.0_num, MIN(u, 1.0_num))
    v = MAX(0.0_num, MIN(v, 1.0_num))

    q11 = vals(i1,   i2)
    q21 = vals(i1+1, i2)
    q12 = vals(i1,   i2+1)
    q22 = vals(i1+1, i2+1)

    interp2d = (1.0_num - u) * (1.0_num - v) * q11 &
             + u * (1.0_num - v) * q21             &
             + (1.0_num - u) * v * q12             &
             + u * v * q22

  END FUNCTION interp2d

  ! Placeholder for 3D spatiotemporal profile injection.
  ! Not yet implemented — the data format for E(y, z, t) is TBD.
  ! Currently prints a warning and returns 1.0 (flat profile).
  REAL(num) FUNCTION custom_laser_profile_3d(laser, pos1, pos2)
    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num), INTENT(IN) :: pos1, pos2

    custom_laser_profile_3d = 1.0_num

  END FUNCTION custom_laser_profile_3d

END MODULE custom_laser
