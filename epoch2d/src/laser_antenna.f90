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

! Interior laser antenna: a Huygens/total-field-scattered-field (TFSF)
! injection plane at interior x = x0, as an alternative to the boundary-
! based custom laser pipeline (laser.f90/custom_laser.f90) for cases that
! need a genuinely vectorial incident field (independent tangential E and
! B, including a real longitudinal field developing self-consistently)
! rather than one derived from the other at a domain edge.
!
! Deliberately NOT built on laser_block/laser.f90's boundary-characteristic
! mechanism -- see the antenna_block comment in shared_data.F90 for why.
! Deliberately does NOT "USE fields": fields.f90 calls into this module
! (apply_antenna_e_correction/apply_antenna_b_correction), so a reverse
! dependency would be circular. The half-step coefficients (hdtx, cnx) are
! cheap to recompute locally from dt/dx/c rather than importing fields'
! cached copies.
MODULE laser_antenna

  USE custom_laser

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_antenna(direction, antenna)

    INTEGER, INTENT(IN) :: direction
    TYPE(antenna_block), INTENT(INOUT) :: antenna

    antenna%direction = direction
    antenna%x0 = 0.0_num
    antenna%amp = 1.0_num
    antenna%t_start = 0.0_num
    antenna%t_end = t_end
    antenna%fields_loaded = .FALSE.
    NULLIFY(antenna%ey_inc_matrix)
    NULLIFY(antenna%ez_inc_matrix)
    NULLIFY(antenna%by_inc_matrix)
    NULLIFY(antenna%bz_inc_matrix)
    NULLIFY(antenna%next)

  END SUBROUTINE init_antenna



  ! Snap the deck-declared x0 to the nearest B-grid seam. By/bz at GLOBAL
  ! index k sit at physical position x_min + k*dx: this follows directly
  ! from update_b_field's curl stencil (by(k,:) is a centred difference of
  ! ez between E-grid indices k and k+1, i.e. physically at the midpoint
  ! x_global(k) + dx/2) combined with setup_grid's x_global(k) = x_min +
  ! dx/2 + (k-1)*dx -- confirmed directly against src/housekeeping/
  ! setup.F90's setup_grid, not assumed from the xb/xb_global arrays (whose
  ! own index convention turned out to be offset from this by a full cell
  ! and would have given a silently-wrong snap).
  SUBROUTINE snap_antenna_plane(antenna)

    TYPE(antenna_block), INTENT(INOUT) :: antenna
    INTEGER :: iseam_global

    iseam_global = NINT((antenna%x0 - x_min) / dx)
    antenna%x0_actual = x_min + REAL(iseam_global, num) * dx
    antenna%i0_global = iseam_global + 1

    IF (rank == 0) THEN
      PRINT *, '>>> Laser antenna plane snapped to x = ', &
          antenna%x0_actual, ' (declared x0 = ', antenna%x0, ') <<<'
    END IF

  END SUBROUTINE snap_antenna_plane



  ! Load the four tangential incident-field files. Mirrors
  ! custom_laser_spatial_setup's grid-declaration check and reuses
  ! load_spatiotemporal_matrix (custom_laser.f90) for the actual read --
  ! same raw-binary/stream convention as the boundary-laser profile/phase
  ! files. Every rank loads the full files (epoch2d has no per-rank
  ! windowing for the boundary laser either; see HANDOFF.md for the
  ! epoch3d follow-on that would need it).
  SUBROUTINE antenna_spatial_setup(antenna)

    TYPE(antenna_block), INTENT(INOUT) :: antenna
    LOGICAL :: ok
    INTEGER :: mpi_err

    IF (antenna%fields_loaded) RETURN

    ok = antenna%n_transverse_points >= 2 &
        .AND. antenna%profile_transverse_max > antenna%profile_transverse_min &
        .AND. antenna%n_t_points >= 2 .AND. antenna%t_end > antenna%t_start &
        .AND. LEN_TRIM(antenna%ey_inc_file) > 0 &
        .AND. LEN_TRIM(antenna%ez_inc_file) > 0 &
        .AND. LEN_TRIM(antenna%by_inc_file) > 0 &
        .AND. LEN_TRIM(antenna%bz_inc_file) > 0
    IF (.NOT. ok) THEN
      IF (rank == 0) THEN
        PRINT *, 'ERROR: begin:laser_antenna requires n_t_points and ', &
            'n_transverse_points (each >= 2), profile_transverse_min < ', &
            'profile_transverse_max, t_start < t_end, and all four of ', &
            'ey_inc_file/ez_inc_file/by_inc_file/bz_inc_file to be set.'
      END IF
      CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
    END IF

    CALL load_spatiotemporal_matrix(antenna%ey_inc_file, &
        antenna%n_transverse_points, antenna%n_t_points, &
        antenna%ey_inc_matrix)
    CALL load_spatiotemporal_matrix(antenna%ez_inc_file, &
        antenna%n_transverse_points, antenna%n_t_points, &
        antenna%ez_inc_matrix)
    CALL load_spatiotemporal_matrix(antenna%by_inc_file, &
        antenna%n_transverse_points, antenna%n_t_points, &
        antenna%by_inc_matrix)
    CALL load_spatiotemporal_matrix(antenna%bz_inc_file, &
        antenna%n_transverse_points, antenna%n_t_points, &
        antenna%bz_inc_matrix)

    antenna%fields_loaded = .TRUE.

    IF (rank == 0) THEN
      PRINT *, '>>> Laser antenna incident fields loaded successfully! <<<'
      PRINT *, '    Grid size: ', antenna%n_transverse_points, &
          ' (transverse) x ', antenna%n_t_points, ' (temporal)'
    END IF

  END SUBROUTINE antenna_spatial_setup



  SUBROUTINE attach_antenna(antenna)

    TYPE(antenna_block), POINTER :: antenna
    TYPE(antenna_block), POINTER :: current

    n_antennas = n_antennas + 1

    IF (ASSOCIATED(antennas)) THEN
      current => antennas
      DO WHILE (ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => antenna
    ELSE
      antennas => antenna
    END IF

    CALL snap_antenna_plane(antenna)
    CALL antenna_spatial_setup(antenna)

  END SUBROUTINE attach_antenna



  SUBROUTINE deallocate_antenna(antenna)

    TYPE(antenna_block), POINTER :: antenna

    IF (ASSOCIATED(antenna%ey_inc_matrix)) DEALLOCATE(antenna%ey_inc_matrix)
    IF (ASSOCIATED(antenna%ez_inc_matrix)) DEALLOCATE(antenna%ez_inc_matrix)
    IF (ASSOCIATED(antenna%by_inc_matrix)) DEALLOCATE(antenna%by_inc_matrix)
    IF (ASSOCIATED(antenna%bz_inc_matrix)) DEALLOCATE(antenna%bz_inc_matrix)
    DEALLOCATE(antenna)

  END SUBROUTINE deallocate_antenna



  SUBROUTINE deallocate_antennas

    TYPE(antenna_block), POINTER :: current, next

    current => antennas
    DO WHILE (ASSOCIATED(current))
      next => current%next
      CALL deallocate_antenna(current)
      current => next
    END DO

  END SUBROUTINE deallocate_antennas



  ! Sample one of the four incident-field matrices at transverse position
  ! y(iy) and time t_sample, scaled by the antenna's amp multiplier.
  REAL(num) FUNCTION sample_inc(antenna, matrix, pos, t_sample)

    TYPE(antenna_block), INTENT(IN) :: antenna
    REAL(num), DIMENSION(:,:), INTENT(IN) :: matrix
    REAL(num), INTENT(IN) :: pos, t_sample

    sample_inc = antenna%amp * sample_spatiotemporal_matrix(matrix, pos, &
        t_sample, antenna%profile_transverse_min, &
        antenna%profile_transverse_max, antenna%n_transverse_points, &
        antenna%t_start, antenna%t_end, antenna%n_t_points)

  END FUNCTION sample_inc



  ! TFSF correction to the B-field update, applied immediately after every
  ! CALL update_b_field in fields.f90 (both update_eb_fields_half and
  ! update_eb_fields_final). t_sample must be the physical time B was just
  ! advanced to -- NOT necessarily fields.f90's ambient 'time' variable,
  ! which is stale by a full half-step at the second call site of a
  ! timestep (see the antenna_t_mid discussion in fields.f90).
  !
  ! Always applied at the single seam B-grid index i0_global-1, regardless
  ! of direction (there is only one seam). The sign flips with direction:
  ! for c_bd_x_max (total field for x > x0) the seam B-value is designated
  ! scattered and the correction removes the total-side E leaking into it;
  ! for c_bd_x_min (total field for x < x0) the same seam B-value is still
  ! designated scattered, but it is now the WEST side that is total, so the
  ! sign of both terms flips. See shared_data.F90's antenna_block comment
  ! and the design report for the full derivation -- this was checked by
  ! hand for both directions, not assumed symmetric.
  SUBROUTINE apply_antenna_b_correction(t_sample)

    REAL(num), INTENT(IN) :: t_sample
    TYPE(antenna_block), POINTER :: current
    INTEGER :: iseam_global, ix_local, iy
    REAL(num) :: hdtx, sgn, ey_inc, ez_inc

    IF (.NOT. ASSOCIATED(antennas)) RETURN

    hdtx = 0.5_num * dt / dx

    current => antennas
    DO WHILE (ASSOCIATED(current))

      iseam_global = current%i0_global - 1
      ix_local = iseam_global - nx_global_min + 1

      IF (ix_local >= 0 .AND. ix_local <= nx) THEN

        sgn = 1.0_num
        IF (current%direction == c_bd_x_min) sgn = -1.0_num

        DO iy = 0, ny
          ey_inc = sample_inc(current, current%ey_inc_matrix, y(iy), &
              t_sample)
          ez_inc = sample_inc(current, current%ez_inc_matrix, y(iy), &
              t_sample)

          by(ix_local, iy) = by(ix_local, iy) - sgn * hdtx * ez_inc
          bz(ix_local, iy) = bz(ix_local, iy) + sgn * hdtx * ey_inc
        END DO

      END IF

      current => current%next
    END DO

  END SUBROUTINE apply_antenna_b_correction



  ! TFSF correction to the E-field update, applied immediately after every
  ! CALL update_e_field in fields.f90. t_sample is the ambient 'time' at
  ! the call site (correct at both call sites, unlike the B-correction).
  !
  ! Unlike the B-correction, the target E-grid index DEPENDS on direction:
  ! c_bd_x_max corrects the first E-grid point strictly east of the seam
  ! (i0_global); c_bd_x_min corrects the seam's own west-adjacent E-grid
  ! point (i0_global-1) instead, because the total/scattered designation of
  ! that point flips with direction while the seam B-value's own
  ! designation (always "scattered") does not. Re-derived by hand from the
  ! base update_e_field stencil for both directions -- do not "simplify" to
  ! a same-index sign flip, that was checked and found wrong.
  SUBROUTINE apply_antenna_e_correction(t_sample)

    REAL(num), INTENT(IN) :: t_sample
    TYPE(antenna_block), POINTER :: current
    INTEGER :: ie_global, ix_local, iy
    REAL(num) :: cnx, sgn, by_inc, bz_inc

    IF (.NOT. ASSOCIATED(antennas)) RETURN

    cnx = 0.5_num * dt / dx * c**2

    current => antennas
    DO WHILE (ASSOCIATED(current))

      sgn = 1.0_num
      ie_global = current%i0_global
      IF (current%direction == c_bd_x_min) THEN
        sgn = -1.0_num
        ie_global = current%i0_global - 1
      END IF

      ix_local = ie_global - nx_global_min + 1

      IF (ix_local >= 0 .AND. ix_local <= nx) THEN

        DO iy = 0, ny
          by_inc = sample_inc(current, current%by_inc_matrix, y(iy), &
              t_sample)
          bz_inc = sample_inc(current, current%bz_inc_matrix, y(iy), &
              t_sample)

          ey(ix_local, iy) = ey(ix_local, iy) + sgn * cnx * bz_inc
          ez(ix_local, iy) = ez(ix_local, iy) - sgn * cnx * by_inc
        END DO

      END IF

      current => current%next
    END DO

  END SUBROUTINE apply_antenna_e_correction

END MODULE laser_antenna
