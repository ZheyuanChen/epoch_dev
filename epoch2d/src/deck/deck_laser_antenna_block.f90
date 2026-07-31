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

! Deck parser for begin:laser_antenna -- an interior Huygens/TFSF
! injection plane at x = x0, distinct from begin:laser (which drives a
! domain-boundary characteristic condition and derives B from E). See
! src/laser_antenna.f90 and the antenna_block comment in shared_data.F90.
MODULE deck_laser_antenna_block

  USE strings_advanced
  USE laser_antenna
  USE fields
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: laser_antenna_deck_initialise, laser_antenna_deck_finalise
  PUBLIC :: laser_antenna_block_start, laser_antenna_block_end
  PUBLIC :: laser_antenna_block_handle_element, laser_antenna_block_check

  TYPE(antenna_block), POINTER :: working_antenna
  LOGICAL :: direction_set = .FALSE.

CONTAINS

  SUBROUTINE laser_antenna_deck_initialise

    n_antennas = 0

  END SUBROUTINE laser_antenna_deck_initialise



  SUBROUTINE laser_antenna_deck_finalise

  END SUBROUTINE laser_antenna_deck_finalise



  SUBROUTINE laser_antenna_block_start

    IF (deck_state == c_ds_first) RETURN

    ALLOCATE(working_antenna)

  END SUBROUTINE laser_antenna_block_start



  SUBROUTINE laser_antenna_block_end

    IF (deck_state == c_ds_first) RETURN

    CALL attach_antenna(working_antenna)
    direction_set = .FALSE.

  END SUBROUTINE laser_antenna_block_end



  FUNCTION laser_antenna_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: direction, io, iu

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'direction')) THEN
      direction = as_boundary_print(value, element, errcode)
      IF (direction /= c_bd_x_min .AND. direction /= c_bd_x_max) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
            WRITE(io,*) '"direction" in block "laser_antenna" must be ', &
                '"x_min" or "x_max" (an interior plane normal to x).'
          END DO
        END IF
        errcode = c_err_bad_value
        RETURN
      END IF
      CALL init_antenna(direction, working_antenna)
      direction_set = .TRUE.
      RETURN
    END IF

    IF (.NOT. direction_set) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Cannot set laser_antenna properties before ', &
              '"direction" is set'
        END DO
      END IF
      extended_error_string = 'direction'
      errcode = c_err_required_element_not_set
      RETURN
    END IF

    IF (str_cmp(element, 'x0')) THEN
      working_antenna%x0 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'amp')) THEN
      working_antenna%amp = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_start')) THEN
      working_antenna%t_start = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_end')) THEN
      working_antenna%t_end = as_time_print(value, element, errcode)
      RETURN
    END IF

    ! Filenames for the four tangential incident-field channels -- raw
    ! binary, access='stream', same convention as the boundary laser's
    ! profile_data_file (see custom_laser.f90). Plain filenames resolve
    ! relative to data_dir; absolute paths are used as-is.
    IF (str_cmp(element, 'ey_inc_file')) THEN
      working_antenna%ey_inc_file = TRIM(ADJUSTL(value))
      RETURN
    END IF

    IF (str_cmp(element, 'ez_inc_file')) THEN
      working_antenna%ez_inc_file = TRIM(ADJUSTL(value))
      RETURN
    END IF

    IF (str_cmp(element, 'by_inc_file')) THEN
      working_antenna%by_inc_file = TRIM(ADJUSTL(value))
      RETURN
    END IF

    IF (str_cmp(element, 'bz_inc_file')) THEN
      working_antenna%bz_inc_file = TRIM(ADJUSTL(value))
      RETURN
    END IF

    ! Shape and bounds of the four incident-field files -- identical
    ! convention to begin:laser's spatiotemporal profile/phase grid
    ! declaration (see deck_laser_block.f90). The temporal extent reuses
    ! t_start/t_end above rather than a separate pair of elements.
    IF (str_cmp(element, 'n_t_points') .OR. str_cmp(element, 'n_t')) THEN
      working_antenna%n_t_points = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'n_transverse_points') &
        .OR. str_cmp(element, 'n_y')) THEN
      working_antenna%n_transverse_points = &
          as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'profile_transverse_min') &
        .OR. str_cmp(element, 'y_min')) THEN
      working_antenna%profile_transverse_min = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'profile_transverse_max') &
        .OR. str_cmp(element, 'y_max')) THEN
      working_antenna%profile_transverse_max = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION laser_antenna_block_handle_element



  FUNCTION laser_antenna_block_check() RESULT(errcode)

    INTEGER :: errcode
    TYPE(antenna_block), POINTER :: current
    INTEGER :: io, iu

    errcode = c_err_none

    current => antennas
    DO WHILE(ASSOCIATED(current))
      IF (current%t_end <= current%t_start) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Must have "t_end" > "t_start" for every ', &
                'laser_antenna.'
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF
      current => current%next
    END DO

    ! The TFSF correction (laser_antenna.f90) was derived from ONE
    ! specific branch of update_e_field/update_b_field: field_order = 2,
    ! maxwell_solver = yee, cpml_boundaries = F (the plain nearest-
    ! neighbour Yee curl, cx1 = hdtx/cnx with no alphax/betaxy/deltax
    ! cross terms). field_order = 4/6 reaches wider neighbours; the Lehe/
    ! Lehe_x/Lehe_y/Pukhov/Cowan solvers and the cpml_boundaries branch
    ! all use different stencil coefficients entirely (confirmed by
    ! re-reading fields.f90's update_e_field/update_b_field branching).
    ! None of these would crash if combined with begin:laser_antenna --
    ! they would silently under-cancel or misapply the correction, since
    ! it does not know about the extra terms. Require the exact
    ! combination the correction assumes, explicitly, whenever any
    ! laser_antenna is present. Errors returned via this function's
    ! errcode (like the t_start/t_end check above) are informational only
    ! -- check_compulsory_blocks' aggregate result is never upgraded to
    ! c_err_terminate the way per-element parse errors are (traced through
    ! handle_deck_element), so a real hard-stop needs a direct abort_code
    ! call here, matching control_deck_finalise's own pattern for its
    ! mandatory "nx"/"ny" checks.
    IF (ASSOCIATED(antennas) .AND. (field_order /= 2 &
        .OR. maxwell_solver /= c_maxwell_solver_yee &
        .OR. cpml_boundaries)) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'begin:laser_antenna requires "field_order = 2", ', &
              '"maxwell_solver = yee" (both defaults) and no CPML ', &
              'boundaries -- the TFSF correction was derived from that ', &
              'one specific Yee curl stencil and does not account for ', &
              'the extra terms any other combination introduces.'
        END DO
      END IF
      CALL abort_code(c_err_bad_value)
    END IF

  END FUNCTION laser_antenna_block_check

END MODULE deck_laser_antenna_block
