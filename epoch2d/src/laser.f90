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


! Some comments from Gemini:

!!!!! The Big Picture: How laser.f90 is Organised

!!! The laser module acts as an object-oriented class (within Fortran limits)
!!! for managing lasers.
! It uses a custom derived data type called laser_block (defined elsewhere,
! likely in a core data module)
! which acts like a struct containing all the properties of a single laser
! (amplitude, frequency, temporal profile,
! spatial profile). The subroutines in this file fall into four main categories:

!!!A. Initialisation and Memory Management

! init_laser: Sets up a new laser_block with safe default values (e.g.,
! -1.0_num)
! and allocates the spatial profile and phase arrays.

! deallocate_laser(s): Cleans up memory when the simulation ends.

! attach_laser: EPOCH allows multiple lasers on the same boundary.
! This routine links a newly created laser into a "linked list" for a specific
! boundary,
! allowing the code to loop through all active lasers later.

! allocate_with_boundary: Allocates the arrays for spatial profiles. Notice the
! 1-ng : ny+ng syntax.
! ng stands for "number of ghost cells". EPOCH needs these arrays to cover not
! just the physical domain, but the overlapping parallel boundary cells too.

!!! B. The Math Evaluators

! populate_pack_from_laser & laser_update_{phase, profile, omega}: EPOCH has a
! built-in mathematical parser.
! When you write profile = gauss(y, 0, 10e-6) in your input.deck, these routines
! translate that text
! into numerical values across the grid at every time step using
! evaluate_with_parameters.

! laser_time_profile: Calculates the temporal envelope (the time-varying
! amplitude of the laser pulse)
! at the current simulation time.

!!! C. Diagnostics
! calc_absorption: Computes the Poynting flux (electromagnetic energy transfer)
! across the boundaries
! to track how much laser energy was injected and how much scattered light left
! the domain.

!!! D. The Core Physics (The Boundary Conditions)

! outflow_bcs_x_min, outflow_bcs_x_max, outflow_bcs_y_min, outflow_bcs_y_max:
! These are the most critical
! subroutines for your project. They perform the actual physics of injecting the
! laser wave into the simulation box.




MODULE laser

  USE custom_laser
  USE evaluator

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_laser(boundary, laser)

    INTEGER, INTENT(IN) :: boundary
    TYPE(laser_block), INTENT(INOUT) :: laser

    laser%boundary = boundary
    laser%id = -1
    laser%use_time_function = .FALSE.
    laser%use_phase_function = .FALSE.
    laser%use_profile_function = .FALSE.
    laser%use_omega_function = .FALSE.
    laser%amp = -1.0_num
    laser%omega = -1.0_num
    laser%pol_angle = 0.0_num
    laser%t_start = 0.0_num
    laser%t_end = t_end
    laser%current_integral_phase = 0.0_num
    NULLIFY(laser%profile)
    NULLIFY(laser%phase)
    NULLIFY(laser%next)

    CALL allocate_with_boundary(laser%profile, boundary)
    CALL allocate_with_boundary(laser%phase, boundary)
    laser%profile = 1.0_num
    laser%phase = 0.0_num

  END SUBROUTINE init_laser



  SUBROUTINE setup_laser_phases(phases)

    REAL(num), DIMENSION(:), INTENT(IN) :: phases
    TYPE(laser_block), POINTER :: laser
    INTEGER :: ilas

    ilas = 1
    laser => lasers
    DO WHILE(ASSOCIATED(laser))
      laser%current_integral_phase = phases(ilas)
      ilas = ilas + 1
      laser => laser%next
    END DO

  END SUBROUTINE setup_laser_phases


! Deallocate all memory at the end of the simulation for a single laser, and for
! the linked list of lasers on each boundary.
  SUBROUTINE deallocate_laser(laser)

    TYPE(laser_block), POINTER :: laser

    IF (ASSOCIATED(laser%profile)) DEALLOCATE(laser%profile)
    IF (ASSOCIATED(laser%phase)) DEALLOCATE(laser%phase)
    IF (laser%use_profile_function) &
        CALL deallocate_stack(laser%profile_function)
    IF (laser%use_phase_function) &
        CALL deallocate_stack(laser%phase_function)
    IF (laser%use_time_function) &
        CALL deallocate_stack(laser%time_function)
    IF (laser%use_omega_function) &
        CALL deallocate_stack(laser%omega_function)
    DEALLOCATE(laser)

  END SUBROUTINE deallocate_laser


! Deallocate the linked list of lasers on each boundary at the end of the
! simulation. This loops through the linked list and calls deallocate_laser for
! each one.
  SUBROUTINE deallocate_lasers

    TYPE(laser_block), POINTER :: current, next

    current => lasers
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_laser(current)
      current => next
    END DO

  END SUBROUTINE deallocate_lasers


! Subroutine to attach a created laser object to the correct boundary
  SUBROUTINE attach_laser(laser)

    TYPE(laser_block), POINTER :: laser
    TYPE(laser_block), POINTER :: current
    INTEGER :: boundary

    boundary = laser%boundary

    n_lasers(boundary) = n_lasers(boundary) + 1

    IF (ASSOCIATED(lasers)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => laser
    ELSE
      lasers => laser
    END IF

    CALL custom_laser_spatial_setup(laser)

  END SUBROUTINE attach_laser



  ! This routine populates the constant elements of a parameter pack
  ! from a laser

  SUBROUTINE populate_pack_from_laser(laser, parameters)

    TYPE(laser_block), POINTER :: laser
    TYPE(parameter_pack), INTENT(INOUT) :: parameters

    parameters%pack_ix = 0
    parameters%pack_iy = 0

    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min)
        parameters%pack_ix = 0
      CASE(c_bd_x_max)
        parameters%pack_ix = nx
      CASE(c_bd_y_min)
        parameters%pack_iy = 0
      CASE(c_bd_y_max)
        parameters%pack_iy = ny
    END SELECT

  END SUBROUTINE populate_pack_from_laser


  ! There is a potential issue of overwriting the time_profile.
  ! If, in addition to using a temporal_spatial customised profile, one also
  ! specifies a t_profile in the laser block, then the overall laser profile
  ! will
  ! be the product of the two. This is not necessarily a problem, but it is
  ! something to be aware of. The logic here is that if a t_profile is
  ! specified, it will be used, otherwise the custom profile will be used.
  ! Therefore, if one wants to use a custom temporal profile, one should not
  ! specify a t_profile in the laser block. This is something that could be
  ! improved in future versions of EPOCH.
  FUNCTION laser_time_profile(laser)

    TYPE(laser_block), POINTER :: laser
    REAL(num) :: laser_time_profile
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    ! use_time_function is set to TRUE if there is a line specifying
    ! t_profile=... in laser block
    IF (laser%use_time_function) THEN
      laser_time_profile = evaluate_with_parameters(laser%time_function, &
          parameters, err)
      RETURN
    END IF

    ! There might be an overwriting issue about time_profile. Need to
    ! recheck the logic.
    laser_time_profile = custom_laser_time_profile(laser)

  END FUNCTION laser_time_profile



  SUBROUTINE laser_update_phase(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: i, err
    TYPE(parameter_pack) :: parameters

    !!! ADD THIS DECLARATION for custom (file-based) laser phase
    REAL(num) :: pos

    err = 0
    CALL populate_pack_from_laser(laser, parameters)

    ! Mirrors laser_update_profile: if the phase comes from a file, interpolate
    ! it from the phase data file (custom_laser_phase) at every grid point;
    ! otherwise evaluate the deck 'phase = ...' expression as before. When
    ! use_phase_from_file is set, any deck 'phase' expression is ignored.
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min, c_bd_x_max)
        DO i = 0,ny
          IF (laser%use_phase_from_file) THEN
            ! Use y(i) (cell centre) to match the analytical evaluator and the
            ! amplitude profile, which resolve the spatial coordinate at y(i).
            pos = y(i)
            laser%phase(i) = custom_laser_phase(laser, pos)
          ELSE
            parameters%pack_iy = i
            laser%phase(i) = &
                evaluate_with_parameters(laser%phase_function, parameters, err)
          END IF
        END DO
      CASE(c_bd_y_min, c_bd_y_max)
        DO i = 0,nx
          IF (laser%use_phase_from_file) THEN
            pos = x(i)
            laser%phase(i) = custom_laser_phase(laser, pos)
          ELSE
            parameters%pack_ix = i
            laser%phase(i) = &
                evaluate_with_parameters(laser%phase_function, parameters, err)
          END IF
        END DO
    END SELECT

  END SUBROUTINE laser_update_phase


  SUBROUTINE laser_update_profile(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: i, err
    TYPE(parameter_pack) :: parameters

    !!!  ADD THIS DECLARATION for custom laser profile
    REAL(num) :: pos

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    SELECT CASE(laser%boundary)
      !!! Changed the logic here to allow for custom laser profiles. If a custom
      !!! profile is specified, it will be used instead of the profile_function.
      !!! This allows for more flexibility in defining the spatial profile of
      !!! the laser.

      ! Note that use_profile_function is set to be FALSE in default (in
      ! laser.f90).
      ! it is set to TRUE in deck_laser_block.f90 if the profile is
      ! time-varying.

      ! So, here the logic is:
      ! If we use an analytical time-independent profile, then
      ! use_profile_function is set to FALSE and we do not need to update the
      ! profile at every time step (it is set at the start of the simulation)
      ! If we use an analytical time-varying profile, then use_profile_function
      ! is set to TRUE and we call evaluate_with_parameters to update the
      ! profile at every time step
      ! If we use a custom spatiotemporal profile, then use_custom_profile and
      ! use_spatiotemporal are both set to TRUE and we call custom_laser_profile
      ! to update the profile at every time step

      ! If we use a custom spatial profile, then use_custom_profile is set to
      ! TRUE and use_spatiotemporal is set to FALSE, and we do not need to
      ! update the profile at every time step (it is set at the start of the
      ! simulation)


      ! Fix proposed by Issue 1: rewrite the condition here, adding
      ! laser_profile_function%init and changing the sequence.
      CASE(c_bd_x_min, c_bd_x_max)
        DO i = 0,ny
          IF (laser%use_custom_profile .AND. laser%use_spatiotemporal) THEN
            ! Use y(i) (cell centre) to match the analytical evaluator,
            ! which resolves the deck variable 'y' at y(pack_iy).
            pos = y(i)
            laser%profile(i) = custom_laser_profile(laser, pos)
          ELSE IF (laser%use_profile_function &
              .OR. laser%profile_function%init) THEN
            parameters%pack_iy = i
            laser%profile(i) = evaluate_with_parameters( &
                laser%profile_function, parameters, err)
          END IF
        END DO

      CASE(c_bd_y_min, c_bd_y_max)
        DO i = 0,nx
          IF (laser%use_custom_profile .AND. laser%use_spatiotemporal) THEN
            pos = x(i)
            laser%profile(i) = custom_laser_profile(laser, pos)
          ELSE IF (laser%use_profile_function &
              .OR. laser%profile_function%init) THEN
            parameters%pack_ix = i
            laser%profile(i) = evaluate_with_parameters( &
                laser%profile_function, parameters, err)
          END IF
        END DO


    END SELECT

  END SUBROUTINE laser_update_profile



  SUBROUTINE laser_update_omega(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    laser%omega = &
        evaluate_with_parameters(laser%omega_function, parameters, err)
    IF (laser%omega_func_type == c_of_freq) &
        laser%omega = 2.0_num * pi * laser%omega
    IF (laser%omega_func_type == c_of_lambda) &
        laser%omega = 2.0_num * pi * c / laser%omega

  END SUBROUTINE laser_update_omega



  SUBROUTINE update_laser_omegas

    TYPE(laser_block), POINTER :: current

    current => lasers
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) THEN
        CALL laser_update_omega(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE
        current%current_integral_phase = current%omega * time
      END IF
      current => current%next
    END DO

  END SUBROUTINE update_laser_omegas



  SUBROUTINE allocate_with_boundary(array, boundary)

    REAL(num), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: boundary

    IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
      ALLOCATE(array(1-ng:ny+ng))
    ELSE IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
      ALLOCATE(array(1-ng:nx+ng))
    END IF

  END SUBROUTINE allocate_with_boundary



  SUBROUTINE set_laser_dt

    REAL(num) :: dt_local
    TYPE(laser_block), POINTER :: current

    dt_laser = HUGE(1.0_num)

    current => lasers
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    END DO

    ! Need at least two iterations per laser period
    ! (Nyquist)
    dt_laser = dt_laser / 2.0_num

  END SUBROUTINE set_laser_dt


  ! In a normal running: outflow_bcs_x_min is part of the final B-boundary
  ! update once per main timestep, plus one startup call during initial field
  ! setup.
  SUBROUTINE outflow_bcs_x_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_min

    laserpos = 1
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_x_min_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:ny))
    ALLOCATE(source2(0:ny))
    source1 = 0.0_num
    source2 = 0.0_num

    bx(laserpos-1, 0:ny) = bx_x_min(0:ny)

    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_x_min) THEN
          ! evaluate the temporal evolution of the laser
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            IF (current%use_phase_function .OR. current%use_phase_from_file) &
                CALL laser_update_phase(current)

            ! ---> TRIGGER FOR BOTH MATH STRINGS AND OUR 2D FILE <---
            !!! Here, the logic is: if we use the normal way to declare a laser
            !!! profile AND we don't declare a t_profile,
            !!! then the first condition is met (I think that if we use a
            !!! time-independent laser profile, EPOCH initialises
            !!! it at the start and then TURN use_profile_function to FALSE, in
            !!! which case we don't need and it doesn't use
            !!! update_laser_profile)
            !!! If we use spatiotemporal profile (in which case both conditions
            !!! in the second line are met), then we need to call
            !!! laser_update_profile
            !!! If we only use a spatial profile, then there is no need to call
            !!! laser_update_profile as the profile is initiated at the start.

            ! For reference, in laser_block_handle_element function in
            ! deck_laser_block.f90, we have the following:
            ! One can see that only if the profile is time varying, then
            ! use_profile_function is declared TRUE.
            !IF (working_laser%profile_function%is_time_varying) THEN
            !working_laser%use_profile_function = .TRUE.
            !ELSE
            !  CALL deallocate_stack(working_laser%profile_function)

            ! Based on this, I think this is the correct way.
            IF (current%use_profile_function .OR. &
                (current%use_custom_profile .AND. current%use_spatiotemporal)) &
                CALL laser_update_profile(current)

            ! SO this line multiplies the profile set by t_profile and the
            ! profile set by either profile = ... or the customised dat file.
            t_env = laser_time_profile(current) * current%amp
            DO i = 0,ny
              base = t_env * current%profile(i) &
                * SIN(current%current_integral_phase + current%phase(i))
              source1(i) = source1(i) + base * COS(current%pol_angle)
              source2(i) = source2(i) + base * SIN(current%pol_angle)
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    bz(laserpos-1, 0:ny) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ey_x_min(0:ny) + c * bz_x_min(0:ny)) &
        - 2.0_num * ey(laserpos, 0:ny) &
        + dt_eps * jy(laserpos, 0:ny) &
        + diff * bz(laserpos, 0:ny))

    by(laserpos-1, 0:ny) = sum * (-4.0_num * source2 &
        - 2.0_num * (ez_x_min(0:ny) - c * by_x_min(0:ny)) &
        + 2.0_num * ez(laserpos, 0:ny) &
        - ly * (bx(laserpos, 0:ny) - bx(laserpos, -1:ny-1)) &
        - dt_eps * jz(laserpos, 0:ny) &
        + diff * by(laserpos, 0:ny))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_x_min, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_x_min)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_x_min



  SUBROUTINE outflow_bcs_x_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_max

    laserpos = nx
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_x_max_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:ny))
    ALLOCATE(source2(0:ny))
    source1 = 0.0_num
    source2 = 0.0_num

    bx(laserpos+1, 0:ny) = bx_x_max(0:ny)

    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_x_max) THEN
          ! evaluate the temporal evolution of the laser
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            IF (current%use_phase_function .OR. current%use_phase_from_file) &
                CALL laser_update_phase(current)
            IF (current%use_profile_function .OR. &
                (current%use_custom_profile .AND. current%use_spatiotemporal)) &
                CALL laser_update_profile(current)
            t_env = laser_time_profile(current) * current%amp
            DO i = 0,ny
              base = t_env * current%profile(i) &
                * SIN(current%current_integral_phase + current%phase(i))
              source1(i) = source1(i) + base * COS(current%pol_angle)
              source2(i) = source2(i) + base * SIN(current%pol_angle)
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    bz(laserpos, 0:ny) = sum * (-4.0_num * source1 &
        - 2.0_num * (ey_x_max(0:ny) - c * bz_x_max(0:ny)) &
        + 2.0_num * ey(laserpos, 0:ny) &
        - dt_eps * jy(laserpos, 0:ny) &
        + diff * bz(laserpos-1, 0:ny))

    by(laserpos, 0:ny) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ez_x_max(0:ny) + c * by_x_max(0:ny)) &
        - 2.0_num * ez(laserpos, 0:ny) &
        + ly * (bx(laserpos, 0:ny) - bx(laserpos, -1:ny-1)) &
        + dt_eps * jz(laserpos, 0:ny) &
        + diff * by(laserpos-1, 0:ny))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_x_max, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_x_max)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_x_max



  SUBROUTINE outflow_bcs_y_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_y_min

    laserpos = 1
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_y_min_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:nx))
    ALLOCATE(source2(0:nx))
    source1 = 0.0_num
    source2 = 0.0_num

    by(0:nx, laserpos-1) = by_y_min(0:nx)

    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_y_min) THEN
          ! evaluate the temporal evolution of the laser
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            IF (current%use_phase_function .OR. current%use_phase_from_file) &
                CALL laser_update_phase(current)
            IF (current%use_profile_function .OR. &
                (current%use_custom_profile .AND. current%use_spatiotemporal)) &
                CALL laser_update_profile(current)
            t_env = laser_time_profile(current) * current%amp
            DO i = 0,nx
              base = t_env * current%profile(i) &
                * SIN(current%current_integral_phase + current%phase(i))
              source1(i) = source1(i) + base * COS(current%pol_angle)
              source2(i) = source2(i) + base * SIN(current%pol_angle)
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    bx(0:nx, laserpos-1) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ez_y_min(0:nx) + c * bx_y_min(0:nx)) &
        - 2.0_num * ez(0:nx, laserpos) &
        - lx * (by(0:nx, laserpos) - by(-1:nx-1, laserpos)) &
        + dt_eps * jz(0:nx, laserpos) &
        + diff * bx(0:nx, laserpos))

    bz(0:nx, laserpos-1) = sum * (-4.0_num * source2 &
        - 2.0_num * (ex_y_min(0:nx) - c * bz_y_min(0:nx)) &
        + 2.0_num * ex(0:nx, laserpos) &
        - dt_eps * jx(0:nx, laserpos) &
        + diff * bz(0:nx, laserpos))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_y_min, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_y_min)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_y_min



  SUBROUTINE outflow_bcs_y_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_y_max

    laserpos = ny
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_y_max_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:nx))
    ALLOCATE(source2(0:nx))
    source1 = 0.0_num
    source2 = 0.0_num

    by(0:nx, laserpos+1) = by_y_max(0:nx)

    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_y_max) THEN
          ! evaluate the temporal evolution of the laser
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            IF (current%use_phase_function .OR. current%use_phase_from_file) &
                CALL laser_update_phase(current)
            IF (current%use_profile_function .OR. &
                (current%use_custom_profile .AND. current%use_spatiotemporal)) &
                CALL laser_update_profile(current)
            t_env = laser_time_profile(current) * current%amp
            DO i = 0,nx
              base = t_env * current%profile(i) &
                * SIN(current%current_integral_phase + current%phase(i))
              source1(i) = source1(i) + base * COS(current%pol_angle)
              source2(i) = source2(i) + base * SIN(current%pol_angle)
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    bx(0:nx, laserpos) = sum * (-4.0_num * source1 &
        - 2.0_num * (ez_y_max(0:nx) - c * bx_y_max(0:nx)) &
        + 2.0_num * ez(0:nx, laserpos) &
        + lx * (by(0:nx, laserpos) - by(-1:nx-1, laserpos)) &
        - dt_eps * jz(0:nx, laserpos) &
        + diff * bx(0:nx, laserpos-1))

    bz(0:nx, laserpos) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ex_y_max(0:nx) + c * bz_y_max(0:nx)) &
        - 2.0_num * ex(0:nx, laserpos) &
        + dt_eps * jx(0:nx, laserpos) &
        + diff * bz(0:nx, laserpos-1))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_y_max, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_y_max)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_y_max



  SUBROUTINE calc_absorption(bd, lasers)

    TYPE(laser_block), POINTER, OPTIONAL :: lasers
    INTEGER, INTENT(IN) :: bd
    TYPE(laser_block), POINTER :: current
    REAL(num) :: t_env, dir, dd, factor, lfactor, laser_inject_sum
    REAL(num), DIMENSION(:), ALLOCATABLE :: e1, e2, b1, b2
    INTEGER :: mm, ibc, icell

    ! Note: ideally e1, e2, b1, b2 should be face-centred. However, this is not
    ! possible with 'open' boundaries since E-fields are not defined in the
    ! ghost cell, so we use the cell-centred quantities in the first cell.

    dir = 1.0_num
    mm = 1

    SELECT CASE(bd)
      CASE(c_bd_x_min, c_bd_x_max)
        dd = dy
        mm = ny
        ALLOCATE(e1(mm), e2(mm), b1(mm), b2(mm))

        ibc = 1
        IF (bd == c_bd_x_max) THEN
          dir = -1.0_num
          ibc = nx
        END IF

        e1 = 0.5_num  * (ey(ibc  , 0:ny-1) + ey(ibc, 1:ny  ))
        e2 = ez(ibc, 1:ny)
        b1 = 0.25_num * (bz(ibc-1, 0:ny-1) + bz(ibc, 0:ny-1) &
                       + bz(ibc-1, 1:ny  ) + bz(ibc, 1:ny  ))
        b2 = 0.5_num  * (by(ibc-1, 1:ny  ) + by(ibc, 1:ny  ))

      CASE(c_bd_y_min, c_bd_y_max)
        dd = dx
        mm = nx
        ALLOCATE(e1(mm), e2(mm), b1(mm), b2(mm))

        ibc = 1
        IF (bd == c_bd_y_max) THEN
          dir = -1.0_num
          ibc = ny
        END IF

        e1 = ez(1:nx, ibc)
        e2 = 0.5_num  * (ex(0:nx-1, ibc  ) + ex(1:nx  , ibc))
        b1 = 0.5_num  * (bx(1:nx  , ibc-1) + bx(1:nx  , ibc))
        b2 = 0.25_num * (bz(0:nx-1, ibc-1) + bz(0:nx-1, ibc) &
                       + bz(1:nx  , ibc-1) + bz(1:nx  , ibc))

      CASE DEFAULT
        dd = 0.0_num
        ALLOCATE(e1(mm), e2(mm), b1(mm), b2(mm))

        e1 = 0.0_num
        e2 = 0.0_num
        b1 = 0.0_num
        b2 = 0.0_num
    END SELECT

    factor = dt * dd * dir
    laser_absorb_local = laser_absorb_local &
        + (factor / mu0) * SUM(e1 * b1 - e2 * b2)

    IF (PRESENT(lasers)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == bd) THEN
          laser_inject_sum = 0.0_num
          DO icell = 1, mm
            laser_inject_sum = laser_inject_sum + current%profile(icell)**2
          END DO
          t_env = laser_time_profile(current)
          lfactor = 0.5_num * epsilon0 * c * ABS(factor) &
              * (t_env * current%amp)**2
          laser_inject_local = laser_inject_local + lfactor * laser_inject_sum
        END IF
        current => current%next
      END DO
    END IF

    DEALLOCATE(e1, e2, b1, b2)

  END SUBROUTINE calc_absorption

END MODULE laser
