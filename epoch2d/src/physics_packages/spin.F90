! Copyright (C) 2009-2023 University of Warwick
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


! This for now is directly copied from Holger's repo.
! Need to update


MODULE spin

#ifdef SPIN

  USE constants
  USE shared_data
  USE random_generator
  USE utilities
  IMPLICIT NONE

  ! a(chi)/a(0) on a uniform log10(chi) grid (anomalous_moment.table);
  ! amm_table_loaded is set once the table has been read (use_qed only)
  LOGICAL :: amm_table_loaded = .FALSE.
  INTEGER :: n_amm
  REAL(num) :: amm_log_chi_min, amm_log_chi_max, amm_idlog
  REAL(num) :: amm_chi2_min, amm_tail_scale
  REAL(num), ALLOCATABLE :: amm_log_f(:)

CONTAINS


  SUBROUTINE setup_particle_spin(part_species)
    TYPE(particle_species), POINTER :: part_species
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart

    partlist => part_species%attached_list
    current => partlist%head
    ipart = 0

    DO WHILE(ipart < partlist%count)
      CALL init_particle_spin(part_species, current)
      current => current%next
      ipart = ipart + 1
    END DO

  END SUBROUTINE setup_particle_spin

  SUBROUTINE init_particle_spin(species, new_particle)
    TYPE(particle_species), INTENT(IN) :: species
    TYPE(particle), POINTER :: new_particle

    SELECT CASE(species%spin_distribution)
      CASE(c_spin_uniform)
        CALL init_particle_spin_uniform(species, new_particle)
      CASE(c_spin_directed)
        CALL init_particle_spin_directed(species, new_particle)
      CASE DEFAULT; new_particle%spin = 0.0_num
    END SELECT

  END SUBROUTINE init_particle_spin


  SUBROUTINE init_particle_spin_directed(species, new_particle)
    TYPE(particle_species), INTENT(IN) :: species
    TYPE(particle), POINTER :: new_particle

    new_particle%spin = species%spin_orientation
  END SUBROUTINE init_particle_spin_directed

  SUBROUTINE init_particle_spin_uniform(species, new_particle)
    TYPE(particle_species), INTENT(IN) :: species
    TYPE(particle), POINTER :: new_particle

    REAL(num) :: theta, phi
    REAL(num) :: sx, sy, sz

    phi = 2.0_num * pi * random()
    theta = acos(1.0_num - 2.0_num * random())
    sx = sin(theta) * cos(phi)
    sy = sin(theta) * sin(phi)
    sz = cos(theta)

    new_particle%spin = (/ sx, sy, sz /)
  END SUBROUTINE init_particle_spin_uniform
  


  ! Read chi-dependent anomalous_moment from 
  ! anomalous_moment.table
  ! rank 0 reads the table then broadcasts
  ! Note that bisection search is not needed because
  ! the table is generated on a uniform grid (log chi)
  ! This is called in photons.f90 so no need to do ifdef here
  SUBROUTINE setup_anomalous_moment_table

    INTEGER :: i
    REAL(num) :: log_chi, dlog, buf(2), y

    IF (rank == 0) THEN
      OPEN(unit=lu, file=TRIM(qed_table_location) &
          // '/anomalous_moment.table', status='OLD')
      READ(lu,*) n_amm, amm_log_chi_min, amm_log_chi_max
      ALLOCATE(amm_log_f(n_amm))
      dlog = (amm_log_chi_max - amm_log_chi_min) / (n_amm - 1)
      DO i = 1, n_amm
        READ(lu,*) log_chi, amm_log_f(i)
        ! The direct-index lookup assumes a uniform grid
        IF (ABS(log_chi - amm_log_chi_min - (i - 1) * dlog) > 1.0e-8_num) &
            CALL abort_code(c_err_bad_value)
      END DO
      CLOSE(unit=lu)
    END IF

    CALL MPI_BCAST(n_amm, 1, MPI_INTEGER, 0, comm, errcode)
    buf = (/ amm_log_chi_min, amm_log_chi_max /)
    CALL MPI_BCAST(buf, 2, mpireal, 0, comm, errcode)
    amm_log_chi_min = buf(1)
    amm_log_chi_max = buf(2)
    IF (rank /= 0) ALLOCATE(amm_log_f(n_amm))
    CALL MPI_BCAST(amm_log_f, n_amm, mpireal, 0, comm, errcode)

    amm_idlog = (n_amm - 1) / (amm_log_chi_max - amm_log_chi_min)
    amm_chi2_min = 10.0_num**(2.0_num * amm_log_chi_min)
    ! Scale the analytic tail so it meets the last table entry
    y = 10.0_num**(-2.0_num * amm_log_chi_max / 3.0_num)
    amm_tail_scale = 10.0_num**amm_log_f(n_amm) &
        / (y * (c_amm_tail0 + c_amm_tail1 * y))
    amm_table_loaded = .TRUE.

  END SUBROUTINE setup_anomalous_moment_table



  SUBROUTINE deallocate_spin_tables
    DEALLOCATE(amm_log_f)
  END SUBROUTINE deallocate_spin_tables


  ! FUnction for using anomalous_moment.table
  FUNCTION anomalous_moment_factor(chi2) RESULT(f)

    ! a(chi)/a(0) for chi^2 = chi2: 1 below the table, log-log linear
    ! interpolation inside it, scaled two-term asymptote above it
    REAL(num), INTENT(IN) :: chi2
    REAL(num) :: f, x, w, y
    INTEGER :: i

    IF (chi2 <= amm_chi2_min) THEN
      f = 1.0_num
      RETURN
    END IF

    ! log10(chi) straight from chi^2: no SQRT needed
    x = 0.5_num * LOG10(chi2)
    ! Upper truncation: use large-chi expansion
    ! Note: if upper truncation is needed, the Ritus-Narozhny 
    ! condition has already broken down and the loop contribution
    ! won't be valid. 
    IF (x >= amm_log_chi_max) THEN
      y = chi2**(-1.0_num / 3.0_num)
      f = amm_tail_scale * y * (c_amm_tail0 + c_amm_tail1 * y)
    ELSE
      w = (x - amm_log_chi_min) * amm_idlog
      i = INT(w)
      w = w - i
      i = i + 1
      f = 10.0_num**((1.0_num - w) * amm_log_f(i) + w * amm_log_f(i+1))
    END IF

  END FUNCTION anomalous_moment_factor



  ! The following four helpers implement Algorithms A, R and P in
  ! SPIN_DOCUMENTATION.md. They do not gather fields or change particles.
  SUBROUTINE compute_local_field_axes(px_photon, py_photon, pz_photon, &
      ex_part, ey_part, ez_part, bx_part, by_part, bz_part, &
      n, e_hat, b_hat, valid, charge_sign)

    REAL(num), INTENT(IN) :: px_photon, py_photon, pz_photon
    REAL(num), INTENT(IN) :: ex_part, ey_part, ez_part
    REAL(num), INTENT(IN) :: bx_part, by_part, bz_part
    REAL(num), INTENT(OUT) :: n(3), e_hat(3), b_hat(3)
    LOGICAL, INTENT(OUT) :: valid
    REAL(num), OPTIONAL, INTENT(IN) :: charge_sign
    REAL(num) :: p(3), electric(3), magnetic(3), crossed_e(3)
    REAL(num) :: scale, sign_field

    valid = .FALSE.
    n = 0.0_num
    e_hat = 0.0_num
    b_hat = 0.0_num
    p = (/ px_photon, py_photon, pz_photon /)
    ! Scaling avoids squaring the very small SI photon momenta.
    scale = MAXVAL(ABS(p))
    IF (scale <= TINY(1.0_num)) RETURN
    n = p / scale
    n = n / SQRT(DOT_PRODUCT(n, n))

    electric = (/ ex_part, ey_part, ez_part /)
    magnetic = c * (/ bx_part, by_part, bz_part /)
    ! crossed_E = E + c n x B - n (n dot E).
    ! Scale before projection; a relative tolerance rejects cancellation
    ! at round-off without imposing an absolute physical field cutoff.
    scale = MAX(MAXVAL(ABS(electric)), MAXVAL(ABS(magnetic)))
    IF (scale <= TINY(1.0_num)) RETURN
    electric = electric / scale
    magnetic = magnetic / scale
    crossed_e = electric - n * DOT_PRODUCT(n, electric) &
        + (/ n(2)*magnetic(3) - n(3)*magnetic(2), &
             n(3)*magnetic(1) - n(1)*magnetic(3), &
             n(1)*magnetic(2) - n(2)*magnetic(1) /)
    crossed_e = crossed_e - n * DOT_PRODUCT(n, crossed_e)
    scale = SQRT(DOT_PRODUCT(crossed_e, crossed_e))
    IF (scale <= 16.0_num * EPSILON(1.0_num)) RETURN

    ! T21 Eq. (51): charge is absorbed into the field axes. Default to
    ! the electron convention; use only the sign, never charge in SI.
    ! Physical spin/helicity signs in event kernels await a separate audit.
    sign_field = -1.0_num
    IF (PRESENT(charge_sign)) THEN
      IF (ABS(charge_sign) <= TINY(1.0_num)) RETURN
      sign_field = SIGN(1.0_num, charge_sign)
    END IF
    e_hat = sign_field * crossed_e / scale
    ! crossed_B = cB - cn(n dot B) - n cross E = crossed_E cross n.
    ! B_hat = E_hat cross n; no independent magnetic normalisation.
    b_hat = (/ e_hat(2)*n(3) - e_hat(3)*n(2), &
               e_hat(3)*n(1) - e_hat(1)*n(3), &
               e_hat(1)*n(2) - e_hat(2)*n(1) /)
    valid = .TRUE.

  END SUBROUTINE compute_local_field_axes


  SUBROUTINE reference_axis(n, r1, r2)

    ! Requires unit n. r1 = z cross n, r2 = n cross r1.
    REAL(num), INTENT(IN) :: n(3)
    REAL(num), INTENT(OUT) :: r1(3), r2(3)
    REAL(num) :: norm

    r1 = (/ -n(2), n(1), 0.0_num /)
    norm = SQRT(DOT_PRODUCT(r1, r1))
    IF (norm < 1.0e-12_num) THEN
      ! The documented x-axis fallback, projected transverse to n.
      ! This also handles directions extremely close to the poles.
      r1 = (/ 1.0_num, 0.0_num, 0.0_num /) - n(1) * n
      norm = SQRT(DOT_PRODUCT(r1, r1))
    END IF
    r1 = r1 / norm
    r2 = (/ n(2)*r1(3) - n(3)*r1(2), &
            n(3)*r1(1) - n(1)*r1(3), &
            n(1)*r1(2) - n(2)*r1(1) /)

  END SUBROUTINE reference_axis


  FUNCTION pol_stored_to_obs(r1, r2, e_hat, xi_stored) RESULT(xi_obs)

    REAL(num), INTENT(IN) :: r1(3), r2(3), e_hat(3), xi_stored(3)
    REAL(num) :: xi_obs(3), cw, sw, c2w, s2w

    ! Unit transverse axes supplied by Algorithms A and R.
    ! cos(2omega) = (E_hat dot r1)^2 - (E_hat dot r2)^2;
    ! sin(2omega) = 2 (E_hat dot r1) (E_hat dot r2).
    cw = DOT_PRODUCT(e_hat, r1)
    sw = DOT_PRODUCT(e_hat, r2)
    c2w = cw*cw - sw*sw
    s2w = 2.0_num * cw * sw
    ! Rotate xi_1 and xi_3 by (cos2omega -sin2omega;
    !                           sin2omega  cos2omega).
    ! xi_2 is circular; mixed-state length is preserved, not normalised.
    xi_obs(1) = c2w*xi_stored(1) - s2w*xi_stored(3)
    xi_obs(2) = xi_stored(2)
    xi_obs(3) = s2w*xi_stored(1) + c2w*xi_stored(3)

  END FUNCTION pol_stored_to_obs


  FUNCTION pol_obs_to_stored(r1, r2, e_hat, xi_obs) RESULT(xi_stored)

    REAL(num), INTENT(IN) :: r1(3), r2(3), e_hat(3), xi_obs(3)
    REAL(num) :: xi_stored(3), cw, sw, c2w, s2w

    ! Same double angle; inverse rotation is the transpose, not another
    ! stored-to-observed rotation. No photon direction argument is needed.
    cw = DOT_PRODUCT(e_hat, r1)
    sw = DOT_PRODUCT(e_hat, r2)
    c2w = cw*cw - sw*sw
    s2w = 2.0_num * cw * sw
    xi_stored(1) = c2w*xi_obs(1) + s2w*xi_obs(3)
    xi_stored(2) = xi_obs(2)
    xi_stored(3) = -s2w*xi_obs(1) + c2w*xi_obs(3)

  END FUNCTION pol_obs_to_stored


  !Get the lepton axes from compute_local_field_axes(..., charge_sign = 1.0_num). 
    !With +1 this returns the physical ε̂ and β̂ = ε̂ × n directly.
  !Take the spin components s·ε̂, s·β̂, s·v̂.
  !Look up R₁ and R₂ at x = 2λ/(3χ(1−λ)).
  !Evaluate w and g from polarised_lcfa.py, with qs = the sign of the charge.
  !Set the final spin, either to the mean g/w or by collapsing along ĝ.
  !Evaluate F/F₀ with that final spin to get ξ_obs.
  !Convert with pol_obs_to_stored(r1, r2, e_hat_photon, xi_obs). 
  !Here r₁, r₂ come from reference_axis(n) and Ê_photon = −ε̂; its sign doesn't matter.
  SUBROUTINE nlcs_final_states(p_lepton, e_field, b_field, charge_sign, &
      lam, chi, spin, xi_stored)
    ! p_lepton: momentum direction before emission; lam = photon energy /
    ! lepton energy; chi = EPOCH's eta. spin is INTENT(INOUT); xi_stored
    ! is the photon's Stokes vector on its reference axes.






#endif

END MODULE spin
