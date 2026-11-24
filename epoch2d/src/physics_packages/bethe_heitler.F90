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
!
!-------------------------------------------------------------------------------
!
! This module contains scripts used to calculate the Bethe-Heitler decay of high
! energy photons into electron-positron pairs. These routines are called from
! bremsstrahlung.F90, and describe the optical depth update, e-/e+ energy
! sampling, and angular deflections of the pair particles.
!
! These algorithms have been adapted from Geant4, and details are provided in
! the Geant4 Physics Reference Manual. The model we are using is valid for
! photon energies between 1.5 MeV and 100 GeV.
!
! To save on particle variable count, the bremsstrahlung optical depth will be
! used for Bethe-Heiter pair production, which was previously unused for photons

MODULE bethe_heitler
#ifdef BREMSSTRAHLUNG

  USE partlist
  USE particles
  USE calc_df
  USE setup
#ifdef SPIN
  USE spin, ONLY: reference_axis
#endif

  IMPLICIT NONE

  ! Pre-calculate constants to be used by these functions/subroutines
  REAL(num), PRIVATE, PARAMETER :: ubarn = 1.0e-34_num        ! 1 microbarn [m²]
  REAL(num), PRIVATE, PARAMETER :: gev_100 = 1.0e11_num * q0  ! 100 GeV [J]
  REAL(num), PRIVATE, PARAMETER :: mev_50 = 50.0e6_num * q0   ! 50 MeV [J]
  REAL(num), PRIVATE, PARAMETER :: mev_20 = 20.0e6_num * q0   ! 20 MeV [J]
  REAL(num), PRIVATE, PARAMETER :: mev_2 = 2.0e6_num * q0     ! 2 MeV [J]
  REAL(num), PRIVATE, PARAMETER :: kev_1500 = 1.5e6_num * q0  ! 1.5 MeV [J]
  REAL(num), PRIVATE, PARAMETER :: two_mc2 = 2.0_num * m0c2
  REAL(num), PRIVATE, PARAMETER :: extrap_fac = 1.0_num/(kev_1500 - two_mc2)

  ! Pre-calculated functions
  REAL(num), PRIVATE :: cdt
  INTEGER, PRIVATE :: z2ind(100)
  REAL(num), ALLOCATABLE, PRIVATE :: screen_fac(:), fz_low(:), fz_high(:)

CONTAINS

  SUBROUTINE init_bethe_heitler(size_z, z_values, z_to_index)

    ! The Bethe-Heitler algorithm uses some functions which are only Z-dependent
    ! so it makes no sense to recalculate these each step if the Z materials
    ! present don't change over a simulation. Hence, we will pre-calculate some
    ! functions and store them in private arrays
    !
    ! Here, z_values is an array containing all the unique z values present in
    ! the simulation, z2index is an array such that z2index(z) returns the index
    ! in z_values corresponding to z, and size_z is the size of z_values
    !
    ! The functions fz_low and fz_high describe screening and Coulomb
    ! corrections

    INTEGER, INTENT(IN) :: z_values(:)
    INTEGER, INTENT(IN) :: z_to_index(100)
    INTEGER, INTENT(IN) :: size_z
    INTEGER :: iz
    REAL(num) :: z_ion, iz13
    REAL(num) :: alpha_z2, alpha_z4, alpha_z6, coulomb_correction

    cdt = c*dt
    z2ind = z_to_index

    ! Allocate pre-calculated arrays, same indexing as z_values
    ALLOCATE(screen_fac(size_z), fz_low(size_z), fz_high(size_z))

    ! Loop over all possible atomic numbers, calculate terms if this number is
    ! present
    DO iz = 1,size_z

        z_ion = REAL(z_values(iz), num)

        ! Useful pre-factors
        iz13 = z_ion**(-1.0_num/3.0_num)
        screen_fac(iz) = 136.0_num * iz13

        ! Screening terms
        fz_low(iz) = 8.0_num/3.0_num * LOG(z_ion)
        alpha_z2 = (alpha*z_ion)**2
        alpha_z4 = alpha_z2**2
        alpha_z6 = alpha_z4*alpha_z2
        coulomb_correction = alpha_z2 * (1.0_num/(1.0_num + alpha_z2) &
            + 0.20206_num - 0.0369_num*alpha_z2 + 0.0083_num*alpha_z4 &
            - 0.002_num*(alpha_z6))
        fz_high(iz) = fz_low(iz) + 8.0_num * coulomb_correction

    END DO

  END SUBROUTINE init_bethe_heitler



  SUBROUTINE bethe_heitler_update_depth

    ! Cycles through all photon species and updates the Bethe-Heitler optical
    ! depths (labelled bremsstrahlung for photons)

    INTEGER :: iz, z_temp, ispecies
    REAL(num), ALLOCATABLE :: grid_num_density_ion(:,:)
    TYPE(particle), POINTER :: photon, next_photon
    REAL(num) :: part_x, part_y
    REAL(num) :: cross_sec, part_ni, delta_opdep

    ALLOCATE(grid_num_density_ion(1-ng:nx+ng, 1-ng:ny+ng))

    ! Calculate the number density of each ion species
    DO iz = 1, n_species
      ! Identify if the charge is greater than 1
      z_temp = species_list(iz)%atomic_no
      IF (z_temp < 1 .OR. z_temp > 100) CYCLE

      CALL calc_number_density(grid_num_density_ion, iz)
      CALL field_bc(grid_num_density_ion, ng)

      ! Update the optical depth for each photon species
      DO ispecies = 1, n_species
        ! Only update optical_depth_bremsstrahlung for the photon species
        IF (species_list(ispecies)%species_type == c_species_id_photon) THEN
          ! Cycle through all photons in this species
          photon => species_list(ispecies)%attached_list%head
          DO WHILE(ASSOCIATED(photon))
            ! Remember the next photon in the list, as the current photon could
            ! be removed
            next_photon => photon%next

            ! Photon energy is calculated on creation, and this energy never
            ! changes in the current version of the code (11/Feb/2021). Check it
            ! is high enough to pair produce.
            IF (photon%particle_energy < two_mc2) THEN
              photon => next_photon
              CYCLE
            END IF

            ! Calculate the cross section at this photon energy
            cross_sec = calc_bh_cross_sec(z_temp, photon%particle_energy)

            ! Get background number density at photon
            part_x = photon%part_pos(1) - x_grid_min_local
            part_y = photon%part_pos(2) - y_grid_min_local
            CALL grid_centred_var_at_particle(part_x, part_y, part_ni, &
                grid_num_density_ion)

            ! Update the photon optical depth
            delta_opdep = part_ni * cross_sec * cdt * boost_pairs
            photon%optical_depth_bremsstrahlung = &
                photon%optical_depth_bremsstrahlung - delta_opdep

            ! If photon optical depth drops below 0, create an e-/e+ pair and
            ! remove the photon from the simulation
            IF (photon%optical_depth_bremsstrahlung < 0.0_num) THEN
              CALL generate_pair(photon, z_temp, &
                  bethe_heitler_electron_species, &
                  bethe_heitler_positron_species, ispecies)
            END IF

            photon => next_photon
          END DO
        END IF
      END DO
    END DO

    DEALLOCATE(grid_num_density_ion)

  END SUBROUTINE bethe_heitler_update_depth



  function calc_bh_cross_sec(z_int, e_gamma)

    ! Calculates the Bethe-Heitler cross section in the presence of an ion
    ! species with atomic number z_int, for a photon of energy e_gamma.
    !
    ! Parametrisation taken from Geant4 source code, which acts as a fit to
    ! cross section data and is estimated to have less than 5% error.
    !
    ! For e_gamma above 100 GeV, the cross section is assumed constant. Below
    ! 1.5 MeV, results are extrapolated down to e_gamma = 2*m0*c².

    INTEGER, INTENT(IN) :: z_int
    REAL(num), INTENT(IN) :: e_gamma
    REAL(num) :: z_ion, e_gamma_calc
    REAL(num) :: x_fit, x2_fit, x3_fit, x4_fit, x5_fit
    REAL(num) :: f1_fit, f2_fit, f3_fit
    LOGICAL :: extrapolate_low
    REAL(num) :: calc_bh_cross_sec

    ! Convert integer atomic number to real
    z_ion = REAL(z_int, num)

    ! High/low energy treatment
    extrapolate_low = .FALSE.
    IF (e_gamma > gev_100) THEN
      e_gamma_calc = gev_100
    ELSE IF (e_gamma < kev_1500) THEN
      e_gamma_calc = kev_1500
      extrapolate_low = .TRUE.
    ELSE
      e_gamma_calc = e_gamma
    END IF

    ! Obtain fit variables
    x_fit = LOG(e_gamma_calc/m0c2)
    x2_fit = x_fit *x_fit
    x3_fit = x2_fit*x_fit
    x4_fit = x3_fit*x_fit
    x5_fit = x4_fit*x_fit

    ! Calculate fit functions
    f1_fit = (878.42_num - 1962.5_num*x_fit + 1294.9_num*x2_fit &
        - 200.28_num*x3_fit + 12.575_num*x4_fit - 0.28333_num*x5_fit)
    f2_fit = (-10.342_num + 17.692_num*x_fit - 8.2381_num*x2_fit &
        + 1.3063_num*x3_fit - 0.090815_num*x4_fit &
        + 2.3586e-3_num*x5_fit)
    f3_fit =  (-452.63_num + 1116.1_num*x_fit - 867.49_num*x2_fit &
        + 217.73_num*x3_fit - 20.467_num*x4_fit + 0.65372_num*x5_fit)

    ! Combine for the total cross section
    calc_bh_cross_sec = (z_ion + 1.0_num) &
        * (f1_fit*z_ion + f2_fit*z_ion**2 + f3_fit) * ubarn

    ! Consider extrapolation if photon energy is below 1.5 MeV
    IF (extrapolate_low) THEN
      calc_bh_cross_sec = calc_bh_cross_sec &
          * ((e_gamma - two_mc2) * extrap_fac)**2
    END IF

  END FUNCTION calc_bh_cross_sec



  SUBROUTINE generate_pair(photon, z_int, ielectron, ipositron, iphoton)

    ! This subroutine is called when "photon" undergoes Bethe-Heitler pair
    ! production. Macro-particles for e- and e+ are created, and the energies
    ! and directions are sampled using Geant4 algorithms. A full description of
    ! these equations may be found in the Geant4 Physics Reference Manual, and
    ! in the references therein. The photon is removed from the simulation here

    TYPE(particle), POINTER :: photon
    INTEGER, INTENT(IN) ::  z_int, ielectron, ipositron, iphoton
    TYPE(particle), POINTER :: new_electron, new_positron
    REAL(num) :: e_frac
    REAL(num) :: energy_1, momentum_1, theta_1, phi_1
    REAL(num) :: energy_2, momentum_2, theta_2, phi_2
    REAl(num) :: norm, dir_x, dir_y, dir_z
    REAL(num) :: p_save, e_save, theta_save, phi_save
    LOGICAL :: use_basis
#ifdef SPIN
    REAL(num) :: psi_1, psi_2, xi(3), n_gamma(3), r1(3), r2(3)
    LOGICAL :: bh_pol
#endif

    ! Create pair at photon position
    CALL create_particle(new_electron)
    CALL create_particle(new_positron)
    new_electron%part_pos = photon%part_pos
    new_positron%part_pos = photon%part_pos

    ! e- and e+ have the same weights as generating photon
    new_electron%weight = photon%weight / boost_pairs
    new_positron%weight = photon%weight / boost_pairs

    ! Calculate fractional energy split
    e_frac = energy_split(photon, z_int)

    ! Calculate momentum magnitude going to each particle
    energy_1 = e_frac * photon%particle_energy
    energy_2 = (1.0_num - e_frac) * photon%particle_energy
    momentum_1 = SQRT(energy_1**2 - m0c2**2)/c
    momentum_2 = SQRT(energy_2**2 - m0c2**2)/c

    ! Calculate scatter angles of each particle
    CALL sample_theta(energy_1, energy_2, theta_1, theta_2)
    phi_1 = random()*2.0_num*pi
    use_basis = .FALSE.
#ifdef SPIN
    ! Pair azimuths: screening functions at the sampled split, then
    ! the pair azimuth, measured from the photon's Stokes reference
    ! axis r1. 
    bh_pol = .FALSE.
    IF (bh_spin_transfer) THEN
      CALL bh_pol_setup(photon, z_int, energy_1, bh_pol, psi_1, psi_2, xi)
      IF (bh_pol .AND. bh_pol_azimuth) &
          CALL bh_pair_azimuth(photon%part_p, e_frac, psi_1, psi_2, xi, &
          phi_1, n_gamma, r1, r2, use_basis)
    END IF
#endif
    phi_2 = phi_1 + pi

    ! Calculate momentum direction of generating photon
    norm  = c / photon%particle_energy
    dir_x = photon%part_p(1) * norm
    dir_y = photon%part_p(2) * norm
    dir_z = photon%part_p(3) * norm

    ! Randomly swap properties of 1 and 2 (could be either e+ or e-)
    IF (random() > 0.5_num) THEN
      p_save = momentum_1
      e_save = energy_1
      theta_save = theta_1
      phi_save = phi_1

      momentum_1 = momentum_2
      energy_1 = energy_2
      theta_1 = theta_2
      phi_1 = phi_2

      momentum_2 = p_save
      energy_2 = e_save
      theta_2 = theta_save
      phi_2 = phi_save
    END IF

    IF (use_basis) THEN
#ifdef SPIN
      ! Azimuths are measured from r1 towards r2 about the photon
      ! direction (rotate_p uses a different azimuthal origin)
      new_electron%part_p(1:3) = momentum_1 &
          * bh_direction(n_gamma, r1, r2, COS(theta_1), phi_1)
      new_positron%part_p(1:3) = momentum_2 &
          * bh_direction(n_gamma, r1, r2, COS(theta_2), phi_2)
#endif
    ELSE
      ! Assign momenta to particles
      new_electron%part_p(1:3) = momentum_1 * (/ dir_x, dir_y, dir_z /)
      new_positron%part_p(1:3) = momentum_2 * (/ dir_x, dir_y, dir_z /)

      ! Rotate momenta for each particle
      CALL rotate_p(new_electron, COS(theta_1), phi_1, momentum_1)
      CALL rotate_p(new_positron, COS(theta_2), phi_2, momentum_2)
    END IF

    ! Set particle energies
    new_electron%particle_energy = energy_1
    new_positron%particle_energy = energy_2

#ifdef SPIN
    IF (bh_pol .AND. spin_bh /= c_spin_process_off) THEN
      IF (spin_bh == c_spin_process_unpolarised) xi(2) = 0.0_num
      CALL bh_pair_spins(energy_1 / photon%particle_energy, psi_1, psi_2, &
          xi(2), new_electron%part_p, new_positron%part_p, &
          new_electron%spin, new_positron%spin)
    END IF
#endif

    ! Save pair to particle list
    CALL add_particle_to_partlist(species_list(ielectron)%attached_list, &
        new_electron)
    CALL add_particle_to_partlist(species_list(ipositron)%attached_list, &
        new_positron)

    ! Remove photon
    IF (random() < 1.0_num / boost_pairs) THEN
      CALL remove_particle_from_partlist(species_list(iphoton)%attached_list, &
          photon)
      CALL destroy_particle(photon)
    ELSE
      photon%optical_depth_bremsstrahlung = &
          photon%optical_depth_bremsstrahlung - LOG(random())
    END IF

  END SUBROUTINE generate_pair


#ifdef SPIN
  SUBROUTINE bh_pol_setup(photon, z_int, e_lepton, valid, psi_1, psi_2, xi)

    ! Common set-up of the BH pair polarisation (bh_spin_transfer): the
    ! validity checks, the Olsen-Maximon screening functions at the
    ! sampled split, and the photon's stored Stokes vector. valid = F for
    ! photons below 20 MeV or a lepton below gamma = 10, where the
    ! high-energy, small-angle kernels do not apply; such pairs get zero
    ! spins and a uniform azimuth. psi_1, psi_2 use the same fits and
    ! f(Z) as energy_split, so the spin-summed kernel is the spectrum
    ! that energy_split samples. photon_polarisation = F (qed block)
    ! treats the photon as unpolarised, as for the other processes.
    TYPE(particle), POINTER, INTENT(IN) :: photon
    INTEGER, INTENT(IN) :: z_int
    REAL(num), INTENT(IN) :: e_lepton
    LOGICAL, INTENT(OUT) :: valid
    REAL(num), INTENT(OUT) :: psi_1, psi_2, xi(3)
    REAL(num) :: e_gamma, x, screen_term, f1, f2, fz

    valid = .FALSE.
    psi_1 = 0.0_num
    psi_2 = 0.0_num
    xi = 0.0_num
    e_gamma = photon%particle_energy
    IF (e_gamma < mev_20) RETURN
    IF (MIN(e_lepton, e_gamma - e_lepton) < 10.0_num * m0c2) RETURN
    valid = .TRUE.

    ! 0 < x < 1; psi is symmetric in x and 1 - x
    x = e_lepton / e_gamma
    screen_term = screen_fac(z2ind(z_int)) * m0c2 &
        / (e_gamma * x * (1.0_num - x))
    CALL get_screening_funcs(photon, screen_term, z_int, f1, f2, fz, &
        psi_1, psi_2)
    ! Clamp each component of the photon polarisation to [-1,1]
    xi = MIN(MAX(photon%spin, -1.0_num), 1.0_num)
#ifdef PHOTONS
    IF (.NOT. photon_polarisation) xi = 0.0_num
#endif

  END SUBROUTINE bh_pol_setup



  SUBROUTINE bh_pair_azimuth(p_photon, x, psi_1, psi_2, xi, phi, n, r1, &
      r2, use_basis)

    ! Azimuth of the pair for a linearly polarised photon. 
    ! Rejection sampling: accept the azimuth angle phi if
    ! r(1+A P_lin) <= 1+ A[xi_1 (cos2phi) - xi_3 sin(2phi)]  
      
    ! OM (Phys. Rev. 114, 887 (1959)) Eq. 10.1 summed over spins,
    ! integrated over one lepton's direction and over the other's polar
    ! angle, gives that lepton's azimuth phi about the photon direction:
    !   p(phi) = [1 + A (xi_1 cos 2phi - xi_3 sin 2phi)] / (2 pi),
    ! with phi measured from the Stokes reference axis r1 towards r2
    ! (reference_axis; (xi_1, xi_3) = P_lin (cos 2psi, -sin 2psi) for a
    ! linear axis at psi from r1, so leptons favour the polarisation
    ! plane) and A from bh_azimuth_asymmetry. A is symmetric in the two
    ! leptons, and EPOCH's pair is coplanar (phi_2 = phi_1 + pi), so both
    ! leptons follow this law. It is the one-lepton law, not the
    ! pair-bisector law of Gros and Bernard; spin-azimuth correlations
    ! are not modelled.
    ! The input phi, already uniform, is the first rejection candidate
    ! (envelope 1 + A P_lin). use_basis = F leaves phi unchanged (no
    ! linear polarisation), so unpolarised photons keep the old draw;
    ! otherwise phi is on (r1, r2) and n, r1, r2 are returned.
    REAL(num), INTENT(IN) :: p_photon(3), x, psi_1, psi_2, xi(3)
    REAL(num), INTENT(INOUT) :: phi
    REAL(num), INTENT(OUT) :: n(3), r1(3), r2(3)
    LOGICAL, INTENT(OUT) :: use_basis
    REAL(num) :: a, p_lin, bound
    LOGICAL :: valid

    use_basis = .FALSE.
    n = 0.0_num
    r1 = 0.0_num
    r2 = 0.0_num
    p_lin = SQRT(xi(1)**2 + xi(3)**2)
    a = bh_azimuth_asymmetry(x, psi_1, psi_2)
    IF (a * p_lin <= 0.0_num) RETURN
    CALL bh_unit_vector(p_photon, n, valid)
    IF (.NOT. valid) RETURN
    CALL reference_axis(n, r1, r2)

    bound = 1.0_num + a * p_lin
    DO
      IF (random() * bound <= 1.0_num + a * (xi(1) * COS(2.0_num * phi) &
          - xi(3) * SIN(2.0_num * phi))) EXIT
      phi = random() * 2.0_num * pi
    END DO
    use_basis = .TRUE.

  END SUBROUTINE bh_pair_azimuth



  PURE FUNCTION bh_azimuth_asymmetry(x, psi_1, psi_2) RESULT(a)
    ! The A function used for sampling the pair azimuth.
    ! x is the electron's energy fraction and psi1 and psi2 are the 
    ! screening functions.

    ! Linear-polarisation azimuthal asymmetry of one BH lepton, from
    ! Olsen and Maximon Eq. 10.1 (spins summed) integrated over the
    ! lepton's polar angle. With xi = 1/(1 + u^2), Eq. 10.1 has
    ! 1/2 (e1^2 + e2^2)(3 + 2 Gamma) + e1 e2 + 8 e1 e2 xi^2 Gamma |u.e|^2,
    ! and xi^2 u^2 = xi (1 - xi). OM Eq. 8.6, psi_1 = 6 + 4 int Gamma and
    ! psi_2 = 6 + 24 int xi (1 - xi) Gamma, turns the integral into
    !   A = 2/3 x y (psi_2 - 6) / [(x^2 + y^2) psi_1 + 2/3 x y psi_2],
    ! whose denominator is the spin-summed spectrum. The absolute
    ! normalisation of psi matters here: Geant4's phi_i - F(Z)/2 matches
    ! OM's (complete screening: psi_1 = 4 ln(183 Z^(-1/3)) - 4 f in both,
    ! up to the fit). A is clamped to [0, 1]; psi_2 < 6 is outside the
    ! model.
    REAL(num), INTENT(IN) :: x, psi_1, psi_2
    REAL(num) :: a, y, d

    a = 0.0_num
    y = 1.0_num - x
    d = (x * x + y * y) * psi_1 + 2.0_num * x * y * psi_2 / 3.0_num
    IF (d <= 0.0_num) RETURN
    a = 2.0_num * x * y * (psi_2 - 6.0_num) / (3.0_num * d)
    a = MIN(MAX(a, 0.0_num), 1.0_num)

  END FUNCTION bh_azimuth_asymmetry



  PURE FUNCTION bh_direction(n, r1, r2, cos_theta, phi) RESULT(dir)

    ! Unit vector at polar angle theta from n and azimuth phi from r1
    ! towards r2. sin(theta) >= 0 as in rotate_p.
    REAL(num), INTENT(IN) :: n(3), r1(3), r2(3), cos_theta, phi
    REAL(num) :: dir(3), sin_theta

    sin_theta = SQRT(MAX(1.0_num - cos_theta**2, 0.0_num))
    dir = cos_theta * n + sin_theta * (COS(phi) * r1 + SIN(phi) * r2)

  END FUNCTION bh_direction



  SUBROUTINE bh_pair_spins(x, psi_1, psi_2, xi2, p_electron, p_positron, &
      spin_electron, spin_positron)

    ! Helicities of a Bethe-Heitler pair from the angle-integrated
    ! Olsen-Maximon kernel (Phys. Rev. 114, 887 (1959), Eq. 10.8); see
    ! bh_spin_kernel. x is the electron's share of the photon energy and
    ! psi_1, psi_2, xi2 come from bh_pol_setup. The spins are rest-frame
    ! Cartesian vectors along each lepton's own momentum (Olsen and
    ! Maximon quantise along the photon direction, which differs by an
    ! angle of order m c^2 / E).
    ! spin_final_states (qed block) selects, as for NBW and LBW pairs,
    !   mean:     s_e = P_e p_e_hat, s_p = P_p p_p_hat; no random
    !             number and no pair correlation.
    !   collapse: s_e = +-p_e_hat with P(+) = (1 + P_e) / 2; given its
    !             sign s, s_p = +-p_p_hat with P(+) = (1 + k) / 2,
    !             k = (P_p + s C_zz) / (1 + s P_e). Two random numbers.
    ! Averaged over the random numbers, collapse gives the mean spins and
    ! the longitudinal correlation C_zz. The transverse correlations of
    ! Eq. 10.8, the only place linear photon polarisation enters the
    ! spins, are not carried: the pair is entangled there, and two spin
    ! vectors cannot represent it. Without -DPHOTONS there is no qed
    ! block, and mean is used.
    ! Test switches: bh_spin_screening = F uses psi_2 = psi_1 (the
    ! unscreened xi2 T(x)); bh_spin_correlation = F collapses the
    ! positron on its own marginal (k = P_p, so C_zz = P_e P_p).
    REAL(num), INTENT(IN) :: x, psi_1, psi_2, xi2
    REAL(num), INTENT(IN) :: p_electron(3), p_positron(3)
    REAL(num), INTENT(INOUT) :: spin_electron(3), spin_positron(3)
    REAL(num) :: psi_2_used, p_e, p_p, c_zz, s_e, s_p, k
    REAL(num) :: n_e(3), n_p(3)
    LOGICAL :: valid_e, valid_p, collapse

    CALL bh_unit_vector(p_electron, n_e, valid_e)
    CALL bh_unit_vector(p_positron, n_p, valid_p)
    IF (.NOT. (valid_e .AND. valid_p)) RETURN

    psi_2_used = psi_2
    IF (.NOT. bh_spin_screening) psi_2_used = psi_1
    CALL bh_spin_kernel(x, psi_1, psi_2_used, xi2, p_e, p_p, c_zz)
    IF (.NOT. bh_spin_correlation) c_zz = p_e * p_p

#ifdef PHOTONS
    collapse = spin_final_states == c_spin_final_collapse
#else
    collapse = .FALSE.
#endif
    IF (collapse) THEN
      s_e = 1.0_num
      IF (random() >= 0.5_num * (1.0_num + p_e)) s_e = -1.0_num
      ! 1 + s_e P_e > 0 for the sign drawn; this only guards round-off
      k = (p_p + s_e * c_zz) / MAX(1.0_num + s_e * p_e, TINY(1.0_num))
      k = MIN(MAX(k, -1.0_num), 1.0_num)
      s_p = 1.0_num
      IF (random() >= 0.5_num * (1.0_num + k)) s_p = -1.0_num
      spin_electron = s_e * n_e
      spin_positron = s_p * n_p
    ELSE
      spin_electron = p_e * n_e
      spin_positron = p_p * n_p
    END IF

  END SUBROUTINE bh_pair_spins



  PURE SUBROUTINE bh_spin_kernel(x, psi_1, psi_2, xi2, p_e, p_p, c_zz)

    ! Longitudinal part of the Olsen-Maximon angle-integrated pair-spin
    ! kernel (Phys. Rev. 114, 887 (1959), Eqs. 10.8 and 10.11). x is the
    ! electron's share of the photon energy, y = 1 - x the positron's,
    ! psi_1 and psi_2 the screening functions with the Coulomb
    ! correction (OM Eq. 8.6; psi_i = phi_i - f(Z) / 2 in the Geant4
    ! notation of get_screening_funcs) and xi2 the photon's circular
    ! Stokes parameter. With psi_d = psi_1 - 2 psi_2 / 3 and
    ! d = (x^2 + y^2) psi_1 + 2 x y psi_2 / 3, the spin-summed spectrum,
    !   P_e  = xi2 (x psi_1 - y psi_d) / d,
    !   P_p  = xi2 (y psi_1 - x psi_d) / d,
    !   C_zz = (2 x y psi_2 / 3 - (x^2 + y^2) psi_d) / d
    ! are the electron and positron helicities and their correlation
    ! <h_e h_p>. C_zz does not depend on the photon polarisation. For
    ! psi_1 = psi_2 (no screening), P_e = xi2 (4x - 1) / (4x^2 - 4x + 3)
    ! (Muller and Muller 2011) and C_zz = (4xy - 1) / (3 - 4xy).
    ! The outputs are clamped to [-1, 1] against fit round-off.
    REAL(num), INTENT(IN) :: x, psi_1, psi_2, xi2
    REAL(num), INTENT(OUT) :: p_e, p_p, c_zz
    REAL(num) :: y, xi, psi_d, d, x2y2

    p_e = 0.0_num
    p_p = 0.0_num
    c_zz = 0.0_num
    y = 1.0_num - x
    x2y2 = x * x + y * y
    psi_d = psi_1 - 2.0_num * psi_2 / 3.0_num
    d = x2y2 * psi_1 + 2.0_num * x * y * psi_2 / 3.0_num
    ! d > 0 wherever energy_split samples; this only guards misuse
    IF (d <= 0.0_num) RETURN

    xi = MIN(MAX(xi2, -1.0_num), 1.0_num)
    p_e = xi * (x * psi_1 - y * psi_d) / d
    p_p = xi * (y * psi_1 - x * psi_d) / d
    c_zz = (2.0_num * x * y * psi_2 / 3.0_num - x2y2 * psi_d) / d
    p_e = MIN(MAX(p_e, -1.0_num), 1.0_num)
    p_p = MIN(MAX(p_p, -1.0_num), 1.0_num)
    c_zz = MIN(MAX(c_zz, -1.0_num), 1.0_num)

  END SUBROUTINE bh_spin_kernel



  SUBROUTINE bh_unit_vector(v, n, valid)

    ! n = v / |v|, scaled first so that tiny momenta do not underflow;
    ! valid = F (and n = 0) for a zero vector
    REAL(num), INTENT(IN) :: v(3)
    REAL(num), INTENT(OUT) :: n(3)
    LOGICAL, INTENT(OUT) :: valid
    REAL(num) :: scale

    n = 0.0_num
    scale = MAXVAL(ABS(v))
    valid = scale > TINY(1.0_num)
    IF (.NOT. valid) RETURN
    n = v / scale
    n = n / SQRT(DOT_PRODUCT(n, n))

  END SUBROUTINE bh_unit_vector
#endif



  FUNCTION energy_split(photon, z_int)

    ! This subroutine calculates the energy fraction given to each particle in
    ! the e-/e+ pair

    TYPE(particle), POINTER, INTENT(IN) :: photon
    INTEGER, INTENT(IN) :: z_int
    REAL(num) :: split_0, split_1, split_min
    REAL(num) :: screen_min, screen_max, screen_term
    REAL(num) :: f10, f20, screen_func
    REAL(num) :: split_n1, split_n2, samp_cut
    LOGICAL :: reject
    REAL(num) :: r1, r2, r3
    REAL(num) :: f1_test, f2_test, fz_test
    REAL(num) :: energy_split

    ! Energy split which corresponds to one particle being at rest
    split_0 = m0c2 / photon%particle_energy

    ! Low energy sampling is uniform between split_0 and 2 MeV
    IF (photon%particle_energy < mev_2) THEN
      energy_split = (0.5_num - split_0)*random() + split_0
      RETURN
    END IF

    ! Parameters describing atomic screening
    screen_min = screen_fac(z2ind(z_int)) * (4.0_num * split_0)
    CALL get_screening_funcs(photon, screen_min, z_int, f10, f20, screen_func)

    ! Energy split which prevents screening from making the cross section go -ve
    screen_max = EXP((42.24_num - screen_func)/8.368_num) - 0.952_num
    split_1 = 0.5_num * (1.0_num - SQRT(1.0_num - screen_min/screen_max))

    ! Calculate the minimum energy fraction which may be transferred to one
    ! particle in the pair
    split_min = MAX(split_0, split_1)

    ! Term for splitting treatment in accept/reject algorithm
    split_n1 = (0.5_num - split_min)**2 * f10
    split_n2 = 1.5_num * f20
    samp_cut = split_n1/(split_n1 + split_n2)

    ! Perform accept/reject algorithm
    reject = .TRUE.
    DO WHILE (reject)
      r1 = random()
      r2 = random()
      r3 = random()

      ! Let r1 choose the sampling method
      IF (r1 < samp_cut) THEN
        energy_split = 0.5_num - (0.5_num - split_min) * r2**(1.0_num/3.0_num)
      ELSE
        energy_split = split_min + (0.5_num - split_min) * r2
      END IF

      ! Calculate corresponding screening functions
      screen_term = screen_fac(z2ind(z_int)) * split_0 &
          / (energy_split * (1.0_num - energy_split))
      CALL get_screening_funcs(photon, screen_term, z_int, f1_test, f2_test, &
          fz_test)

      ! Test for accept/reject condition
      IF (r1 < samp_cut) THEN
        IF (f1_test/f10 > r3) reject = .FALSE.
      ELSE
        IF (f2_test/f20 > r3) reject = .FALSE.
      END IF
    END DO

  END FUNCTION energy_split



  SUBROUTINE get_screening_funcs(photon, screen_term, iz, f1, f2, fz, &
      psi_1, psi_2)

    ! The Geant4 physics reference manual describes three functions, F1, F2 and
    ! f(Z) which are used to sample the energy split in pair creation. These
    ! are calculated for each sampled split in the accept/reject algorithm, and
    ! so are moved to this subroutine for conciceness.
    !
    ! The optional psi_1, psi_2 = phi_1,2 - f(Z)/2 are the screening
    ! functions of the spin-summed spectrum (Olsen and Maximon's psi_1,2),
    ! used by the BH spin kernel.

    TYPE(particle), POINTER, INTENT(IN) :: photon
    REAL(num), INTENT(IN) :: screen_term
    REAL(num) :: phi_1, phi_2
    REAL(num), INTENT(OUT) :: f1, f2, fz
    REAL(num), INTENT(OUT), OPTIONAL :: psi_1, psi_2
    INTEGER :: iz

    ! Auxiliary phi functions
    IF (screen_term <= 1.0_num) THEN
      phi_1 = 20.867_num - 3.242_num*screen_term + 0.625_num*screen_term**2
      phi_2 = 20.209_num - 1.930_num*screen_term - 0.086_num*screen_term**2
    ELSE
      phi_1 = 21.12_num - 4.184_num*LOG(screen_term + 0.952_num)
      phi_2 = phi_1
    END IF

    IF (photon%particle_energy < mev_50) THEN
      fz = fz_low(z2ind(iz))
    ELSE
      fz = fz_high(z2ind(iz))
    END IF

    ! Auxiliary screening functions. The spectrum
    ! [e^2 + (1-e)^2] (phi_1 - fz/2) + 2/3 e (1-e) (phi_2 - fz/2)
    ! equals 2/3 (1/2 - e)^2 f1 + f2/3 only with +0.5 phi_2 in f2, as in
    ! Geant4 (G4BetheHeitlerModel ScreenFunction2 = 1.5 PHI1 + 0.5 PHI2).
    ! Upstream EPOCH had -0.5 phi_2, which oversamples asymmetric splits.
    f1 = 3.0_num * phi_1 - phi_2 - fz
    f2 = 1.5_num * phi_1 + 0.5_num * phi_2 - fz

    IF (PRESENT(psi_1)) psi_1 = phi_1 - 0.5_num * fz
    IF (PRESENT(psi_2)) psi_2 = phi_2 - 0.5_num * fz

  END SUBROUTINE get_screening_funcs



  SUBROUTINE sample_theta(energy_1, energy_2, theta_1, theta_2)

    ! Let the polar angle be parallel to the photon momentum direction. This
    ! subroutine samples the angles theta_1 and theta_2 for the scattered
    ! electron/positron pair, using an algorithm described by the Geant4 physics
    ! reference manual

    REAL(num), INTENT(IN) :: energy_1, energy_2
    REAL(num), INTENT(OUT) :: theta_1, theta_2
    REAL(num) :: r1, r2, r3, u_samp

    r1 = random()
    r2 = random()
    r3 = random()

    IF (r1 < 0.25_num) THEN
      u_samp = -1.6_num * LOG(r2*r3)
    ELSE
      u_samp = -8.0_num * LOG(r2*r3) / 15.0_num
    END IF

    theta_1 = u_samp * m0c2 / energy_1
    theta_2 = u_samp * m0c2 / energy_2

  END SUBROUTINE sample_theta



  ! Calculates the value of a grid-centred variable part_var stored in the grid
  ! grid_var, averaged over the particle shape for a particle at position
  ! (part_x, part_y)

  SUBROUTINE grid_centred_var_at_particle(part_x, part_y, part_var, &
      grid_var)

    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(IN) :: grid_var(1-ng:,1-ng:)
    REAL(num), INTENT(OUT) :: part_var
    INTEGER :: cell_x1, cell_y1
    REAL(num) :: cell_x_r, cell_y_r
    REAL(num) :: cell_frac_x, cell_frac_y
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = (1.0_num)**c_ndims
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    ! The following method is lifted from photons.F90 (field_at_particle), for
    ! the cell-centered fields, taking into account the various particle shapes
#ifdef PARTICLE_SHAPE_TOPHAT
    cell_x_r = part_x / dx - 0.5_num
    cell_y_r = part_y / dy - 0.5_num
#else
    cell_x_r = part_x / dx
    cell_y_r = part_y / dy
#endif
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1
    cell_y1 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y1, num) - cell_y_r
    cell_y1 = cell_y1 + 1

#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#include "bspline3/part_var.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#include "tophat/part_var.inc"
#else
#include "triangle/gx.inc"
#include "triangle/part_var.inc"
#endif
    part_var = fac * part_var

  END SUBROUTINE grid_centred_var_at_particle

#endif
END MODULE bethe_heitler
