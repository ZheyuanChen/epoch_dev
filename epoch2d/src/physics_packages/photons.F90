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

MODULE photons
#ifdef PHOTONS

  USE collisions
  USE partlist
  USE utilities
  USE spin
  IMPLICIT NONE

  REAL(num), PRIVATE :: sig2cdt_dV_lbw
  REAL(num), PRIVATE :: cdt_dV
  REAL(num), PRIVATE :: i_LBW_amp_factor
#ifdef SPIN
  ! For Linear Breit-Wheeler usage.
  ! A candidate pair of photons i, j.
  ! In their centre-of-momentum (CM) frame photon i moves along axis q
  ! and photon j along -q; t1, t2 are the transverse axes of
  ! get_orthonormal, from which the pair azimuth is measured. Collision
  ! axes: (t1, t2, q) for photon i, (t1, -t2, -q) for photon j, both
  ! right-handed. ri1, ri2 (rj1, rj2) are the reference axes of the
  ! stored Stokes vector of i (j) in the CM frame (photon_cm_axes) and
  ! xi_i, xi_j the Stokes vectors on the collision axes. sig_0, sig_l,
  ! sig_c and sig_pol are sigma_0, sigma_L, sigma_C and the polarised
  ! cross-section. polarised = .FALSE. when both Stokes vectors are
  ! zero; then only sig_pol (= sigma_0) is set.
  TYPE lbw_pol_pair
    LOGICAL :: polarised
    REAL(num), DIMENSION(3) :: q, t1, t2, ri1, ri2, rj1, rj2, xi_i, xi_j
    REAL(num) :: sig_0, sig_l, sig_c, sig_pol
  END TYPE lbw_pol_pair
#endif
  SAVE
  REAL(num) :: part_pos_global, gamma_global, eta_global

CONTAINS

  SUBROUTINE setup_qed_module

    INTEGER :: ispecies, iu, io
    LOGICAL :: found

    ! Sanity check
    IF (produce_photons .AND. .NOT.ic_from_restart) THEN
      ! Identify if there exists any *populated* electron/positron species
      found = .FALSE.
      DO ispecies = 1, n_species
        IF (species_list(ispecies)%count <= 0) CYCLE

        IF (species_list(ispecies)%species_type == c_species_id_electron &
            .OR. species_list(ispecies)%species_type == c_species_id_positron &
            .OR. ispecies == photon_species ) THEN
          found = .TRUE.
          EXIT
        END IF
      END DO

      IF (.NOT.found) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Electron, positron and photon species are either ', &
                'unspecified or contain no'
            WRITE(io,*) 'particles. QED routines will do nothing.'
          END DO
        END IF
      END IF
    END IF

    ! Load the tables for the QED routines
    CALL setup_tables_qed
    IF (.NOT.ic_from_restart) THEN
      DO ispecies = 1, n_species
        CALL initialise_optical_depth(species_list(ispecies))
      END DO
    END IF

    ! The polarised cross-section has a larger maximum
    IF (LBW_polarisation) THEN
      sig2cdt_dV_lbw = 2.0_num * sigma_lBW_pol_max * c * dt / dx / dy &
          * LBW_amp_factor
    ELSE
      sig2cdt_dV_lbw = 2.0_num * sigma_lBW_max * c * dt / dx / dy &
          * LBW_amp_factor
    END IF
    cdt_dV = c * dt / dx / dy * LBW_amp_factor
    i_LBW_amp_factor = 1.0_num / LBW_amp_factor

  END SUBROUTINE setup_qed_module



  SUBROUTINE shutdown_qed_module

    CALL deallocate_tables_qed

  END SUBROUTINE shutdown_qed_module



  FUNCTION check_qed_variables()

    INTEGER :: check_qed_variables
    INTEGER :: io, iu, ispecies
    INTEGER :: first_electron = -1, first_positron = -1

    check_qed_variables = c_err_none

    IF (.NOT.use_qed) RETURN

    ! Linear Breit-Wheeler needs its pair species, with or without
    ! photon emission
    IF (use_LBW) THEN
      check_qed_variables = check_lbw_species()
      IF (check_qed_variables /= c_err_none) RETURN
    END IF

    ! If you're only doing radiation reaction force then don't need any special
    ! species, so don't do any checking here
    IF (.NOT.produce_photons) RETURN

    ! Identify if there exists any electron species
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type == c_species_id_electron &
          .AND. first_electron == -1) THEN
        first_electron = ispecies
      ELSE IF (species_list(ispecies)%species_type == c_species_id_positron &
          .AND. first_positron == -1) THEN
        first_positron = ispecies
      END IF
    END DO

    IF (first_electron < 0 .AND. first_positron < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No electron or positron species specified.'
          WRITE(io,*) 'Specify using "identify:electron" or "identify:positron"'
          WRITE(io,*) 'QED routines require at least one species of ', &
              'electrons or positrons.'
        END DO
      END IF
      check_qed_variables = c_err_missing_elements
      RETURN
    END IF

    IF (photon_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No photon species specified. Specify using ', &
              '"identify:photon"'
        END DO
      END IF
      check_qed_variables = c_err_missing_elements
      RETURN
    END IF

    ! If you're not producing pairs then you don't have to designate special
    ! electron or positron species so just return
    IF (.NOT.produce_pairs) RETURN

    IF (first_electron < 0 .OR. first_positron < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'To use pair production routines need at least one ', &
              'positron species and one ', 'electron species. Specify ', &
              'using "identify:electron" or "identify:positron"'
        END DO
      END IF
      check_qed_variables = c_err_missing_elements
      RETURN
    END IF

    IF (breit_wheeler_positron_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No Breit-Wheeler positron species specified.'
          WRITE(io,*) 'Specify using "identify:breit_wheeler_positron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_positron)%name), ' instead.'
        END DO
      END IF
      breit_wheeler_positron_species = first_positron
    END IF

#ifdef TRIDENT_PHOTONS
    IF (trident_positron_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No trident positron species specified.'
          WRITE(io,*) 'Specify using "identify:trident_positron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_positron)%name), ' instead.'
        END DO
      END IF
      trident_positron_species = first_positron
    END IF
#endif

    IF (breit_wheeler_electron_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No Breit-Wheeler electron species specified.'
          WRITE(io,*) 'Specify using "identify:breit_wheeler_electron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_electron)%name), ' instead.'
        END DO
      END IF
      breit_wheeler_electron_species = first_electron
    END IF

#ifdef TRIDENT_PHOTONS
    IF (trident_electron_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No trident electron species specified.'
          WRITE(io,*) 'Specify using "identify:trident_electron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_electron)%name), ' instead.'
        END DO
      END IF
      trident_electron_species = first_electron
    END IF
#endif
  END FUNCTION check_qed_variables



  FUNCTION check_lbw_species() RESULT(errcode)

    ! Pair species of linear Breit-Wheeler: the first electron and
    ! positron species unless identify:LBW_electron and
    ! identify:LBW_positron name them.
    INTEGER :: errcode
    INTEGER :: io, iu, ispecies
    INTEGER :: first_electron, first_positron

    errcode = c_err_none
    first_electron = -1
    first_positron = -1
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type == c_species_id_electron &
          .AND. first_electron == -1) THEN
        first_electron = ispecies
      ELSE IF (species_list(ispecies)%species_type &
          == c_species_id_positron .AND. first_positron == -1) THEN
        first_positron = ispecies
      END IF
    END DO

    IF ((lbw_electron_species < 0 .AND. first_electron < 0) &
        .OR. (lbw_positron_species < 0 .AND. first_positron < 0)) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'linear_breit_wheeler=T needs an electron and a ', &
              'positron species. Specify using "identify:LBW_electron" ', &
              'and "identify:LBW_positron".'
        END DO
      END IF
      errcode = c_err_missing_elements
      RETURN
    END IF

    IF (lbw_electron_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No linear Breit-Wheeler electron species specified.'
          WRITE(io,*) 'Specify using "identify:LBW_electron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_electron)%name), ' instead.'
        END DO
      END IF
      lbw_electron_species = first_electron
    END IF

    IF (lbw_positron_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No linear Breit-Wheeler positron species specified.'
          WRITE(io,*) 'Specify using "identify:LBW_positron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_positron)%name), ' instead.'
        END DO
      END IF
      lbw_positron_species = first_positron
    END IF

  END FUNCTION check_lbw_species



  SUBROUTINE setup_tables_qed

    ! Reads files epsilon.table, log_chi.table, energy_split.table
    ! and sets up appropriate tables

    REAL(num) :: etalog_min = 0.0_num, etalog_max = 0.0_num
    REAL(num) :: etalog_dx, chi_min, chi_dx
    REAL(num), ALLOCATABLE :: realbuf(:)
    INTEGER :: i, n, ichi2, iepsilon, ieta, ichi, bufsize, intbuf(6)

    IF (rank == 0) THEN
      OPEN(unit=lu, file=TRIM(qed_table_location)//'/hsokolov.table', &
          status='OLD')
      READ(lu,*) n_sample_h
      ALLOCATE(log_hsokolov(n_sample_h,2))
      DO ieta = 1, n_sample_h
        READ(lu,*) log_hsokolov(ieta,1), log_hsokolov(ieta,2)
      END DO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/pairprod.table', &
          status='OLD')
      READ(lu,*) n_sample_t
      ALLOCATE(log_tpair(n_sample_t,2))
      ALLOCATE(log_omegahat(n_sample_t,2))
      DO ichi = 1, n_sample_t
        READ(lu,*) log_tpair(ichi,1), log_omegahat(ichi,2), log_tpair(ichi,2)
      END DO
      CLOSE(unit=lu)
#ifdef SPIN
      ! pair_tables = lcfa: the exact LCFA T(chi) in place of Erber's
      ! approximation (the trident column is kept)
      IF (pair_tables == c_pair_tables_lcfa) &
          CALL read_pairprod_lcfa(n_sample_t, log_tpair(:,1), log_tpair(:,2))
#endif

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/ksi_sokolov.table', &
          status='OLD')
      READ(lu,*) n_sample_eta, n_sample_chi, etalog_min, etalog_max
      ALLOCATE(p_photon_energy(n_sample_eta,n_sample_chi))
      DO ieta = 1, n_sample_eta
        READ(lu,*) (p_photon_energy(ieta,ichi), ichi=1,n_sample_chi)
      END DO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/chimin.table', &
          status='OLD')
      ALLOCATE(chimin_table(n_sample_eta))
      READ(lu,*) (chimin_table(ieta), ieta=1,n_sample_eta)
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/log_chi2.table', &
          status='OLD')
      READ(lu,*) n_sample_chi2
      ALLOCATE(log_chi2(n_sample_chi2))
      DO ichi2 = 1, n_sample_chi2
        READ(lu,*) log_chi2(ichi2)
      END DO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/epsilon.table', &
          status='OLD')
      READ(lu,*) n_sample_epsilon
      ALLOCATE(epsilon_split(n_sample_epsilon))
      DO iepsilon = 1, n_sample_epsilon
        READ(lu,*) epsilon_split(iepsilon)
      END DO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/energy_split.table', &
          status='OLD')
      ALLOCATE(p_energy(n_sample_chi2,n_sample_epsilon))
      DO ichi2 = 1, n_sample_chi2
        DO iepsilon = 1, n_sample_epsilon
          READ(lu,*) p_energy(ichi2,iepsilon)
        END DO
      END DO
      CLOSE(unit=lu)

      ! Pack data for broadcasting to other processes

      intbuf(1) = n_sample_h
      intbuf(2) = n_sample_t
      intbuf(3) = n_sample_eta
      intbuf(4) = n_sample_chi
      intbuf(5) = n_sample_chi2
      intbuf(6) = n_sample_epsilon

      CALL MPI_BCAST(intbuf, 6, MPI_INTEGER, 0, comm, errcode)

      bufsize = n_sample_h * 2
      bufsize = bufsize + n_sample_t * 3
      bufsize = bufsize + 2 + n_sample_eta * n_sample_chi
      bufsize = bufsize + n_sample_eta
      bufsize = bufsize + n_sample_chi2
      bufsize = bufsize + n_sample_epsilon
      bufsize = bufsize + n_sample_chi2 * n_sample_epsilon

      ALLOCATE(realbuf(bufsize))
      n = 1

      DO i = 1, 2
        DO ieta = 1, n_sample_h
          realbuf(n) = log_hsokolov(ieta,i)
          n = n + 1
        END DO
      END DO

      DO ichi = 1, n_sample_t
        realbuf(n) = log_tpair(ichi,1)
        n = n + 1
        realbuf(n) = log_tpair(ichi,2)
        n = n + 1
        realbuf(n) = log_omegahat(ichi,2)
        n = n + 1
      END DO

      realbuf(n) = etalog_min
      n = n + 1
      realbuf(n) = etalog_max
      n = n + 1

      DO ichi = 1, n_sample_chi
        DO ieta = 1, n_sample_eta
          realbuf(n) = p_photon_energy(ieta,ichi)
          n = n + 1
        END DO
      END DO

      DO ieta = 1, n_sample_eta
        realbuf(n) = chimin_table(ieta)
        n = n + 1
      END DO

      DO ichi2 = 1, n_sample_chi2
        realbuf(n) = log_chi2(ichi2)
        n = n + 1
      END DO

      DO iepsilon = 1, n_sample_epsilon
        realbuf(n) = epsilon_split(iepsilon)
        n = n + 1
      END DO

      DO iepsilon = 1, n_sample_epsilon
        DO ichi2 = 1, n_sample_chi2
          realbuf(n) = p_energy(ichi2,iepsilon)
          n = n + 1
        END DO
      END DO

      CALL MPI_BCAST(realbuf, bufsize, mpireal, 0, comm, errcode)

      DEALLOCATE(realbuf)
    ELSE
      ! Unpack data from rank zero

      CALL MPI_BCAST(intbuf, 6, MPI_INTEGER, 0, comm, errcode)

      n_sample_h       = intbuf(1)
      n_sample_t       = intbuf(2)
      n_sample_eta     = intbuf(3)
      n_sample_chi     = intbuf(4)
      n_sample_chi2    = intbuf(5)
      n_sample_epsilon = intbuf(6)

      bufsize = n_sample_h * 2
      bufsize = bufsize + n_sample_t * 3
      bufsize = bufsize + 2 + n_sample_eta * n_sample_chi
      bufsize = bufsize + n_sample_eta
      bufsize = bufsize + n_sample_chi2
      bufsize = bufsize + n_sample_epsilon
      bufsize = bufsize + n_sample_chi2 * n_sample_epsilon

      ALLOCATE(realbuf(bufsize))
      n = 1

      CALL MPI_BCAST(realbuf, bufsize, mpireal, 0, comm, errcode)

      ALLOCATE(log_hsokolov(n_sample_h,2))
      DO i = 1, 2
        DO ieta = 1, n_sample_h
          log_hsokolov(ieta,i) = realbuf(n)
          n = n + 1
        END DO
      END DO

      ALLOCATE(log_tpair(n_sample_t,2))
      ALLOCATE(log_omegahat(n_sample_t,2))
      DO ichi = 1, n_sample_t
        log_tpair(ichi,1)    = realbuf(n)
        n = n + 1
        log_tpair(ichi,2)    = realbuf(n)
        n = n + 1
        log_omegahat(ichi,2) = realbuf(n)
        n = n + 1
      END DO

      etalog_min = realbuf(n)
      n = n + 1
      etalog_max = realbuf(n)
      n = n + 1

      ALLOCATE(p_photon_energy(n_sample_eta,n_sample_chi))
      DO ichi = 1, n_sample_chi
        DO ieta = 1, n_sample_eta
          p_photon_energy(ieta,ichi) = realbuf(n)
          n = n + 1
        END DO
      END DO

      ALLOCATE(chimin_table(n_sample_eta))
      DO ieta = 1, n_sample_eta
        chimin_table(ieta) = realbuf(n)
        n = n + 1
      END DO

      ALLOCATE(log_chi2(n_sample_chi2))
      DO ichi2 = 1, n_sample_chi2
        log_chi2(ichi2) = realbuf(n)
        n = n + 1
      END DO

      ALLOCATE(epsilon_split(n_sample_epsilon))
      DO iepsilon = 1, n_sample_epsilon
        epsilon_split(iepsilon) = realbuf(n)
        n = n + 1
      END DO

      ALLOCATE(p_energy(n_sample_chi2,n_sample_epsilon))
      DO iepsilon = 1, n_sample_epsilon
        DO ichi2 = 1, n_sample_chi2
          p_energy(ichi2,iepsilon) = realbuf(n)
          n = n + 1
        END DO
      END DO

      DEALLOCATE(realbuf)
    END IF
#ifdef SPIN
    ! Table for a(chi); without it leptons keep their constant moment
    IF (use_chi_dependent_moment) CALL setup_anomalous_moment_table
    ! Spin asymmetry A(eta) of the NLCS rate
    IF (use_spin_qed) CALL setup_nlcs_asymmetry_table
    ! K_1/3/K_2/3 and IntK_1/3/K_2/3
    IF (use_spin_qed) CALL setup_bessel_k_table
    ! Polarisation asymmetry A(chi) and energy split of NBW
    IF (use_spin_qed .OR. pair_tables == c_pair_tables_lcfa) THEN
      CALL setup_nbw_asymmetry_table
      CALL setup_nbw_split_tables
    END IF
#endif
    log_omegahat(:,1) = log_tpair(:,1)

    ALLOCATE(log_eta(n_sample_eta))
    ALLOCATE(log_chi(n_sample_eta,n_sample_chi))

    etalog_dx = (etalog_max - etalog_min) / REAL(n_sample_eta-1,num)
    DO ieta = 1, n_sample_eta
      log_eta(ieta) = etalog_min + REAL(ieta-1,num) * etalog_dx
      chi_min = LOG10(chimin_table(ieta))
      chi_dx  = (log_eta(ieta) - LOG10(2.0_num) - chi_min) &
          / REAL(n_sample_chi-1,num)
      DO ichi = 1, n_sample_chi
        log_chi(ieta,ichi) = chi_min + REAL(ichi-1,num) * chi_dx
      END DO
    END DO

    DEALLOCATE(chimin_table)

#ifdef SPIN
    ! Spin-resolved photon energies: quantile tables, or the legacy CDFs
    ! on the grid just built
    IF (use_spin_qed) THEN
      IF (use_nlcs_quantile_tables) THEN
        CALL setup_nlcs_quantile_tables
      ELSE
        CALL setup_nlcs_spin_cdf_tables
      END IF
    END IF
#endif

  END SUBROUTINE setup_tables_qed



  SUBROUTINE deallocate_tables_qed

    DEALLOCATE(log_chi2)
    DEALLOCATE(epsilon_split)
    DEALLOCATE(p_energy)
    DEALLOCATE(log_hsokolov)
    DEALLOCATE(log_eta)
    DEALLOCATE(log_chi)
    DEALLOCATE(p_photon_energy)
    DEALLOCATE(log_tpair)
    DEALLOCATE(log_omegahat)
#ifdef SPIN
    CALL deallocate_spin_tables
#endif

  END SUBROUTINE deallocate_tables_qed



  SUBROUTINE initialise_optical_depth(current_species)

    ! Resets optical depth (to random number) of all particles
    TYPE(particle_species) :: current_species
    TYPE(particle), POINTER :: current
    REAL(num) :: p_tau

    current => current_species%attached_list%head
    DO WHILE(ASSOCIATED(current))
      p_tau = random()
      current%optical_depth = -LOG(1.0_num - p_tau)

#ifdef TRIDENT_PHOTONS
      p_tau = random()
      current%optical_depth_tri = -LOG(1.0_num - p_tau)
#endif
      current => current%next
    END DO

  END SUBROUTINE initialise_optical_depth



  FUNCTION reset_optical_depth()

    ! Resets optical depth of particle
    REAL(num) :: reset_optical_depth
    REAL(num) :: p_tau

    p_tau = random()
    reset_optical_depth = -LOG(1.0_num - p_tau)

  END FUNCTION reset_optical_depth


  ! A 14/10/26 update enables sub-timestep sampling
  SUBROUTINE qed_update_optical_depth

    ! Updates the optical depth for electrons and photons
    INTEGER :: ispecies
    TYPE(particle), POINTER :: current, next_pt

    REAL(num) :: part_x, part_y
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: dir_x, dir_y, dir_z
    REAL(num) :: eta, chi_val, part_e, gamma_rel, norm
    REAL(num) :: dtau, f_ll(3)
    LOGICAL :: legacy_step
#ifdef SPIN
    REAL(num) :: dtau_avg, f_event, f_rest
    REAL(num) :: e_at_part(3), b_at_part(3), zeta, n_hat(3)
    REAL(num) :: a_spin, x_spin
    REAL(num) :: e_hat(3), ax1(3), ax2(3), xi_obs(3)
    LOGICAL :: spin_valid, pol_valid
#endif

    ! Emission handled here, one event at most per step (legacy), or
    ! inside lepton_emission_substeps
    legacy_step = use_continuous_emission .OR. .NOT. use_emission_substeps

    DO ispecies = 1, n_species

      ! First consider electrons and positrons
      IF (species_list(ispecies)%species_type == c_species_id_electron &
          .OR. species_list(ispecies)%species_type == c_species_id_positron) &
          THEN
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          ! Find eta at particle position
          part_x  = current%part_pos(1) - x_grid_min_local
          part_y  = current%part_pos(2) - y_grid_min_local
          part_ux = current%part_p(1) / mc0
          part_uy = current%part_p(2) / mc0
          part_uz = current%part_p(3) / mc0
          gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)

#ifdef SPIN
          eta = calculate_eta(part_x, part_y, part_ux, part_uy, &
              part_uz, gamma_rel, e_at_part, b_at_part)
#else
          eta = calculate_eta(part_x, part_y, part_ux, part_uy, &
              part_uz, gamma_rel)
#endif

          ! If using continuous emission emit a photon macroparticle at every timestep
          IF (use_continuous_emission) THEN
            ! The other Landau-Lifshitz terms (field-gradient and
            ! quadratic fields), and the midpoint kick's correction, from
            ! the momentum before this step's radiation reaction
            IF (rr_correction) CALL rr_correction_force(part_x, part_y, &
                part_ux, part_uy, part_uz, gamma_rel, &
                species_list(ispecies)%charge, eta, f_ll)
            CALL generate_photon(current, photon_species, eta)
            IF (rr_correction) &
                current%part_p = current%part_p + dt * f_ll
          ELSE IF (use_emission_substeps) THEN
            ! The whole step, event by event (emissions and spin)
#ifdef SPIN
            CALL lepton_emission_substeps(current, ispecies, part_x, &
                part_y, eta, gamma_rel, e_at_part, b_at_part)
#else
            CALL lepton_emission_substeps(current, part_x, part_y, eta, &
                gamma_rel)
#endif
          ELSE
            ! Legacy step: decrement here, at most one emission below
            ! Maybe a better logic is needed here to avoid double-setting
            ! dtaut
            dtau = delta_optical_depth(eta, gamma_rel)
#ifdef SPIN
            ! If spin is used with full-sampling (legacy) method
            ! the only change here is to use spin-dependent rather than
            ! independent NLCS rate.
            ! Spin-dependent rate W_bar (1 + A zeta), zeta = qs s.beta_hat
            IF (use_spin_qed) THEN
              CALL lepton_spin_projection(part_ux, part_uy, part_uz, &
                  gamma_rel, e_at_part, b_at_part, current%spin, &
                  species_list(ispecies)%charge, zeta, n_hat, spin_valid)
              a_spin = nlcs_asymmetry(eta)
              x_spin = a_spin * dtau
              ! spin_nlcs (tests): unpolarised and photon_only keep
              ! W_bar, off emits nothing; the survivor update keeps
              ! x_spin either way (photon_only skips it)
              IF (spin_nlcs == c_spin_process_resolved) THEN
                dtau = spin_delta_optical_depth(dtau, a_spin, zeta, &
                    spin_second_order_depth)
              ELSE IF (spin_nlcs == c_spin_process_off) THEN
                dtau = 0.0_num
              END IF
            END IF
#endif
            current%optical_depth = current%optical_depth - dtau
          ENDIF

#ifdef TRIDENT_PHOTONS
          ! Trident pair production is part of produce_pairs, exactly
          ! like the Breit-Wheeler photon branch below -- skip it
          ! entirely when produce_pairs = F, matching
          ! check_qed_variables, which already skips assigning default
          ! trident_electron_species/trident_positron_species in that
          ! case. Without this guard those stay at their unset -1
          ! default and generate_pair_tri segfaults on
          ! species_list(-1).
          IF (produce_pairs) THEN
            current%optical_depth_tri = current%optical_depth_tri &
                - delta_optical_depth_tri(eta, gamma_rel)
          END IF
#endif
          ! Legacy step only (use_emission_substeps = F): if the
          ! optical depth dropped below zero generate photon...
          IF (legacy_step) THEN
            IF (current%optical_depth < 0.0_num) THEN
#ifdef SPIN 
              ! A few flags for testing purposes.
              IF (use_spin_qed .AND. spin_nlcs == c_spin_process_resolved) &
                  THEN
                CALL generate_photon(current, photon_species, eta, zeta, &
                    a_spin, e_at_part, b_at_part, &
                    species_list(ispecies)%charge)
              ELSE IF (use_spin_qed &
                  .AND. spin_nlcs == c_spin_process_photon_only) THEN
                ! Spin-averaged spectrum (zeta = 0); generate_photon
                ! leaves the lepton spin unchanged in this mode
                CALL generate_photon(current, photon_species, eta, &
                    0.0_num, a_spin, e_at_part, b_at_part, &
                    species_list(ispecies)%charge)
              ELSE
                CALL generate_photon(current, photon_species, eta)
              END IF
#else
              CALL generate_photon(current, photon_species, eta)
#endif
              ! ... and reset optical depth
              current%optical_depth = reset_optical_depth()
#ifdef SPIN
            ELSE IF (use_spin_qed .AND. spin_no_transition &
                .AND. spin_nlcs /= c_spin_process_photon_only &
                .AND. .NOT. use_continuous_emission) THEN
              ! No emission this step: survivor (non-transition) spin
              IF (spin_valid) CALL lepton_no_emission_update(current%spin, &
                  n_hat, x_spin)
#endif
            END IF
          END IF

#ifdef TRIDENT_PHOTONS
          IF (produce_pairs .AND. current%optical_depth_tri < 0.0_num) THEN
            CALL generate_pair_tri(current, trident_electron_species, &
                trident_positron_species)
            ! ... and reset optical depth
            current%optical_depth_tri = &
                current%optical_depth_tri + reset_optical_depth()
          END IF
#endif
          current => current%next
        END DO

      ! and finally photons
      ELSE IF (species_list(ispecies)%species_type == c_species_id_photon &
          .AND. produce_pairs) THEN
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          ! Current may be deleted
          next_pt => current%next
          part_x  = current%part_pos(1) - x_grid_min_local
          part_y  = current%part_pos(2) - y_grid_min_local
          norm  = c / current%particle_energy
          dir_x = current%part_p(1) * norm
          dir_y = current%part_p(2) * norm
          dir_z = current%part_p(3) * norm
          part_e  = current%particle_energy / m0 / c**2
          
          ! We modified the calculate_chi function such that it could also 
          ! return the field at particle values should we need (and we 
          ! need them for spin NBW/NLCS)
#ifdef SPIN
          chi_val = calculate_chi(part_x, part_y, dir_x, dir_y, &
              dir_z, part_e, e_at_part, b_at_part)
#else
          chi_val = calculate_chi(part_x, part_y, dir_x, dir_y, &
              dir_z, part_e)
#endif
          ! Reuse the variable dtau to update photon optical depth
          dtau = delta_optical_depth_photon(chi_val, part_e)
#ifdef SPIN
          ! Polarisation-dependent rate W_bar (1 - A xi1_obs)
          IF (use_spin_qed) THEN
            CALL photon_pol_frame(current%part_p, e_at_part, b_at_part, &
                n_hat, e_hat, ax1, ax2, pol_valid)
            ! photon_polarisation = F: unpolarised rate, split and pair
            ! spins, and no survivor update
            IF (.NOT. photon_polarisation) pol_valid = .FALSE.
            xi_obs = 0.0_num
            IF (pol_valid) xi_obs = pol_stored_to_obs(ax1, ax2, e_hat, &
                current%spin)
            a_spin = nbw_asymmetry(chi_val)
            ! The survivor update uses the physical (unboosted) decrement
            x_spin = a_spin * dtau / boost_pairs
            dtau_avg = dtau
            ! spin_nbw (tests): unpolarised keeps W_bar, off decays
            ! nothing; the survivor update keeps x_spin either way
            IF (spin_nbw == c_spin_process_resolved) THEN
              dtau = spin_delta_optical_depth(dtau, a_spin, -xi_obs(1), &
                  spin_second_order_depth)
            ELSE IF (spin_nbw == c_spin_process_off) THEN
              dtau = 0.0_num
            END IF
          END IF
#endif
          current%optical_depth = current%optical_depth - dtau
          ! If optical depth dropped below zero generate pair...
          ! (the photon is removed from its own species, ispecies, which
          ! need not be the QED photon species)

          ! Need to think about boost_pairs
          IF (current%optical_depth < 0.0_num) THEN
#ifdef SPIN
            ! use_emission_substeps: the photon is a survivor until it
            ! decays, at the fraction f of the step
            f_rest = 0.0_num
            IF (use_spin_qed .AND. use_emission_substeps &
                .AND. spin_no_transition .AND. pol_valid) THEN
              f_event = 1.0_num + current%optical_depth / dtau
              IF (spin_nbw == c_spin_process_resolved &
                  .AND. spin_second_order_depth) &
                  f_event = spin_event_fraction(current%optical_depth &
                  + dtau, dtau_avg, a_spin, -xi_obs(1))
              xi_obs = photon_no_decay_update(xi_obs, TANH(f_event * x_spin))
              current%spin = pol_obs_to_stored(ax1, ax2, e_hat, xi_obs)
              f_rest = 1.0_num - f_event
            END IF
            IF (use_spin_qed .AND. spin_nbw == c_spin_process_resolved) THEN
              CALL generate_pair(current, chi_val, ispecies, &
                  breit_wheeler_electron_species, &
                  breit_wheeler_positron_species, xi_obs, a_spin, &
                  e_at_part, b_at_part)
            ELSE
              CALL generate_pair(current, chi_val, ispecies, &
                  breit_wheeler_electron_species, &
                  breit_wheeler_positron_species)
            END IF
            ! A photon kept by boost_pairs survives the rest of the step
            ! Need to decide whether the polarisation needs to be updated
            ! As a result of partial decay.
            IF (ASSOCIATED(current) .AND. f_rest > 0.0_num) &
                current%spin = pol_obs_to_stored(ax1, ax2, e_hat, &
                photon_no_decay_update(xi_obs, TANH(f_rest * x_spin)))
#else
            CALL generate_pair(current, chi_val, ispecies, &
                breit_wheeler_electron_species, breit_wheeler_positron_species)
#endif
#ifdef SPIN
          ELSE IF (use_spin_qed .AND. spin_no_transition .AND. pol_valid) &
              THEN
            ! No decay this step: survivor (non-transition) polarisation,
            ! dichroism only (no birefringence table yet)
            current%spin = pol_obs_to_stored(ax1, ax2, e_hat, &
                photon_no_decay_update(xi_obs, TANH(x_spin)))
#endif
          END IF
          current => next_pt
        END DO
      END IF
    END DO

  END SUBROUTINE qed_update_optical_depth



#ifdef SPIN
  SUBROUTINE lepton_emission_substeps(current, ispecies, part_x, part_y, &
      eta_start, gamma_start, e_at_part, b_at_part)
#else
  SUBROUTINE lepton_emission_substeps(current, part_x, part_y, &
      eta_start, gamma_start)
#endif

    ! This is a substantial change to both spin and QED (NLCS).
    ! The crucial problem is that EPOCH previously emits only at full 
    ! time steps where the optical depth is updated and checked against
    ! 0. However, the quantum non-transition event is (previously) either
    ! a yes or no. For the case where an emission happens, there's no quantum
    ! non-transition update, even if the optical depth drops to 0 (in theory)
    ! at, say, 0.1*time_step. This introduces noticeable but not-so-significant
    ! errors. Also, this error depends on the resolution of time-steps--a
    ! shorter timestep introduces smaller error. 

    ! Nevertheless, we can fix this probelm quite easily.

    ! Discrete photon emission (NLCS) of one lepton over one step, with
    ! the fields at the lepton frozen for the step (qed block
    ! use_emission_substeps = T, the default). The step is followed
    ! event by event. If the optical depth runs out at a fraction f of
    ! the remaining step, the spin first gets the survivor
    ! (non-transition) update for that fraction, the photon is then
    ! emitted with the spin of that moment, a new optical depth is
    ! drawn, and the rest of the step, 1 - f, is treated the same way.
    ! A lepton can therefore emit more than once per step. The legacy
    ! step (one emission at most, the overshoot discarded, no survivor
    ! update in an emission step) emits at W (1 - W dt / 2) instead of
    ! W, and applies the survivor update at a rate W (1 - W dt).
    ! f = tau / dtau for the first-order decrement; with
    ! spin_second_order_depth the decrement is not linear in time and f
    ! comes from spin_event_fraction. After a recoil, eta and gamma are
    ! recomputed from the new momentum in the same fields.
    TYPE(particle), POINTER :: current
#ifdef SPIN
    INTEGER, INTENT(IN) :: ispecies
#endif
    REAL(num), INTENT(IN) :: part_x, part_y, eta_start, gamma_start
#ifdef SPIN
    REAL(num), INTENT(IN) :: e_at_part(3), b_at_part(3)
#endif
    ! Safety net: more events than this in one step means a time
    ! step far too large for the emission rate
    INTEGER, PARAMETER :: max_events = 10000
    INTEGER :: n_event
    REAL(num) :: eta, gamma_rel, frac, dtau0, dtau, f
    REAL(num) :: part_ux, part_uy, part_uz
#ifdef SPIN
    REAL(num) :: zeta, n_hat(3), a_spin, x_spin, charge
    LOGICAL :: spin_valid, survivor, resolved
#endif

    eta = eta_start
    gamma_rel = gamma_start
    frac = 1.0_num
#ifdef SPIN
    charge = species_list(ispecies)%charge
    resolved = use_spin_qed .AND. spin_nlcs == c_spin_process_resolved
    survivor = use_spin_qed .AND. spin_no_transition &
        .AND. spin_nlcs /= c_spin_process_photon_only
    zeta = 0.0_num
    x_spin = 0.0_num
    a_spin = 0.0_num
    spin_valid = .FALSE.
#endif

    DO n_event = 1, max_events
      ! Decrement over the remaining fraction of the step
      dtau0 = frac * delta_optical_depth(eta, gamma_rel)
      dtau = dtau0
#ifdef SPIN
      ! Spin-dependent rate W_bar (1 + A zeta), zeta = qs s.beta_hat;
      ! spin_nlcs (tests): unpolarised and photon_only keep W_bar; off
      ! emits nothing; the survivor update keeps x_spin either way
      IF (use_spin_qed) THEN
        part_ux = current%part_p(1) / mc0
        part_uy = current%part_p(2) / mc0
        part_uz = current%part_p(3) / mc0
        CALL lepton_spin_projection(part_ux, part_uy, part_uz, &
            gamma_rel, e_at_part, b_at_part, current%spin, charge, zeta, &
            n_hat, spin_valid)
        a_spin = nlcs_asymmetry(eta)
        x_spin = a_spin * dtau0
        IF (resolved) THEN
          dtau = spin_delta_optical_depth(dtau0, a_spin, zeta, &
              spin_second_order_depth)
        ELSE IF (spin_nlcs == c_spin_process_off) THEN
          dtau = 0.0_num
        END IF
      END IF
#endif

      ! No (further) emission in this step
      IF (current%optical_depth >= dtau) THEN
        current%optical_depth = current%optical_depth - dtau
#ifdef SPIN
        IF (survivor .AND. spin_valid) &
            CALL lepton_no_emission_update(current%spin, n_hat, x_spin)
#endif
        RETURN
      END IF

      ! An emission at the fraction f of the remaining step
      f = current%optical_depth / dtau
#ifdef SPIN
      IF (resolved .AND. spin_second_order_depth) &
          f = spin_event_fraction(current%optical_depth, dtau0, a_spin, &
          zeta)
      ! The spin evolves as a survivor's until the emission
      IF (survivor .AND. spin_valid) THEN
        CALL lepton_no_emission_update(current%spin, n_hat, f * x_spin)
        zeta = DOT_PRODUCT(current%spin, n_hat)
      END IF
      IF (resolved) THEN
        CALL generate_photon(current, photon_species, eta, zeta, a_spin, &
            e_at_part, b_at_part, charge)
      ELSE IF (use_spin_qed &
          .AND. spin_nlcs == c_spin_process_photon_only) THEN
        ! Spin-averaged spectrum (zeta = 0); generate_photon leaves the
        ! lepton spin unchanged in this mode
        CALL generate_photon(current, photon_species, eta, 0.0_num, &
            a_spin, e_at_part, b_at_part, charge)
      ELSE
        CALL generate_photon(current, photon_species, eta)
      END IF
#else
      CALL generate_photon(current, photon_species, eta)
#endif
      current%optical_depth = reset_optical_depth()

      ! The rest of the step, with the recoiled momentum
      frac = frac * (1.0_num - f)
      IF (frac <= 0.0_num) RETURN
      IF (use_radiation_reaction) THEN
        part_ux = current%part_p(1) / mc0
        part_uy = current%part_p(2) / mc0
        part_uz = current%part_p(3) / mc0
        gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
        eta = calculate_eta(part_x, part_y, part_ux, part_uy, part_uz, &
            gamma_rel)
      END IF
    END DO

  END SUBROUTINE lepton_emission_substeps



  FUNCTION delta_optical_depth(eta, gamma_rel)

    ! Function that calcualtes the change to the optical depth
    ! c.f. Ridgers et al 2014 eqn (2)
    REAL(num) :: delta_optical_depth
    REAL(num), INTENT(IN) :: eta, gamma_rel
    REAL(num) :: hsokolov

    hsokolov = find_value_from_table_1d(eta, n_sample_h, log_hsokolov(:,1), &
        log_hsokolov(:,2))

    delta_optical_depth = dt * eta * alpha_f * SQRT(3.0_num) * hsokolov &
        / (2.0_num * pi * tau_c * gamma_rel)

  END FUNCTION delta_optical_depth



  FUNCTION delta_optical_depth_tri(eta, gamma_rel)

    ! Function that calcualtes the change to the optical depth
    REAL(num) :: delta_optical_depth_tri
    REAL(num), INTENT(IN) :: eta, gamma_rel
    REAL(num) :: omegahat

    omegahat = find_value_from_table_1d(eta, n_sample_t, log_omegahat(:,1), &
        log_omegahat(:,2))

    delta_optical_depth_tri = dt * eta * alpha_f**2 * 0.64_num * omegahat &
        / (2.0_num * pi * tau_c * gamma_rel) * boost_pairs

  END FUNCTION delta_optical_depth_tri



  FUNCTION delta_optical_depth_photon(chi_val, part_e)

    ! Function that calcualtes the change to the optical depth
    REAL(num) :: delta_optical_depth_photon
    REAL(num), INTENT(IN) :: chi_val, part_e
    REAL(num) :: tpair

    tpair = find_value_from_table_1d(chi_val, n_sample_t, log_tpair(:,1), &
        log_tpair(:,2))

    delta_optical_depth_photon = dt / tau_c * alpha_f / part_e * chi_val &
        * tpair * boost_pairs

  END FUNCTION delta_optical_depth_photon



  FUNCTION calculate_eta(part_x, part_y, part_ux, part_uy, part_uz, &
      gamma_rel, e_out, b_out)

    REAL(num) :: calculate_eta
    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(IN) :: part_ux, part_uy, part_uz, gamma_rel
    REAL(num) :: e_at_part(3), b_at_part(3)
    REAL(num) :: beta_x, beta_y, beta_z
    REAL(num) :: flperp(3), i_e, tau0, roland_eta
    REAL(num) :: lambdac, coeff_eta, moduclip, moduclip2, u_dot_e
    ! To be used by the spin routine
    REAL(num), INTENT(OUT), OPTIONAL :: e_out(3), b_out(3)    

    CALL field_at_particle(part_x, part_y, e_at_part, b_at_part)

    moduclip2 = MAX(part_ux**2 + part_uy**2 + part_uz**2, c_tiny)
    moduclip = SQRT(moduclip2)

    beta_x = part_ux / gamma_rel
    beta_y = part_uy / gamma_rel
    beta_z = part_uz / gamma_rel

    lambdac = h_bar / mc0

    coeff_eta = SQRT(3.0_num * lambdac / (2.0_num * alpha_f * m0 * c**3))

    u_dot_e = (part_ux * e_at_part(1) + part_uy * e_at_part(2) &
        + part_uz * e_at_part(3)) / moduclip2

    flperp(1) = q0 * (e_at_part(1) - u_dot_e * part_ux &
        + c * (beta_y * b_at_part(3) - beta_z * b_at_part(2)))

    flperp(2) = q0 * (e_at_part(2) - u_dot_e * part_uy &
        + c * (beta_z * b_at_part(1) - beta_x * b_at_part(3)))

    flperp(3) = q0 * (e_at_part(3) - u_dot_e * part_uz &
        + c * (beta_x * b_at_part(2) - beta_y * b_at_part(1)))

    ! Dipole emission intensity

    tau0 = q0**2 / (6.0_num * pi * epsilon0 * m0 * c**3)

    i_e = tau0 * gamma_rel**2 * (flperp(1)**2 + flperp(2)**2 + flperp(3)**2 &
        + (q0 * (beta_x * e_at_part(1) + beta_y * e_at_part(2) &
        + beta_z * e_at_part(3)) / moduclip)**2) / m0

    roland_eta = coeff_eta * SQRT(i_e)

    ! Determine eta from fields
    calculate_eta = roland_eta

    ! Used by the spin routine
    IF (PRESENT(e_out)) e_out = e_at_part
    IF (PRESENT(b_out)) b_out = b_at_part
  END FUNCTION calculate_eta



  FUNCTION calculate_chi(part_x, part_y, dir_x, dir_y, dir_z, part_e, &
      e_out, b_out)

    REAL(num) :: calculate_chi
    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(IN) :: dir_x, dir_y, dir_z, part_e
    REAL(num) :: e_at_part(3), b_at_part(3), q(3)
    REAL(num) :: e_dot_dir, calculate_chi_roland
    ! Used by the spin routine
    REAL(num), INTENT(OUT), OPTIONAL :: e_out(3), b_out(3)
    CALL field_at_particle(part_x, part_y, e_at_part, b_at_part)

    e_dot_dir = e_at_part(1) * dir_x + e_at_part(2) * dir_y &
        + e_at_part(3) * dir_z

    q(1) = e_at_part(1) - e_dot_dir * dir_x &
        + c * (dir_y * b_at_part(3) - dir_z * b_at_part(2))
    q(2) = e_at_part(2) - e_dot_dir * dir_y &
        + c * (dir_z * b_at_part(1) - dir_x * b_at_part(3))
    q(3) = e_at_part(3) - e_dot_dir * dir_z &
        + c * (dir_x * b_at_part(2) - dir_y * b_at_part(1))

    calculate_chi_roland = 0.5_num * SQRT(q(1)**2 + q(2)**2 + q(3)**2) &
        * part_e / e_s

    ! Determine chi from fields
    calculate_chi = calculate_chi_roland
    IF (PRESENT(e_out)) e_out = e_at_part
    IF (PRESENT(b_out)) b_out = b_at_part
  END FUNCTION calculate_chi



  SUBROUTINE field_at_particle(part_x, part_y, e_at_part, b_at_part)

    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(OUT) :: e_at_part(3), b_at_part(3)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x1, cell_x2, cell_y1, cell_y2

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r

    ! The fraction of a cell between the particle position and cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min:sf_max) :: hx, hy
    ! Temporary variables
    INTEGER :: dcellx, dcelly
    REAL(num) :: ex_part, ey_part, ez_part
    REAL(num) :: bx_part, by_part, bz_part

    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = (1.0_num)**c_ndims
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    ! Grid cell position as a fraction.
#ifdef PARTICLE_SHAPE_TOPHAT
    cell_x_r = part_x / dx - 0.5_num
    cell_y_r = part_y / dy - 0.5_num
#else
    cell_x_r = part_x / dx
    cell_y_r = part_y / dy
#endif
    ! Round cell position to nearest cell
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    ! Calculate fraction of cell between nearest cell boundary and particle
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1

    cell_y1 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y1, num) - cell_y_r
    cell_y1 = cell_y1 + 1

    ! Particle weight factors as described in the manual, page25
    ! These weight grid properties onto particles
    ! Also used to weight particle properties onto grid, used later
    ! to calculate J
    ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#else
#include "triangle/gx.inc"
#endif

    ! Now redo shifted by half a cell due to grid stagger.
    ! Use shifted version for ex in X, ey in Y, ez in Z
    ! And in Y&Z for bx, X&Z for by, X&Y for bz
    cell_x2 = FLOOR(cell_x_r)
    cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
    cell_x2 = cell_x2 + 1

    cell_y2 = FLOOR(cell_y_r)
    cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
    cell_y2 = cell_y2 + 1

    dcellx = 0
    dcelly = 0
    ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

    ! These are the electric and magnetic fields interpolated to the
    ! particle position. They have been checked and are correct.
    ! Actually checking this is messy.
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/e_part.inc"
#include "bspline3/b_part.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/e_part.inc"
#include "tophat/b_part.inc"
#else
#include "triangle/e_part.inc"
#include "triangle/b_part.inc"
#endif

    ! update particle momenta using weighted fields
    ! ex_part etc are NOT fields at particle, but fac times
    ! field

    e_at_part(1) = fac * ex_part
    e_at_part(2) = fac * ey_part
    e_at_part(3) = fac * ez_part

    b_at_part(1) = fac * bx_part
    b_at_part(2) = fac * by_part
    b_at_part(3) = fac * bz_part

    ! Uniform external fields (fields block, for tests)
    IF (use_external_fields) THEN
      e_at_part = e_at_part + external_e
      b_at_part = b_at_part + external_b
    END IF

  END SUBROUTINE field_at_particle



  SUBROUTINE generate_photon(generating_electron, iphoton, eta, zeta, &
      a_spin, e_field, b_field, charge)
    ! SPIN change:
    ! Generates a photon moving in same direction as electron
    ! With zeta and a_spin (SPIN, use_spin_qed), the photon energy 
    ! comes from the spin-resolved spectrum of the emitting lepton; 
    !  with the fields at the lepton and  its charge as well, 
    ! the lepton's final spin and the photon's
    ! polarisation are set by nlcs_final_states.

    TYPE(particle), POINTER :: generating_electron
    INTEGER, INTENT(IN) :: iphoton
    REAL(num), INTENT(IN) :: eta
    REAL(num), INTENT(IN), OPTIONAL :: zeta, a_spin
    REAL(num), INTENT(IN), OPTIONAL :: e_field(3), b_field(3), charge
    REAL(num) :: dir_x, dir_y, dir_z, mag_p, generating_gamma
    REAL(num) :: rand_temp, photon_energy
    REAL(num) :: g_eta, sync_force, sync_rate, taubar_c, beta
    REAL(num) :: rad_limit = 1.0_num
    LOGICAL :: samp_photon, make_photon
    TYPE(particle), POINTER :: new_photon
#ifdef SPIN
    REAL(num) :: u_mag, xi_stored(3), spin_avg(3), chi_photon
    LOGICAL :: spin_event, spin_down
#endif

    mag_p = MAX(SQRT(generating_electron%part_p(1)**2 &
        + generating_electron%part_p(2)**2 &
        + generating_electron%part_p(3)**2), c_tiny)

    dir_x = generating_electron%part_p(1) / mag_p
    dir_y = generating_electron%part_p(2) / mag_p
    dir_z = generating_electron%part_p(3) / mag_p

    generating_gamma = SQRT(1.0_num + (mag_p / m0 / c)**2)
#ifdef SPIN
    ! The final states use the lepton before recoil
    u_mag = mag_p / (m0 * c)
#endif

    ! Determine photon energy

    rand_temp = random()
#ifdef SPIN
    ! Use Int F^-1_{up,down} to sample momentum transfer
    IF (PRESENT(zeta) .AND. PRESENT(a_spin)) THEN
      ! W = W_up (1 + zeta)/2 + W_down (1 - zeta)/2: pick the spectrum
      ! by its weight at this step (exact, by linearity in zeta)
      spin_down = random() < nlcs_down_probability(zeta, a_spin)
      IF (use_nlcs_quantile_tables) THEN
        IF (spin_down) THEN
          chi_photon = nlcs_quantile_chi(rand_temp, eta, nlcs_q_down)
        ELSE
          chi_photon = nlcs_quantile_chi(rand_temp, eta, nlcs_q_up)
        END IF
        ! As in calculate_photon_energy
        photon_energy = (2.0_num * chi_photon / eta) * generating_gamma &
            * m0 * c**2
      ELSE IF (spin_down) THEN
        photon_energy = calculate_photon_energy(rand_temp, eta, &
            generating_gamma, nlcs_cdf_down)
      ELSE
        photon_energy = calculate_photon_energy(rand_temp, eta, &
            generating_gamma, nlcs_cdf_up)
      END IF
    ELSE ! If something goes wrong or use_spin_qed not activated
          ! use the old way
      photon_energy = calculate_photon_energy(rand_temp, eta, &
          generating_gamma)
    END IF
#else
    photon_energy = calculate_photon_energy(rand_temp, eta, generating_gamma)
#endif

    IF (use_radiation_reaction) THEN
      IF (use_continuous_emission) THEN
	! Calculate the momentum loss from the Landau-Lifshitz force
        ! As gaunt_g, written out here: with -ffast-math, calling
        ! gaunt_g instead changes the rounding of this force
        IF (use_classical_emission) THEN
          g_eta = 1.0_num
        ELSE IF (gaunt_factor == c_gaunt_esirkepov2015) THEN
          ! Esirkepov et al. (2015)
          g_eta = (1.0_num + eta * (18.0_num + eta * (69.0_num &
              + eta * (73.0_num + 5.806_num * eta))))**(-1.0_num / 3.0_num)
        ELSE
          ! Ridgers et al. (2014) -- the legacy default way
          g_eta = (1.0_num+4.8_num*(1.0_num+eta)*&
                log(1.0_num+1.7_num*eta)+2.44_num*eta*eta)**(-2.0_num/3.0_num)
        END IF
        taubar_c = h_bar/m0/c/c
        beta = SQRT(1.0_num - 1.0_num/generating_gamma/generating_gamma)
        sync_force = 2.0_num*alpha_f*m0*c/(3.0_num*taubar_c) &
                      * beta*eta*eta*g_eta

        ! Ensure particle can't give away more momentum than it starts with
        IF (mag_p > dt*sync_force) THEN
          mag_p = mag_p - dt*sync_force
        ELSE
          rad_limit = mag_p / (dt*sync_force)
          mag_p = 0.0_num
        END IF
      ELSE
        ! Calculate electron recoil from photon energy
        mag_p = mag_p - photon_energy / c
      END IF

      generating_electron%part_p(1) = dir_x * mag_p
      generating_electron%part_p(2) = dir_y * mag_p
      generating_electron%part_p(3) = dir_z * mag_p
    END IF

    samp_photon = .TRUE.
    IF (photon_sample_fraction < 1.0_num) THEN
      samp_photon = random() < photon_sample_fraction
    END IF

    ! This will only create photons that have energies above a user specified
    ! cutoff and if photon generation is turned on. E+/- recoil is always
    ! considered
    make_photon = photon_energy > photon_energy_min .AND. produce_photons &
        .AND. samp_photon

#ifdef SPIN
    ! Final lepton spin (always) and photon polarisation (if the photon
    ! is created); lam = photon energy / nitial lepton energy = 2 chi / eta
    spin_event = PRESENT(zeta) .AND. PRESENT(a_spin) &
        .AND. PRESENT(e_field) .AND. PRESENT(b_field) .AND. PRESENT(charge)
    IF (spin_event .AND. spin_nlcs == c_spin_process_photon_only) THEN
      ! Photon from an unpolarised lepton, summed over its final spin;
      ! the lepton's own spin is left as it is
      spin_avg = 0.0_num
      CALL nlcs_final_states(dir_x * u_mag, dir_y * u_mag, &
          dir_z * u_mag, generating_gamma, e_field, b_field, charge, &
          photon_energy / (generating_gamma * m0 * c**2), eta, make_photon, &
          spin_avg, xi_stored)
    ELSE IF (spin_event) THEN
      ! radiative_spin_flip = F: the final states are drawn as usual
      ! The lepton retains its initial spin.
      ! The photon's Stokes vector, conditional on the drawn final spin,
      ! averages to the one summed over the final spin.
      spin_avg = generating_electron%spin
      CALL nlcs_final_states(dir_x * u_mag, dir_y * u_mag, &
          dir_z * u_mag, generating_gamma, e_field, b_field, charge, &
          photon_energy / (generating_gamma * m0 * c**2), eta, make_photon, &
          generating_electron%spin, xi_stored)
      IF (.NOT. radiative_spin_flip) generating_electron%spin = spin_avg
    END IF
    ! photon_polarisation = F: every photon is unpolarised
    IF (.NOT. photon_polarisation) xi_stored = 0.0_num
#endif

    IF (make_photon) THEN
      IF (photon_energy < c_tiny) photon_energy = c_tiny

      CALL create_particle(new_photon)
      new_photon%part_pos = generating_electron%part_pos

      new_photon%part_p(1) = dir_x * photon_energy / c
      new_photon%part_p(2) = dir_y * photon_energy / c
      new_photon%part_p(3) = dir_z * photon_energy / c

      new_photon%optical_depth = reset_optical_depth()
      new_photon%particle_energy = photon_energy
#ifdef SPIN
      ! Stokes vector on the photon's reference axes (zero otherwise)
      IF (spin_event) new_photon%spin = xi_stored
#endif

      IF (use_continuous_emission) THEN
        ! Calculate photon weight from synchrotron emission rate
        sync_rate = delta_optical_depth(eta, generating_gamma) ! This already has *dt included
        new_photon%weight = generating_electron%weight * sync_rate &
                                                       / photon_sample_fraction
      ELSE
        new_photon%weight = generating_electron%weight / photon_sample_fraction
      ENDIF

      CALL add_particle_to_partlist(species_list(iphoton)%attached_list, &
          new_photon)
    END IF

  END SUBROUTINE generate_photon



  FUNCTION gaunt_g(eta)

    ! Gaunt factor g(eta) of the radiation reaction force: 1 for
    ! classical emission, a fit (qed block gaunt_factor) otherwise.
    ! generate_photon writes the same out (keep the two in step)
    REAL(num), INTENT(IN) :: eta
    REAL(num) :: gaunt_g

    IF (use_classical_emission) THEN
      gaunt_g = 1.0_num
    ELSE IF (gaunt_factor == c_gaunt_esirkepov2015) THEN
      ! Esirkepov et al. (2015)
      gaunt_g = (1.0_num + eta * (18.0_num + eta * (69.0_num &
          + eta * (73.0_num + 5.806_num * eta))))**(-1.0_num / 3.0_num)
    ELSE
      ! Ridgers et al. (2014) -- the legacy default way
      gaunt_g = (1.0_num+4.8_num*(1.0_num+eta)*&
            log(1.0_num+1.7_num*eta)+2.44_num*eta*eta)**(-2.0_num/3.0_num)
    END IF

  END FUNCTION gaunt_g



  SUBROUTINE ll_gather(part_x, part_y, e_at_part, b_at_part, grad)

    ! Fields at the particle and, for the gradient term, their space
    ! derivatives: grad(:,1:4) = dE/dx, dB/dx, dE/dy, dB/dy. These are
    ! central differences of the interpolated fields over one cell
    ! (d/dz = 0 in 2D); zero without the gradient term
    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(OUT) :: e_at_part(3), b_at_part(3), grad(3,4)
    REAL(num) :: e_p(3), b_p(3), e_m(3), b_m(3)

    CALL field_at_particle(part_x, part_y, e_at_part, b_at_part)
    grad = 0.0_num
    IF (.NOT. (use_full_landau_lifshitz .AND. ll_gradient_term)) RETURN

    CALL field_at_particle(part_x + 0.5_num * dx, part_y, e_p, b_p)
    CALL field_at_particle(part_x - 0.5_num * dx, part_y, e_m, b_m)
    grad(:,1) = (e_p - e_m) / dx
    grad(:,2) = (b_p - b_m) / dx
    CALL field_at_particle(part_x, part_y + 0.5_num * dy, e_p, b_p)
    CALL field_at_particle(part_x, part_y - 0.5_num * dy, e_m, b_m)
    grad(:,3) = (e_p - e_m) / dy
    grad(:,4) = (b_p - b_m) / dy

  END SUBROUTINE ll_gather



  SUBROUTINE ll_terms(v, gamma_rel, charge, g_eta, e_at_part, b_at_part, &
      grad, force)

    ! The field-gradient and F F u terms of the Landau-Lifshitz force,
    ! dp/dt, for velocity v (Landau and Lifshitz, Classical Theory of
    ! Fields, section 76), in SI with tau0 = q0^2 / (6 pi epsilon0 m0 c^3):
    !   gradient:  tau0 q gamma [D E/Dt + v x D B/Dt]
    !   quadratic: g tau0 (q^2 / m0) [E x B + B x (B x v) + E (v.E) / c^2]
    ! with D/Dt = d/dt + v.grad. The quadratic term is the space part of
    ! the term in F F u, divided by gamma. It carries the same g as the
    ! term in (F u)^2 u, so that the two together stay orthogonal to u
    ! (the force keeps u.u = 1); the gradient term, orthogonal to u by
    ! itself, has no g. The time derivatives assume source-free fields
    ! at the particle (dE/dt = c^2 curl B, dB/dt = -curl E). The uniform
    ! external fields enter the quadratic term only.
    REAL(num), INTENT(IN) :: v(3), gamma_rel, charge, g_eta
    REAL(num), INTENT(IN) :: e_at_part(3), b_at_part(3), grad(3,4)
    REAL(num), INTENT(OUT) :: force(3)
    REAL(num) :: de(3), db(3), tau0

    tau0 = q0**2 / (6.0_num * pi * epsilon0 * m0 * c**3)
    force = 0.0_num

    IF (ll_gradient_term) THEN
      ! d/dt from the source-free Maxwell equations ...
      de(1) = c**2 * grad(3,4)
      de(2) = -c**2 * grad(3,2)
      de(3) = c**2 * (grad(2,2) - grad(1,4))
      db(1) = -grad(3,3)
      db(2) = grad(3,1)
      db(3) = grad(1,3) - grad(2,1)
      ! ... plus v.grad
      de = de + v(1) * grad(:,1) + v(2) * grad(:,3)
      db = db + v(1) * grad(:,2) + v(2) * grad(:,4)

      force = tau0 * charge * gamma_rel * (de + ll_cross(v, db))
    END IF

    IF (ll_quadratic_term) THEN
      force = force + g_eta * tau0 * (charge**2 / m0) &
          * (ll_cross(e_at_part, b_at_part) &
          + ll_cross(b_at_part, ll_cross(b_at_part, v)) &
          + e_at_part * DOT_PRODUCT(v, e_at_part) / c**2)
    END IF

  CONTAINS

    FUNCTION ll_cross(a, b)

      REAL(num), INTENT(IN) :: a(3), b(3)
      REAL(num) :: ll_cross(3)

      ll_cross(1) = a(2) * b(3) - a(3) * b(2)
      ll_cross(2) = a(3) * b(1) - a(1) * b(3)
      ll_cross(3) = a(1) * b(2) - a(2) * b(1)

    END FUNCTION ll_cross

  END SUBROUTINE ll_terms



  SUBROUTINE rr_correction_force(part_x, part_y, part_ux, part_uy, &
      part_uz, gamma_rel, charge, eta, force)

    ! What qed_update_optical_depth adds, as dp/dt over the step, to the
    ! radiation reaction of generate_photon (classical and semiclassical
    ! emission). generate_photon applies the term in (F u)^2 u as an
    ! Euler kick from the momentum at the start of the step, p_start.

    ! rr_midpoint = F: the force is ll_terms at p_start (needs
    ! use_full_landau_lifshitz). rr_midpoint = T: the force turns that
    ! kick into a midpoint kick in the fields at the start of the
    ! step,
    !   p_half = p_start + (dt / 2) F(p_start)
    !   p_mid  = p_start + dt F(p_half)
    !   force  = (p_mid - p_euler) / dt
    ! F being the leading term plus, with use_full_landau_lifshitz,
    ! ll_terms; p_euler is generate_photon's kick, worked out here. (A
    ! kick applied in place of generate_photon's would change the code
    ! around that call, which changed the default results in the last
    ! bit under -ffast-math.) If the leading term would take more than
    ! |p_start| in one step (dt far too large), the Euler kick is kept,
    ! with ll_terms at p_start.

    ! To put it simple, F is the original first order solver (Euler kick)
    ! while T is an attempt to make RR 2nd order.
    ! That is, for F the error decreases linearly with Delta t 
    ! where delta t is the time step and for T quadratically  
    ! But note that atm T is not a true second-order. For this to be 
    ! true second order we need to place the RR routine in the middle of 
    ! the particle pusher. Note that the photon (QED) routine takes
    ! place at full time steps. To make it true second-order requires
    ! invasive modifications.
    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(IN) :: part_ux, part_uy, part_uz, gamma_rel, charge
    REAL(num), INTENT(IN) :: eta
    REAL(num), INTENT(OUT) :: force(3)
    REAL(num) :: e_at_part(3), b_at_part(3), grad(3,4)
    REAL(num) :: p_start(3), p_half(3), p_euler(3), v(3)
    REAL(num) :: f_start(3), f_ll_start(3), f_half(3), f_ll_half(3)
    REAL(num) :: mag_p, sync_start, sync_half, sync_euler, beta

    CALL ll_gather(part_x, part_y, e_at_part, b_at_part, grad)

    IF (.NOT. rr_midpoint) THEN
      v(1) = c * part_ux / gamma_rel
      v(2) = c * part_uy / gamma_rel
      v(3) = c * part_uz / gamma_rel
      CALL ll_terms(v, gamma_rel, charge, gaunt_g(eta), e_at_part, &
          b_at_part, grad, force)
      RETURN
    END IF

    p_start(1) = part_ux * mc0
    p_start(2) = part_uy * mc0
    p_start(3) = part_uz * mc0
    mag_p = MAX(SQRT(SUM(p_start**2)), c_tiny)

    CALL rr_force(p_start, f_start, f_ll_start, sync_start)
    p_half = p_start + 0.5_num * dt * f_start
    CALL rr_force(p_half, f_half, f_ll_half, sync_half)

    ! generate_photon's kick, with its eta (from calculate_eta)
    beta = SQRT(1.0_num - 1.0_num / (1.0_num + (mag_p / mc0)**2))
    sync_euler = 2.0_num * alpha_f * m0 * c &
        / (3.0_num * h_bar / m0 / c / c) * beta * eta**2 * gaunt_g(eta)

    IF (dt * MAX(sync_start, sync_half, sync_euler) < mag_p) THEN
      p_euler = p_start * (1.0_num - dt * sync_euler / mag_p)
      force = (p_start + dt * f_half - p_euler) / dt
    ELSE
      force = f_ll_start
    END IF

  CONTAINS

    SUBROUTINE rr_force(p, force, f_ll, sync_force)

      ! Radiation reaction force (all terms) and its extra
      ! Landau-Lifshitz part f_ll at momentum p, in the gathered fields;
      ! sync_force is the magnitude of the leading term
      REAL(num), INTENT(IN) :: p(3)
      REAL(num), INTENT(OUT) :: force(3), f_ll(3), sync_force
      REAL(num) :: u(3), v(3), gamma_p, mag, eta_p, g_p, f_perp(3)

      u = p / mc0
      gamma_p = SQRT(1.0_num + SUM(u**2))
      v = c * u / gamma_p
      mag = MAX(SQRT(SUM(p**2)), c_tiny)

      ! eta as calculate_eta, from the fields:
      ! (eta e_s)^2 = gamma^2 [(E + v x B)^2 - (v.E / c)^2]
      f_perp(1) = e_at_part(1) + v(2) * b_at_part(3) - v(3) * b_at_part(2)
      f_perp(2) = e_at_part(2) + v(3) * b_at_part(1) - v(1) * b_at_part(3)
      f_perp(3) = e_at_part(3) + v(1) * b_at_part(2) - v(2) * b_at_part(1)
      eta_p = gamma_p * SQRT(MAX(SUM(f_perp**2) &
          - (DOT_PRODUCT(v, e_at_part) / c)**2, 0.0_num)) / e_s
      g_p = gaunt_g(eta_p)

      sync_force = 2.0_num * alpha_f * m0 * c &
          / (3.0_num * h_bar / m0 / c / c) * (mag / mc0 / gamma_p) &
          * eta_p**2 * g_p
      force = -sync_force * p / mag

      f_ll = 0.0_num
      IF (use_full_landau_lifshitz) THEN
        CALL ll_terms(v, gamma_p, charge, g_p, e_at_part, b_at_part, &
            grad, f_ll)
        force = force + f_ll
      END IF

    END SUBROUTINE rr_force

  END SUBROUTINE rr_correction_force



  FUNCTION calculate_photon_energy(rand_seed, eta, generating_gamma, &
      p_table)

    ! p_table: a photon-energy CDF on the grid of p_photon_energy
    ! (default p_photon_energy, the spin-averaged ksi_sokolov.table)
    REAL(num) :: calculate_photon_energy
    REAL(num), INTENT(IN) :: rand_seed, eta, generating_gamma
    REAL(num), INTENT(IN), OPTIONAL :: p_table(:,:)
    REAL(num) :: chi_final

    ! Check whether SPIN is used. 
    IF (PRESENT(p_table)) THEN
      chi_final = sample_chi(p_table)
    ELSE
      chi_final = sample_chi(p_photon_energy)
    END IF

    calculate_photon_energy = (2.0_num * chi_final / eta) * generating_gamma &
        * m0 * c**2

  CONTAINS

    FUNCTION sample_chi(table)

      REAL(num) :: sample_chi
      REAL(num), INTENT(IN) :: table(n_sample_eta,n_sample_chi)
      REAL(num) :: eta_min, chi_tmp

      eta_min = 10.0_num**MINVAL(log_eta)
      ! In the classical case, always use the spectrum at minimum eta
      ! and extrapolate downwards with chi \propto eta^2
      IF (use_classical_emission .OR. (eta < eta_min)) THEN
        chi_tmp = find_value_from_table_alt(eta_min, rand_seed, &
            n_sample_eta, n_sample_chi, log_eta, log_chi, table)
        sample_chi = chi_tmp * (eta / eta_min)**2
      ELSE
        sample_chi = find_value_from_table_alt(eta, rand_seed, &
            n_sample_eta, n_sample_chi, log_eta, log_chi, table)
      END IF

    END FUNCTION sample_chi

  END FUNCTION calculate_photon_energy



  SUBROUTINE generate_pair(generating_photon, chi_val, iphoton, ielectron, &
      ipositron, xi_obs, a_nbw, e_field, b_field)

    ! Generates a pair moving in same direction as photon. With xi_obs
    ! and a_nbw (SPIN, use_spin_qed) the energy split comes from the
    ! polarisation-resolved spectrum; with the fields at the photon as
    ! well, the pair spins are set by nbw_final_states.
    ! The implementation is non-invasive
    TYPE(particle), POINTER :: generating_photon
    REAL(num), INTENT(IN) :: chi_val
    INTEGER, INTENT(IN) :: iphoton, ielectron, ipositron
    REAL(num), INTENT(IN), OPTIONAL :: xi_obs(3), a_nbw
    REAL(num), INTENT(IN), OPTIONAL :: e_field(3), b_field(3)
    REAL(num) :: dir_x, dir_y, dir_z, mag_p
    REAL(num) :: probability_split, epsilon_frac, norm
    TYPE(particle), POINTER :: new_electron, new_positron
#ifdef SPIN
    LOGICAL :: spin_event
#endif

    CALL create_particle(new_electron)
    CALL create_particle(new_positron)

    new_electron%part_pos = generating_photon%part_pos
    new_positron%part_pos = generating_photon%part_pos

    norm  = c / generating_photon%particle_energy
    dir_x = generating_photon%part_p(1) * norm
    dir_y = generating_photon%part_p(2) * norm
    dir_z = generating_photon%part_p(3) * norm

    ! Determine how to split the energy amoung e-/e+
    ! IS CHI HERE SAME AS ROLAND'S? DEFINED BSinT/B_s

    probability_split = random()

#ifdef SPIN
    spin_event = PRESENT(xi_obs) .AND. PRESENT(a_nbw) &
        .AND. PRESENT(e_field) .AND. PRESENT(b_field)
    IF (spin_event) THEN
      ! W = W_par (1 + xi1)/2 + W_perp (1 - xi1)/2: pick the spectrum by
      ! its weight (exact, by linearity in xi1); epsilon_frac is the
      ! electron fraction (the spectrum is symmetric)
      IF (random() < nbw_par_probability(xi_obs(1), a_nbw)) THEN
        epsilon_frac = find_value_from_table(chi_val, probability_split, &
            n_nbw_chi, n_nbw_split, nbw_log_chi, nbw_split, nbw_cdf_par)
      ELSE
        epsilon_frac = find_value_from_table(chi_val, probability_split, &
            n_nbw_chi, n_nbw_split, nbw_log_chi, nbw_split, nbw_cdf_perp)
      END IF
    ELSE IF (pair_tables == c_pair_tables_lcfa) THEN
      ! Unpolarised photon (spin off): by linearity in xi1 the mixture
      ! with par weight (1 - A)/2 is exactly the LCFA spectrum
      IF (random() < nbw_par_probability(0.0_num, &
          nbw_asymmetry(chi_val))) THEN
        epsilon_frac = find_value_from_table(chi_val, probability_split, &
            n_nbw_chi, n_nbw_split, nbw_log_chi, nbw_split, nbw_cdf_par)
      ELSE
        epsilon_frac = find_value_from_table(chi_val, probability_split, &
            n_nbw_chi, n_nbw_split, nbw_log_chi, nbw_split, nbw_cdf_perp)
      END IF
    ELSE
      epsilon_frac = find_value_from_table(chi_val, probability_split, &
          n_sample_chi2, n_sample_epsilon, log_chi2, epsilon_split, &
          p_energy)
    END IF
#else
    epsilon_frac = find_value_from_table(chi_val, probability_split, &
        n_sample_chi2, n_sample_epsilon, log_chi2, epsilon_split, p_energy)
#endif

    mag_p = MAX(generating_photon%particle_energy / c, c_tiny)

    new_electron%part_p(1) = epsilon_frac * mag_p * dir_x
    new_electron%part_p(2) = epsilon_frac * mag_p * dir_y
    new_electron%part_p(3) = epsilon_frac * mag_p * dir_z

    new_positron%part_p(1) = (1.0_num - epsilon_frac) * mag_p * dir_x
    new_positron%part_p(2) = (1.0_num - epsilon_frac) * mag_p * dir_y
    new_positron%part_p(3) = (1.0_num - epsilon_frac) * mag_p * dir_z

#ifdef SPIN
    ! Pair spins (zero, i.e. unpolarised, otherwise)
    IF (spin_event) CALL nbw_final_states(generating_photon%part_p, &
        e_field, b_field, xi_obs, chi_val, epsilon_frac, new_electron%spin, &
        new_positron%spin)
#endif

    new_electron%optical_depth = reset_optical_depth()
    new_positron%optical_depth = reset_optical_depth()

#ifdef TRIDENT_PHOTONS
    new_electron%optical_depth_tri = reset_optical_depth()
    new_positron%optical_depth_tri = reset_optical_depth()
#endif

    new_electron%weight = generating_photon%weight / boost_pairs
    new_positron%weight = generating_photon%weight / boost_pairs

    CALL add_particle_to_partlist(species_list(ielectron)%attached_list, &
        new_electron)
    CALL add_particle_to_partlist(species_list(ipositron)%attached_list, &
        new_positron)

    ! Remove photon
    IF (random() < 1.0_num / boost_pairs) THEN
      CALL remove_particle_from_partlist(species_list(iphoton)%attached_list, &
          generating_photon)
      DEALLOCATE(generating_photon)
    ELSE
      generating_photon%optical_depth = &
          generating_photon%optical_depth + reset_optical_depth()
    END IF

  END SUBROUTINE generate_pair



  SUBROUTINE generate_pair_tri(generating_electron, ielectron, ipositron)

    ! Generates a pair moving in same direction as photon
    TYPE(particle), POINTER :: generating_electron
    INTEGER, INTENT(IN) :: ielectron, ipositron
    TYPE(particle), POINTER :: new_electron, new_positron

    CALL create_particle(new_electron)
    CALL create_particle(new_positron)

    new_electron%part_pos = generating_electron%part_pos
    new_positron%part_pos = generating_electron%part_pos

    new_electron%part_p = 0.0_num
    new_positron%part_p = 0.0_num

    new_electron%optical_depth = reset_optical_depth()
    new_positron%optical_depth = reset_optical_depth()

#ifdef TRIDENT_PHOTONS
    new_electron%optical_depth_tri = reset_optical_depth()
    new_positron%optical_depth_tri = reset_optical_depth()
#endif

    new_electron%weight = generating_electron%weight / boost_pairs
    new_positron%weight = generating_electron%weight / boost_pairs

    CALL add_particle_to_partlist(species_list(ielectron)%attached_list, &
        new_electron)
    CALL add_particle_to_partlist(species_list(ipositron)%attached_list, &
        new_positron)

  END SUBROUTINE generate_pair_tri



  FUNCTION find_value_from_table_1d(x_in, nx, x, values)

    REAL(num) :: find_value_from_table_1d
    REAL(num), INTENT(IN) :: x_in
    INTEGER, INTENT(IN) :: nx
    REAL(num), INTENT(IN) :: x(nx), values(nx)
    REAL(num) :: fx, x_value, value_interp, xdif1, xdif2, xdifm
    INTEGER :: i1, i2, im
    LOGICAL, SAVE :: warning = .TRUE.

    x_value = LOG10(MAX(x_in,c_tiny))

    i1 = 1
    i2 = nx
    xdif1 = x(i1) - x_value
    xdif2 = x(i2) - x_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = x(im) - x_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fx = (x_value - x(i1)) / (x(i2) - x(i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_1d" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      END IF
    END IF

    value_interp = (1.0_num - fx) * values(i1) + fx * values(i2)

    find_value_from_table_1d = 10.0_num**value_interp

  END FUNCTION find_value_from_table_1d



  FUNCTION find_value_from_table_alt(x_in, p_value, nx, ny, x, y, p_table)

    REAL(num) :: find_value_from_table_alt
    REAL(num), INTENT(IN) :: x_in, p_value
    INTEGER, INTENT(IN) :: nx, ny
    REAL(num), INTENT(IN) :: x(nx), y(nx,ny), p_table(nx,ny)
    INTEGER :: ix, index_lt, index_gt, i1, i2, im
    REAL(num) :: fx, fp, y_lt, y_gt, x_value, y_interp, xdif1, xdif2, xdifm
    LOGICAL, SAVE :: warning = .TRUE.

    x_value = LOG10(x_in) ! The interpolation is performed in log10-space

    ! Scan through x to find correct row of table
    i1 = 1
    i2 = nx
    xdif1 = x(i1) - x_value
    xdif2 = x(i2) - x_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = x(im) - x_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fx = (x_value - x(i1)) / (x(i2) - x(i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      END IF
    END IF

    index_lt = i1
    index_gt = i2

    ix = index_lt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
    END IF

    y_lt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ix = index_gt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
    END IF

    y_gt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ! Interpolate in x

    y_interp = (1.0_num - fx) * y_lt + fx * y_gt

    find_value_from_table_alt = 10.0_num**y_interp

  END FUNCTION find_value_from_table_alt



  FUNCTION find_value_from_table(x_in, p_value, nx, ny, x, y, p_table)

    REAL(num) :: find_value_from_table
    REAL(num), INTENT(IN) :: x_in, p_value
    INTEGER, INTENT(IN) :: nx, ny
    REAL(num), INTENT(IN) :: x(nx), y(ny), p_table(nx,ny)
    INTEGER :: ix, index_lt, index_gt, i1, i2, im
    REAL(num) :: fx, fp, y_lt, y_gt, x_value, y_interp, xdif1, xdif2, xdifm
    LOGICAL, SAVE :: warning = .TRUE.

    x_value = LOG10(x_in)

    ! Scan through x to find correct row of table
    i1 = 1
    i2 = nx
    xdif1 = x(i1) - x_value
    xdif2 = x(i2) - x_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = x(im) - x_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fx = (x_value - x(i1)) / (x(i2) - x(i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      END IF
    END IF

    index_lt = i1
    index_gt = i2

    ix = index_lt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
    END IF

    y_lt = (1.0_num - fp) * y(i1) + fp * y(i2)

    ix = index_gt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
    END IF

    y_gt = (1.0_num - fp) * y(i1) + fp * y(i2)

    ! Interpolate in x

    y_interp = (1.0_num - fx) * y_lt + fx * y_gt

    find_value_from_table = y_interp

  END FUNCTION find_value_from_table



  SUBROUTINE do_binary_collisions

    INTEGER :: ispecies, jspecies
    INTEGER(i8) :: ix, iy
    TYPE(particle_list), POINTER :: p_list1, p_list2
    TYPE(particle_list) :: new_lbw_electrons, new_lbw_positrons

    !!! Linear Breit-Wheeler
    IF (use_LBW) THEN
      DO ispecies = 1, n_species
        IF (species_list(ispecies)%species_type /= c_species_id_photon) CYCLE
        DO jspecies = ispecies, n_species
          IF (species_list(jspecies)%species_type /= c_species_id_photon) CYCLE
          DO ix = 1, nx
            DO iy = 1, ny
              CALL create_empty_partlist(new_lbw_electrons)
              CALL create_empty_partlist(new_lbw_positrons)

              IF (ispecies == jspecies) THEN
                p_list1 => species_list(ispecies)%secondary_list(ix,iy)
                CALL linear_Breit_Wheeler_intra( &
                    p_list1, ispecies, ix, iy, &
                    new_lbw_electrons, new_lbw_positrons)
              ELSE
                p_list1 => species_list(ispecies)%secondary_list(ix,iy)
                p_list2 => species_list(jspecies)%secondary_list(ix,iy)
                CALL linear_Breit_Wheeler_inter( &
                    p_list1, p_list2, ispecies, jspecies, ix, iy, &
                    new_lbw_electrons, new_lbw_positrons)
              END IF

              IF (new_lbw_electrons%count > 0) THEN
                CALL append_partlist(species_list(lbw_electron_species &
                    )%secondary_list(ix,iy), new_lbw_electrons)
              END IF
              IF (new_lbw_positrons%count > 0) THEN
                CALL append_partlist(species_list(lbw_positron_species &
                    )%secondary_list(ix,iy), new_lbw_positrons)
              END IF
            END DO ! do iy = 1, ny
          END DO ! do ix = 1, nx
        END DO ! jspecies
      END DO ! ispecies
    END IF ! if use_LBW

  END SUBROUTINE do_binary_collisions



  SUBROUTINE linear_Breit_Wheeler_intra(p_list_i, ispe, &
    ixx, iyy, lbw_elec_list, lbw_posi_list)

    TYPE(particle_list), INTENT(INOUT) :: p_list_i
    INTEGER, INTENT(IN) :: ispe
    INTEGER(i8), INTENT(IN) :: ixx, iyy
    TYPE(particle_list), INTENT(INOUT) :: lbw_elec_list, lbw_posi_list
    INTEGER(i8) :: icount
    REAL(num) :: P_max, N_max, cost, sint, rand_phi
    INTEGER :: N_coll, iu, io
    TYPE(particle), POINTER :: current_i, current_j
    REAL(num) :: i_Pmax
    INTEGER :: N_scattered
    REAL(num), DIMENSION(3) :: moment_i, moment_j
    REAL(num) :: pi_abs, pj_abs, pipj, pipj_cos, kappa, en_com2
    REAL(num) :: weight_i, weight_j
    REAL(num) :: com_beta, sigma_lbw, sigma_coll, P_coll, bv2
#ifdef SPIN
    TYPE(lbw_pol_pair) :: pol
#endif
    REAL(num) :: Y_ij
    REAL(num), DIMENSION(2) :: lepton_pos
    REAL(num), DIMENSION(3) :: beta_v, n_v, p_phot_com
    REAL(num) :: gamma_v
    REAL(num), DIMENSION(3) :: e1, e2, e3
    REAL(num) :: p_lep_com_mag
    REAL(num), DIMENSION(3) :: p_lep_com, p_elec_lab, p_posi_lab
    TYPE(particle), POINTER :: new_electron, new_positron
    TYPE(particle), POINTER :: next_i, next_j

    ! If there aren't enough particles to collide, then don't bother
    icount = p_list_i%count
    IF (icount < 2) RETURN

    !!! Determine how many macro-macro pairings are up to be checked
    ! maximal possible collisional probability P_max
    P_max = sig2cdt_dV_lbw * max_weight(p_list_i)

    ! determine how many macro-particle pairs are up to collide (N_coll)
    N_max = P_max * icount * (icount - 1.0_num) / 2.0_num

    If (random()< (N_max-FLOOR(N_max))) THEN
      N_coll = CEILING(N_max)
    ELSE
      N_coll = FLOOR(N_max)
    END IF

    ! Check the collision of these N_coll pairings of macro-photon
    IF (N_coll <= 0) RETURN

    IF (N_coll > (icount*0.5_num)) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Too many linear Breit-Wheeler collisions'
        END DO
      END IF
      CALL abort_code(c_err_qed)
    END IF

    ! Shuffle particle list only if there are pairs to check
    CALL shuffle_particle_list_random(p_list_i)
    current_i => p_list_i%head
    current_j => current_i%next
    i_Pmax = 1.0_num/P_max
    N_scattered = 0

    DO WHILE (N_scattered < N_coll)
      !!! check if four-momentum conservation

      moment_i = current_i%part_p
      moment_j = current_j%part_p

      ! |p_i| and |p_j|
      pi_abs = current_i%particle_energy * inv_c
      pj_abs = current_j%particle_energy * inv_c

      ! |p_i|*|p_j|
      pipj = pi_abs * pj_abs

      ! |p_i|*|p_j|*(1-cos(theta))
      pipj_cos = pipj - DOT_PRODUCT(moment_i,moment_j)

      ! dimensionless photon energy squared in c.o.m.
      en_com2 = pipj_cos * 0.5_num * inv_mc0_sq

      IF (en_com2 < 1.0_num) THEN
        current_i => current_j%next
        IF (ASSOCIATED(current_i)) current_j => current_i%next
        N_scattered = N_scattered +1
        CYCLE
      END IF

      ! kinematic factor (1-cos(theta)), after the threshold test,
      ! which skips photons of zero energy (pipj = 0)
      kappa = pipj_cos / pipj

      !!! calculate collisional probability P_ij/P_max
      weight_i = current_i%weight
      weight_j = current_j%weight

      ! normalized lepton velocity in c.o.m. frame
      com_beta = SQRT(1.0_num - 1.0_num/en_com2)

      ! cross section for the linear BW (without kinematic factor)
      ! Notice this is a Lorentz invariant.
      ! Later we will re-use it in the c.o.m frame.
      sigma_lbw = lbw_cross_sec(com_beta)
      sigma_coll = sigma_lbw
#ifdef SPIN
      ! LBW_polarisation: the cross-section of the two Stokes vectors
      ! (pol is used again for the pair angles, the pair spins and the
      ! survivors). LBW_pol_cross_section = F (test switch) keeps the
      ! unpolarised cross-section.
      IF (LBW_polarisation) THEN
        CALL lbw_pol_setup(moment_i, moment_j, pi_abs, pj_abs, &
            current_i%spin, current_j%spin, com_beta, sigma_lbw, pol)
        IF (LBW_pol_cross_section) sigma_coll = pol%sig_pol
      END IF
#endif

      ! collisional probability after modification by P_max
      P_coll = sigma_coll * cdt_dV * MAX(weight_i, weight_j) * kappa * i_Pmax

      IF (random() > P_coll) THEN
        ! These two macro-photons do not collide (due to collisional probability)
        ! Now move pointer to next particle
        current_i => current_j%next
        IF (ASSOCIATED(current_i)) current_j=>current_i%next
        ! Increment counter
        N_scattered = N_scattered + 1
        CYCLE
      END IF

      ! Now, collide these two macro-photons.
      ! Pair yield
      Y_ij = MIN(weight_i, weight_j) * i_LBW_amp_factor

      ! Pair position
      lepton_pos = (current_i%part_pos + current_j%part_pos)*0.5_num

      ! Lorentz transformation velocity vector from lab to c.o.m. (norm. by c)
      beta_v = (moment_i + moment_j) / (pi_abs + pj_abs)

      ! Unit vector of this Lorentz transformation velocity; there is
      ! none for opposite momenta of equal size (the lab is the c.o.m.)
      bv2 = DOT_PRODUCT(beta_v, beta_v)
      n_v = 0.0_num
      IF (bv2 > 0.0_num) n_v = beta_v / SQRT(bv2)

      ! Gamma factor of this Lorentz transformatin
      gamma_v = 1.0_num/SQRT(1.0_num-DOT_PRODUCT(beta_v,beta_v))

      ! Magnitude of lepton momentum in c.o.m. (in SI)
      p_lep_com_mag = SQRT(en_com2 - 1.0_num)*mc0

      IF (use_LBW_diff) THEN
        ! Photon momentum in c.o.m. (in SI)
        p_phot_com = moment_i + (gamma_v-1.0_num) &
            * DOT_PRODUCT(moment_i,n_v)*n_v - gamma_v*pi_abs*beta_v

        ! Orthonormal basis (e1, e2, e3) s.t. e1//photon momentum in c.o.m.
        CALL get_orthonormal(p_phot_com, e1, e2, e3)

        DO
          ! Random azimuthal angle in c.o.m.
          rand_phi = 2.0_num * pi * random()

          ! Random polar angle in c.o.m. (cosine of this random angle)
          cost = random_polar_lbw(com_beta, sigma_lbw)
#ifdef SPIN
          ! LBW_polarisation: keep these (unpolarised) angles with the
          ! ratio of polarised to unpolarised distributions (rejection).
          ! LBW_pol_angles = F (test switch) keeps the first draw.
          IF (.NOT. LBW_polarisation .OR. .NOT. LBW_pol_angles) EXIT
          IF (lbw_pol_keep_angles(pol, com_beta, cost, rand_phi)) EXIT
#else
          EXIT
#endif
        END DO
        sint = SQRT(1.0_num - cost**2)

        ! Lepton momentum in c.o.m. (in SI)
        p_lep_com = p_lep_com_mag &
            * (e1*cost + e2*sint*COS(rand_phi) + e3*sint*SIN(rand_phi))
      ELSE
        ! Uniform distribution on sphere surface

        ! Random azimuthal angle in c.o.m.
        rand_phi = 2.0_num * pi * random()

        ! Random polar angle
        cost = 2.0_num * random() - 1.0_num
        sint = SQRT(1.0_num - cost**2)

        ! Lepton momentum in c.o.m. (in SI)
        p_lep_com(1) = p_lep_com_mag * sint * COS(rand_phi)
        p_lep_com(2) = p_lep_com_mag * sint * SIN(rand_phi)
        p_lep_com(3) = p_lep_com_mag * cost
      END IF

      ! Electron momentum in lab
      p_elec_lab = p_lep_com + (gamma_v-1.0_num) * DOT_PRODUCT( &
          p_lep_com,n_v)*n_v + gamma_v*SQRT(en_com2)*mc0*beta_v

      ! Positron momentum in lab
      p_posi_lab = - p_lep_com - (gamma_v-1.0_num) * DOT_PRODUCT( &
          p_lep_com,n_v)*n_v + gamma_v*SQRT(en_com2)*mc0*beta_v

      ! Create linear Breit-Wheeler electron and positron
      CALL create_particle(new_electron)
      new_electron%weight   = Y_ij
      new_electron%part_pos = lepton_pos
      new_electron%part_p   = p_elec_lab
      new_electron%optical_depth = reset_optical_depth()
      CALL add_particle_to_partlist(lbw_elec_list, new_electron)

      CALL create_particle(new_positron)
      new_positron%weight   = Y_ij
      new_positron%part_pos = lepton_pos
      new_positron%part_p   = p_posi_lab
      new_positron%optical_depth = reset_optical_depth()
      CALL add_particle_to_partlist(lbw_posi_list, new_positron)

#ifdef SPIN
      ! LBW_pair_spin: spins of the pair at these angles, from both
      ! Stokes vectors (otherwise unpolarised leptons, spin 0)
      IF (LBW_polarisation .AND. LBW_pair_spin) THEN
        CALL lbw_pol_pair_spins(pol, current_i%part_p, &
            current_j%part_p, pi_abs, pj_abs, com_beta, cost, rand_phi, &
            new_electron%part_p, new_positron%part_p, &
            new_electron%spin, new_positron%spin)
      END IF

      ! LBW_survivor_update: Stokes vectors of the photons left in each
      ! macro-photon (f = 1, no update, for one that is deleted below).
      ! With LBW_pol_cross_section = F the photons are absorbed whatever
      ! their polarisation, so the survivors keep their Stokes vector.
      IF (LBW_polarisation .AND. LBW_survivor_update &
          .AND. LBW_pol_cross_section) THEN
        CALL lbw_pol_survivor(pol, .TRUE., Y_ij / weight_i, current_i%spin)
        CALL lbw_pol_survivor(pol, .FALSE., Y_ij / weight_j, current_j%spin)
      END IF
#endif

      ! Annihilate photons
      ! First, store pointers, as photons may be deleted
      next_i => current_j%next
      IF (ASSOCIATED(next_i)) next_j => next_i%next

      ! now, annihilate photons by decreasing weight, or delete macro-photon
      IF (lbw_amp_factor > 1.0_num) THEN
        ! In this case, no macro-photons are deleted
        current_i%weight = weight_i - Y_ij
        current_j%weight = weight_j - Y_ij
      ELSE
        ! no amplification of the cross-section.
        ! In this case, need to delete the smaller macro-photon,
        ! or delete both, if they have same weight
        IF ((weight_i/weight_j) >1.000001_num) THEN ! photon i is larger
          current_i%weight = weight_i - Y_ij
          CALL remove_particle_from_partlist(species_list(ispe &
              )%secondary_list(ixx,iyy), current_j)
          DEALLOCATE(current_j)
        ELSE IF ((weight_j/weight_i) >1.000001_num) THEN ! photon j is larger
          current_j%weight = weight_j - Y_ij
          CALL remove_particle_from_partlist(species_list(ispe &
              )%secondary_list(ixx,iyy), current_i)
          DEALLOCATE(current_i)
        ELSE ! photon i and j have same weight
          CALL remove_particle_from_partlist(species_list(ispe &
              )%secondary_list(ixx,iyy), current_i)
          DEALLOCATE(current_i)
            CALL remove_particle_from_partlist(species_list(ispe &
                )%secondary_list(ixx,iyy), current_j)
          DEALLOCATE(current_j)
        END IF ! comparing weights
      END IF ! if amplification_factor > 1.0_num

      ! Now, collision of these two macro photons is done,
      ! Move pointer to next particle
      current_i => next_i
      current_j => next_j

      ! Now, collision done, and pointer moved
      ! Increment counter
      N_scattered = N_scattered + 1
    END DO ! While more to scatter

  END SUBROUTINE linear_Breit_Wheeler_intra



  SUBROUTINE linear_Breit_Wheeler_inter(p_list_i, p_list_j, ispe, jspe, &
    ixx, iyy, lbw_elec_list, lbw_posi_list)

    TYPE(particle_list), INTENT(INOUT) :: p_list_i, p_list_j
    INTEGER, INTENT(IN) :: ispe, jspe
    INTEGER(i8), INTENT(IN) :: ixx, iyy
    TYPE(particle_list), INTENT(INOUT) :: lbw_elec_list, lbw_posi_list
    INTEGER(i8) :: icount, jcount
    REAL(num) :: q_i, q_j, P_max, N_max, cost, sint, rand_phi
    INTEGER :: N_coll, iu, io
    TYPE(particle), POINTER :: current_i, current_j
    REAL(num) :: i_Pmax
    INTEGER :: N_scattered
    REAL(num), DIMENSION(3) :: moment_i, moment_j
    REAL(num) :: pi_abs, pj_abs, pipj, pipj_cos, kappa, en_com2
    REAL(num) :: weight_i, weight_j
    REAL(num) :: com_beta, sigma_lbw, sigma_coll, P_coll, bv2
#ifdef SPIN
    TYPE(lbw_pol_pair) :: pol
#endif
    REAL(num) :: Y_ij
    REAL(num), DIMENSION(2) :: lepton_pos
    REAL(num), DIMENSION(3) :: beta_v, n_v, p_phot_com
    REAL(num) :: gamma_v
    REAL(num), DIMENSION(3) :: e1, e2, e3
    REAL(num) :: p_lep_com_mag
    REAL(num), DIMENSION(3) :: p_lep_com, p_elec_lab, p_posi_lab
    TYPE(particle), POINTER :: new_electron, new_positron
    TYPE(particle), POINTER :: next_i, next_j

    ! If there aren't enough particles to collide, then don't bother
    icount = p_list_i%count
    jcount = p_list_j%count
    IF (icount < 1 .OR. jcount < 1) RETURN

    !!! Determine how many macro-macro pairings are up to be checked
    ! maximal possible collisional probability P_max
    q_i   = max_weight(p_list_i)
    q_j   = max_weight(p_list_j)
    P_max = sig2cdt_dV_lbw * MAX(q_i, q_j)

    ! determine how many macro-particle pairs are up to collide (N_coll)
    N_max = P_max * icount * jcount

    If (random()< (N_max-FLOOR(N_max))) THEN
      N_coll = CEILING(N_max)
    ELSE
      N_coll = FLOOR(N_max)
    END IF

    ! Check the collision of these N_coll pairings of macro-photon
    IF (N_coll <= 0) RETURN

    IF (N_coll > MIN(icount, jcount)) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Too many linear Breit-Wheeler collisions'
        END DO
      END IF
      CALL abort_code(c_err_qed)
    END IF

    ! Shuffle particle list only if there are pairs to check
    CALL shuffle_particle_list_random(p_list_i)
    CALL shuffle_particle_list_random(p_list_j)
    current_i => p_list_i%head
    current_j => p_list_j%head
    i_Pmax = 1.0_num/P_max
    N_scattered = 0

    DO WHILE (N_scattered < N_coll)
      !!! check if four-momentum conservation

      moment_i = current_i%part_p
      moment_j = current_j%part_p

      ! |p_i| and |p_j|
      pi_abs = current_i%particle_energy * inv_c
      pj_abs = current_j%particle_energy * inv_c

      ! |p_i|*|p_j|
      pipj = pi_abs * pj_abs

      ! |p_i|*|p_j|*(1-cos(theta))
      pipj_cos = pipj - DOT_PRODUCT(moment_i,moment_j)

      ! dimensionless photon energy squared in c.o.m.
      en_com2 = pipj_cos * 0.5_num * inv_mc0_sq

      IF (en_com2 < 1.0_num) THEN
        current_i => current_i%next
        current_j => current_j%next
        N_scattered = N_scattered +1
        CYCLE
      END IF

      ! kinematic factor (1-cos(theta)), after the threshold test,
      ! which skips photons of zero energy (pipj = 0)
      kappa = pipj_cos / pipj

      !!! calculate collisional probability P_ij/P_max
      weight_i = current_i%weight
      weight_j = current_j%weight

      ! normalized lepton velocity in c.o.m. frame
      com_beta = SQRT(1.0_num - 1.0_num/en_com2)

      ! cross section for the linear BW (without kinematic factor)
      ! Notice this is a Lorentz invariant.
      ! Later we will re-use it in the c.o.m frame.
      sigma_lbw = lbw_cross_sec(com_beta)
      sigma_coll = sigma_lbw
#ifdef SPIN
      ! LBW_polarisation: the cross-section of the two Stokes vectors
      ! (pol is used again for the pair angles, the pair spins and the
      ! survivors). LBW_pol_cross_section = F (test switch) keeps the
      ! unpolarised cross-section.
      IF (LBW_polarisation) THEN
        CALL lbw_pol_setup(moment_i, moment_j, pi_abs, pj_abs, &
            current_i%spin, current_j%spin, com_beta, sigma_lbw, pol)
        IF (LBW_pol_cross_section) sigma_coll = pol%sig_pol
      END IF
#endif

      ! collisional probability after modification by P_max
      P_coll = sigma_coll * cdt_dV * MAX(weight_i, weight_j) * kappa * i_Pmax

      IF (random() > P_coll) THEN
        ! These two macro-photons do not collide (due to collisional probability)
        ! Now move pointer to next particle
        current_i => current_i%next
        current_j => current_j%next
        ! Increment counter
        N_scattered = N_scattered + 1
        CYCLE
      END IF

      ! Now, collide these two macro-photons.
      ! Pair yield
      Y_ij = MIN(weight_i, weight_j) * i_LBW_amp_factor

      ! Pair position
      lepton_pos = (current_i%part_pos + current_j%part_pos)*0.5_num

      ! Lorentz transformation velocity vector from lab to c.o.m. (norm. by c)
      beta_v = (moment_i + moment_j) / (pi_abs + pj_abs)

      ! Unit vector of this Lorentz transformation velocity; there is
      ! none for opposite momenta of equal size (the lab is the c.o.m.)
      bv2 = DOT_PRODUCT(beta_v, beta_v)
      n_v = 0.0_num
      IF (bv2 > 0.0_num) n_v = beta_v / SQRT(bv2)

      ! Gamma factor of this Lorentz transformatin
      gamma_v = 1.0_num/SQRT(1.0_num-DOT_PRODUCT(beta_v,beta_v))

      ! Magnitude of lepton momentum in c.o.m. (in SI)
      p_lep_com_mag = SQRT(en_com2 - 1.0_num)*mc0

      IF (use_LBW_diff) THEN
        ! Photon momentum in c.o.m. (in SI)
        p_phot_com = moment_i + (gamma_v-1.0_num) &
            * DOT_PRODUCT(moment_i,n_v)*n_v - gamma_v*pi_abs*beta_v

        ! Orthonormal basis (e1, e2, e3) s.t. e1//photon momentum in c.o.m.
        CALL get_orthonormal(p_phot_com, e1, e2, e3)

        DO
          ! Random azimuthal angle in c.o.m.
          rand_phi = 2.0_num * pi * random()

          ! Random polar angle in c.o.m. (cosine of this random angle)
          cost = random_polar_lbw(com_beta, sigma_lbw)
#ifdef SPIN
          ! LBW_polarisation: keep these (unpolarised) angles with the
          ! ratio of polarised to unpolarised distributions (rejection).
          ! LBW_pol_angles = F (test switch) keeps the first draw.
          IF (.NOT. LBW_polarisation .OR. .NOT. LBW_pol_angles) EXIT
          IF (lbw_pol_keep_angles(pol, com_beta, cost, rand_phi)) EXIT
#else
          EXIT
#endif
        END DO
        sint = SQRT(1.0_num - cost**2)

        ! Lepton momentum in c.o.m. (in SI)
        p_lep_com = p_lep_com_mag &
            * (e1*cost + e2*sint*COS(rand_phi) + e3*sint*SIN(rand_phi))
      ELSE
        ! Uniform distribution on sphere surface

        ! Random azimuthal angle in c.o.m.
        rand_phi = 2.0_num * pi * random()

        ! Random polar angle
        cost = 2.0_num * random() - 1.0_num
        sint = SQRT(1.0_num - cost**2)

        ! Lepton momentum in c.o.m. (in SI)
        p_lep_com(1) = p_lep_com_mag * sint * COS(rand_phi)
        p_lep_com(2) = p_lep_com_mag * sint * SIN(rand_phi)
        p_lep_com(3) = p_lep_com_mag * cost
      END IF

      ! Electron momentum in lab
      p_elec_lab = p_lep_com + (gamma_v-1.0_num) * DOT_PRODUCT( &
          p_lep_com,n_v)*n_v + gamma_v*SQRT(en_com2)*mc0*beta_v

      ! Positron momentum in lab
      p_posi_lab = - p_lep_com - (gamma_v-1.0_num) * DOT_PRODUCT( &
          p_lep_com,n_v)*n_v + gamma_v*SQRT(en_com2)*mc0*beta_v

      ! Create linear Breit-Wheeler electron and positron
      CALL create_particle(new_electron)
      new_electron%weight   = Y_ij
      new_electron%part_pos = lepton_pos
      new_electron%part_p   = p_elec_lab
      new_electron%optical_depth = reset_optical_depth()
      CALL add_particle_to_partlist(lbw_elec_list, new_electron)

      CALL create_particle(new_positron)
      new_positron%weight   = Y_ij
      new_positron%part_pos = lepton_pos
      new_positron%part_p   = p_posi_lab
      new_positron%optical_depth = reset_optical_depth()
      CALL add_particle_to_partlist(lbw_posi_list, new_positron)

#ifdef SPIN
      ! LBW_pair_spin: spins of the pair at these angles, from both
      ! Stokes vectors (otherwise unpolarised leptons, spin 0)
      IF (LBW_polarisation .AND. LBW_pair_spin) THEN
        CALL lbw_pol_pair_spins(pol, current_i%part_p, &
            current_j%part_p, pi_abs, pj_abs, com_beta, cost, rand_phi, &
            new_electron%part_p, new_positron%part_p, &
            new_electron%spin, new_positron%spin)
      END IF

      ! LBW_survivor_update: Stokes vectors of the photons left in each
      ! macro-photon (f = 1, no update, for one that is deleted below).
      ! With LBW_pol_cross_section = F the photons are absorbed whatever
      ! their polarisation, so the survivors keep their Stokes vector.
      IF (LBW_polarisation .AND. LBW_survivor_update &
          .AND. LBW_pol_cross_section) THEN
        CALL lbw_pol_survivor(pol, .TRUE., Y_ij / weight_i, current_i%spin)
        CALL lbw_pol_survivor(pol, .FALSE., Y_ij / weight_j, current_j%spin)
      END IF
#endif

      ! Annihilate photons
      ! First, store pointers, as photons may be deleted
      next_i => current_i%next
      next_j => current_j%next

      ! Now, annihilate photons by decreasing weight, or delete macro-photon
      IF (lbw_amp_factor > 1.0_num) THEN
        ! In this case, no macro-photons are deleted
        current_i%weight = weight_i - Y_ij
        current_j%weight = weight_j - Y_ij
      ELSE
        ! no amplification of the cross-section.
        ! In this case, need to delete the smaller macro-photon,
        ! or delete both, if they have same weight
        IF ((weight_i/weight_j) >1.000001_num) THEN ! photon i is larger
          current_i%weight = weight_i - Y_ij
          CALL remove_particle_from_partlist(species_list(jspe &
              )%secondary_list(ixx,iyy), current_j)
          DEALLOCATE(current_j)
        ELSE IF ((weight_j/weight_i) >1.000001_num) THEN ! photon j is larger
          current_j%weight = weight_j - Y_ij
          CALL remove_particle_from_partlist(species_list(ispe &
              )%secondary_list(ixx,iyy), current_i)
          DEALLOCATE(current_i)
        ELSE ! photon i and j have same weight, delete both
          CALL remove_particle_from_partlist(species_list(ispe &
              )%secondary_list(ixx,iyy), current_i)
          DEALLOCATE(current_i)
          CALL remove_particle_from_partlist(species_list(jspe &
              )%secondary_list(ixx,iyy), current_j)
          DEALLOCATE(current_j)
        END IF ! comparing weights
      END IF ! if lbw_amp_factor > 1.0

      ! Now, collision of these two macro photons is done,
      ! Move pointer to next particle
      current_i => next_i
      current_j => next_j

      ! Now, collision done, and pointer moved
      ! Increment counter
      N_scattered = N_scattered + 1
    END DO ! While more to scatter

  END SUBROUTINE linear_Breit_Wheeler_inter



  FUNCTION max_weight(p_list)

    TYPE(particle_list), INTENT(IN) :: p_list
    TYPE(particle), POINTER :: current
    REAL(num) :: max_weight

    max_weight = 0.0_num

    current => p_list%head

    DO WHILE (ASSOCIATED(current))
      max_weight = MAX(max_weight, current%weight)
      current => current%next
    END DO

  END FUNCTION max_weight



  SUBROUTINE get_orthonormal(vec, e1, e2, e3)

    ! output a set of 3D orthonormal basis (e1,e2,e3) s.t. e1 // vec
    REAL(num), DIMENSION(3), INTENT(IN) :: vec
    REAL(num), DIMENSION(3), INTENT(OUT) :: e1, e2, e3

    e1 = vec/ SQRT(DOT_PRODUCT(vec, vec))

    IF (ABS(e1(2))<=c_tiny .AND. ABS(e1(3))<=c_tiny) THEN
      ! if e1 // (1,0,0), cross product e1 with (0,1,0) to get e2
      e2(1) = -e1(3)
      e2(2) = 0
      e2(3) = e1(1)
    ELSE
      ! cross product e1 with (1,0,0) to get e2
      e2(1) = 0
      e2(2) = e1(3)
      e2(3) = -e1(2)
    END IF

    e2 = e2 / SQRT(DOT_PRODUCT(e2, e2))

    e3(1) = e1(2)*e2(3) - e1(3)*e2(2)
    e3(2) = e1(3)*e2(1) - e1(1)*e2(3)
    e3(3) = e1(1)*e2(2) - e1(2)*e2(1)

  END SUBROUTINE get_orthonormal



  FUNCTION lbw_cross_sec(ze)

    REAL(num), INTENT(IN) :: ze
    REAL(num) :: lbw_cross_sec
    ! input ze is dimensionless lepton velocity in c.o.m.

    lbw_cross_sec = half_pire2 * (1.0_num - ze**2) &
      * (-2.0_num*ze*(2.0_num-ze**2) + (3.0_num-ze**4) &
      * LOG((1.0_num+ze)/(1.0_num-ze)))

  END FUNCTION lbw_cross_sec



  FUNCTION random_polar_lbw(v, sigma)

  ! This function generates (the cosine of) a random polar angle
  ! following the distribution of the differential cross section
  ! of the lbw process in the c.o.m. frame, by inverse transform
  ! sampling method using bisection method.
  REAL(num), INTENT(IN) :: v, sigma
  REAL(num) :: rnd_cdf
  REAL(num) :: lb, ub, mid_point
  REAL(num) :: cdf_err_lb, cdf_err_mp
  REAL(num) :: random_polar_lbw

  rnd_cdf = random()

  ! upper bound, lower bound, and mid point of bisection method.
  lb = -1.0_num
  ub = 1.0_num
  mid_point = 0.0_num

  cdf_err_lb = lbw_polar_cdf_err(lb,        rnd_cdf, v, sigma)
  cdf_err_mp = lbw_polar_cdf_err(mid_point, rnd_cdf, v, sigma)

  DO WHILE (ABS(cdf_err_mp) > tolerance_cdf .AND. &
      (ABS(ub-lb)>tolerance_cos_angle))
    IF (ABS(SIGN(1.0_num, cdf_err_lb) - SIGN(1.0_num, cdf_err_mp)) &
        < 0.5_num) THEN
      ! lb and mp have same sign
      ! meaning correct value is between mp and ub
      lb = mid_point
      cdf_err_lb = lbw_polar_cdf_err(lb, rnd_cdf, v, sigma)
    ELSE
      ! lb and mp have different sign
      ! meaning correct value is between lb and mp
        ub = mid_point
    END IF
    mid_point = (lb+ub)*0.5_num
    cdf_err_mp = lbw_polar_cdf_err(mid_point, rnd_cdf, v, sigma)
  END DO

  random_polar_lbw = mid_point

  END FUNCTION random_polar_lbw


  FUNCTION lbw_polar_cdf_err(mu, rnd, b, sig)

    ! error relative to the random cummulative distribution function
    ! for the polar angle of lbw differential cross-section
    REAL(num), INTENT(IN) :: mu, rnd, b, sig
    REAL(num) :: temp, big_bracket
    REAL(num) :: term1, term2, term3, term4
    REAL(num) :: lbw_polar_cdf_err

    temp  = 1.0_num - b**2
    term1 = 1.0_num + b * mu
    term2 = 1.0_num - b * mu
    term3 = 1.0_num + b
    term4 = 1.0_num - b

    big_bracket = 2*b*(1.0_num-mu) + (3.0_num-b**4) &
      * LOG(term1*term4/term2/term3) + temp**2 &
      * (1.0_num/term1 - 1.0_num/term2 - 1.0_num/term3 + 1.0_num/term4)

    ! the minus sign come from the definition for solid angle Omega in
    ! eq. (1) of ppcf. 60. 104001 (2018) by Ribeyre et. al.
    lbw_polar_cdf_err = - quarter_pire2 * temp * big_bracket / sig - rnd

  END FUNCTION lbw_polar_cdf_err



#ifdef SPIN
  SUBROUTINE lbw_pol_setup(p_i, p_j, pabs_i, pabs_j, spin_i, spin_j, ze, &
      sigma0, pol)

    ! LBW_polarisation: the cross-section of a candidate pair of photons
    ! (lab momenta p_i, p_j of size pabs_i, pabs_j; stored Stokes
    ! vectors spin_i, spin_j) and the frame of lbw_pol_pair, with the
    ! boost of the collision routines. ze is the lepton velocity in the
    ! CM frame and sigma0 = lbw_cross_sec(ze). With xi_i, xi_j on the
    ! collision axes, the pair spins summed,
    !   sigma = sigma_0 + sigma_L (xi_i1 xi_j1 - xi_i3 xi_j3)
    !         + sigma_C xi_i2 xi_j2,
    !   sigma_L = -(pi r_e^2 / 2) (1 - ze^2)^2 (ze + Theta (1 - ze^2)),
    !   sigma_C = pi r_e^2 (1 - ze^2) (3 ze - 2 Theta),
    ! where 1 - ze^2 = 1 / C^2 and Theta = atanh(ze). Photons linearly
    ! polarised along one axis give sigma_0 + sigma_L, along
    ! perpendicular axes sigma_0 - sigma_L; equal helicities give
    ! sigma_0 + sigma_C, opposite ones sigma_0 - sigma_C (Breit and
    ! Wheeler, Phys. Rev. 46, 1087 (1934); Zhao et al., Phys. Rev. D
    ! 107, 096013 (2023)). sigma does not depend on the choice of t1.
    REAL(num), INTENT(IN) :: p_i(3), p_j(3), pabs_i, pabs_j
    REAL(num), INTENT(IN) :: spin_i(3), spin_j(3), ze, sigma0
    TYPE(lbw_pol_pair), INTENT(OUT) :: pol
    REAL(num) :: beta_v(3), n_v(3), p_com(3), bv2, gamma_v, ic2, theta

    pol%polarised = MAXVAL(ABS(spin_i)) > 0.0_num &
        .OR. MAXVAL(ABS(spin_j)) > 0.0_num
    pol%sig_pol = sigma0
    IF (.NOT. pol%polarised) RETURN

    ! Boost to the CM frame (none for opposite momenta of equal size)
    ! and photon i's direction there
    beta_v = (p_i + p_j) / (pabs_i + pabs_j)
    bv2 = DOT_PRODUCT(beta_v, beta_v)
    n_v = 0.0_num
    IF (bv2 > 0.0_num) n_v = beta_v / SQRT(bv2)
    gamma_v = 1.0_num / SQRT(1.0_num - bv2)
    p_com = p_i + (gamma_v - 1.0_num) * DOT_PRODUCT(p_i, n_v) * n_v &
        - gamma_v * pabs_i * beta_v
    CALL get_orthonormal(p_com, pol%q, pol%t1, pol%t2)

    ! The stored Stokes vectors on the collision axes
    CALL photon_cm_axes(p_i, beta_v, gamma_v, pol%q, pol%ri1, pol%ri2)
    CALL photon_cm_axes(p_j, beta_v, gamma_v, -pol%q, pol%rj1, pol%rj2)
    pol%xi_i = pol_stored_to_obs(pol%ri1, pol%ri2, pol%t1, spin_i)
    pol%xi_j = pol_stored_to_obs(pol%rj1, pol%rj2, pol%t1, spin_j)

    ic2 = 1.0_num - ze**2
    theta = 0.5_num * LOG((1.0_num + ze) / (1.0_num - ze))
    pol%sig_0 = sigma0
    pol%sig_l = -half_pire2 * ic2**2 * (ze + theta * ic2)
    pol%sig_c = 2.0_num * half_pire2 * ic2 &
        * (3.0_num * ze - 2.0_num * theta)
    pol%sig_pol = sigma0 + pol%sig_l * (pol%xi_i(1) * pol%xi_j(1) &
        - pol%xi_i(3) * pol%xi_j(3)) + pol%sig_c * pol%xi_i(2) * pol%xi_j(2)
    ! sigma >= 0; this only removes round-off near the threshold
    pol%sig_pol = MAX(pol%sig_pol, 0.0_num)

  END SUBROUTINE lbw_pol_setup



  SUBROUTINE lbw_pol_kernel(ze, mu, xa, xb, b0, b)

    ! Angular distribution of the LBW pair in the CM frame, summed over
    ! the pair spins, for the electron at cos(theta) = mu from photon i.
    ! xa, xb: the Stokes vectors of photons i and j on the
    ! scattering-plane axes (u, q x u, q) and (u, -q x u, -q), where u
    ! is the direction of the electron normal to q. Then
    !   d sigma / d Omega = r_e^2 ze / (8 C^2) b,
    !   b = K00 + K10 (xa1 + xb1) + K11 xa1 xb1 + K33 xa3 xb3
    !     + K22 xa2 xb2,
    ! with s = 1 - mu^2, D = 1 - ze^2 mu^2, g = 2 ze^2 s / D,
    !   K00 = 2 (1 + ze^2 mu^2) / D + 1 - (1 - g)^2,
    !   K10 = g (2 - g),  K11 = -1 - (1 - g)^2,  K33 = 2 (1 - g),
    !   K22 = 2 (1 - ze^4 - ze^2 s (2 - ze^2 s)) / D^2.
    ! b0 = K00 is the unpolarised distribution of random_polar_lbw. The
    ! linear terms are the polarised Klein-Nishina result crossed to
    ! pair production, K22 follows from the helicity cross-sections;
    ! all agree with direct Dirac traces (tests/spin/lbw_polarisation).
    ! The other Stokes products vanish (parity; real tree amplitudes).
    ! Over all angles K10 drops out, K11 - K33 integrates to 2 sigma_L
    ! and K22 to sigma_C (lbw_pol_setup).
    REAL(num), INTENT(IN) :: ze, mu, xa(3), xb(3)
    REAL(num), INTENT(OUT) :: b0, b
    REAL(num) :: z2, s, d, g, h

    z2 = ze**2
    s = 1.0_num - mu**2
    d = 1.0_num - z2 * mu**2
    g = 2.0_num * z2 * s / d
    h = (1.0_num - g)**2
    b0 = 2.0_num * (1.0_num + z2 * mu**2) / d + 1.0_num - h
    b = b0 + g * (2.0_num - g) * (xa(1) + xb(1)) &
        - (1.0_num + h) * xa(1) * xb(1) &
        + 2.0_num * (1.0_num - g) * xa(3) * xb(3) &
        + 2.0_num * (1.0_num - z2**2 - z2 * s * (2.0_num - z2 * s)) &
        / d**2 * xa(2) * xb(2)

  END SUBROUTINE lbw_pol_kernel



  FUNCTION lbw_pol_keep_angles(pol, ze, mu, phi) RESULT(keep)

    ! Rejection step of the polarised pair angles (LBW_polarisation):
    ! mu and the azimuth phi (from t1 towards t2) were drawn from the
    ! unpolarised distribution; keep them with probability
    ! b / (env b0) (lbw_pol_kernel), env = (1 + |xi_i|) (1 + |xi_j|).
    ! env bounds b / b0 everywhere: each photon's density matrix is at
    ! most (1 + |xi|) times the unpolarised one, and the reaction
    ! operator is positive. The mean number of draws is
    ! env sigma_0 / sigma. Unpolarised pairs are kept without a random
    ! number.
    ! KNOWN (review, 4 November 2026; section 13 of the report): the
    ! caller's loop has no bound; near threshold a pure pair needs
    ! about 1.2 / ze^2 draws (parallel linear) or 1.5 / ze^2 (opposite
    ! helicities), though at most 3 per candidate test on average.
    TYPE(lbw_pol_pair), INTENT(IN) :: pol
    REAL(num), INTENT(IN) :: ze, mu, phi
    LOGICAL :: keep
    REAL(num) :: u(3), xa(3), xb(3), env, b0, b

    keep = .TRUE.
    IF (.NOT. pol%polarised) RETURN

    ! Stokes vectors from the collision axes to the scattering plane:
    ! u is at phi from t1, about q for photon i and about -q for j
    u = COS(phi) * pol%t1 + SIN(phi) * pol%t2
    xa = pol_stored_to_obs(pol%t1, pol%t2, u, pol%xi_i)
    xb = pol_stored_to_obs(pol%t1, -pol%t2, u, pol%xi_j)
    CALL lbw_pol_kernel(ze, mu, xa, xb, b0, b)
    env = (1.0_num + SQRT(DOT_PRODUCT(pol%xi_i, pol%xi_i))) &
        * (1.0_num + SQRT(DOT_PRODUCT(pol%xi_j, pol%xi_j)))
    keep = random() * env * b0 <= b

  END FUNCTION lbw_pol_keep_angles



  SUBROUTINE lbw_pol_survivor(pol, photon_i, f, spin)

    ! LBW_survivor_update: Stokes vector of the photons left in
    ! macro-photon i (photon_i = .TRUE.) or j after a collision removed
    ! the fraction f < 1 of them. The cross-section of lbw_pol_setup is
    ! sigma_0 (1 + a . xi) for photon i, with
    ! sigma_0 a = (sigma_L xi_j1, sigma_C xi_j2, -sigma_L xi_j3) on its
    ! collision axes, and symmetrically for j. Under this absorption
    ! N xi changes at the rate of N (1 + a . xi), times the Stokes
    ! vector xi_abs = (xi + a) / (1 + a . xi), so the removed photons
    ! carry xi_abs and the survivors keep
    !   xi' = (xi - f xi_abs) / (1 - f)
    !       = xi - f / (1 - f) (a - (a . xi) xi) / (1 + a . xi).
    ! Collisions happen at a rate proportional to sigma, each removing
    ! the fraction f, so the mean eigen-populations along a follow the
    ! continuous dichroic absorption by the partner photons for any
    ! f < 1. |xi'| is clamped to 1. The result is put back on the
    ! stored (lab) axes.
    ! To be fixed:
    ! - nothing is updated at f = 1, nor are photons that are tested and
    !   not absorbed: with amplify_LBW_factor = 1 the lighter (or
    !   equal-weight) macro-photon of a pair never changes;
    ! - a pure state with a component normal to a gets |xi'| > 1 for
    !   every f > 0, so the clamp biases the mean by O(f).
    TYPE(lbw_pol_pair), INTENT(IN) :: pol
    LOGICAL, INTENT(IN) :: photon_i
    REAL(num), INTENT(IN) :: f
    REAL(num), INTENT(INOUT) :: spin(3)
    REAL(num) :: xi(3), xo(3), sa(3), xi_abs(3), norm

    ! A photon that is used up (f = 1) is deleted
    IF (.NOT. pol%polarised .OR. pol%sig_pol <= 0.0_num &
        .OR. f >= 1.0_num) RETURN
    IF (photon_i) THEN
      xi = pol%xi_i
      xo = pol%xi_j
    ELSE
      xi = pol%xi_j
      xo = pol%xi_i
    END IF
    sa = (/ pol%sig_l * xo(1), pol%sig_c * xo(2), -pol%sig_l * xo(3) /)
    ! An unpolarised partner absorbs every polarisation equally
    IF (MAXVAL(ABS(sa)) <= 0.0_num) RETURN
    xi_abs = (pol%sig_0 * xi + sa) / pol%sig_pol
    xi = (xi - f * xi_abs) / (1.0_num - f)
    norm = SQRT(DOT_PRODUCT(xi, xi))
    IF (norm > 1.0_num) xi = xi / norm
    IF (photon_i) THEN
      spin = pol_obs_to_stored(pol%ri1, pol%ri2, pol%t1, xi)
    ELSE
      spin = pol_obs_to_stored(pol%rj1, pol%rj2, pol%t1, xi)
    END IF

  END SUBROUTINE lbw_pol_survivor



  SUBROUTINE lbw_pol_spin_kernel(ze, mu, xa, xb, b, gm, gp, cc)

    ! Angular distribution of the LBW pair in the CM frame with the
    ! pair spins kept (report EPOCH2D_LBW_POLARISATION.md, section 16).
    ! ze, mu, xa and xb are as in lbw_pol_kernel. For an electron and a
    ! positron found with rest-frame spins along zm and zp,
    !   d sigma / d Omega = r_e^2 ze / (8 C^2) b(zm, zp),
    !   b(zm, zp) = (b + gm . zm + gp . zp + zm . cc . zp) / 4,
    ! where b is the spin-summed kernel of lbw_pol_kernel. gm / b and
    ! gp / b are the polarisation vectors of the electron and of the
    ! positron, and cc / b is their spin correlation (first index: the
    ! electron). Both spins are on the common axes (T, N, L): L along
    ! the electron in the CM frame, N = q x u normal to the scattering
    ! plane and T = N x L. Each is the spin in the lepton's rest frame
    ! reached from the CM frame by a pure boost.
    ! The amplitude for the photon basis states m (photon i: 1 along u,
    ! 2 along q x u) and n (photon j: 1 along u, 2 along -q x u) is
    ! 2 chi_-^+ N_mn eta_+ with the 2 x 2 matrices
    ! N_mn = n0(m,n) + nv(:,m,n) . sigma,
    !   N_11 = t sigma_T + (l_+ + l_-) sigma_L,
    !   N_22 = t sigma_T + (l_+ - l_-) sigma_L,
    !   N_12 = -i k_0 - k_N sigma_N,   N_21 = -i k_0 + k_N sigma_N,
    ! with nu = sqrt(1 - mu^2), D = 1 - ze^2 mu^2 and
    !   t = ze mu nu / D,  l_+ = ze nu^2 / (C D),  l_- = ze / (C D),
    !   k_0 = 1 / (C D),   k_N = ze nu / D.
    ! With the photons' density matrices ra, rb (EPOCH's Stokes
    ! convention) the kernel is 4 sum ra(m,m') rb(n,n') times the Pauli
    ! trace Tr[rho_- N_mn rho~_+ N_m'n'^+], written out below for
    ! N = a0 + a . sigma and N' = b0 + b . sigma (section 16.3).
    ! Only circular components (xa2, xb2) give gm and gp; cc does not
    ! vanish even for unpolarised photons.
    REAL(num), INTENT(IN) :: ze, mu, xa(3), xb(3)
    REAL(num), INTENT(OUT) :: b, gm(3), gp(3), cc(3,3)
    COMPLEX(num), PARAMETER :: iu = (0.0_num, 1.0_num)
    COMPLEX(num) :: n0(2,2), nv(3,2,2), ra(2,2), rb(2,2)
    COMPLEX(num) :: w, a0, b0, av(3), bv(3), ab(3), s0, sv, e(3)
    COMPLEX(num) :: zb, zgm(3), zgp(3), zcc(3,3)
    REAL(num) :: ic, nu, d, t, lp, lm, k0, kn
    INTEGER :: m, n, mp, np, i, j

    ic = SQRT(1.0_num - ze**2)
    nu = SQRT(MAX(1.0_num - mu**2, 0.0_num))
    d = 1.0_num - ze**2 * mu**2
    t = ze * mu * nu / d
    lp = ze * nu**2 * ic / d
    lm = ze * ic / d
    k0 = ic / d
    kn = ze * nu / d

    n0 = 0.0_num
    nv = 0.0_num
    nv(1,1,1) = t
    nv(3,1,1) = lp + lm
    nv(1,2,2) = t
    nv(3,2,2) = lp - lm
    n0(1,2) = -iu * k0
    nv(2,1,2) = -kn
    n0(2,1) = -iu * k0
    nv(2,2,1) = kn

    ra(1,1) = 0.5_num * (1.0_num + xa(1))
    ra(2,2) = 0.5_num * (1.0_num - xa(1))
    ra(1,2) = -0.5_num * (xa(3) + iu * xa(2))
    ra(2,1) = CONJG(ra(1,2))
    rb(1,1) = 0.5_num * (1.0_num + xb(1))
    rb(2,2) = 0.5_num * (1.0_num - xb(1))
    rb(1,2) = -0.5_num * (xb(3) + iu * xb(2))
    rb(2,1) = CONJG(rb(1,2))

    zb = 0.0_num
    zgm = 0.0_num
    zgp = 0.0_num
    zcc = 0.0_num
    DO np = 1, 2
      DO mp = 1, 2
        DO n = 1, 2
          DO m = 1, 2
            w = 4.0_num * ra(m,mp) * rb(n,np)
            a0 = n0(m,n)
            av = nv(:,m,n)
            b0 = CONJG(n0(mp,np))
            bv = CONJG(nv(:,mp,np))
            ! a . b* and a x b*
            sv = SUM(av * bv)
            ab = (/ av(2) * bv(3) - av(3) * bv(2), &
                av(3) * bv(1) - av(1) * bv(3), &
                av(1) * bv(2) - av(2) * bv(1) /)
            s0 = a0 * b0
            zb = zb + w * (s0 + sv)
            zgm = zgm + w * (a0 * bv + b0 * av + iu * ab)
            zgp = zgp - w * (a0 * bv + b0 * av - iu * ab)
            ! -[(a0 b0* - a . b*) delta_ij + a_i b*_j + b*_i a_j
            !   + i eps_ijl (a0 b*_l - b0* a_l)]
            e = iu * (a0 * bv - b0 * av)
            DO j = 1, 3
              DO i = 1, 3
                zcc(i,j) = zcc(i,j) - w * (av(i) * bv(j) + bv(i) * av(j))
              END DO
              zcc(j,j) = zcc(j,j) - w * (s0 - sv)
            END DO
            zcc(1,2) = zcc(1,2) - w * e(3)
            zcc(2,1) = zcc(2,1) + w * e(3)
            zcc(2,3) = zcc(2,3) - w * e(1)
            zcc(3,2) = zcc(3,2) + w * e(1)
            zcc(3,1) = zcc(3,1) - w * e(2)
            zcc(1,3) = zcc(1,3) + w * e(2)
          END DO
        END DO
      END DO
    END DO

    ! The sums are real: (m'n'; mn) is the conjugate of (mn; m'n')
    b = REAL(zb, num)
    gm = REAL(zgm, num)
    gp = REAL(zgp, num)
    cc = REAL(zcc, num)

  END SUBROUTINE lbw_pol_spin_kernel



  FUNCTION lbw_spin_cm_to_lab(zeta, p_cm, beta_v, n_v, gamma_v, p_lab) &
      RESULT(spin)

    ! Spin of a new LBW lepton as EPOCH stores it. zeta is the spin in
    ! the rest frame reached from the CM frame by a pure boost (the
    ! frame of lbw_pol_spin_kernel), given as a vector in the CM frame;
    ! p_cm and p_lab are the lepton's momenta there and in the lab,
    ! which sees the CM frame move with velocity beta_v c (unit vector
    ! n_v, Lorentz factor gamma_v). EPOCH stores the spin in the rest
    ! frame reached from the lab by a pure boost (spin_rest_to_lab).
    ! The two rest frames differ by a Wigner rotation. It is applied
    ! here through the spin 4-vector:
    !   1. CM frame: a'^0 = u' . zeta,
    !      a' = zeta + (u' . zeta) u' / (gamma' + 1),  u' = p_cm / m c;
    !   2. boost to the lab, as for the momenta in the collision
    !      routines:
    !      a^0 = gamma_v (a'^0 + beta_v . a'),
    !      a = a' + (gamma_v - 1) (a' . n_v) n_v + gamma_v a'^0 beta_v;
    !   3. the lab's rest-frame vector: spin = a - a^0 u / (gamma + 1),
    !      u = p_lab / m c.
    ! A rotation keeps the length; step 3 subtracts terms of order
    ! gamma, so the length of zeta is restored against round-off.
    REAL(num), INTENT(IN) :: zeta(3), p_cm(3), beta_v(3), n_v(3), gamma_v
    REAL(num), INTENT(IN) :: p_lab(3)
    REAL(num) :: spin(3), u(3), a(3), a0, a0_lab, len2, new2

    spin = 0.0_num
    len2 = DOT_PRODUCT(zeta, zeta)
    IF (len2 <= 0.0_num) RETURN

    u = p_cm / mc0
    a0 = DOT_PRODUCT(u, zeta)
    a = zeta + a0 / (SQRT(1.0_num + DOT_PRODUCT(u, u)) + 1.0_num) * u

    a0_lab = gamma_v * (a0 + DOT_PRODUCT(beta_v, a))
    a = a + (gamma_v - 1.0_num) * DOT_PRODUCT(a, n_v) * n_v &
        + gamma_v * a0 * beta_v

    u = p_lab / mc0
    spin = a - a0_lab / (SQRT(1.0_num + DOT_PRODUCT(u, u)) + 1.0_num) * u
    new2 = DOT_PRODUCT(spin, spin)
    IF (new2 > 0.0_num) spin = spin * SQRT(len2 / new2)

  END FUNCTION lbw_spin_cm_to_lab



  SUBROUTINE lbw_pol_pair_spins(pol, p_i, p_j, pabs_i, pabs_j, ze, cost, &
      phi, p_elec_lab, p_posi_lab, spin_electron, spin_positron)

    ! LBW_pair_spin: spins of the pair created by a collision
    ! (LBW_polarisation), at the angles already drawn. pol is the
    ! candidate's lbw_pol_pair, p_i and p_j the photons' lab momenta
    ! (sizes pabs_i, pabs_j), ze the lepton velocity in the CM frame,
    ! cost and phi the pair angles of the collision routines, and
    ! p_elec_lab, p_posi_lab the leptons' lab momenta. The outputs are
    ! the spins as EPOCH stores them (lbw_spin_cm_to_lab).
    !  1. The boost and the electron's CM direction L, as in the
    !     collision routines (about q with LBW_differential_cross,
    !     otherwise about the lab z axis). Scattering-plane axes:
    !     mu = L . q, u along the part of L normal to q, N = q x u,
    !     T = N x L; the Stokes vectors xa, xb on them as in
    !     lbw_pol_keep_angles.
    !  2. b, gm, gp, cc from lbw_pol_spin_kernel.
    !  3. spin_final_states (qed block), as in nbw_final_states:
    !     mean:     se = gm / b, sp = gp / b; no random number. Only
    !               circular photon polarisation gives a mean spin, so
    !               linearly polarised and unpolarised photons give
    !               unpolarised leptons.
    !     collapse: se = +-gm_hat with P(+) = (1 + |gm| / b) / 2 (the N
    !               axis if gm = 0). Given se, the positron has the
    !               weight b + se . gm and the vector k = gp + se . cc:
    !               sp = +-k_hat with P(+) = (1 + |k| / (b + se . gm))
    !               / 2. Two random numbers.
    ! Averaged over the random numbers, collapse gives the mean spins.
    ! It also keeps the pair correlation along the electron's axis,
    ! <se_c sp_d> = gm_hat_c (gm_hat . cc)_d / b; the rest of cc / b
    ! cannot be carried by two separate spin states. Mean mode keeps no
    ! correlation.
    ! Unpolarised photons (pol%polarised = F): gm = gp = 0, so the mean
    ! spins are zero, and collapse still draws correlated spins.
    ! The arguments are the stored particle data and the drawn angles,
    ! not the collision routines' intermediate vectors, and the boost
    ! is worked out again here: handing those vectors to this routine
    ! lets the compiler regroup the routines' arithmetic (-ffast-math),
    ! and the pair momenta then differ in the last bit from a run
    ! without pair spins.
    TYPE(lbw_pol_pair), INTENT(IN) :: pol
    REAL(num), INTENT(IN) :: p_i(3), p_j(3), pabs_i, pabs_j, ze, cost, phi
    REAL(num), INTENT(IN) :: p_elec_lab(3), p_posi_lab(3)
    REAL(num), INTENT(OUT) :: spin_electron(3), spin_positron(3)
    REAL(num) :: beta_v(3), n_v(3), gamma_v, bv2, p_lep_com(3), sint
    REAL(num) :: q(3), t1(3), t2(3), p_com(3), l_hat(3), u(3)
    REAL(num) :: t_hat(3), n_hat(3), xa(3), xb(3)
    REAL(num) :: b, gm(3), gp(3), cc(3,3), se(3), sp(3), k(3)
    REAL(num) :: mu, un, gn, kn, wk, p_plus

    spin_electron = 0.0_num
    spin_positron = 0.0_num
    IF (.NOT. pol%polarised &
        .AND. spin_final_states /= c_spin_final_collapse) RETURN

    ! Boost to the CM frame (none for opposite momenta of equal size)
    beta_v = (p_i + p_j) / (pabs_i + pabs_j)
    bv2 = DOT_PRODUCT(beta_v, beta_v)
    n_v = 0.0_num
    IF (bv2 > 0.0_num) n_v = beta_v / SQRT(bv2)
    gamma_v = 1.0_num / SQRT(1.0_num - bv2)

    ! Collision axes; lbw_pol_setup does not set them for unpolarised
    ! photons
    IF (pol%polarised) THEN
      q = pol%q
      t1 = pol%t1
      t2 = pol%t2
    ELSE
      p_com = p_i + (gamma_v - 1.0_num) * DOT_PRODUCT(p_i, n_v) * n_v &
          - gamma_v * pabs_i * beta_v
      CALL get_orthonormal(p_com, q, t1, t2)
    END IF

    ! The electron's direction and momentum in the CM frame
    sint = SQRT(MAX(1.0_num - cost**2, 0.0_num))
    IF (use_LBW_diff) THEN
      l_hat = q * cost + t1 * sint * COS(phi) + t2 * sint * SIN(phi)
    ELSE
      l_hat = (/ sint * COS(phi), sint * SIN(phi), cost /)
    END IF
    p_lep_com = mc0 * ze / SQRT(1.0_num - ze**2) * l_hat

    ! Spin axes (T, N, L). At threshold, or for a lepton along q, the
    ! kernel does not depend on the choice of u.
    mu = MIN(MAX(DOT_PRODUCT(l_hat, q), -1.0_num), 1.0_num)
    u = l_hat - DOT_PRODUCT(l_hat, q) * q
    un = SQRT(DOT_PRODUCT(u, u))
    IF (un > 1.0e-12_num) THEN
      u = u / un
    ELSE
      u = t1
      l_hat = SIGN(1.0_num, mu) * q
    END IF
    n_hat = (/ q(2) * u(3) - q(3) * u(2), q(3) * u(1) - q(1) * u(3), &
        q(1) * u(2) - q(2) * u(1) /)
    t_hat = (/ n_hat(2) * l_hat(3) - n_hat(3) * l_hat(2), &
        n_hat(3) * l_hat(1) - n_hat(1) * l_hat(3), &
        n_hat(1) * l_hat(2) - n_hat(2) * l_hat(1) /)

    xa = 0.0_num
    xb = 0.0_num
    IF (pol%polarised) THEN
      xa = pol_stored_to_obs(t1, t2, u, pol%xi_i)
      xb = pol_stored_to_obs(t1, -t2, u, pol%xi_j)
    END IF
    CALL lbw_pol_spin_kernel(ze, mu, xa, xb, b, gm, gp, cc)
    ! b > 0 for the angles that were drawn; this only guards round-off
    IF (b <= 0.0_num) RETURN

    gn = SQRT(DOT_PRODUCT(gm, gm))
    IF (spin_final_states == c_spin_final_collapse) THEN
      ! Electron along the axis of its polarisation vector
      IF (gn > TINY(1.0_num)) THEN
        se = gm / gn
      ELSE
        ! Unpolarised: any axis, with equal probabilities
        se = (/ 0.0_num, 1.0_num, 0.0_num /)
      END IF
      IF (random() >= 0.5_num * (1.0_num + MIN(gn / b, 1.0_num))) se = -se
      ! Positron given the electron
      k = gp + MATMUL(se, cc)
      wk = b + DOT_PRODUCT(se, gm)
      kn = SQRT(DOT_PRODUCT(k, k))
      IF (kn > TINY(1.0_num)) THEN
        sp = k / kn
      ELSE
        sp = (/ 0.0_num, 1.0_num, 0.0_num /)
      END IF
      ! wk >= |k| > 0 up to round-off; P(+) = 1 at the bound
      p_plus = 1.0_num
      IF (kn < wk) p_plus = 0.5_num * (1.0_num + kn / wk)
      IF (random() >= p_plus) sp = -sp
    ELSE
      se = gm / MAX(b, gn)
      sp = gp / MAX(b, SQRT(DOT_PRODUCT(gp, gp)))
    END IF

    ! CM-frame vectors, then the stored (lab) spins
    spin_electron = lbw_spin_cm_to_lab( &
        se(1) * t_hat + se(2) * n_hat + se(3) * l_hat, p_lep_com, &
        beta_v, n_v, gamma_v, p_elec_lab)
    spin_positron = lbw_spin_cm_to_lab( &
        sp(1) * t_hat + sp(2) * n_hat + sp(3) * l_hat, -p_lep_com, &
        beta_v, n_v, gamma_v, p_posi_lab)

  END SUBROUTINE lbw_pol_pair_spins
#endif
#endif
END MODULE photons
