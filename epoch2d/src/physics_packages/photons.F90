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

    sig2cdt_dV_lbw = 2.0_num * sigma_lBW_max * c * dt / dx / dy &
        * LBW_amp_factor
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
    ! Spin-resolved photon-energy CDFs, on the grid just built
    IF (use_spin_qed) CALL setup_nlcs_spin_cdf_tables
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
    REAL(num) :: dtau
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
            CALL generate_photon(current, photon_species, eta)
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
    REAL(num) :: u_mag, xi_stored(3), spin_avg(3)
    LOGICAL :: spin_event
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
      IF (random() < nlcs_down_probability(zeta, a_spin)) THEN
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
        IF (use_classical_emission) THEN
          g_eta = 1.0_num
        ELSE
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
      CALL nlcs_final_states(dir_x * u_mag, dir_y * u_mag, &
          dir_z * u_mag, generating_gamma, e_field, b_field, charge, &
          photon_energy / (generating_gamma * m0 * c**2), eta, make_photon, &
          generating_electron%spin, xi_stored)
    END IF
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
    REAL(num) :: com_beta, sigma_lbw, P_coll
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

      ! kinematic factor (1-cos(theta))
      kappa = pipj_cos / pipj

      ! dimensionless photon energy squared in c.o.m.
      en_com2 = pipj_cos * 0.5_num * inv_mc0_sq

      IF (en_com2 < 1.0_num) THEN
        current_i => current_j%next
        IF (ASSOCIATED(current_i)) current_j => current_i%next
        N_scattered = N_scattered +1
        CYCLE
      END IF

      !!! calculate collisional probability P_ij/P_max
      weight_i = current_i%weight
      weight_j = current_j%weight

      ! normalized lepton velocity in c.o.m. frame
      com_beta = SQRT(1.0_num - 1.0_num/en_com2)

      ! cross section for the linear BW (without kinematic factor)
      ! Notice this is a Lorentz invariant.
      ! Later we will re-use it in the c.o.m frame.
      sigma_lbw = lbw_cross_sec(com_beta)

      ! collisional probability after modification by P_max
      P_coll = sigma_lbw * cdt_dV * MAX(weight_i, weight_j) * kappa * i_Pmax

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

      ! Unit vector of this Lorentz transformation velocity
      n_v = beta_v /SQRT(DOT_PRODUCT(beta_v, beta_v))

      ! Gamma factor of this Lorentz transformatin
      gamma_v = 1.0_num/SQRT(1.0_num-DOT_PRODUCT(beta_v,beta_v))

      ! Magnitude of lepton momentum in c.o.m. (in SI)
      p_lep_com_mag = SQRT(en_com2 - 1.0_num)*mc0

      IF (use_LBW_diff) THEN
        ! Random azimuthal angle in c.o.m.
        rand_phi = 2.0_num * pi * random()

        ! Random polar angle in c.o.m. (cosine of this random angle)
        cost = random_polar_lbw(com_beta, sigma_lbw)
        sint = SQRT(1.0_num - cost**2)

        ! Photon momentum in c.o.m. (in SI)
        p_phot_com = moment_i + (gamma_v-1.0_num) &
            * DOT_PRODUCT(moment_i,n_v)*n_v - gamma_v*pi_abs*beta_v

        ! Orthonormal basis (e1, e2, e3) s.t. e1//photon momentum in c.o.m.
        CALL get_orthonormal(p_phot_com, e1, e2, e3)

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
    REAL(num) :: com_beta, sigma_lbw, P_coll
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

      ! kinematic factor (1-cos(theta))
      kappa = pipj_cos / pipj

      ! dimensionless photon energy squared in c.o.m.
      en_com2 = pipj_cos * 0.5_num * inv_mc0_sq

      IF (en_com2 < 1.0_num) THEN
        current_i => current_i%next
        current_j => current_j%next
        N_scattered = N_scattered +1
        CYCLE
      END IF

      !!! calculate collisional probability P_ij/P_max
      weight_i = current_i%weight
      weight_j = current_j%weight

      ! normalized lepton velocity in c.o.m. frame
      com_beta = SQRT(1.0_num - 1.0_num/en_com2)

      ! cross section for the linear BW (without kinematic factor)
      ! Notice this is a Lorentz invariant.
      ! Later we will re-use it in the c.o.m frame.
      sigma_lbw = lbw_cross_sec(com_beta)

      ! collisional probability after modification by P_max
      P_coll = sigma_lbw * cdt_dV * MAX(weight_i, weight_j) * kappa * i_Pmax

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

      ! Unit vector of this Lorentz transformation velocity
      n_v = beta_v /SQRT(DOT_PRODUCT(beta_v, beta_v))

      ! Gamma factor of this Lorentz transformatin
      gamma_v = 1.0_num/SQRT(1.0_num-DOT_PRODUCT(beta_v,beta_v))

      ! Magnitude of lepton momentum in c.o.m. (in SI)
      p_lep_com_mag = SQRT(en_com2 - 1.0_num)*mc0

      IF (use_LBW_diff) THEN
        ! Random azimuthal angle in c.o.m.
        rand_phi = 2.0_num * pi * random()

        ! Random polar angle in c.o.m. (cosine of this random angle)
        cost = random_polar_lbw(com_beta, sigma_lbw)
        sint = SQRT(1.0_num - cost**2)

        ! Photon momentum in c.o.m. (in SI)
        p_phot_com = moment_i + (gamma_v-1.0_num) &
            * DOT_PRODUCT(moment_i,n_v)*n_v - gamma_v*pi_abs*beta_v

        ! Orthonormal basis (e1, e2, e3) s.t. e1//photon momentum in c.o.m.
        CALL get_orthonormal(p_phot_com, e1, e2, e3)

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
#endif
END MODULE photons
