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

MODULE deck_bremsstrahlung_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: bremsstrahlung_deck_initialise, bremsstrahlung_deck_finalise
  PUBLIC :: bremsstrahlung_block_start, bremsstrahlung_block_end
  PUBLIC :: bremsstrahlung_block_handle_element, bremsstrahlung_block_check

CONTAINS

  SUBROUTINE bremsstrahlung_deck_initialise

#ifdef BREMSSTRAHLUNG
    IF (deck_state /= c_ds_first) RETURN
    bremsstrahlung_table_location = 'src/physics_packages/TABLES/br'
    use_bremsstrahlung_recoil = .TRUE.
    use_bremsstrahlung = .FALSE.
    bh_spin_transfer = .FALSE.
    brem_spin_polarisation = .FALSE.
    brem_spin_energy_min = 1.0e6_num * q0
#ifdef SPIN
    spin_bh = c_spin_process_resolved
    bh_spin_screening = .TRUE.
    bh_spin_correlation = .TRUE.
    bh_pol_azimuth = .TRUE.
#endif
    bremsstrahlung_photon_species = -1
#ifndef PHOTONS
    photon_species = -1
#endif
    photon_energy_min_bremsstrahlung = EPSILON(1.0_num)
    bremsstrahlung_start_time = 0.0_num
    photon_weight = 1.0_num
    produce_bremsstrahlung_photons = .FALSE.
    bremsstrahlung_photon_dynamics = .FALSE.
    use_plasma_screening = .FALSE.
#endif

  END SUBROUTINE bremsstrahlung_deck_initialise



  SUBROUTINE bremsstrahlung_deck_finalise

    INTEGER :: io, iu
#ifdef BREMSSTRAHLUNG
    LOGICAL :: exists

    IF (deck_state == c_ds_first) RETURN

    IF (rank == 0 .AND. use_bremsstrahlung) THEN
      INQUIRE(file=TRIM(bremsstrahlung_table_location) // '/br1', exist=exists)
      IF (.NOT.exists) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to find bremsstrahlung tables in the ', &
              'directory "' // TRIM(bremsstrahlung_table_location) // '"'
        END DO
        CALL abort_code(c_err_io_error)
      END IF
    END IF

    IF (use_bremsstrahlung) need_random_state = .TRUE.

#ifndef BREM_TRIDENT
    IF (use_brem_trident) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to set "use_brem_trident=T" in the ', &
              '"bremsstrahlung" block.'
          WRITE(io,*) 'Please recompile with the -DBREM_TRIDENT ', &
              'preprocessor flag.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif

#ifndef BREM_MUON
    IF (use_brem_muon) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to set "use_brem_muon=T" in the ', &
              '"bremsstrahlung" block.'
          WRITE(io,*) 'Please recompile with the -DBREM_MUON ', &
              'preprocessor flag.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif
#else
    IF (use_bremsstrahlung) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to set "use_bremsstrahlung=T" in the ', &
              '"bremsstrahlung" block.'
          WRITE(io,*) 'Please recompile with the -DBREMSSTRAHLUNG ', &
              'preprocessor flag.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif

  END SUBROUTINE bremsstrahlung_deck_finalise



  SUBROUTINE bremsstrahlung_block_start

#ifdef BREMSSTRAHLUNG
    ! Enable by default if any bremsstrahlung block is defined
    use_bremsstrahlung = .TRUE.
#endif

  END SUBROUTINE bremsstrahlung_block_start



  SUBROUTINE bremsstrahlung_block_end

  END SUBROUTINE bremsstrahlung_block_end



  FUNCTION bremsstrahlung_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'use_bremsstrahlung') &
        .OR. str_cmp(element, 'bremsstrahlung') &
        .OR. str_cmp(element, 'enable')) THEN
      use_bremsstrahlung = as_logical_print(value, element, errcode)
      RETURN
    END IF

#ifdef BREMSSTRAHLUNG
    IF (str_cmp(element, 'bremsstrahlung_start_time') &
        .OR. str_cmp(element, 'start_time')) THEN
      bremsstrahlung_start_time = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'produce_bremsstrahlung_photons') &
        .OR. str_cmp(element, 'produce_photons')) THEN
      produce_bremsstrahlung_photons = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_bremsstrahlung_recoil') &
        .OR. str_cmp(element, 'use_radiation_reaction')) THEN
      use_bremsstrahlung_recoil = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'photon_energy_min') &
        .OR. str_cmp(element, 'min_photon_energy') &
        .OR. str_cmp(element, 'photon_energy_min_bremsstrahlung')) THEN
      photon_energy_min_bremsstrahlung = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'photon_weight') &
        .OR. str_cmp(element, 'photon_weight_multiplier')) THEN
      photon_weight = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bremsstrahlung_table_location') &
        .OR. str_cmp(element, 'table_location')) THEN
      bremsstrahlung_table_location = TRIM(ADJUSTL(value))
      RETURN
    END IF

    IF (str_cmp(element, 'photon_dynamics') &
        .OR. str_cmp(element, 'bremsstrahlung_photon_dynamics')) THEN
      bremsstrahlung_photon_dynamics = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_plasma_screening')) THEN
      use_plasma_screening = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_brem_scatter') &
        .OR. str_cmp(element, 'use_bremsstrahlung_scatter')) THEN
      use_brem_scatter = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_bethe_heitler')) THEN
      use_bethe_heitler = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bh_spin_transfer')) THEN
      bh_spin_transfer = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'brem_spin_polarisation')) THEN
      brem_spin_polarisation = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'brem_spin_energy_min')) THEN
      brem_spin_energy_min = as_real_print(value, element, errcode)
      RETURN
    END IF

#ifdef SPIN
    ! BH pair polarisation switches (for tests; need bh_spin_transfer = T)
    IF (str_cmp(element, 'spin_bh')) THEN
      IF (str_cmp(value, 'resolved')) THEN
        spin_bh = c_spin_process_resolved
      ELSE IF (str_cmp(value, 'unpolarised') &
          .OR. str_cmp(value, 'unpolarized')) THEN
        spin_bh = c_spin_process_unpolarised
      ELSE IF (str_cmp(value, 'off')) THEN
        spin_bh = c_spin_process_off
      ELSE
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'bh_spin_screening')) THEN
      bh_spin_screening = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bh_spin_correlation')) THEN
      bh_spin_correlation = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bh_pol_azimuth')) THEN
      bh_pol_azimuth = as_logical_print(value, element, errcode)
      RETURN
    END IF
#endif

    IF (str_cmp(element, 'use_brem_trident')) THEN
      use_brem_trident = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_brem_muon')) THEN
      use_brem_muon = as_logical_print(value, element, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'use_positron_brem') .OR. &
        str_cmp(element, 'positron_brem')) THEN
      positron_brem = as_logical_print(value, element, errcode)
      RETURN
    END IF

    errcode = c_err_unknown_element
#endif

  END FUNCTION bremsstrahlung_block_handle_element



  FUNCTION bremsstrahlung_block_check() RESULT(errcode)

    INTEGER :: errcode
#ifdef BREMSSTRAHLUNG
    INTEGER :: io, iu
#endif

    errcode = c_err_none

#ifdef BREMSSTRAHLUNG
    IF (brem_spin_polarisation) THEN
#if !defined(SPIN) || !defined(PHOTONS)
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'brem_spin_polarisation needs -DSPIN -DPHOTONS.'
        END DO
      END IF
      errcode = c_err_pp_options_missing + c_err_terminate
#else
      IF (.NOT. use_bremsstrahlung .OR. &
          .NOT. produce_bremsstrahlung_photons .OR. &
          .NOT. use_bremsstrahlung_recoil .OR. &
          .NOT. use_brem_scatter .OR. use_plasma_screening .OR. &
          positron_brem .OR. .NOT. photon_polarisation .OR. &
          ABS(photon_weight - 1.0_num) > 1.0e-12_num .OR. &
          brem_spin_energy_min < 1.0e3_num * q0) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'brem_spin_polarisation requires photons,', &
                ' recoil, scatter, photon_weight=1, and', &
                ' brem_spin_energy_min >= 1 keV; plasma screening', &
                ' and positron bremsstrahlung must be off.'
          END DO
        END IF
        errcode = c_err_bad_value + c_err_terminate
      ELSE IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) 'BH bremsstrahlung spin mode: Olsen-Maximon', &
              ' screened electron-nucleus kernel for gamma_i=50..1e8', &
              ' and gamma_f>=10; below/above this range the old', &
              ' unpolarised bremsstrahlung path is used.'
          WRITE(io,*) 'The physical photon event cut is', &
              brem_spin_energy_min / (1.0e6_num * q0), ' MeV.'
        END DO
      END IF
#endif
    END IF

    IF (bh_spin_transfer) THEN
#ifndef SPIN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'bh_spin_transfer requires -DSPIN.'
        END DO
      END IF
      errcode = c_err_pp_options_missing + c_err_terminate
#else
      IF (.NOT. use_bremsstrahlung .OR. .NOT. use_bethe_heitler) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'bh_spin_transfer requires use_bethe_heitler=T', &
                ' and use_bremsstrahlung=T.'
          END DO
        END IF
        errcode = c_err_bad_value + c_err_terminate
      ELSE IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'BH spin transfer uses a high-energy, ', &
              'angle-integrated approximation for photons at least 20 MeV', &
              ' with both lepton Lorentz factors at least 10.'
          WRITE(io,*) 'Outside this range, newborn pair spins are zero.'
          WRITE(io,*) 'Use high-Z, stationary targets: the existing BH', &
              ' total rate includes an atomic-electron contribution.'
          WRITE(io,*) 'Helicities and their correlation follow', &
              ' qed spin_final_states; the pair azimuth follows', &
              ' the linear photon polarisation (one-lepton law).'
          WRITE(io,*) 'Transverse pair-spin correlations and', &
              ' spin-angle correlations are not included.'
        END DO
      END IF
#endif
    END IF

#ifdef SPIN
    ! The BH pair polarisation switches act only with bh_spin_transfer
    IF (.NOT. bh_spin_transfer .AND. (spin_bh /= c_spin_process_resolved &
        .OR. .NOT. bh_spin_screening .OR. .NOT. bh_spin_correlation &
        .OR. .NOT. bh_pol_azimuth)) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'spin_bh, bh_spin_screening, bh_spin_correlation', &
              ' and bh_pol_azimuth need bh_spin_transfer=T.'
        END DO
      END IF
      errcode = c_err_bad_value + c_err_terminate
    END IF
#endif

    IF (photon_weight <= 0.0_num) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'You cannot set the photon_weight to less than or ', &
              'equal to zero. To prevent bremsstrahlung photons from being ', &
              'emitted, set use_bremsstrahlung = F.'
          WRITE(io,*) 'Code will terminate.'
        END DO
      END IF
      errcode = c_err_bad_value + c_err_terminate
    END IF

    IF (photon_weight > 1.0_num) THEN
      photon_weight = 1.0_num
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'You cannot set photon_weight > 1.0. This variable ', &
              'has been truncated to 1.0.'
        END DO
      END IF
    END IF
#endif

  END FUNCTION bremsstrahlung_block_check

END MODULE deck_bremsstrahlung_block
