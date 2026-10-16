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

MODULE deck_qed_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: qed_deck_initialise, qed_deck_finalise
  PUBLIC :: qed_block_start, qed_block_end
  PUBLIC :: qed_block_handle_element, qed_block_check

CONTAINS

  SUBROUTINE qed_deck_initialise

#ifdef PHOTONS
    IF (deck_state == c_ds_first) THEN
      qed_table_location = 'src/physics_packages/TABLES'
      use_radiation_reaction = .TRUE.
      use_qed = .FALSE.
      photon_species = -1
      trident_electron_species = -1
      breit_wheeler_electron_species = -1
      trident_positron_species = -1
      breit_wheeler_positron_species = -1
      photon_energy_min = EPSILON(1.0_num)
      qed_start_time = 0.0_num
      produce_pairs = .FALSE.
      use_radiation_reaction = .TRUE.
      produce_photons = .FALSE.
      photon_dynamics = .FALSE.
      use_emission_substeps = .TRUE.
      use_binary_collisions = .FALSE.
      use_LBW = .FALSE.
      use_LBW_diff = .TRUE.
      LBW_amp_factor = 1.0_num
#ifdef SPIN
      use_spin_qed = .FALSE.
      spin_second_order_depth = .FALSE.
      spin_final_states = c_spin_final_collapse
      pair_tables = c_pair_tables_unset
      spin_no_transition = .TRUE.
      spin_nlcs = c_spin_process_resolved
      spin_nbw = c_spin_process_resolved
      use_chi_dependent_moment = .TRUE.
#endif
    END IF
#endif



  END SUBROUTINE qed_deck_initialise



  SUBROUTINE qed_deck_finalise

    INTEGER :: io, iu
#ifdef PHOTONS
    LOGICAL :: exists
    INTEGER :: j

    IF (deck_state == c_ds_first) RETURN

    IF (rank == 0 .AND. use_qed) THEN
      INQUIRE(file=TRIM(qed_table_location)//'/hsokolov.table', exist=exists)
      IF (.NOT.exists) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to find QED tables in the ', &
              'directory "' // TRIM(qed_table_location) // '"'
        END DO
        CALL abort_code(c_err_io_error)
      END IF
    END IF

    ! Check for chi-dependent anomalous moment table
#if defined(SPIN) && defined(PHOTONS)
    IF (rank == 0 .AND. use_qed .AND. use_chi_dependent_moment) THEN
      INQUIRE(file=TRIM(qed_table_location)//'/anomalous_moment.table'&
          , exist=exists)
      IF (.NOT.exists) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to find SPIN tables in the ', &
              'directory "' // TRIM(qed_table_location) // '"'
        END DO
        CALL abort_code(c_err_io_error)
      END IF
    END IF

    ! NBW tables: the exact LCFA ones by default with use_spin_qed,
    ! EPOCH's legacy ones otherwise
    IF (pair_tables == c_pair_tables_unset) THEN
      IF (use_spin_qed) THEN
        pair_tables = c_pair_tables_lcfa
      ELSE
        pair_tables = c_pair_tables_legacy
      END IF
    END IF
    IF (rank == 0 .AND. use_qed .AND. pair_tables == c_pair_tables_lcfa) &
        THEN
      CALL check_spin_table('pairprod_lcfa.table')
      CALL check_spin_table('nbw_polarisation_asymmetry.table')
      CALL check_spin_table('nbw_energy_split_par.table')
      CALL check_spin_table('nbw_energy_split_perp.table')
    END IF

    ! Spin-dependent NLCS rate and photon-energy spectra
    IF (rank == 0 .AND. use_qed .AND. use_spin_qed) THEN
      CALL check_spin_table('nlcs_spin_asymmetry.table')
      CALL check_spin_table('ksi_sokolov_up.table')
      CALL check_spin_table('ksi_sokolov_down.table')
      CALL check_spin_table('bessel_k.table')
      CALL check_spin_table('nbw_polarisation_asymmetry.table')
      CALL check_spin_table('nbw_energy_split_par.table')
      CALL check_spin_table('nbw_energy_split_perp.table')
    END IF
#endif

    IF (use_qed) need_random_state = .TRUE.

    use_binary_collisions = use_LBW

    IF (use_binary_collisions) THEN
      DO j = 1, n_species
        IF (species_list(j)%species_type == c_species_id_photon) THEN
          species_list(j)%make_secondary_list = .TRUE.
        END IF
        IF (species_list(j)%species_type == c_species_id_electron) THEN
          species_list(j)%make_secondary_list = .TRUE.
        END IF
        IF (species_list(j)%species_type == c_species_id_positron) THEN
          species_list(j)%make_secondary_list = .TRUE.
        END IF
    END DO
  END IF

  lbw_amp_factor = MAX(lbw_amp_factor, 1.0_num)
#else
    IF (use_qed) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to set "use_qed=T" in the "qed" block.'
          WRITE(io,*) 'Please recompile with the -DPHOTONS preprocessor flag.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif
#if defined(SPIN) && defined(PHOTONS)
  CONTAINS

    SUBROUTINE check_spin_table(name)

      CHARACTER(LEN=*), INTENT(IN) :: name
      LOGICAL :: found

      INQUIRE(file=TRIM(qed_table_location) // '/' // name, exist=found)
      IF (.NOT.found) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to find ', name, ' in the ', &
              'directory "' // TRIM(qed_table_location) // '"'
        END DO
        CALL abort_code(c_err_io_error)
      END IF

    END SUBROUTINE check_spin_table
#endif

  END SUBROUTINE qed_deck_finalise



  SUBROUTINE qed_block_start

  END SUBROUTINE qed_block_start



  SUBROUTINE qed_block_end

  END SUBROUTINE qed_block_end



  FUNCTION qed_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'use_qed') .OR. str_cmp(element, 'qed')) THEN
      use_qed = as_logical_print(value, element, errcode)
      RETURN
    END IF

#ifdef PHOTONS
    IF (str_cmp(element, 'qed_start_time')) THEN
      qed_start_time = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'produce_photons')) THEN
      produce_photons = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_radiation_reaction')) THEN
      use_radiation_reaction = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'photon_energy_min') &
        .OR. str_cmp(element, 'min_photon_energy')) THEN
      photon_energy_min = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'produce_pairs')) THEN
      produce_pairs = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'qed_table_location')) THEN
      qed_table_location = TRIM(ADJUSTL(value))
      RETURN
    END IF

    IF (str_cmp(element, 'photon_dynamics')) THEN
      photon_dynamics = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF(str_cmp(element, 'semiclassical_photon_emission')) THEN
       use_continuous_emission = as_logical_print(value, element, errcode)
       RETURN
    ENDIF

    IF(str_cmp(element, 'classical_photon_emission')) THEN
      use_classical_emission = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'photon_sample_fraction') &
        .OR. str_cmp(element, 'photon_downsampling')) THEN
      photon_sample_fraction = as_real_print(value, element, errcode)
      IF (photon_sample_fraction > 1.0) photon_sample_fraction = 1.0
      RETURN
    END IF

    IF(str_cmp(element, 'linear_breit_wheeler')) THEN
      use_LBW = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF(str_cmp(element, 'LBW_differential_cross')) THEN
      use_LBW_diff = as_logical_print(value, element, errcode)
      RETURN
    END IF

    ! Emission followed event by event within a step (default), or
    ! EPOCH's legacy step: at most one emission, overshoot discarded
    IF (str_cmp(element, 'use_emission_substeps')) THEN
      use_emission_substeps = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'amplify_LBW_factor')) THEN
      LBW_amp_factor = as_real_print(value, element, errcode)
      RETURN
    END IF

#ifdef SPIN
    IF (str_cmp(element, 'use_spin_qed') &
        .OR. str_cmp(element, 'spin_qed')) THEN
      use_spin_qed = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'spin_second_order_depth')) THEN
      spin_second_order_depth = as_logical_print(value, element, errcode)
      RETURN
    END IF

    ! NBW rate and energy split: legacy (pairprod.table's Erber
    ! approximation, energy_split.table) or lcfa (exact LCFA). The
    ! default, set when the block is finalised, follows use_spin_qed.
    IF (str_cmp(element, 'pair_tables')) THEN
      IF (str_cmp(value, 'legacy')) THEN
        pair_tables = c_pair_tables_legacy
      ELSE IF (str_cmp(value, 'lcfa')) THEN
        pair_tables = c_pair_tables_lcfa
      ELSE
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    ! Final spin and photon polarisation of an event: collapse onto the
    ! quantisation axes (default) or keep the mean
    IF (str_cmp(element, 'spin_final_states')) THEN
      IF (str_cmp(value, 'collapse')) THEN
        spin_final_states = c_spin_final_collapse
      ELSE IF (str_cmp(value, 'mean')) THEN
        spin_final_states = c_spin_final_mean
      ELSE
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    ! Switches for tests: the survivor (non-transition) updates, and
    ! each process resolved, unpolarised (as with use_spin_qed = F) or
    ! off (no events; the survivor updates still use its rate)
    IF (str_cmp(element, 'spin_no_transition')) THEN
      spin_no_transition = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'spin_nlcs')) THEN
      spin_nlcs = as_spin_process(value, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'spin_nbw')) THEN
      spin_nbw = as_spin_process(value, errcode)
      ! photon_only exists for NLCS only
      IF (spin_nbw == c_spin_process_photon_only) &
          errcode = IOR(errcode, c_err_bad_value)
      RETURN
    END IF

    ! a(chi) for leptons in T-BMT with use_qed = T (default), or the
    ! species' constant moment (for comparison with constant-a theory)
    IF (str_cmp(element, 'use_chi_dependent_moment') &
        .OR. str_cmp(element, 'chi_dependent_moment')) THEN
      use_chi_dependent_moment = as_logical_print(value, element, errcode)
      RETURN
    END IF
#endif

    errcode = c_err_unknown_element
#endif

  END FUNCTION qed_block_handle_element



  FUNCTION qed_block_check() RESULT(errcode)

    INTEGER :: errcode
#ifdef PHOTONS
    INTEGER :: io, iu
#endif

    errcode = c_err_none

#ifdef PHOTONS
    IF (produce_pairs .AND. .NOT. photon_dynamics) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'You cannot set photon_dynamics=F when ', &
            'produce_pairs=T. Without ', 'photon motion, pair ', &
            'creation will be incorrect.'
          WRITE(io,*) 'Code will terminate.'
        END DO
      END IF
      errcode = c_err_bad_value + c_err_terminate
    END IF

    ! Handling for different combinations of semiclassical and classical emission
    IF (use_continuous_emission .AND. use_classical_emission) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'You cannot have both classical_photon_emission=T ', &
            'and semiclassical_photon_emission=T. Assuming only ', &
            'classical emission=T.'
        END DO
      END IF
    END IF
    use_continuous_emission = (use_continuous_emission &
                               .OR. use_classical_emission)

#ifdef SPIN
    ! Continuous emission has no discrete events to resolve spin with
    IF (use_qed .AND. use_spin_qed .AND. use_continuous_emission) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'use_spin_qed=T cannot be combined with ', &
              'classical or semiclassical photon emission.'
          WRITE(io,*) 'Code will terminate.'
        END DO
      END IF
      errcode = c_err_bad_value + c_err_terminate
    END IF

    ! The test switches act on the spin-resolved QED only
    IF (use_qed .AND. .NOT. use_spin_qed .AND. (.NOT. spin_no_transition &
        .OR. spin_nlcs /= c_spin_process_resolved &
        .OR. spin_nbw /= c_spin_process_resolved)) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'spin_no_transition, spin_nlcs and spin_nbw ', &
              'need use_spin_qed=T.'
          WRITE(io,*) 'Code will terminate.'
        END DO
      END IF
      errcode = c_err_bad_value + c_err_terminate
    END IF
#endif
#endif

  END FUNCTION qed_block_check



#if defined(SPIN) && defined(PHOTONS)
  FUNCTION as_spin_process(str_in, err) RESULT(res)

    ! Value of spin_nlcs and spin_nbw: resolved, unpolarised, off or
    ! (spin_nlcs only) photon_only
    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: res

    res = c_spin_process_resolved
    IF (str_cmp(str_in, 'resolved')) THEN
      res = c_spin_process_resolved
    ELSE IF (str_cmp(str_in, 'unpolarised') &
        .OR. str_cmp(str_in, 'unpolarized')) THEN
      res = c_spin_process_unpolarised
    ELSE IF (str_cmp(str_in, 'off')) THEN
      res = c_spin_process_off
    ELSE IF (str_cmp(str_in, 'photon_only')) THEN
      res = c_spin_process_photon_only
    ELSE
      err = IOR(err, c_err_bad_value)
    END IF

  END FUNCTION as_spin_process
#endif

END MODULE deck_qed_block
