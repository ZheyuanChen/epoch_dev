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


#endif

END MODULE spin