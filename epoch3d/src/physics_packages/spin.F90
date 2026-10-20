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

  ! NLCS rate asymmetry A(eta) on a uniform log10(eta) grid
  ! (nlcs_spin_asymmetry.table): W(s) = W_bar (1 + A zeta)
  INTEGER :: n_asym
  REAL(num) :: asym_log_eta_min, asym_log_eta_max, asym_idlog
  REAL(num), ALLOCATABLE :: asym_log_a(:)

  ! Photon-energy CDFs of NLCS for zeta = +1 (up) and -1 (down), on
  ! the grid of ksi_sokolov.table (log_eta, log_chi); sampled with
  ! find_value_from_table_alt like p_photon_energy
  REAL(num), ALLOCATABLE :: nlcs_cdf_up(:,:), nlcs_cdf_down(:,:)

  ! The same photon energies as quantiles (nlcs_quantile_up/down.table,
  ! use_nlcs_quantile_tables): log10 chi at the probabilities nlcs_q_p,
  ! stored as nlcs_q_up(j, i) for row i, uniform in log10(eta) from
  ! nlcs_q_log_eta_min in steps nlcs_q_dlog; sampled by
  ! nlcs_quantile_chi
  INTEGER :: nlcs_q_n_eta, nlcs_q_n_p
  REAL(num) :: nlcs_q_log_eta_min, nlcs_q_log_eta_max, nlcs_q_dlog
  REAL(num), ALLOCATABLE :: nlcs_q_p(:), nlcs_q_up(:,:), nlcs_q_down(:,:)

  ! R1 = K_1/3(x)/K_2/3(x) and R2 = IntK_1/3(x)/K_2/3(x) as log10 R on
  ! a uniform log10(x) grid (bessel_k.table)
  INTEGER :: n_bessel
  REAL(num) :: bessel_log_x_min, bessel_log_x_max, bessel_idlog
  REAL(num), ALLOCATABLE :: bessel_log_r1(:), bessel_log_r2(:)

  ! NBW rate asymmetry A(chi) on a uniform log10(chi) grid
  ! (nbw_polarisation_asymmetry.table): W(xi) = W_bar (1 - A xi1_obs)
  INTEGER :: n_nbw_asym
  REAL(num) :: nbw_asym_log_chi_min, nbw_asym_log_chi_max, nbw_asym_idlog
  REAL(num), ALLOCATABLE :: nbw_asym_log_a(:)

  ! Electron-fraction CDFs of NBW for xi1_obs = +1 (par) and -1 (perp)
  ! (nbw_energy_split_par/perp.table) on the grid nbw_log_chi (log10
  ! chi) x nbw_split (0..1); sampled with find_value_from_table
  INTEGER :: n_nbw_chi, n_nbw_split
  REAL(num), ALLOCATABLE :: nbw_log_chi(:), nbw_split(:)
  REAL(num), ALLOCATABLE :: nbw_cdf_par(:,:), nbw_cdf_perp(:,:)

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

  ! nlcs_spin_asymmetry.table (A(chi))
  ! The grid is uniform in log10(eta)
  ! Called only with use_spin_qed.
  SUBROUTINE setup_nlcs_asymmetry_table

    INTEGER :: i
    REAL(num) :: log_eta, dlog, buf(2)

    IF (rank == 0) THEN
      OPEN(unit=lu, file=TRIM(qed_table_location) &
          // '/nlcs_spin_asymmetry.table', status='OLD')
      READ(lu,*) n_asym, asym_log_eta_min, asym_log_eta_max
      ALLOCATE(asym_log_a(n_asym))
      dlog = (asym_log_eta_max - asym_log_eta_min) / (n_asym - 1)
      DO i = 1, n_asym
        READ(lu,*) log_eta, asym_log_a(i)
        IF (ABS(log_eta - asym_log_eta_min - (i - 1) * dlog) &
            > 1.0e-8_num) CALL abort_code(c_err_bad_value)
      END DO
      CLOSE(unit=lu)
    END IF

    CALL MPI_BCAST(n_asym, 1, MPI_INTEGER, 0, comm, errcode)
    buf = (/ asym_log_eta_min, asym_log_eta_max /)
    CALL MPI_BCAST(buf, 2, mpireal, 0, comm, errcode)
    asym_log_eta_min = buf(1)
    asym_log_eta_max = buf(2)
    IF (rank /= 0) ALLOCATE(asym_log_a(n_asym))
    CALL MPI_BCAST(asym_log_a, n_asym, mpireal, 0, comm, errcode)
    asym_idlog = (n_asym - 1) / (asym_log_eta_max - asym_log_eta_min)

  END SUBROUTINE setup_nlcs_asymmetry_table



  ! Read ksi_sokolov_up.table and ksi_sokolov_down.table
  ! , which share the grid of ksi_sokolov.table. 
  !Call after setup_tables_qed has built log_eta.
  ! Called only with use_spin_qed.
  SUBROUTINE setup_nlcs_spin_cdf_tables

    CALL read_cdf_table('ksi_sokolov_up.table', nlcs_cdf_up)
    CALL read_cdf_table('ksi_sokolov_down.table', nlcs_cdf_down)

  CONTAINS

    SUBROUTINE read_cdf_table(name, cdf)

      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL(num), ALLOCATABLE, INTENT(INOUT) :: cdf(:,:)
      INTEGER :: n_eta, n_chi, ieta, ichi
      REAL(num) :: lmin, lmax

      ALLOCATE(cdf(n_sample_eta, n_sample_chi))
      IF (rank == 0) THEN
        OPEN(unit=lu, file=TRIM(qed_table_location) // '/' // name, &
            status='OLD')
        READ(lu,*) n_eta, n_chi, lmin, lmax
        ! The table must share the grid of ksi_sokolov.table
        IF (n_eta /= n_sample_eta .OR. n_chi /= n_sample_chi &
            .OR. ABS(lmin - log_eta(1)) > 1.0e-8_num &
            .OR. ABS(lmax - log_eta(n_sample_eta)) > 1.0e-8_num) THEN
          PRINT*, '*** ERROR ***'
          PRINT*, TRIM(name), ' is not on the grid of ksi_sokolov.table'
          CALL abort_code(c_err_bad_value)
        END IF
        DO ieta = 1, n_eta
          READ(lu,*) (cdf(ieta,ichi), ichi = 1, n_chi)
        END DO
        CLOSE(unit=lu)
      END IF
      CALL MPI_BCAST(cdf, n_sample_eta * n_sample_chi, mpireal, 0, comm, &
          errcode)

    END SUBROUTINE read_cdf_table

  END SUBROUTINE setup_nlcs_spin_cdf_tables

  ! This table has better resolution than the old one, making the NLCS
  ! momentum sampling more accurate but increasing the computation cost
  ! Read nlcs_quantile_up.table and nlcs_quantile_down.table.
  ! which must share their eta rows and probabilities. 
  ! Called only with use_spin_qed and use_nlcs_quantile_tables.
  SUBROUTINE setup_nlcs_quantile_tables

    INTEGER :: ibuf(2)
    REAL(num) :: buf(2)

    IF (rank == 0) THEN
      CALL read_quantile_table('nlcs_quantile_up.table', nlcs_q_up, .TRUE.)
      CALL read_quantile_table('nlcs_quantile_down.table', nlcs_q_down, &
          .FALSE.)
    END IF

    ibuf = (/ nlcs_q_n_eta, nlcs_q_n_p /)
    CALL MPI_BCAST(ibuf, 2, MPI_INTEGER, 0, comm, errcode)
    nlcs_q_n_eta = ibuf(1)
    nlcs_q_n_p = ibuf(2)
    buf = (/ nlcs_q_log_eta_min, nlcs_q_log_eta_max /)
    CALL MPI_BCAST(buf, 2, mpireal, 0, comm, errcode)
    nlcs_q_log_eta_min = buf(1)
    nlcs_q_log_eta_max = buf(2)
    IF (rank /= 0) THEN
      ALLOCATE(nlcs_q_p(nlcs_q_n_p))
      ALLOCATE(nlcs_q_up(nlcs_q_n_p, nlcs_q_n_eta))
      ALLOCATE(nlcs_q_down(nlcs_q_n_p, nlcs_q_n_eta))
    END IF
    CALL MPI_BCAST(nlcs_q_p, nlcs_q_n_p, mpireal, 0, comm, errcode)
    CALL MPI_BCAST(nlcs_q_up, nlcs_q_n_p * nlcs_q_n_eta, mpireal, 0, comm, &
        errcode)
    CALL MPI_BCAST(nlcs_q_down, nlcs_q_n_p * nlcs_q_n_eta, mpireal, 0, &
        comm, errcode)
    nlcs_q_dlog = (nlcs_q_log_eta_max - nlcs_q_log_eta_min) &
        / REAL(nlcs_q_n_eta - 1, num)

  CONTAINS

    SUBROUTINE read_quantile_table(name, q, first)

      ! Header: n_eta, n_p, log10 eta_min, log10 eta_max; then the n_p
      ! probabilities, then one row of n_p values log10 chi per eta.
      ! The first table sets the grid, the second must match it.
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL(num), ALLOCATABLE, INTENT(INOUT) :: q(:,:)
      LOGICAL, INTENT(IN) :: first
      INTEGER :: n_eta, n_p, ieta, j
      REAL(num) :: lmin, lmax
      REAL(num), ALLOCATABLE :: p(:)
      LOGICAL :: bad

      OPEN(unit=lu, file=TRIM(qed_table_location) // '/' // name, &
          status='OLD')
      READ(lu,*) n_eta, n_p, lmin, lmax
      ALLOCATE(p(n_p))
      READ(lu,*) (p(j), j = 1, n_p)
      IF (first) THEN
        nlcs_q_n_eta = n_eta
        nlcs_q_n_p = n_p
        nlcs_q_log_eta_min = lmin
        nlcs_q_log_eta_max = lmax
        ALLOCATE(nlcs_q_p(n_p))
        nlcs_q_p = p
        bad = n_eta < 2 .OR. n_p < 2 .OR. lmax <= lmin
        ! The lookup needs p from 0 to 1, strictly increasing
        IF (.NOT. bad) bad = p(1) /= 0.0_num .OR. p(n_p) /= 1.0_num &
            .OR. ANY(p(2:n_p) <= p(1:n_p-1))
      ELSE
        bad = n_eta /= nlcs_q_n_eta .OR. n_p /= nlcs_q_n_p &
            .OR. ABS(lmin - nlcs_q_log_eta_min) > 1.0e-8_num &
            .OR. ABS(lmax - nlcs_q_log_eta_max) > 1.0e-8_num
        IF (.NOT. bad) bad = ANY(p /= nlcs_q_p)
      END IF
      IF (bad) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, TRIM(name), ' has a bad header or probability grid, or ', &
            'does not match nlcs_quantile_up.table'
        CALL abort_code(c_err_bad_value)
      END IF
      ALLOCATE(q(n_p, n_eta))
      DO ieta = 1, n_eta
        READ(lu,*) (q(j,ieta), j = 1, n_p)
      END DO
      CLOSE(unit=lu)
      DEALLOCATE(p)

    END SUBROUTINE read_quantile_table

  END SUBROUTINE setup_nlcs_quantile_tables



  ! chi (chi_gamma / 2) of an NLCS photon from
  ! the quantile table q (nlcs_q_up or nlcs_q_down) at the random
  ! number p in [0, 1). The rule is that of find_value_from_table_alt
  ! with the probabilities common to all rows: log10 chi is linear in p
  ! between the bracketing probabilities, in each of the two rows
  ! around eta, and the two are blended linearly in log10 eta. Below the
  ! table chi scales as eta^2 from its first row, as for
  ! ksi_sokolov.table; above it the last row is used.
  FUNCTION nlcs_quantile_chi(p, eta, q) RESULT(chi)

    REAL(num), INTENT(IN) :: p, eta
    REAL(num), INTENT(IN) :: q(:,:)
    REAL(num) :: chi
    REAL(num) :: x, fx, fp, y_lt, y_gt, scale
    INTEGER :: i, j1, j2, jm
    LOGICAL, SAVE :: warning = .TRUE.

    x = LOG10(eta)
    scale = 1.0_num
    IF (x < nlcs_q_log_eta_min) THEN
      scale = (eta / 10.0_num**nlcs_q_log_eta_min)**2
      x = nlcs_q_log_eta_min
    ELSE IF (x > nlcs_q_log_eta_max) THEN
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'eta above the range of the NLCS quantile tables.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      x = nlcs_q_log_eta_max
    END IF

    ! Row i and i + 1 around x, 1 <= i <= n_eta - 1, and 0 <= fx <= 1
    fx = (x - nlcs_q_log_eta_min) / nlcs_q_dlog
    i = MIN(INT(fx) + 1, nlcs_q_n_eta - 1)
    fx = fx - REAL(i - 1, num)

    ! Bisection: nlcs_q_p(j1) <= p < nlcs_q_p(j2) = nlcs_q_p(j1 + 1);
    ! nlcs_q_p runs from 0 to 1
    j1 = 1
    j2 = nlcs_q_n_p
    DO WHILE (j2 - j1 > 1)
      jm = (j1 + j2) / 2
      IF (nlcs_q_p(jm) <= p) THEN
        j1 = jm
      ELSE
        j2 = jm
      END IF
    END DO
    fp = (p - nlcs_q_p(j1)) / (nlcs_q_p(j2) - nlcs_q_p(j1))

    y_lt = (1.0_num - fp) * q(j1,i) + fp * q(j2,i)
    y_gt = (1.0_num - fp) * q(j1,i+1) + fp * q(j2,i+1)
    chi = scale * 10.0_num**((1.0_num - fx) * y_lt + fx * y_gt)

  END FUNCTION nlcs_quantile_chi



  ! Read the R1 and R2 columns of bessel_k.table
  ! the absolute K columns that follow on each row are skipped.
  ! The grid is uniform in log10(x).
  ! Called only with use_spin_qed.
  SUBROUTINE setup_bessel_k_table

    INTEGER :: i
    REAL(num) :: log_x, dlog, buf(2)

    IF (rank == 0) THEN
      OPEN(unit=lu, file=TRIM(qed_table_location) // '/bessel_k.table', &
          status='OLD')
      READ(lu,*) n_bessel, bessel_log_x_min, bessel_log_x_max
      ALLOCATE(bessel_log_r1(n_bessel), bessel_log_r2(n_bessel))
      dlog = (bessel_log_x_max - bessel_log_x_min) / (n_bessel - 1)
      DO i = 1, n_bessel
        READ(lu,*) log_x, bessel_log_r1(i), bessel_log_r2(i)
        IF (ABS(log_x - bessel_log_x_min - (i - 1) * dlog) &
            > 1.0e-8_num) CALL abort_code(c_err_bad_value)
      END DO
      CLOSE(unit=lu)
    END IF

    CALL MPI_BCAST(n_bessel, 1, MPI_INTEGER, 0, comm, errcode)
    buf = (/ bessel_log_x_min, bessel_log_x_max /)
    CALL MPI_BCAST(buf, 2, mpireal, 0, comm, errcode)
    bessel_log_x_min = buf(1)
    bessel_log_x_max = buf(2)
    IF (rank /= 0) ALLOCATE(bessel_log_r1(n_bessel), bessel_log_r2(n_bessel))
    CALL MPI_BCAST(bessel_log_r1, n_bessel, mpireal, 0, comm, errcode)
    CALL MPI_BCAST(bessel_log_r2, n_bessel, mpireal, 0, comm, errcode)
    bessel_idlog = (n_bessel - 1) / (bessel_log_x_max - bessel_log_x_min)

  END SUBROUTINE setup_bessel_k_table



  ! nbw_polarisation_asymmetry.table (A(chi), EPOCH's photon chi)
  ! The grid is uniform in log10(chi)
  ! Called only with use_spin_qed.
  SUBROUTINE setup_nbw_asymmetry_table

    INTEGER :: i
    REAL(num) :: log_chi, dlog, buf(2)

    IF (rank == 0) THEN
      OPEN(unit=lu, file=TRIM(qed_table_location) &
          // '/nbw_polarisation_asymmetry.table', status='OLD')
      READ(lu,*) n_nbw_asym, nbw_asym_log_chi_min, nbw_asym_log_chi_max
      ALLOCATE(nbw_asym_log_a(n_nbw_asym))
      dlog = (nbw_asym_log_chi_max - nbw_asym_log_chi_min) &
          / (n_nbw_asym - 1)
      DO i = 1, n_nbw_asym
        READ(lu,*) log_chi, nbw_asym_log_a(i)
        IF (ABS(log_chi - nbw_asym_log_chi_min - (i - 1) * dlog) &
            > 1.0e-8_num) CALL abort_code(c_err_bad_value)
      END DO
      CLOSE(unit=lu)
    END IF

    CALL MPI_BCAST(n_nbw_asym, 1, MPI_INTEGER, 0, comm, errcode)
    buf = (/ nbw_asym_log_chi_min, nbw_asym_log_chi_max /)
    CALL MPI_BCAST(buf, 2, mpireal, 0, comm, errcode)
    nbw_asym_log_chi_min = buf(1)
    nbw_asym_log_chi_max = buf(2)
    IF (rank /= 0) ALLOCATE(nbw_asym_log_a(n_nbw_asym))
    CALL MPI_BCAST(nbw_asym_log_a, n_nbw_asym, mpireal, 0, comm, errcode)
    nbw_asym_idlog = (n_nbw_asym - 1) &
        / (nbw_asym_log_chi_max - nbw_asym_log_chi_min)

  END SUBROUTINE setup_nbw_asymmetry_table



  ! Read nbw_energy_split_par.table and nbw_energy_split_perp.table
  ! which share one grid: uniform in log10(chi), electron fraction
  ! uniform in [0, 1]. Called only with use_spin_qed.
  SUBROUTINE setup_nbw_split_tables

    INTEGER :: i, n_chi, n_split
    REAL(num) :: buf(2), lmin, lmax

    IF (rank == 0) THEN
      OPEN(unit=lu, file=TRIM(qed_table_location) &
          // '/nbw_energy_split_par.table', status='OLD')
      READ(lu,*) n_nbw_chi, n_nbw_split, buf(1), buf(2)
      CLOSE(unit=lu)
    END IF
    CALL MPI_BCAST(n_nbw_chi, 1, MPI_INTEGER, 0, comm, errcode)
    CALL MPI_BCAST(n_nbw_split, 1, MPI_INTEGER, 0, comm, errcode)
    CALL MPI_BCAST(buf, 2, mpireal, 0, comm, errcode)

    ALLOCATE(nbw_log_chi(n_nbw_chi), nbw_split(n_nbw_split))
    DO i = 1, n_nbw_chi
      nbw_log_chi(i) = buf(1) + (buf(2) - buf(1)) * REAL(i - 1, num) &
          / REAL(n_nbw_chi - 1, num)
    END DO
    DO i = 1, n_nbw_split
      nbw_split(i) = REAL(i - 1, num) / REAL(n_nbw_split - 1, num)
    END DO
    CALL read_split_table('nbw_energy_split_par.table', nbw_cdf_par)
    CALL read_split_table('nbw_energy_split_perp.table', nbw_cdf_perp)

  CONTAINS

    SUBROUTINE read_split_table(name, cdf)

      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL(num), ALLOCATABLE, INTENT(INOUT) :: cdf(:,:)
      INTEGER :: ichi, j

      ALLOCATE(cdf(n_nbw_chi, n_nbw_split))
      IF (rank == 0) THEN
        OPEN(unit=lu, file=TRIM(qed_table_location) // '/' // name, &
            status='OLD')
        READ(lu,*) n_chi, n_split, lmin, lmax
        ! Both tables must share one grid
        IF (n_chi /= n_nbw_chi .OR. n_split /= n_nbw_split &
            .OR. ABS(lmin - buf(1)) > 1.0e-8_num &
            .OR. ABS(lmax - buf(2)) > 1.0e-8_num) THEN
          PRINT*, '*** ERROR ***'
          PRINT*, TRIM(name), ' is not on the grid of ', &
              'nbw_energy_split_par.table'
          CALL abort_code(c_err_bad_value)
        END IF
        DO ichi = 1, n_nbw_chi
          READ(lu,*) (cdf(ichi,j), j = 1, n_nbw_split)
        END DO
        CLOSE(unit=lu)
      END IF
      CALL MPI_BCAST(cdf, n_nbw_chi * n_nbw_split, mpireal, 0, comm, &
          errcode)

    END SUBROUTINE read_split_table

  END SUBROUTINE setup_nbw_split_tables



  ! Replace the T(chi) column of pairprod.table (Erber's approximation
  ! 3.3d) with the exact LCFA one from pairprod_lcfa.table
  ! (Erber 3.3b), which must share pairprod.table's
  ! grid. pairprod.table's trident column is untouched.
  SUBROUTINE read_pairprod_lcfa(n, log_chi, log_t)

    INTEGER, INTENT(IN) :: n
    REAL(num), INTENT(IN) :: log_chi(n)
    REAL(num), INTENT(INOUT) :: log_t(n)
    INTEGER :: i, n_file
    REAL(num) :: lmin, lmax, lc

    OPEN(unit=lu, file=TRIM(qed_table_location) // '/pairprod_lcfa.table', &
        status='OLD')
    READ(lu,*) n_file, lmin, lmax
    IF (n_file /= n) THEN
      PRINT*, '*** ERROR ***'
      PRINT*, 'pairprod_lcfa.table is not on the grid of pairprod.table'
      CALL abort_code(c_err_bad_value)
    END IF
    DO i = 1, n
      READ(lu,*) lc, log_t(i)
      IF (ABS(lc - log_chi(i)) > 1.0e-8_num) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'pairprod_lcfa.table is not on the grid of pairprod.table'
        CALL abort_code(c_err_bad_value)
      END IF
    END DO
    CLOSE(unit=lu)

  END SUBROUTINE read_pairprod_lcfa



  SUBROUTINE deallocate_spin_tables
    IF (ALLOCATED(amm_log_f)) DEALLOCATE(amm_log_f)
    IF (ALLOCATED(asym_log_a)) DEALLOCATE(asym_log_a)
    IF (ALLOCATED(nlcs_cdf_up)) DEALLOCATE(nlcs_cdf_up)
    IF (ALLOCATED(nlcs_cdf_down)) DEALLOCATE(nlcs_cdf_down)
    IF (ALLOCATED(nlcs_q_p)) DEALLOCATE(nlcs_q_p)
    IF (ALLOCATED(nlcs_q_up)) DEALLOCATE(nlcs_q_up)
    IF (ALLOCATED(nlcs_q_down)) DEALLOCATE(nlcs_q_down)
    IF (ALLOCATED(bessel_log_r1)) DEALLOCATE(bessel_log_r1)
    IF (ALLOCATED(bessel_log_r2)) DEALLOCATE(bessel_log_r2)
    IF (ALLOCATED(nbw_asym_log_a)) DEALLOCATE(nbw_asym_log_a)
    IF (ALLOCATED(nbw_log_chi)) DEALLOCATE(nbw_log_chi)
    IF (ALLOCATED(nbw_split)) DEALLOCATE(nbw_split)
    IF (ALLOCATED(nbw_cdf_par)) DEALLOCATE(nbw_cdf_par)
    IF (ALLOCATED(nbw_cdf_perp)) DEALLOCATE(nbw_cdf_perp)
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



  FUNCTION nlcs_asymmetry(eta) RESULT(a)

    ! A(eta) = (W_up - W_down) / (W_up + W_down): log-log linear
    ! interpolation inside the table, the leading order A ~ eta below
    ! it (A -> 3 eta / 10), clamped above it
    REAL(num), INTENT(IN) :: eta
    REAL(num) :: a, x, w
    INTEGER :: i

    IF (eta <= 0.0_num) THEN
      a = 0.0_num
      RETURN
    END IF

    x = LOG10(eta)
    IF (x <= asym_log_eta_min) THEN
      a = 10.0_num**(asym_log_a(1) + x - asym_log_eta_min)
    ELSE IF (x >= asym_log_eta_max) THEN
      a = 10.0_num**asym_log_a(n_asym)
    ELSE
      w = (x - asym_log_eta_min) * asym_idlog
      i = INT(w)
      w = w - i
      i = i + 1
      a = 10.0_num**((1.0_num - w) * asym_log_a(i) + w * asym_log_a(i+1))
    END IF

  END FUNCTION nlcs_asymmetry



  FUNCTION nlcs_down_probability(zeta, a) RESULT(p_down)

    ! At an emission, the photon energy is drawn from the zeta = -1
    ! (down) spectrum with this probability and from the zeta = +1 (up)
    ! spectrum otherwise. With W = W_up (1 + zeta)/2 + W_down (1 - zeta)/2
    ! and W_up/down = W_bar (1 +- A), the weight of the down part is
    ! (1 - zeta)(1 - A) / (2 (1 + A zeta)). Exact, since the spectrum
    ! is linear in zeta (Qian's eqn. 3.3).
    REAL(num), INTENT(IN) :: zeta, a
    REAL(num) :: p_down

    p_down = 0.5_num * (1.0_num - zeta) * (1.0_num - a) &
        / (1.0_num + a * zeta)

  END FUNCTION nlcs_down_probability



  FUNCTION nbw_asymmetry(chi) RESULT(a)

    ! A(chi) = (W_perp - W_par) / (W_perp + W_par) of NBW, chi beingEPOCH's
    ! photon parameter: log-log linear interpolation inside the table,
    ! clamped outside it (A -> 1/3 as chi -> 0, A -> 1/5 as chi -> Inf)
    REAL(num), INTENT(IN) :: chi
    REAL(num) :: a, x, w
    INTEGER :: i

    IF (chi <= 0.0_num) THEN
      a = 10.0_num**nbw_asym_log_a(1)
      RETURN
    END IF

    x = LOG10(chi)
    IF (x <= nbw_asym_log_chi_min) THEN
      a = 10.0_num**nbw_asym_log_a(1)
    ELSE IF (x >= nbw_asym_log_chi_max) THEN
      a = 10.0_num**nbw_asym_log_a(n_nbw_asym)
    ELSE
      w = (x - nbw_asym_log_chi_min) * nbw_asym_idlog
      i = INT(w)
      w = w - i
      i = i + 1
      a = 10.0_num**((1.0_num - w) * nbw_asym_log_a(i) &
          + w * nbw_asym_log_a(i+1))
    END IF

  END FUNCTION nbw_asymmetry



  FUNCTION nbw_par_probability(xi1, a) RESULT(p_par)

    ! At a decay, the electron fraction is drawn from the xi1 = +1 (par)
    ! spectrum with this probability and from the xi1 = -1 (perp)
    ! spectrum otherwise. With W = W_par (1 + xi1)/2 + W_perp (1 - xi1)/2
    ! and W_par/perp = W_bar (1 -+ A), the weight of the par part is
    ! (1 + xi1)(1 - A) / (2 (1 - A xi1)). Exact, since the spectrum is
    ! linear in xi1. (Qian eqn. 3.11)
    REAL(num), INTENT(IN) :: xi1, a
    REAL(num) :: p_par

    p_par = 0.5_num * (1.0_num + xi1) * (1.0_num - a) &
        / (1.0_num - a * xi1)

  END FUNCTION nbw_par_probability



  SUBROUTINE bessel_ratios(x, r1, r2)

    ! R1 = K_1/3(x)/K_2/3(x) and R2 = IntK_1/3(x)/K_2/3(x), both on
    ! [0, 1]: linear interpolation of log10 R in log10 x inside
    ! bessel_k.table, the two-term expansions of gen_bessel_k_table.py
    ! outside it (they join the table to ~1e-7 without rescaling)
    !   x < x_min:  R1 = e1 x^(1/3) - 3x/2,
    !               R2 = d1 x^(2/3) - 3 e1 x^(4/3)/2,
    !   x > x_max:  R1 = 1 - 1/(6x) + 7/(72x^2),  R2 = 1 - 2/(3x) + 1/x^2.
    ! x <= 0 (lam = 0) gives 0 and x = +Inf (lam = 1) gives 1.
    REAL(num), INTENT(IN) :: x
    REAL(num), INTENT(OUT) :: r1, r2
    REAL(num) :: lx, w, y, x13
    INTEGER :: i

    IF (x <= 0.0_num) THEN
      r1 = 0.0_num
      r2 = 0.0_num
      RETURN
    END IF

    lx = LOG10(x)
    IF (lx <= bessel_log_x_min) THEN
      x13 = x**(1.0_num / 3.0_num)
      r1 = c_bessel_e1 * x13 - 1.5_num * x
      r2 = x13 * (c_bessel_d1 * x13 - 1.5_num * c_bessel_e1 * x)
    ELSE IF (lx >= bessel_log_x_max) THEN
      y = 1.0_num / x
      r1 = 1.0_num - y * (1.0_num / 6.0_num - y * 7.0_num / 72.0_num)
      r2 = 1.0_num - y * (2.0_num / 3.0_num - y)
    ELSE
      w = (lx - bessel_log_x_min) * bessel_idlog
      i = INT(w)
      w = w - i
      i = i + 1
      r1 = 10.0_num**((1.0_num - w) * bessel_log_r1(i) &
          + w * bessel_log_r1(i+1))
      r2 = 10.0_num**((1.0_num - w) * bessel_log_r2(i) &
          + w * bessel_log_r2(i+1))
    END IF

  END SUBROUTINE bessel_ratios



  ! nlcs_w_g and nlcs_photon_f are the NLCS final-state functions 
  ! (derived from Torgrimsson 2021. c.f. Qian's Appendix A.1. Note
  ! that the basis in Qian's is slightly different from EPOCH's).
  ! Inputs: lam = photon energy / lepton energy, 0 <= lam < 1;
  ! r1, r2 from bessel_ratios at x = 2 lam / (3 chi (1 - lam)), chi 
  ! (EPOCH's eta); qs the sign of the emitter's charge;
  ! si, sf the initial and final spins as components on the lepton axes
  ! (eps_hat, beta_hat, v_hat) of lepton_spin_axes. 
  ! For a photon of energy lam, summed over final spins and
  ! photon polarisations, the emission density is proportional to w;
  ! resolved in sf and the photon Stokes vector xi_obs it is
  ! proportional to F0 + xi_obs . F, with F0 = w + sf . g.
  ! For reference see Qian's 3.2.1 and appendix A.1, 
  ! and Torgrimsson 2021 (42)-(50). 
  SUBROUTINE nlcs_w_g(lam, r1, r2, qs, si, w, g)

    ! w and g: the final spin, summed over photon polarisations, has
    ! the Stokes-like vector g / w (|g| <= w)
    REAL(num), INTENT(IN) :: lam, r1, r2, qs, si(3)
    REAL(num), INTENT(OUT) :: w, g(3)
    REAL(num) :: s1, kap, a, c2

    s1 = 1.0_num - lam
    kap = s1 + 1.0_num / s1
    a = lam * qs * r1
    c2 = 2.0_num - r2

    w = kap - r2 + a * si(2)
    g(1) = c2 * si(1)
    g(2) = c2 * si(2) + a / s1
    g(3) = (kap - r2 * (lam + s1 * s1) / s1) * si(3)

  END SUBROUTINE nlcs_w_g



  SUBROUTINE nlcs_photon_f(lam, r1, r2, qs, si, sf, f0, f)

    ! F0 and F: given the final spin sf, the photon has the Stokes vector
    ! xi_obs = F / F0 (|F| <= F0) on the axes (xi_1 linear along E_hat,
    ! xi_2 circular, xi_3 diagonal). With sf = 0 this returns F0 = w and
    ! F = F_bar, the photon polarisation summed over final spins.
    REAL(num), INTENT(IN) :: lam, r1, r2, qs, si(3), sf(3)
    REAL(num), INTENT(OUT) :: f0, f(3)
    REAL(num) :: s1, kap, a, c2, h, d, p

    s1 = 1.0_num - lam
    kap = s1 + 1.0_num / s1
    a = lam * qs * r1
    c2 = 2.0_num - r2
    h = 0.5_num * kap
    d = 0.5_num * lam * lam * r2 / s1
    p = lam * (2.0_num - lam) / s1

    f0 = kap - r2 + a * (si(2) + sf(2) / s1) &
        + c2 * (sf(1) * si(1) + sf(2) * si(2)) &
        + (kap - r2 * (lam + s1 * s1) / s1) * sf(3) * si(3)
    f(1) = 1.0_num + a * (sf(2) + si(2) / s1) + (h + d) * sf(2) * si(2) &
        + (h - d) * sf(1) * si(1) + sf(3) * si(3)
    f(2) = a * (sf(2) * si(3) / s1 - sf(3) * si(2)) &
        + (p - lam * r2 / s1) * sf(3) + (p - lam * r2) * si(3)
    f(3) = -a * (sf(1) + si(1) / s1) - (d + 0.5_num * p) * sf(2) * si(1) &
        + (0.5_num * p - d) * sf(1) * si(2)

  END SUBROUTINE nlcs_photon_f



  SUBROUTINE nbw_w0_j_m(s3, r1, r2, xi, w0, je, jp, m)

    ! NBW pair-spin functions 
    ! (derived from Torgrimsson 2021. c.f. Qian's Appendix A.2. Note
    ! that the basis in Qian's is slightly different from EPOCH's)
    ! Inputs: s3 = positron energy / photon energy,
    ! 0 < s3 < 1 (the electron has s2 = 1 - s3); r1, r2 from
    ! bessel_ratios at x = 1 / (3 chi s2 s3), chi EPOCH's photon
    ! parameter; xi the photon's observed Stokes vector (xi_obs). Spins
    ! are components on (eps_hat, beta_hat, v_hat), eps_hat the
    ! normalised E_perp + v x B with v along the photon (E_hat of
    ! xi_obs is -eps_hat; its sign does not matter). For the electron
    ! fraction s2, the probability of the pair spins (se, sp) is
    ! proportional to
    !   w0 + se . je + sp . jp + sp . (m se),
    ! with m(a,c) multiplying sp_a se_c (positron index first; T21's
    ! positron label vector is -sp).
    ! Summed over the pair spins it is 4 w0, so w0 is the spectrum
    ! (rho - xi1) K_2/3 + IntK_1/3 in units of K_2/3. All outputs share
    ! the positive factor K_2/3 / (pi sqrt 3).
    REAL(num), INTENT(IN) :: s3, r1, r2, xi(3)
    REAL(num), INTENT(OUT) :: w0, je(3), jp(3), m(3,3)
    REAL(num) :: s2, p, rho, dd

    s2 = 1.0_num - s3
    p = s2 * s3
    rho = (s2 * s2 + s3 * s3) / p
    dd = (s3 - s2) / p

    w0 = r2 - xi(1) + rho
    je(1) = -r1 * xi(3) / s3
    je(2) = r1 * (xi(1) / s3 - 1.0_num / s2)
    je(3) = xi(2) * (r2 / s2 - dd)
    jp(1) = r1 * xi(3) / s2
    jp(2) = r1 * (1.0_num / s3 - xi(1) / s2)
    jp(3) = xi(2) * (r2 / s3 + dd)

    m(1,1) = r2 - 2.0_num + 0.5_num * xi(1) * (rho - r2 / p)
    m(1,2) = -0.5_num * xi(3) * (r2 / p + dd)
    m(1,3) = 0.0_num
    m(2,1) = 0.5_num * xi(3) * (dd - r2 / p)
    m(2,2) = r2 - 2.0_num + 0.5_num * xi(1) * (rho + r2 / p)
    m(2,3) = r1 * xi(2) / s3
    m(3,1) = 0.0_num
    m(3,2) = -r1 * xi(2) / s2
    ! r2 (1 - p) / p is r2 (s2 + s3^2) / (s2 s3) with s2 + s3 = 1
    m(3,3) = r2 * (1.0_num - p) / p + xi(1) - rho

  END SUBROUTINE nbw_w0_j_m



  SUBROUTINE photon_pol_frame(p_photon, e_field, b_field, n_hat, e_hat, &
      ax1, ax2, valid)

    ! Axes to convert a photon's Stokes vector between stored and
    ! observed (pol_stored_to_obs, pol_obs_to_stored): the direction
    ! n_hat, the local field axis e_hat (Algorithm A, electron
    ! convention: e_hat = -eps_hat, the E_hat of NBW) and the reference
    ! axes ax1, ax2 (Algorithm R). valid = .FALSE. without momentum or
    ! effective field; callers must then skip the conversion.
    REAL(num), INTENT(IN) :: p_photon(3), e_field(3), b_field(3)
    REAL(num), INTENT(OUT) :: n_hat(3), e_hat(3), ax1(3), ax2(3)
    LOGICAL, INTENT(OUT) :: valid
    REAL(num) :: b_hat(3)

    ax1 = 0.0_num
    ax2 = 0.0_num
    CALL compute_local_field_axes(p_photon(1), p_photon(2), p_photon(3), &
        e_field(1), e_field(2), e_field(3), b_field(1), b_field(2), &
        b_field(3), n_hat, e_hat, b_hat, valid)
    IF (valid) CALL reference_axis(n_hat, ax1, ax2)

  END SUBROUTINE photon_pol_frame



  SUBROUTINE lepton_spin_axes(ux, uy, uz, gamma_rel, e_field, b_field, &
      eps_hat, beta_hat, v_hat, valid)

    ! Physical spin axes of a lepton
    ! eps_hat = normalised E_perp + v x B (E_perp transverse to v),
    ! beta_hat = eps_hat x v_hat, v_hat. Unlike compute_local_field_axes
    ! this uses the actual velocity, so it also holds at moderate gamma.
    ! Inputs: u = p / (m c), the Lorentz factor and SI fields.
    REAL(num), INTENT(IN) :: ux, uy, uz, gamma_rel
    REAL(num), INTENT(IN) :: e_field(3), b_field(3)
    REAL(num), INTENT(OUT) :: eps_hat(3), beta_hat(3), v_hat(3)
    LOGICAL, INTENT(OUT) :: valid
    REAL(num) :: u(3), beta(3), electric(3), magnetic(3), crossed_e(3)
    REAL(num) :: umag, scale

    valid = .FALSE.
    eps_hat = 0.0_num
    beta_hat = 0.0_num
    v_hat = 0.0_num
    u = (/ ux, uy, uz /)
    umag = SQRT(DOT_PRODUCT(u, u))
    IF (umag <= TINY(1.0_num)) RETURN
    v_hat = u / umag
    beta = u / gamma_rel

    electric = e_field
    magnetic = c * b_field
    ! Scale first; a relative tolerance rejects round-off cancellation
    scale = MAX(MAXVAL(ABS(electric)), MAXVAL(ABS(magnetic)))
    IF (scale <= TINY(1.0_num)) RETURN
    electric = electric / scale
    magnetic = magnetic / scale
    ! E - v_hat (v_hat . E) + (v / c) x (c B)
    crossed_e = electric - v_hat * DOT_PRODUCT(v_hat, electric) &
        + (/ beta(2)*magnetic(3) - beta(3)*magnetic(2), &
             beta(3)*magnetic(1) - beta(1)*magnetic(3), &
             beta(1)*magnetic(2) - beta(2)*magnetic(1) /)
    scale = SQRT(DOT_PRODUCT(crossed_e, crossed_e))
    IF (scale <= 16.0_num * EPSILON(1.0_num)) RETURN

    eps_hat = crossed_e / scale
    beta_hat = (/ eps_hat(2)*v_hat(3) - eps_hat(3)*v_hat(2), &
                  eps_hat(3)*v_hat(1) - eps_hat(1)*v_hat(3), &
                  eps_hat(1)*v_hat(2) - eps_hat(2)*v_hat(1) /)
    valid = .TRUE.

  END SUBROUTINE lepton_spin_axes



  SUBROUTINE lepton_spin_projection(ux, uy, uz, gamma_rel, e_field, &
      b_field, spin, charge, zeta, n_hat, valid)

    ! n_hat = qs beta_hat (qs the sign of the charge) and zeta = s . n_hat:
    ! the only spin dependence of the NLCS rate, W = W_bar (1 + A zeta),
    ! and the axis of the non-emission update. zeta = 0, n_hat = 0 and
    ! valid = .FALSE. when the axes are undefined (no momentum or no
    ! effective field).
    REAL(num), INTENT(IN) :: ux, uy, uz, gamma_rel
    REAL(num), INTENT(IN) :: e_field(3), b_field(3), spin(3), charge
    REAL(num), INTENT(OUT) :: zeta, n_hat(3)
    LOGICAL, INTENT(OUT) :: valid
    REAL(num) :: eps_hat(3), beta_hat(3), v_hat(3)

    zeta = 0.0_num
    n_hat = 0.0_num
    CALL lepton_spin_axes(ux, uy, uz, gamma_rel, e_field, b_field, &
        eps_hat, beta_hat, v_hat, valid)
    IF (.NOT. valid) RETURN
    n_hat = SIGN(1.0_num, charge) * beta_hat
    zeta = DOT_PRODUCT(spin, n_hat)

  END SUBROUTINE lepton_spin_projection



  FUNCTION spin_delta_optical_depth(dtau0, a, zeta, second_order) &
      RESULT(dtau)

    ! Optical-depth decrement given the spin-averaged one, dtau0 =
    ! W_bar dt. 
    ! First order: dtau = dtau0 (1 + A zeta). 
    !  Second order: the exact two-state survival probability,
    ! dtau = dtau0 - ln(cosh x - zeta sinh x), x = A dtau0.
    ! The two agree for zeta = +-1 and differ at O(x^2) otherwise.
    ! See notes for theory
    REAL(num), INTENT(IN) :: dtau0, a, zeta
    LOGICAL, INTENT(IN) :: second_order
    REAL(num) :: dtau, x, y

    IF (second_order) THEN
      x = a * dtau0
      ! y = cosh x - zeta sinh x - 1, written without cancellation
      y = 2.0_num * SINH(0.5_num * x)**2 - zeta * SINH(x)
      ! ln(1 + y); the series keeps full precision for small y
      IF (ABS(y) < 1.0e-4_num) THEN
        dtau = dtau0 - y * (1.0_num - y * (0.5_num - y / 3.0_num))
      ELSE
        dtau = dtau0 - LOG(1.0_num + y)
      END IF
    ELSE
      dtau = dtau0 * (1.0_num + a * zeta)
    END IF

  END FUNCTION spin_delta_optical_depth



  FUNCTION spin_event_fraction(tau, dtau0, a, zeta) RESULT(g)

    ! Fraction g of a (sub)step at which an optical depth tau runs out
    ! for the exact two-state decrement (spin_second_order_depth),
    ! D(g) = g dtau0 - ln(cosh(g x) - zeta sinh(g x)),  x = a dtau0,
    ! for 0 < tau < D(1). D increases: D'(g) = dtau0 (1 + a zeta_g),
    ! zeta_g = (zeta - t) / (1 - t zeta), t = tanh(g x), the survivor's
    ! zeta at g. 
    REAL(num), INTENT(IN) :: tau, dtau0, a, zeta
    REAL(num) :: g, d, t, zeta_g
    INTEGER :: i

    g = tau / spin_delta_optical_depth(dtau0, a, zeta, .TRUE.)
    DO i = 1, 4
      d = spin_delta_optical_depth(g * dtau0, a, zeta, .TRUE.) - tau
      t = TANH(g * a * dtau0)
      zeta_g = (zeta - t) / (1.0_num - t * zeta)
      g = g - d / (dtau0 * (1.0_num + a * zeta_g))
      g = MIN(MAX(g, 0.0_num), 1.0_num)
    END DO

  END FUNCTION spin_event_fraction



  SUBROUTINE lepton_no_emission_update(spin, n_hat, x)

    ! Spin of a lepton that did not emit in the step: Qian's quantum
    ! non-transition update (section 3.3), T21 Eqs. (133)-(134).
    ! A brief explnation: Leptons whose spin
    ! state emits less (zeta = -1) are more likely to survive, so
    !  s' = [(zeta - t) n_hat + sqrt(1 - t^2) s_perp] / (1 - t zeta),
    ! where zeta = s . n_hat, s_perp = s - zeta n_hat, t = tanh x, with
    ! x = A W_bar dt (A the nlcs_asymmetry, W_bar dt the spin-averaged
    ! optical-depth decrement) and n_hat = qs beta_hat. This is the
    ! normalised survivor state whose norm is the survival probability
    ! used by spin_delta_optical_depth. Pure states stay pure, |s| <= 1
    ! is kept and states along +-n_hat do not change. Steps compose:
    ! n steps of x equal one step of n x. The rotation part of the loop
    ! is the a(chi) term of T-BMT and is not added here.
    ! See notes for moire theory
    REAL(num), INTENT(INOUT) :: spin(3)
    REAL(num), INTENT(IN) :: n_hat(3), x
    REAL(num) :: t, zeta, s_perp(3)

    t = TANH(x)
    zeta = DOT_PRODUCT(spin, n_hat)
    s_perp = spin - zeta * n_hat
    spin = ((zeta - t) * n_hat + SQRT(1.0_num - t * t) * s_perp) &
        / (1.0_num - t * zeta)

  END SUBROUTINE lepton_no_emission_update



  FUNCTION photon_no_decay_update(xi_obs, t, phi) RESULT(xi_new)

    ! The photon equivalence of the above
    ! Stokes vector of a photon that did not decay (NBW) in the step,
    ! in the observation basis: T21 Eqs. (147)-(150)
    ! The same algebra as lepton_no_emission_update with
    ! the Stokes axis xi_1 in place of n_hat: photons polarised parallel
    ! to E_hat (xi_1 = +1) decay less, so with t = tanh nu,
    ! nu = (W_perp - W_par) dt / 2 (dichroism),
    ! D = 1 + t xi_1,  xi_1' = (xi_1 + t) / D,
    ! and (xi_3, xi_2) rotate by the birefringence phase phi (T21 Eq.
    ! 149) and shrink by sqrt(1 - t^2) / D. Pure states stay pure;
    ! |xi| <= 1 is kept. Convert stored -> observed before and observed
    ! -> stored after (pol_stored_to_obs, pol_obs_to_stored).
    REAL(num), INTENT(IN) :: xi_obs(3), t
    REAL(num), INTENT(IN), OPTIONAL :: phi
    REAL(num) :: xi_new(3), d, amp, cphi, sphi

    cphi = 1.0_num
    sphi = 0.0_num
    IF (PRESENT(phi)) THEN
      cphi = COS(phi)
      sphi = SIN(phi)
    END IF
    d = 1.0_num + t * xi_obs(1)
    amp = SQRT(1.0_num - t * t) / d
    xi_new(1) = (xi_obs(1) + t) / d
    xi_new(3) = amp * (xi_obs(3) * cphi - xi_obs(2) * sphi)
    xi_new(2) = amp * (xi_obs(3) * sphi + xi_obs(2) * cphi)

  END FUNCTION photon_no_decay_update



  ! Four helpers
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


  SUBROUTINE nlcs_final_states(ux, uy, uz, gamma_rel, e_field, b_field, &
      charge, lam, eta, want_photon, spin, xi_stored)

    ! Final lepton spin and photon Stokes vector of an NLCS emission
    ! (Qian's thesis section 3, Appendix A.1, in EPOCH's bases).
    ! Inputs: the lepton before recoil (u = p / (m c) and its Lorentz
    ! factor), the fields at it, the charge (only the sign is used),
    ! the photon energy fraction lam and EPOCH's eta (aka chi_e).
    ! spin is always updated, so the lepton changes even when the photon
    ! is below photon_energy_min or not sampled. xi_stored is computed
    ! only if want_photon; otherwise, or if the axes are undefined, it is
    ! zero (and spin is unchanged in the latter case).
    !  1. Lepton axes (eps_hat, beta_hat, v_hat) from lepton_spin_axes;
    !     si = (s.eps_hat, s.beta_hat, s.v_hat).
    !  2. R1, R2 at x = 2 lam / (3 eta (1 - lam)); w, g from nlcs_w_g.
    !  3. Final spin. spin_final_states (qed block) selects
    !     collapse: s_f = +-g_hat (the SQA) with P(+) = (1 + |g|/w)/2,
    !               one more random number is used;
    !     mean:     s_f = g/w.
    !     The collapse (basically Qian's method) vs mean is still under
    !     debate. The default is collapse method. Note that the collapse
    !     method correlates the polarisation with the final lepton spin.
    !     Personally, I think the mean method is better in the sense that
    !     it preserves the Bloch vector of the density matrix and carries
    !     it in the following simulation. That is, we do not perform a 
    !     quantum measurement and do not collapse the spin/pol to a 
    !     given axis. Note that in the mean method, the vector length is
    !     not in general kept at unity, in which case the physical meaning
    !     of it becomes the degree of polarisation instead of physical
    !     spin/pol, which is ok but caution is needed.
    !  4. Photon (nlcs_photon_f). 
    !      collapse: F, F0 given the s_f just drawn,  
    !     xi_obs = +-F_hat (the PQA) with P(+) = (1 + |F|/F0)/2,
    !     a second random number is used 
    !     mean: xi_obs = F_bar/w (summed over s_f; not F/F0 at s_f = g/w).
    !     note that 2 more random numbers are used in the collapse method
    !     which slows the program down.
    !  5. The photon moves along v_hat with E_hat = qs eps_hat (its sign
    !     does not matter): xi_stored from reference_axis and
    !     pol_obs_to_stored.
    ! Averaged over the random numbers, collapse gives exactly the mean
    ! spin and polarisation; it also keeps their correlation per event.
    ! Mean-mode vectors are clamped to |s_f|, |xi| <= 1 against table
    ! round-off (pure states sit at the bound).
    REAL(num), INTENT(IN) :: ux, uy, uz, gamma_rel
    REAL(num), INTENT(IN) :: e_field(3), b_field(3), charge, lam, eta
    LOGICAL, INTENT(IN) :: want_photon
    REAL(num), INTENT(INOUT) :: spin(3)
    REAL(num), INTENT(OUT) :: xi_stored(3)
    REAL(num), PARAMETER :: zero3(3) = 0.0_num
    REAL(num) :: eps_hat(3), beta_hat(3), v_hat(3), ax1(3), ax2(3)
    REAL(num) :: si(3), sf(3), g(3), f(3), xi_obs(3)
    REAL(num) :: qs, lam_c, r1, r2, w, f0, gn, fn, p_plus
    LOGICAL :: valid, collapse

    xi_stored = 0.0_num
    IF (eta <= 0.0_num) RETURN
    CALL lepton_spin_axes(ux, uy, uz, gamma_rel, e_field, b_field, &
        eps_hat, beta_hat, v_hat, valid)
    IF (.NOT. valid) RETURN

    collapse = spin_final_states == c_spin_final_collapse
    qs = SIGN(1.0_num, charge)
    si = (/ DOT_PRODUCT(spin, eps_hat), DOT_PRODUCT(spin, beta_hat), &
        DOT_PRODUCT(spin, v_hat) /)
    ! lam = 1 is the end point of the spectrum; keep 1 - lam > 0
    lam_c = MIN(MAX(lam, 0.0_num), 1.0_num - EPSILON(1.0_num))
    CALL bessel_ratios(2.0_num * lam_c / (3.0_num * eta * (1.0_num - lam_c)), &
        r1, r2)

    ! Final lepton spin
    CALL nlcs_w_g(lam_c, r1, r2, qs, si, w, g)
    gn = SQRT(DOT_PRODUCT(g, g))
    IF (collapse) THEN
      IF (gn > TINY(1.0_num)) THEN
        sf = g / gn
      ELSE
        ! Unpolarised: any axis, with equal probabilities
        sf = (/ 0.0_num, 1.0_num, 0.0_num /)
      END IF
      IF (random() >= 0.5_num * (1.0_num + MIN(gn / w, 1.0_num))) sf = -sf
    ELSE
      sf = g / MAX(w, gn)
    END IF

    ! Photon polarisation
    IF (want_photon) THEN
      IF (collapse) THEN
        CALL nlcs_photon_f(lam_c, r1, r2, qs, si, sf, f0, f)
        fn = SQRT(DOT_PRODUCT(f, f))
        IF (fn > TINY(1.0_num)) THEN
          xi_obs = f / fn
        ELSE
          xi_obs = (/ 1.0_num, 0.0_num, 0.0_num /)
        END IF
        ! F0 >= |F| > 0 up to round-off; P(+) = 1 at the bound
        p_plus = 1.0_num
        IF (fn < f0) p_plus = 0.5_num * (1.0_num + fn / f0)
        IF (random() >= p_plus) xi_obs = -xi_obs
      ELSE
        CALL nlcs_photon_f(lam_c, r1, r2, qs, si, zero3, f0, f)
        fn = SQRT(DOT_PRODUCT(f, f))
        xi_obs = f / MAX(f0, fn)
      END IF
      CALL reference_axis(v_hat, ax1, ax2)
      xi_stored = pol_obs_to_stored(ax1, ax2, eps_hat, xi_obs)
    END IF

    ! Back to lab axes; the final lepton moves along v_hat too
    spin = sf(1) * eps_hat + sf(2) * beta_hat + sf(3) * v_hat

  END SUBROUTINE nlcs_final_states



  SUBROUTINE nbw_final_states(p_photon, e_field, b_field, xi_obs, chi, &
      s_electron, spin_electron, spin_positron)

    ! Spins of the pair created by an NBW decay (Qian's thesis section
    ! 3, Appendix A.2, in EPOCH's bases). 
    ! Inputs: the photon momentum, the fields at it, its observed 
    !  Stokes vector xi_obs (before this step's survivor update), 
    !  EPOCH's photon chi (note the 1/2 prefactor) and the electron
    !  energy fraction s_electron. Outputs are lab vectors; both are zero
    !  (unpolarised) if the axes are undefined.
    !  1. Axes (eps_hat, beta_hat, v_hat) along the photon;
    !     R1, R2 at x = 1 / (3 chi s2 s3); w0, je, jp, m (nbw_w0_j_m).
    !  2. Electron. spin_final_states (qed block) selects
    !     collapse: se = +-je_hat (the SQA) with P(+) = (1 + |je|/w0)/2,
    !               one random number;
    !     mean:     se = je/w0.
    !  3. Positron. collapse: given se, the weight is w0 + se . je and
    !     the vector k = jp + m se; sp = +-k_hat with P(+) = (1 + |k| /
    !     (w0 + se . je))/2, a second random number; mean: sp = jp/w0.
    ! Averaged over the random numbers, collapse gives exactly the mean
    ! spins. It also keeps the pair correlation along the electron's
    ! SQA, <sp_a se_c> = (m je_hat)_a je_hat_c / w0; correlations across
    ! it (the rest of m / w0) cannot be carried by two separate states.
    ! Mean mode keeps no correlation. Its vectors are clamped to <= 1.
    REAL(num), INTENT(IN) :: p_photon(3), e_field(3), b_field(3)
    REAL(num), INTENT(IN) :: xi_obs(3), chi, s_electron
    REAL(num), INTENT(OUT) :: spin_electron(3), spin_positron(3)
    REAL(num) :: n_hat(3), eps_hat(3), beta_hat(3)
    REAL(num) :: je(3), jp(3), m(3,3), se(3), sp(3), k(3)
    REAL(num) :: s2, s3, r1, r2, w0, jn, kn, wk, p_plus
    LOGICAL :: valid

    spin_electron = 0.0_num
    spin_positron = 0.0_num
    IF (chi <= 0.0_num) RETURN
    ! charge_sign = +1 gives eps_hat itself and beta_hat = eps_hat x n
    CALL compute_local_field_axes(p_photon(1), p_photon(2), p_photon(3), &
        e_field(1), e_field(2), e_field(3), b_field(1), b_field(2), &
        b_field(3), n_hat, eps_hat, beta_hat, valid, 1.0_num)
    IF (.NOT. valid) RETURN

    ! The sampled fraction may reach 0 or 1; keep s2 s3 > 0
    s2 = MIN(MAX(s_electron, EPSILON(1.0_num)), 1.0_num - EPSILON(1.0_num))
    s3 = 1.0_num - s2
    CALL bessel_ratios(1.0_num / (3.0_num * chi * s2 * s3), r1, r2)
    CALL nbw_w0_j_m(s3, r1, r2, xi_obs, w0, je, jp, m)

    jn = SQRT(DOT_PRODUCT(je, je))
    IF (spin_final_states == c_spin_final_collapse) THEN
      ! Electron along its SQA
      IF (jn > TINY(1.0_num)) THEN
        se = je / jn
      ELSE
        ! Unpolarised: any axis, with equal probabilities
        se = (/ 0.0_num, 1.0_num, 0.0_num /)
      END IF
      IF (random() >= 0.5_num * (1.0_num + MIN(jn / w0, 1.0_num))) se = -se
      ! Positron given the electron
      k = jp + MATMUL(m, se)
      wk = w0 + DOT_PRODUCT(se, je)
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
      se = je / MAX(w0, jn)
      sp = jp / MAX(w0, SQRT(DOT_PRODUCT(jp, jp)))
    END IF

    spin_electron = se(1) * eps_hat + se(2) * beta_hat + se(3) * n_hat
    spin_positron = sp(1) * eps_hat + sp(2) * beta_hat + sp(3) * n_hat

  END SUBROUTINE nbw_final_states





#endif

END MODULE spin
