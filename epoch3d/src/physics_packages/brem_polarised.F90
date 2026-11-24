! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.

MODULE brem_polarised

#ifdef SPIN
  USE constants
  USE random_generator, ONLY: random
  USE spin, ONLY: reference_axis, pol_obs_to_stored

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bh_brem_setup, bh_brem_active, bh_brem_rate
  PUBLIC :: bh_brem_sample, bh_brem_transfer

  ! Olsen and Maximon, Phys. Rev. 114, 887 (1959), Eqs. 5.7,
  ! 6.29-6.32 and 7.2. The angular measure is xi**2 d(u**2) = d(xi),
  ! which reproduces Tsai's complete-screening spectrum. The upper u
  ! cut (u_max = gamma_min, so sin(theta) <= 1) loses at most 2e-3 of
  ! the rate for Z <= 79 and gamma <= 1e5. The lower photon energy is
  ! a physical event cut, distinct from the output cut.
  INTEGER, PARAMETER :: n_rate = 256, n_quad = 32
  REAL(num), PARAMETER :: gamma_min = 50.0_num
  REAL(num), PARAMETER :: gamma_max = 1.0e8_num
  REAL(num), PARAMETER :: gamma_final_min = 10.0_num
  REAL(num), PARAMETER :: u_max = 50.0_num
  REAL(num), PARAMETER :: xi_min = 1.0_num / (1.0_num + u_max**2)
  REAL(num), PARAMETER :: re2 = (alpha * h_bar / (m0 * c))**2
  REAL(num), PARAMETER :: mol_alpha(3) = &
      (/ 0.10_num, 0.55_num, 0.35_num /)
  REAL(num), PARAMETER :: mol_b(3) = &
      (/ 6.0_num, 1.20_num, 0.30_num /)

  REAL(num), SAVE :: rate_table(100,n_rate), quad_x(n_quad)
  REAL(num), SAVE :: quad_w(n_quad), energy_cut
  LOGICAL, SAVE :: z_present(100) = .FALSE.

CONTAINS

  SUBROUTINE bh_brem_setup(z_values, cut)

    INTEGER, INTENT(IN) :: z_values(:)
    REAL(num), INTENT(IN) :: cut
    INTEGER :: iz, i, z
    REAL(num) :: log_step, gamma_i

    energy_cut = cut
    z_present = .FALSE.
    rate_table = 0.0_num
    CALL bh_brem_quadrature
    log_step = LOG(gamma_max / gamma_min) / REAL(n_rate - 1, num)
    DO iz = 1, SIZE(z_values)
      z = z_values(iz)
      IF (z < 1 .OR. z > 100) CYCLE
      z_present(z) = .TRUE.
      DO i = 1, n_rate
        gamma_i = gamma_min * EXP(REAL(i - 1, num) * log_step)
        rate_table(z,i) = bh_brem_integral(gamma_i, z)
      END DO
    END DO

  END SUBROUTINE bh_brem_setup



  PURE LOGICAL FUNCTION bh_brem_active(energy)

    REAL(num), INTENT(IN) :: energy
    REAL(num) :: gamma_i

    gamma_i = energy / m0c2
    bh_brem_active = gamma_i >= gamma_min .AND. gamma_i <= gamma_max

  END FUNCTION bh_brem_active



  FUNCTION bh_brem_rate(energy, z) RESULT(sigma)

    REAL(num), INTENT(IN) :: energy
    INTEGER, INTENT(IN) :: z
    REAL(num) :: sigma, x, t
    INTEGER :: i

    sigma = 0.0_num
    IF (z < 1 .OR. z > 100) RETURN
    IF (.NOT. z_present(z)) RETURN
    x = LOG(energy / (gamma_min * m0c2)) &
        / LOG(gamma_max / gamma_min) * REAL(n_rate - 1, num)
    i = MIN(MAX(INT(x) + 1, 1), n_rate - 1)
    t = MIN(MAX(x - REAL(i - 1, num), 0.0_num), 1.0_num)
    sigma = (1.0_num - t) * rate_table(z,i) &
        + t * rate_table(z,i+1)

  END FUNCTION bh_brem_rate



  SUBROUTINE bh_brem_sample(energy, z, y, u, screen_gamma)

    REAL(num), INTENT(IN) :: energy
    INTEGER, INTENT(IN) :: z
    REAL(num), INTENT(OUT) :: y, u, screen_gamma
    REAL(num) :: gamma_i, y_min, y_max, xi, delta_min
    REAL(num) :: upper, shape

    gamma_i = energy / m0c2
    y_min = energy_cut / energy
    y_max = 1.0_num - gamma_final_min / gamma_i
    ! Only called after a positive optical-depth decrement.
    delta_min = y_min / (2.0_num * gamma_i * (1.0_num - y_min))
    ! Eq. 6.29 has a non-positive screening correction. Combined with
    ! Eq. 7.2 this gives a global rejection bound over y and u.
    upper = 2.0_num * (3.0_num + 2.0_num * MAX(0.0_num, &
        LOG(1.0_num / delta_min) - 2.0_num - bh_brem_coulomb(z)))
    DO
      y = y_min * EXP(random() * LOG(y_max / y_min))
      ! The measure is uniform in xi; the target density is D.
      xi = xi_min + random() * (1.0_num - xi_min)
      u = SQRT(1.0_num / xi - 1.0_num)
      screen_gamma = bh_brem_gamma(gamma_i, y, xi, z)
      shape = bh_brem_shape(y, u, xi, screen_gamma)
      IF (random() * upper <= shape) EXIT
    END DO

  END SUBROUTINE bh_brem_sample



  SUBROUTINE bh_brem_transfer(p_in, p_ph, y, u, screen_gamma, &
      spin_in, spin_out, xi_stored)

    REAL(num), INTENT(IN) :: p_in(3), p_ph(3), y, u, screen_gamma
    REAL(num), INTENT(IN) :: spin_in(3)
    REAL(num), INTENT(OUT) :: spin_out(3), xi_stored(3)
    REAL(num) :: n(3), pin_hat(3), a1(3), a2(3), r1(3), r2(3)
    REAL(num) :: xi, g, r, d, db, tb, lb, mb, eb, fb, pb
    REAL(num) :: sx, sy, sz, norm, xi_event(3)

    n = p_ph / SQRT(DOT_PRODUCT(p_ph, p_ph))
    pin_hat = p_in / SQRT(DOT_PRODUCT(p_in, p_in))
    CALL reference_axis(n, r1, r2)
    a1 = pin_hat - DOT_PRODUCT(pin_hat, n) * n
    norm = SQRT(DOT_PRODUCT(a1, a1))
    IF (norm > 1.0e-12_num) THEN
      a1 = a1 / norm
    ELSE
      a1 = r1
    END IF
    a2 = (/ n(2)*a1(3) - n(3)*a1(2), &
            n(3)*a1(1) - n(1)*a1(3), &
            n(1)*a1(2) - n(2)*a1(1) /)

    sx = DOT_PRODUCT(spin_in, a1)
    sy = DOT_PRODUCT(spin_in, a2)
    sz = DOT_PRODUCT(spin_in, n)
    xi = 1.0_num / (1.0_num + u*u)
    g = screen_gamma
    r = 1.0_num - y
    d = bh_brem_shape(y, u, xi, g)
    ! Eq. 5.7 traced over the other final particle. These are the
    ! angle-conditioned one-particle means, never unit-normalised.
    db = 8.0_num*r*u*u*xi*xi*g / d
    tb = 4.0_num*y*r*xi*u*(2.0_num*xi - 1.0_num)*g / d
    lb = y*((1.0_num + r)*(3.0_num + 2.0_num*g) &
        - 2.0_num*r*(1.0_num + 4.0_num*u*u*xi*xi*g)) / d
    mb = 4.0_num*r*(1.0_num + g - 2.0_num*u*u*xi*xi*g) / d
    eb = 4.0_num*y*xi*u*(2.0_num*xi - 1.0_num)*g / d
    fb = 4.0_num*r*y*xi*u*(1.0_num - 2.0_num*xi)*g / d
    pb = y*y*(1.0_num + 8.0_num*g*(xi - 0.5_num)**2) / d

    spin_out = (mb*sx + eb*sz)*a1 + mb*sy*a2 &
        + (fb*sx + (mb + pb)*sz)*n
    ! a1 is in the event plane; the unpolarised electron's photon
    ! is linearly polarised perpendicular to it (Q = -db).
    xi_event = (/ -db, tb*sx + lb*sz, 0.0_num /)
    xi_stored = pol_obs_to_stored(r1, r2, a1, xi_event)

  END SUBROUTINE bh_brem_transfer



  PURE FUNCTION bh_brem_shape(y, u, xi, g) RESULT(d)

    REAL(num), INTENT(IN) :: y, u, xi, g
    REAL(num) :: d, r

    r = 1.0_num - y
    d = (1.0_num + r*r)*(3.0_num + 2.0_num*g) &
        - 2.0_num*r*(1.0_num + 4.0_num*u*u*xi*xi*g)

  END FUNCTION bh_brem_shape



  PURE FUNCTION bh_brem_coulomb(z) RESULT(f)

    INTEGER, INTENT(IN) :: z
    REAL(num) :: f, a

    ! The Bethe-Maximon f(Z) approximation used in EPOCH's BH pair code.
    a = (alpha * REAL(z, num))**2
    f = a * (1.0_num / (1.0_num + a) + 0.20206_num &
        - 0.0369_num*a + 0.0083_num*a*a - 0.002_num*a**3)

  END FUNCTION bh_brem_coulomb



  PURE FUNCTION bh_brem_gamma(gamma_i, y, xi, z) RESULT(g)

    REAL(num), INTENT(IN) :: gamma_i, y, xi
    INTEGER, INTENT(IN) :: z
    REAL(num) :: g, delta, beta(3), big_b(3), screening
    INTEGER :: i, j

    delta = y / (2.0_num * gamma_i * (1.0_num - y))
    beta = REAL(z, num)**(1.0_num/3.0_num) * mol_b / 121.0_num
    big_b = (beta * xi / delta)**2
    screening = 0.0_num
    DO i = 1, 3
      screening = screening - 0.5_num * mol_alpha(i)**2 &
          * LOG(1.0_num + big_b(i))
      DO j = 1, 3
        IF (i == j) CYCLE
        screening = screening + mol_alpha(i)*mol_alpha(j) &
            * ((1.0_num + big_b(j)) &
            / (big_b(i) - big_b(j)) &
            * LOG(1.0_num + big_b(j)) + 0.5_num)
      END DO
    END DO
    g = LOG(1.0_num / delta) - 2.0_num &
        - bh_brem_coulomb(z) + screening

  END FUNCTION bh_brem_gamma



  FUNCTION bh_brem_integral(gamma_i, z) RESULT(sigma)

    REAL(num), INTENT(IN) :: gamma_i
    INTEGER, INTENT(IN) :: z
    REAL(num) :: sigma, y_min, y_max, log_span, log_xi_min
    REAL(num) :: y, xi, u, g, d, total
    INTEGER :: i, j

    sigma = 0.0_num
    y_min = energy_cut / (gamma_i * m0c2)
    y_max = 1.0_num - gamma_final_min / gamma_i
    IF (y_min >= y_max) RETURN
    log_span = LOG(y_max / y_min)
    ! Integrate d(xi) = xi d(log xi): the integrand is then smooth.
    log_xi_min = LOG(xi_min)
    total = 0.0_num
    DO i = 1, n_quad
      y = y_min * EXP(quad_x(i) * log_span)
      DO j = 1, n_quad
        xi = EXP((1.0_num - quad_x(j)) * log_xi_min)
        u = SQRT(1.0_num / xi - 1.0_num)
        g = bh_brem_gamma(gamma_i, y, xi, z)
        d = bh_brem_shape(y, u, xi, g)
        total = total + quad_w(i) * quad_w(j) * d * xi
      END DO
    END DO
    sigma = 2.0_num * REAL(z, num)**2 * alpha * re2 &
        * log_span * (-log_xi_min) * total

  END FUNCTION bh_brem_integral



  SUBROUTINE bh_brem_quadrature

    ! Gauss-Legendre nodes and weights mapped to [0,1].
    INTEGER :: i, j
    REAL(num) :: z, z_old, p0, p1, p2, dp

    DO i = 1, n_quad / 2
      z = COS(pi * (REAL(i, num) - 0.25_num) &
          / (REAL(n_quad, num) + 0.5_num))
      DO
        p0 = 1.0_num
        p1 = z
        DO j = 2, n_quad
          p2 = ((2.0_num*j - 1.0_num)*z*p1 &
              - (j - 1.0_num)*p0) / REAL(j, num)
          p0 = p1
          p1 = p2
        END DO
        dp = n_quad * (z*p1 - p0) / (z*z - 1.0_num)
        z_old = z
        z = z_old - p1 / dp
        IF (ABS(z - z_old) < 1.0e-15_num) EXIT
      END DO
      quad_x(i) = (1.0_num - z) / 2.0_num
      quad_x(n_quad+1-i) = (1.0_num + z) / 2.0_num
      quad_w(i) = 1.0_num / ((1.0_num - z*z)*dp*dp)
      quad_w(n_quad+1-i) = quad_w(i)
    END DO

  END SUBROUTINE bh_brem_quadrature
#endif

END MODULE brem_polarised
