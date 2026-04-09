!======================================================================
! main.f90 — Modern Fortran skeleton for the Time-Dependent
!            Fokker-Planck solver (Hamilton, Lu & Petrosian 1990).
!
! Operator-splitting scheme per time step:
!   1. advect_space    — spatial transport along the loop  (tau)
!   2. energy_loss     — collisional + synchrotron energy losses
!   3. pitch_diffusion — Coulomb pitch-angle scattering   (Crank-Nicholson + TDMA)
!
! Reference: Hamilton, Lu & Petrosian, ApJ 354, 726 (1990)
!            Leach & Petrosian, ApJ 251, 781 (1981)
!======================================================================

module fp_parameters
  !! Compile-time grid dimensions (mirrors original fkrplk.h).
  implicit none
  integer, parameter :: nemax   = 120   ! max energy grid points
  integer, parameter :: nmumax  =  60   ! max pitch-angle grid points (each half)
  integer, parameter :: ntaumax = 120   ! max spatial grid points
  real,    parameter :: mc2     = 511.0 ! electron rest-mass energy [keV]
end module fp_parameters


module fp_grids
  !! Runtime grid arrays and physical constants.
  use fp_parameters
  implicit none

  ! Active grid sizes (read from input or set in init)
  integer :: ne, nmu, ntau

  ! --- Energy grid (quasi-logarithmic) ---
  real, allocatable :: e(:)          ! dimensionless kinetic energy E/mc2  [-1 : ne+1]
  real, allocatable :: de(:)         ! energy step widths                  [1  : ne]
  real, allocatable :: beta(:)       ! relativistic velocity v/c           [0  : ne]
  real, allocatable :: beta2(:)      ! beta at half-grid points            [0  : ne]
  real, allocatable :: rb3g2(:)      ! 1 / (beta^3 * gamma^2)             [0  : ne]

  ! --- Pitch-angle grid (linear in mu = cos theta) ---
  real, allocatable :: amu(:)        ! pitch-angle cosines                 [-nmu : nmu]
  real, allocatable :: dmu(:)        ! pitch-angle step widths             [-nmu+1 : nmu]
  real, allocatable :: onmu2(:)      ! 1 - mu^2                           [-nmu : nmu]
  real, allocatable :: onmu22(:)     ! 1 - 0.25*(mu_l + mu_{l-1})^2      [-nmu+1 : nmu]

  ! --- Spatial grid (column depth tau) ---
  real, allocatable :: dtau(:)       ! spatial step widths                 [1 : ntau]
  real, allocatable :: tau(:)        ! column depth values                 [0 : ntau]
  real, allocatable :: rlmbda(:)     ! normalised density rho/rho_max     [0 : ntau]
  real, allocatable :: rlmbda2(:)    ! density at half-grid points        [1 : ntau]
  real, allocatable :: bconv(:)      ! (a/(lnLam*rho_max)) * dlnB/ds     [0 : ntau]
  real, allocatable :: bm(:)         ! magnetic field strength B(s)       [0 : ntau]

  ! --- Distribution function ---
  ! phi(0:ne, -nmu:nmu, 0:ntau, 0:1)
  ! Last index: alternating time buffers (in / out)
  real, allocatable :: phi(:,:,:,:)

  ! --- Time stepping ---
  real :: dy       ! dimensionless time step
  real :: dt       ! physical time step [s]
  real :: tmax     ! maximum simulation time [s]

end module fp_grids


module fp_physics
  !! Physical model subroutines: density, magnetic field.
  implicit none
contains

  subroutine dands(t_col, rho, s)
    !! Map column depth t_col -> density rho [cm^-3] and position s [cm].
    !! Uniform-density model: rho = rho0 = const.
    real, intent(in)  :: t_col
    real, intent(out) :: rho, s
    real, parameter :: xlmbda = 20.0,  a_norm = 1.0e24, rho0 = 1.0e11
    real :: c_conv
    c_conv = xlmbda * rho0 / a_norm
    rho = rho0
    s   = t_col / c_conv
  end subroutine dands

  subroutine bfield(s, b, dbds)
    !! Magnetic field B(s) = B0 * [1 + (rm-1)*s^2/Lc^2]
    !! Returns B and d(lnB)/ds.
    real, intent(in)  :: s
    real, intent(out) :: b, dbds
    real, parameter :: b0 = 100.0, xc = 1.0e9, rm = 2.0
    b    = b0 * (1.0 + (rm - 1.0) * s**2 / xc**2)
    dbds = 2.0 * s * (rm - 1.0) / (xc**2 + (rm - 1.0) * s**2)
  end subroutine bfield

end module fp_physics


module fp_solver
  !! Core numerical routines: TDMA, operator-split steps.
  use fp_parameters
  use fp_grids
  use fp_physics
  implicit none
contains

  !--------------------------------------------------------------------
  ! Thomas Algorithm (TDMA) — tridiagonal matrix solver
  !
  ! Solves the system   A x = r   where A is tridiagonal:
  !
  !   | b(-n)  c(-n)     0   ...          0      |   | x(-n)   |   | r(-n)   |
  !   | a(-n+1) b(-n+1) c(-n+1) ...       0      |   | x(-n+1) | = | r(-n+1) |
  !   |   0     a(-n+2) b(-n+2) c(-n+2)   :      |   |  ...    |   |  ...    |
  !   |   :                        ...   c(n-1)   |   | x(n-1)  |   | r(n-1)  |
  !   |   0      ...    0     a(n)  b(n)          |   | x(n)    |   | r(n)    |
  !
  ! Inputs:   a, b, c  — sub-, main-, super-diagonal  [-n : n]
  !           rhs      — right-hand side vector        [-n : n]
  !           n        — half-bandwidth (system size = 2n+1)
  ! Output:   x        — solution vector               [-n : n]
  !
  ! Adapted from Numerical Recipes (Press et al., 1986, p. 40).
  !--------------------------------------------------------------------
  subroutine tdma(a, b, c, rhs, n, x)
    integer, intent(in) :: n
    real, intent(in)    :: a(-n:n), b(-n:n), c(-n:n), rhs(-n:n)
    real, intent(out)   :: x(-n:n)

    real :: gam(-nmumax:nmumax)
    real :: bet
    integer :: j

    ! --- Forward elimination ---
    bet  = b(-n)
    x(-n) = rhs(-n) / bet

    do j = -n + 1, n
       gam(j) = c(j - 1) / bet
       bet    = b(j) - a(j) * gam(j)
       x(j)   = (rhs(j) - a(j) * x(j - 1)) / bet
    end do

    ! --- Back substitution ---
    do j = n - 1, -n, -1
       x(j) = x(j) - gam(j + 1) * x(j + 1)
    end do

  end subroutine tdma


  !--------------------------------------------------------------------
  ! advect_space — spatial transport along the magnetic loop.
  !
  ! Uses Barton's monotonic scheme (upwind + anti-diffusion limiter).
  ! Separate sweeps for positive (mu>0) and negative (mu<0) velocities.
  !--------------------------------------------------------------------
  subroutine advect_space(in_buf, out_buf, symloop)
    integer, intent(in) :: in_buf, out_buf
    logical, intent(in) :: symloop

    real :: tphi(-nmumax:nmumax, 0:ntaumax)
    real :: tphi1(-nmumax:nmumax, 0:ntaumax)
    real :: a0, a1, a_vel, a1pa, a1ma
    real :: phic, phin, phip, dphi_val, dphip, delt, delto
    integer :: k, l, m

    a0 = dy / dtau(1)

    do k = 0, ne
       a1 = a0 * beta(k)

       ! Copy current time layer
       do m = 0, ntau
          do l = -nmu, nmu
             tphi(l, m)  = phi(k, l, m, in_buf)
             tphi1(l, m) = tphi(l, m)
          end do
       end do

       ! --- Negative velocities (mu < 0): sweep from ntau toward 0 ---
       do l = -nmu, -1
          a_vel = a1 * amu(l)          ! negative
          a1pa  = a_vel * (1.0 + a_vel)

          ! Boundary at m = ntau: phi = 0 beyond ntau
          phic = tphi(l, ntau)
          phin = tphi(l, ntau - 1)
          dphi_val = phic - phin
          delt = -dphi_val * phic
          if (delt > 0.0) then
             delt = -delt / phin
             tphi1(l, ntau) = phic * (1.0 + a_vel) - a1pa * delt
          else
             delt = 0.0
             tphi1(l, ntau) = phic * (1.0 + a_vel)
          end if

          ! Interior points
          do m = ntau - 1, 1, -1
             phip = phic
             phic = phin
             phin = tphi(l, m - 1)
             dphip = dphi_val
             dphi_val = phic - phin
             delto = delt      ! reuse name from original: deltp
             delt = dphi_val * dphip
             if (delt > 0.0) then
                delt = delt / (phip - phin)
             else
                delt = 0.0
             end if
             tphi1(l, m) = phic - a_vel * dphip + a1pa * (delto - delt)
          end do

          ! Boundary at m = 0
          phip = phic
          phic = phin
          dphip = dphi_val
          delto = delt
          if (symloop) then
             dphi_val = phic - tphi(-l, 1)
             delt = dphi_val * dphip
             if (delt > 0.0) then
                delt = delt / (phip - tphi(-l, 1))
             else
                delt = 0.0
             end if
          else
             delt = tphi(l, 0) * dphip
             if (delt > 0.0) then
                delt = delt / tphi(l, 1)
             else
                delt = 0.0
             end if
          end if
          tphi1(l, 0) = phic - a_vel * dphip + a1pa * (delto - delt)
       end do

       ! --- Positive velocities (mu > 0): sweep from 0 toward ntau ---
       do l = 1, nmu
          a_vel = a1 * amu(l)          ! positive
          a1ma  = a_vel * (1.0 - a_vel)

          ! Boundary at m = 0
          phic = tphi(l, 0)
          phin = tphi(l, 1)
          if (symloop) then
             dphi_val = phic - tphi(-l, 1)
             dphip    = phin - phic
             delt  = dphi_val * dphip
             if (delt > 0.0) then
                delt = delt / (phin - tphi(-l, 1))
             else
                delt = 0.0
             end if
             delto = dphi_val * (tphi(-l, 1) - tphi(-l, 2))
             if (delto > 0.0) then
                delto = delto / (phic - tphi(-l, 2))
             else
                delto = 0.0
             end if
             tphi1(l, 0) = phic - a_vel * dphi_val - a1ma * (delt - delto)
          else
             dphip = phin - phic
             delt  = phic * dphip
             if (delt > 0.0) then
                delt = delt / phin
             else
                delt = 0.0
             end if
             tphi1(l, 0) = phic * (1.0 - a_vel) - a1ma * delt
          end if

          ! Interior points
          do m = 1, ntau - 1
             phip = phic
             phic = phin
             phin = tphi(l, m + 1)
             delto = delt
             dphi_val = dphip
             dphip = phin - phic
             delt  = dphi_val * dphip
             if (delt > 0.0) then
                delt = delt / (phin - phip)
             else
                delt = 0.0
             end if
             tphi1(l, m) = phic - a_vel * dphi_val - a1ma * (delt - delto)
          end do

          ! Boundary at m = ntau: phi = 0 beyond ntau
          phip = phic
          phic = phin
          delto = delt
          dphi_val = dphip
          delt = -phic * dphi_val
          if (delt > 0.0) then
             delt = -delt / phip
          else
             delt = 0.0
          end if
          tphi1(l, ntau) = phic - a_vel * dphi_val - a1ma * (delt - delto)
       end do

       ! Write back
       do m = 0, ntau
          do l = -nmu, nmu
             phi(k, l, m, out_buf) = tphi1(l, m)
          end do
       end do
    end do

  end subroutine advect_space


  !--------------------------------------------------------------------
  ! energy_loss — collisional + radiative energy transport.
  !
  ! Barton's monotonic upwind scheme sweeping from high to low energy.
  !--------------------------------------------------------------------
  subroutine energy_loss(in_buf, out_buf)
    integer, intent(in) :: in_buf, out_buf

    real :: tphi(0:nemax), tphi1(0:nemax)
    real :: a_coeff, d1, d2, d3, dk, delt, f1, f1p
    integer :: k, l, m, kp1, km1

    do m = 0, ntau
       a_coeff = 2.0 * dy * rlmbda(m)

       do l = -nmu, nmu
          do k = 0, ne
             tphi(k) = phi(k, l, m, in_buf)
          end do

          ! Top of energy grid (k = ne)
          k   = ne
          km1 = k - 1
          d3  = tphi(k)
          d2  = 0.5 * (tphi(k) + tphi(km1))
          if (tphi(k) <= tphi(km1)) then
             dk = max(d3, d2)
          else
             dk = min(d3, d2)
          end if
          f1 = dk / beta2(k)
          tphi1(k) = d3 - a_coeff * f1 / (e(ne + 1) - e(km1))

          ! Interior points (k = ne-1 down to 1)
          do k = ne - 1, 1, -1
             kp1 = k + 1
             km1 = k - 1
             f1p = f1
             d3  = tphi(k)
             d2  = 0.5 * (d3 + tphi(km1))
             delt = 0.5 * (e(k) - e(km1)) / de(kp1)
             d1   = d3 * (1.0 + delt) - tphi(kp1) * delt
             if (d3 <= tphi(km1)) then
                dk = max(d3, min(d1, d2))
             else
                dk = min(d3, max(d1, d2))
             end if
             f1 = dk / beta2(k)
             tphi1(k) = d3 + a_coeff * (f1p - f1) / (e(kp1) - e(km1))
          end do

          ! Bottom of energy grid (k = 0)
          k   = 0
          kp1 = 1
          km1 = -1           ! ghost index handled by e(-1)
          f1p = f1
          d3  = tphi(0)
          delt = 0.5 * (e(0) - e(-1)) / de(1)
          d1   = d3 * (1.0 + delt) - tphi(1) * delt
          dk   = min(d3, d1)
          f1   = dk / beta2(0)
          tphi1(0) = d3 + a_coeff * (f1p - f1) / (e(1) - e(-1))

          do k = 0, ne
             phi(k, l, m, out_buf) = tphi1(k)
          end do
       end do
    end do

  end subroutine energy_loss


  !--------------------------------------------------------------------
  ! pitch_diffusion — Coulomb pitch-angle scattering.
  !
  ! Crank-Nicholson implicit scheme solved by the Thomas algorithm.
  ! The diffusion coefficient is D_mu = lambda * (1 - mu^2) / (beta^3 * gamma^2).
  !--------------------------------------------------------------------
  subroutine pitch_diffusion(in_buf, out_buf)
    integer, intent(in) :: in_buf, out_buf

    real :: a_sub(-nmumax:nmumax), b_diag(-nmumax:nmumax)
    real :: c_sup(-nmumax:nmumax), rhs(-nmumax:nmumax)
    real :: tphi(-nmumax:nmumax), sol(-nmumax:nmumax)
    real :: b1, b4, b5, aa, bb1, dydsum
    integer :: k, l, m, lp1, lm1

    do k = 0, ne
       b1 = beta(k) * 0.25

       do m = 0, ntau
          b4 = rlmbda(m) * rb3g2(k)

          do l = -nmu, nmu
             tphi(l) = phi(k, l, m, in_buf)
          end do

          ! Interior points: Crank-Nicholson discretisation
          do l = -nmu + 1, nmu - 1
             lp1 = l + 1
             lm1 = l - 1
             dydsum = dy / (dmu(lp1) + dmu(l))
             b5  = b4 * onmu2(l)
             aa  = b5  * dydsum
             bb1 = (-b4 * amu(l)) * dydsum

             ! Tridiagonal matrix coefficients (implicit side)
             a_sub(l) = bb1 - aa / dmu(l)
             b_diag(l) = 1.0 + b5 * dy / (dmu(l) * dmu(lp1))
             c_sup(l) = -bb1 - aa / dmu(lp1)

             ! Right-hand side (explicit side)
             rhs(l) = aa * ((tphi(lp1) - tphi(l)) / dmu(lp1) &
                          - (tphi(l) - tphi(lm1)) / dmu(l))  &
                    + bb1 * (tphi(lp1) - tphi(lm1))          &
                    + tphi(l)
          end do

          ! Boundary at mu = -1: forward difference
          c_sup(-nmu) = -b4 * dy / dmu(-nmu + 1)
          b_diag(-nmu) = 1.0 - c_sup(-nmu)
          rhs(-nmu) = -c_sup(-nmu) * (tphi(-nmu + 1) - tphi(-nmu)) &
                    + tphi(-nmu)

          ! Boundary at mu = +1: backward difference
          a_sub(nmu) = -b4 * dy / dmu(nmu)
          b_diag(nmu) = 1.0 - a_sub(nmu)
          rhs(nmu) = a_sub(nmu) * (tphi(nmu) - tphi(nmu - 1)) &
                   + tphi(nmu)

          ! Solve tridiagonal system
          call tdma(a_sub, b_diag, c_sup, rhs, nmu, sol)

          do l = -nmu, nmu
             phi(k, l, m, out_buf) = sol(l)
          end do
       end do
    end do

  end subroutine pitch_diffusion

end module fp_solver


!======================================================================
! Main program
!======================================================================
program fokker_planck
  use fp_parameters
  use fp_grids
  use fp_physics
  use fp_solver
  implicit none

  ! --- Local variables ---
  real    :: emin, emax, taumin, taumax
  real    :: rhomax, ytot, ttoy
  real    :: t, dtnew, dyold
  integer :: i, ii, jj, itr, ntr, k, l, m
  logical :: symloop

  ! Output reporting times (from LP88 configuration)
  integer, parameter :: ntr_max = 99
  real :: tr(0:ntr_max)

  ! ------------------------------------------------------------------
  ! 1. Physical parameters (matching original configuration)
  ! ------------------------------------------------------------------
  emin   = 10.0       ! minimum kinetic energy [keV]
  emax   = 10000.0    ! maximum kinetic energy [keV]
  ne     = 120        ! energy grid points
  nmu    = 60         ! pitch-angle grid points per hemisphere
  ntau   = 120        ! spatial grid points
  taumin = 1.0e-5     ! minimum column depth
  taumax = 2.0e-3     ! maximum column depth

  symloop = .true.    ! symmetric coronal loop

  ! Reporting times (uniform spacing to tmax ~ 10.5 s)
  ntr = ntr_max
  do i = 0, ntr
     tr(i) = real(i) * 10.51 / real(ntr)
  end do
  tmax = tr(ntr)

  ! ------------------------------------------------------------------
  ! 2. Allocate arrays
  ! ------------------------------------------------------------------
  call allocate_grids()

  ! ------------------------------------------------------------------
  ! 3. Initialise grids
  ! ------------------------------------------------------------------
  call init_energy_grid(emin, emax)
  call init_mu_grid()
  call init_tau_grid(taumax, rhomax)
  call compute_time_step()

  ! Physical time conversion: ytot = 1.667e12 / rhomax
  ytot = 1.6666667e12 / rhomax
  ttoy = 1.0 / ytot
  dt   = ytot * dy
  t    = 0.0

  ! Initialise distribution to zero
  phi = 0.0

  ! ------------------------------------------------------------------
  ! 4. Write header to output file
  ! ------------------------------------------------------------------
  open(unit=10, file='fkrplk.out', status='replace')
  write(10, '(1P6E13.5)') dt, tmax, dy, rhomax
  write(10, '(15I5)')      ne, nmu, ntau, ntr
  write(10, '(1P6E13.5)') (mc2 * e(k), k = 0, ne)
  write(10, '(1P6E13.5)') (amu(l), l = -nmu, nmu)
  write(10, '(1P6E13.5)') (dtau(m), m = 1, ntau)
  write(10, '(1P6E13.5)') (tau(m), m = 1, ntau)
  write(10, '(1P6E13.5)') (bm(m), m = 0, ntau)

  ! ------------------------------------------------------------------
  ! 5. Main time loop — Operator Splitting
  ! ------------------------------------------------------------------
  ii  = 1          ! "in"  buffer index
  jj  = 0          ! "out" buffer index
  itr = 0
  i   = 0

  write(*, '(A,ES10.3,A)') ' Starting integration, tmax = ', tmax, ' s'

  time_loop: do
     ! --- (a) Inject particles ---
     ! call inject(ii, jj, t)           ! <-- plug in your source term

     i = i + 1
     t = real(i) * dt

     ! --- (b) Check for output time ---
     if (t >= tr(itr)) then
        ! Adjust step to land exactly on the reporting time
        dtnew = tr(itr) - real(i - 1) * dt
        dyold = dy
        dy    = ttoy * dtnew

        ! Advance with reduced step
        call advect_space(jj, ii, symloop)
        call energy_loss(ii, jj)
        call pitch_diffusion(jj, ii)

        ! --- Write snapshot ---
        write(10, *) tr(itr), i
        call write_snapshot(10, ii)

        ! Advance remainder of original step
        dtnew = real(i) * dt - tr(itr)
        dy    = ttoy * dtnew

        call advect_space(ii, jj, symloop)
        call energy_loss(jj, ii)
        call pitch_diffusion(ii, jj)

        itr = itr + 1
        dy  = dyold

        write(*, '(A,I3,A,F8.3,A,I7,A)') &
             '  snapshot ', itr, '  t = ', tr(itr - 1), ' s  (step ', i, ')'
     else
        ! Standard full step
        call advect_space(jj, ii, symloop)
        call energy_loss(ii, jj)
        call pitch_diffusion(jj, ii)
     end if

     ! Swap buffers
     ii = 1 - ii
     jj = 1 - jj

     if (t >= tmax) exit time_loop
  end do time_loop

  close(10)
  write(*, '(A,I8,A)') ' Done — ', i, ' time steps completed.'

  call deallocate_grids()

contains

  ! ================================================================
  ! Allocation helpers
  ! ================================================================
  subroutine allocate_grids()
    allocate(e(-1:ne+1))
    allocate(de(1:ne))
    allocate(beta(0:ne))
    allocate(beta2(0:ne))
    allocate(rb3g2(0:ne))

    allocate(amu(-nmu:nmu))
    allocate(dmu(-nmu+1:nmu))
    allocate(onmu2(-nmu:nmu))
    allocate(onmu22(-nmu+1:nmu))

    allocate(dtau(1:ntau))
    allocate(tau(0:ntau))
    allocate(rlmbda(0:ntau))
    allocate(rlmbda2(1:ntau))
    allocate(bconv(0:ntau))
    allocate(bm(0:ntau))

    allocate(phi(0:ne, -nmu:nmu, 0:ntau, 0:1))
  end subroutine allocate_grids

  subroutine deallocate_grids()
    deallocate(e, de, beta, beta2, rb3g2)
    deallocate(amu, dmu, onmu2, onmu22)
    deallocate(dtau, tau, rlmbda, rlmbda2, bconv, bm)
    deallocate(phi)
  end subroutine deallocate_grids

  ! ================================================================
  ! Grid initialisation
  ! ================================================================

  subroutine init_energy_grid(emin_kev, emax_kev)
    !! Build the quasi-logarithmic energy grid.
    !!
    !! Energy steps increase as  de(k) = de0 + (E_{k-1} - E_0) * slope,
    !! where 'slope' is iterated until the grid exactly spans [emin, emax].
    !! When the uniform step is small enough, a simple linear grid is used.
    real, intent(in) :: emin_kev, emax_kev

    real :: g, des, deo, slope, r_fac, delt, s_sum
    integer :: k

    ! First energy point (dimensionless: E / mc2)
    e(0) = emin_kev / mc2
    deo  = e(0) / 4.0

    g       = e(0) + 1.0
    beta(0) = sqrt(1.0 - 1.0 / g**2)
    rb3g2(0) = 1.0 / (beta(0)**3 * g**2)

    e(-1) = e(0) - deo
    g     = 0.5 * (e(0) + e(-1)) + 1.0
    beta2(0) = sqrt(1.0 - 1.0 / g**2)

    des = (emax_kev - emin_kev) / (mc2 * real(ne))

    if (des <= deo) then
       ! --- Uniform spacing ---
       do k = 1, ne
          e(k)  = e(0) + real(k) * des
          de(k) = e(k) - e(k - 1)
          g     = e(k) + 1.0
          beta(k)  = sqrt(1.0 - 1.0 / g**2)
          rb3g2(k) = 1.0 / (beta(k)**3 * g**2)
          g        = 0.5 * (e(k) + e(k - 1)) + 1.0
          beta2(k) = sqrt(1.0 - 1.0 / g**2)
       end do
       e(ne + 1) = e(ne) + des
    else
       ! --- Quasi-logarithmic: adaptive slope ---
       slope = 1.0e-6
       r_fac = 10.0
       delt  = emax_kev / mc2 - e(0) - deo * real(ne)

       do while (mc2 * delt / emax_kev > 1.0e-4)
          slope = slope * r_fac
          s_sum = 1.0
          do k = 1, ne - 1
             s_sum = s_sum + (slope + 1.0)**k
          end do
          delt = emax_kev / mc2 - e(0) - deo * s_sum
          if (delt < 0.0) then
             slope = slope / r_fac
             r_fac = sqrt(r_fac)
             delt  = e(0)      ! reset to enter while again
          end if
       end do

       do k = 1, ne
          de(k) = deo + (e(k - 1) - e(0)) * slope
          e(k)  = e(k - 1) + de(k)
          g     = e(k) + 1.0
          beta(k)  = sqrt(1.0 - 1.0 / g**2)
          rb3g2(k) = 1.0 / (beta(k)**3 * g**2)
          g        = 0.5 * (e(k) + e(k - 1)) + 1.0
          beta2(k) = sqrt(1.0 - 1.0 / g**2)
       end do
       e(ne + 1) = e(ne) + de(ne)
    end if

    write(*, '(A,F8.2,A,F10.1,A)') &
         ' Energy grid: ', emin_kev, ' – ', emax_kev, ' keV'
    write(*, '(A,I4,A,ES10.3,A,ES10.3)') &
         '   ne = ', ne, '   de(1) = ', de(1) * mc2, ' keV   de(ne) = ', de(ne) * mc2

  end subroutine init_energy_grid


  subroutine init_mu_grid()
    !! Build the linear pitch-angle grid: mu = cos(theta) in [-1, +1].
    integer :: l

    amu(-nmu) = -1.0
    amu(nmu)  =  1.0
    onmu2(-nmu) = 0.0
    onmu2(nmu)  = 0.0

    do l = -nmu + 1, nmu - 1
       amu(l)   = real(l) / real(nmu)
       onmu2(l) = 1.0 - amu(l)**2
       onmu22(l) = 1.0 - 0.25 * (amu(l) + amu(l - 1))**2
       dmu(l)   = amu(l) - amu(l - 1)
    end do

    dmu(nmu)   = 1.0 - amu(nmu - 1)
    onmu22(nmu) = 1.0 - 0.25 * (amu(nmu) + amu(nmu - 1))**2
    amu(0)     = 0.0
    onmu2(0)   = 1.0

    write(*, '(A,I4,A,F6.4)') '   nmu = ', nmu, '   dmu = ', dmu(nmu)
  end subroutine init_mu_grid


  subroutine init_tau_grid(taumax_in, rhomax_out)
    !! Build the column-depth grid (linear spacing).
    real, intent(in)  :: taumax_in
    real, intent(out) :: rhomax_out

    real, parameter :: xlmbda = 20.0, a_norm = 1.0e24
    real :: taus, tn, to, rho, s, b, dbds, bfvar
    integer :: m

    taus = taumax_in / real(ntau)
    to   = 0.0
    tau(0) = 0.0

    call dands(taumax_in, rhomax_out, s)
    call dands(to, rho, s)

    bfvar     = a_norm / (xlmbda * rhomax_out)
    rlmbda(0) = rho / rhomax_out

    call bfield(s, b, dbds)   ! s = 0 at first call
    bconv(0) = bfvar * dbds
    bm(0)    = b

    do m = 1, ntau
       tn = taus * real(m)
       dtau(m) = tn - to

       call dands(tn, rho, s)
       tau(m)    = tn
       rlmbda(m) = rho / rhomax_out

       call bfield(s, b, dbds)
       bconv(m) = bfvar * dbds
       bm(m)    = b

       call dands(0.5 * (tn + to), rho, s)
       rlmbda2(m) = rho / rhomax_out

       to = tn
    end do

    write(*, '(A,I4,A,ES10.3)') '   ntau = ', ntau, '   taumax = ', taumax_in
  end subroutine init_tau_grid


  subroutine compute_time_step()
    !! Compute the maximum stable dimensionless time step dy.
    !! CFL-like conditions from collisions, advection, and mirroring.

    real :: dycoul, dytau, dymir, dymiro
    integer :: m, l

    ! Collision terms
    dycoul = 0.9 * min(0.5 * dmu(nmu) / rb3g2(0), beta(0) * de(1))

    ! Advective term
    dytau = dtau(1) / rlmbda2(1)
    do m = 1, ntau
       dytau = min(dytau, dtau(m) / rlmbda2(m))
    end do

    ! Mirroring term
    dymir = 1.0e20
    if (abs(bm(ntau) - bm(0)) > 1.0e-30) then
       do l = 0, nmu - 1
          dymir = min(dymir, dmu(l + 1) / onmu2(l))
       end do
       dymiro = dymir
       if (abs(bconv(ntau)) > 1.0e-30) dymir = dymiro / bconv(ntau)
       do m = 0, ntau - 1
          if (abs(bconv(m)) > 1.0e-30) dymir = min(dymir, dymiro / bconv(m))
       end do
       dymir = 0.5 * dymir / beta(ne)
    end if

    dy = min(dymir, dytau)
    dy = min(dy, dycoul)

    write(*, '(A,ES10.3)') ' Dimensionless time step dy = ', dy
  end subroutine compute_time_step


  subroutine write_snapshot(unit_out, buf)
    !! Write one time snapshot of phi to the output file.
    integer, intent(in) :: unit_out, buf
    integer :: m, l

    do m = 0, ntau
       do l = -nmu, nmu
          write(unit_out, '(1P6E13.5)') (phi(k, l, m, buf), k = 0, ne)
       end do
    end do
  end subroutine write_snapshot

end program fokker_planck
