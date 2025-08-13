program generate_pm_lightcurve
  use lunamod
  implicit none

  ! ---------- constants ----------
  real(8), parameter :: pi = 3.14159265358979323846d0
  real(8), parameter :: twopi = 2.0d0*pi
  real(8), parameter :: G = 6.67430d-11          ! m^3 kg^-1 s^-2
  real(8), parameter :: Rsun = 6.957e8          ! meters
  real(8), parameter :: day2s = 86400.d0        ! seconds per day

  ! ---------- user-set params (converted from Python) ----------
  ! Planet parameters
  real(8) :: p_in    = 0.124107       ! Rp/R*
  real(8) :: rho_star = 2720.77       ! stellar density [kg/m^3]? Check units
  real(8) :: bb      = 0.274943       ! impact parameter
  real(8) :: Pb      = 1071.23        ! orbital period [days]
  real(8) :: T_b     = 2455253.3      ! mid-transit BJD
  real(8) :: eb      = 0.0
  real(8) :: wb      = 0.0
  real(8) :: ab      ! a_B / R*

  ! Limb-darkening parameters (quadratic, converted from q1,q2)
  real(8) :: q1 = 0.428801, q2 = 0.483676
  real(8) :: u1, u2

  ! Moon parameters
  real(8) :: Ps     = 0.630702       ! moon orbital period [days]
  real(8) :: phi_s  = 3.447          ! moon phase in radians
  real(8) :: s_in   = 0.05417        ! Rmoon/R* = Rp/R* * Rmoon/Rplanet
  real(8) :: as     = 1.059          ! moon distance in stellar radii
  real(8) :: es     = 0.0
  real(8) :: ws     = 0.0
  real(8) :: is     = 0.0            ! inclination (radians)
  real(8) :: Os     = -0.02753       ! longitude of ascending node (radians)
  real(8) :: Msp_in = 0.00018904     ! moon/planet mass ratio

  ! Reference true anomaly for barycenter
  real(8) :: fplan

  ! Time grid
  integer, parameter :: nz = 2000
  real(8) :: tspan = 1.0d0
  real(8) :: t0, dt
  real(8), dimension(nz) :: t, ItotL
  integer :: i, ierr

  ! LUNA flags
  logical :: reverse
  integer :: show, animate, limb

  ! ---------- compute derived quantities ----------
  ! Convert Kipping q1,q2 to quadratic limb darkening
  u1 = 2.0d0*sqrt(q1)*q2
  u2 = sqrt(q1)*(1.0d0 - 2.0d0*q2)

  ! Compute a_B / R* from stellar density and orbital period
  ! LUNA expects a_B / R*, using Kepler's law: a/R* = ((G * rho_star)/(3*pi))^(1/3) * (Pb*day2s)^(2/3)
  ab = ((G * rho_star)/(3.0d0*pi))**(1.0d0/3.0d0) * (Pb*day2s)**(2.0d0/3.0d0)

  ! Set reference true anomaly
  fplan = (0.5d0*pi - wb)
  if (fplan < 0.d0) fplan = fplan + twopi
  if (fplan >= twopi) fplan = fplan - twopi

  ! ---------- build time array ----------
  dt = tspan / real(nz-1, 8)
  t0 = T_b - 0.5d0*tspan
  do i = 1, nz
     t(i) = t0 + dt*real(i-1,8)
  end do

  ! ---------- LUNA options ----------
  reverse = .false.
  show    = 0
  animate = 0
  limb    = 1

  ! ---------- call LUNA ----------
  call luna(t, Pb, T_b, p_in, ab, eb, wb, bb, &
            Ps, phi_s, s_in, as, es, ws, is, Os, Msp_in, &
            u1, u2, fplan, reverse, nz, show, animate, limb, ItotL)

  ! ---------- write output ----------
  open(unit=77, file='lightcurve_pm.txt', status='replace', iostat=ierr)
  if (ierr /= 0) then
     print *, 'Error opening lightcurve_pm.txt'
     stop
  end if
  do i = 1, nz
     write(77,'(F18.8,1X,ES24.16)') t(i), ItotL(i)
  end do
  close(77)

  print *, 'Wrote planet+moon light curve to lightcurve_pm.txt'

end program generate_pm_lightcurve
