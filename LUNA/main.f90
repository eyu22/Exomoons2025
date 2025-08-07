program run_luna
  use lunamod
  implicit none

  integer, parameter :: nz = 1000
  integer :: i

  real(8), dimension(nz) :: t       ! time array (e.g. days)
  real(8), dimension(nz) :: ItotL   ! output light curve

  ! Orbital and physical parameters (typical posterior values from Kipping et al.)
  real(8) :: Pb      ! Planet orbital period (days)
  real(8) :: T_b     ! Time of planet transit center (days)
  real(8) :: p_in    ! Planet radius ratio Rp/Rs
  real(8) :: ab      ! Semi-major axis scaled by stellar radius a/Rs (planet)
  real(8) :: eb      ! Planet eccentricity
  real(8) :: wb      ! Planet argument of periastron (radians)
  real(8) :: bb      ! Impact parameter of planet

  real(8) :: Ps      ! Moon orbital period (days)
  real(8) :: phi_s   ! Phase offset for moon orbit (radians)
  real(8) :: s_in    ! Moon radius ratio Rs/Rs (moon/star radius ratio)
  real(8) :: as      ! Moon semi-major axis scaled by stellar radius a/Rs (moon)
  real(8) :: es      ! Moon eccentricity
  real(8) :: ws      ! Moon argument of periastron (radians)
  real(8) :: is      ! Moon inclination (radians)
  real(8) :: Os      ! Moon longitude of ascending node (radians)
  real(8) :: Msp_in  ! Mass ratio moon to planet

  real(8) :: u1, u2  ! Quadratic limb darkening coefficients
  real(8) :: fplan   ! Some flag (usually 0 or 1)
  integer :: show, animate, limb
  logical :: reverse

  ! Assign typical posterior values (example values, replace with your actual posteriors)
  Pb = 10.0d0                ! planet orbital period in days
  T_b = 0.0d0                ! planet mid-transit time
  p_in = 0.1d0               ! planet/star radius ratio ~0.1
  ab = 15.0d0                ! a/Rs for planet
  eb = 0.0d0                 ! planet eccentricity
  wb = 0.0d0                 ! argument of periastron in radians
  bb = 0.0d0                 ! impact parameter

  Ps = 1.5d0                 ! moon orbital period (days)
  phi_s = 0.0d0              ! moon phase offset
  s_in = 0.03d0              ! moon/star radius ratio ~0.03
  as = 5.0d0                 ! moon semi-major axis scaled by stellar radius
  es = 0.0d0                 ! moon eccentricity
  ws = 0.0d0                 ! moon argument of periastron
  is = 0.0d0                 ! moon inclination
  Os = 0.0d0                 ! moon longitude of ascending node
  Msp_in = 0.01d0            ! moon-to-planet mass ratio

  u1 = 0.3d0                 ! limb darkening quadratic coeff 1
  u2 = 0.2d0                 ! limb darkening quadratic coeff 2
  fplan = 0.0d0              ! flag related to planet rotation (0 or 1)
  show = 0
  animate = 0
  limb = 1                   ! use limb darkening
  reverse = .false.

  ! Generate time array centered around transit (e.g., from -0.1 to +0.1 days)
  do i = 1, nz
    t(i) = T_b - 0.1d0 + 0.2d0 * (i-1) / (nz-1)
  end do

  ! Call the LUNA subroutine
  call luna(t, Pb, T_b, p_in, ab, eb, wb, bb, Ps, phi_s, s_in, as, es, ws, is, Os, Msp_in, &
            u1, u2, fplan, reverse, nz, show, animate, limb, ItotL)

  ! Output the light curve (for example, print first 10 points)
  print *, 'Time (days)   Flux (normalized)'
  do i = 1, 10
    print '(F10.6, 2X, F10.8)', t(i), ItotL(i)
  end do

end program run_luna
