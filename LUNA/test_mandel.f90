program test_mandel
  use mandelmod
  implicit none

  integer, parameter :: n = 1000
  real(8) :: z(n), flux(n)
  real(8) :: p, c1, c2
  real(8) :: dz
  integer :: i
  open(unit=10, file='lightcurve.txt', status='replace')

  ! ---- Transit parameters ----
  p = 0.1d0        ! Planet-to-star radius ratio
  c1 = 0.3d0       ! Limb darkening coeff 1
  c2 = 0.2d0       ! Limb darkening coeff 2

  ! ---- Generate z values (impact parameter) ----
  dz = 2.0d0 / n
  do i = 1, n
     z(i) = -1.0d0 + dz * i
     call occultquad_wrapper(abs(z(i)), p, c1, c2, flux(i))
     write(10, *) z(i), flux(i)
  end do

  close(10)
  print *, "Light curve written to lightcurve.txt"
end program test_mandel
