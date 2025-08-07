program test_mandel
  use mandelmod
  implicit none

  integer, parameter :: n = 100
  real(8) :: z(n), flux(n)
  real(8) :: p, c1, c2
  real(8) :: dz
  integer :: i
  integer :: ierr


  ! Your parameters:
  p = 0.124107d0       ! planet-to-star radius ratio
  c1 = 0.633746d0      ! limb darkening coefficient 1
  c2 = -0.030103d0     ! limb darkening coefficient 2

  ! Generate z values (impact parameter) around your b=0.274943
  ! Let's sample +/- 0.5 around b to see transit shape
  dz = 1.0d0 / n
  open(unit=10, file='lightcurve.txt', status='new', iostat=ierr)
  do i = 1, n
    z(i) = 0.274943d0 - 0.5d0 + dz * i
    call occultquad_wrapper(abs(z(i)), p, c1, c2, flux(i))
    print *, "z =", z(i), " flux =", flux(i)
    WRITE(10, *) 'This is a line of text.'

    write(10, *) z(i), flux(i)
  end do
  close(10)

  open(unit=10, file='lightcurve.txt', status='unknown', iostat=ierr)
  if (ierr /= 0) then
   print *, "Error opening lightcurve.txt, IOSTAT=", ierr
   stop
  end if
  do i = 1, n
    write(10, *) z(i), flux(i)
  end do
  close(10)

  print *, "Light curve written to lightcurve.txt"
end program test_mandel
