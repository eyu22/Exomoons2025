module mandelmod
  implicit none
contains

  subroutine occultquad_wrapper(z, p, c1, c2, mu)
    implicit none
    real(8), intent(in) :: z, p, c1, c2
    real(8), intent(out) :: mu

    integer :: nz
    real(8), allocatable :: z0(:), muo1(:), mu0(:)

    nz = 1
    allocate(z0(nz), muo1(nz), mu0(nz))

    z0(1) = z
    call occultquad(z0, c1, c2, p, muo1, mu0, nz)

    mu = muo1(1)  ! or mu0(1) depending on which output you want

    deallocate(z0, muo1, mu0)
  end subroutine occultquad_wrapper

end module mandelmod
