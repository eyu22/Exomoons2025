module mandelmod
  implicit none
  private
  public :: occultquad_wrapper

contains

  subroutine occultquad_wrapper(z, p, c1, c2, mu)
    implicit none
    real(8), intent(in)  :: z, p, c1, c2
    real(8), intent(out) :: mu

    call occultquad(z, p, c1, c2, mu)
  end subroutine occultquad_wrapper

end module mandelmod
