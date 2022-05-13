
!======================
module fitting_m
  !======================
  !
  use fitting_base_m
  use fitting_genetic_m
  implicit none
  private

  interface
     function pole_fitting_pywrap(nx,xdata,ydata,npoles_fit,params) result(ierr)
        use kinds, only : dp
        integer    :: nx,npoles_fit
        real(dp)   :: xdata(nx)
        real(dp)   :: ydata(nx)
        real(dp)   :: params(*)
        integer    :: ierr
     end function
  end interface

  public :: fitting_verbosity

  public :: dlsei_wrap 
  public :: nnls_wrap 
  public :: pole_fitting_pywrap
  public :: pole_fitting_LSQ
  public :: pole_fitting_LSQ_cmplx
  public :: pole_fitting_LSQ_cmplx_lsd
  public :: pole_fitting_LSQ_NL
  public :: pole_fitting_LSQ_NL_cmplx
  public :: pole_fitting_LSQ_NL_cmplx_c
  public :: pole_fitting_genetic_cmplx
  public :: pole_fitting_pade_cmplx
  !
  public :: pole_interp_pade_cmplx
  public :: pole_interp_pade_cmplx_dn
  public :: pole_interp_pade_cmplx_noLS
  public :: pole_rational_interp_cmplx
  !
  public :: rational_interp_A
  public :: rational_interp_B

end module fitting_m

