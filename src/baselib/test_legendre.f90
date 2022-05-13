
!==============
program test
!==============
use kinds
use constants
use gauss_module
implicit none
  !
  integer, parameter :: stdout=6
  !
  integer :: ierr, narg
  integer :: pol_order
  real(dp), allocatable :: pol_coeff(:)
  real(dp), allocatable :: pol_roots(:)
  character(64) :: str
  !character(64) :: subname="test_legendre"

  !
  ! input
  !
  narg = command_argument_count()
  if (narg<1) then
    write(stdout,"(' USAGE: test_legendre.x  <order>')")
    write(stdout,"('     <order> = order of the Legendre polynomial to compute ')")
    stop
  endif
  !
  call get_command_argument(1,str) 
  read(str,*,iostat=ierr) pol_order
  if (ierr/=0) stop "reading <order>"

  allocate(pol_coeff(0:pol_order))
  allocate(pol_roots(pol_order))
  !
  call legendre_pol_coeff(pol_coeff,pol_order)

  write(stdout,"(a,i5)") "Legendre coeff: order = ", pol_order
  call write_poly(stdout, pol_coeff, pol_order)
  !
  call legendre_pol_roots(pol_order, pol_roots)
  write(stdout,"(/,a,i5)") "Roots of Pn(x): "
  write(stdout,"(3x,f25.16)") pol_roots(:)
  ! 
  deallocate(pol_coeff) 
  deallocate(pol_roots) 
  !
end program test


subroutine write_poly(iun,coeff,n)
use kinds
implicit none
  !
  integer   :: iun, n
  real(dbl) :: coeff(0:n)
  integer :: i
  character(10) :: str
  !
  if (n<0) return
  write(iun,"(4x,'P =',f12.6,'     +')") coeff(0)
  !
  do i = 1,n
     write(str,"(i4)") i
     str=trim(adjustl(str))
     if (i<n) str=trim(str)//" +"
     write(iun,"(7x,f12.6,' x^',a)") coeff(i),trim(str)
  enddo
  !
end subroutine
