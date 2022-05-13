
!==============
program test
!==============
use kinds
use constants
use bspline_module
use gauss_module
implicit none
  !
  integer, parameter :: stdout=6
  integer, parameter :: iun=10
  !
  integer :: i,j,ib
  integer :: nbasis
  integer :: bspline_order
  integer :: gauleg_order
  integer :: nx
  integer :: nxp
  real(dp), allocatable :: xgrid(:)
  real(dp), allocatable :: xpgrid(:)
  real(dp), allocatable :: bplot(:,:)
  !
  real(dp) :: xmin, xpmin
  real(dp) :: xmax, xpmax
  !
  type(bspline_type), allocatable :: bspline(:)
  !
  !character(64) :: subname="test_bspline"

  !
  ! params 
  !  
  bspline_order=3
  gauleg_order=8
  nbasis = 15
  !
  nx   = nbasis
  nxp  = 3000
  xmin =   0.0
  xmax =  10.0
  xpmin = -1.0
  xpmax = 13.0
  
  !
  ! workspace
  !
  allocate(xgrid(-bspline_order:nx+bspline_order))
  allocate(xpgrid(nxp))
  allocate(bplot(nxp,nbasis))
  allocate(bspline(nbasis))

  !
  ! grid init
  !
  do i = -bspline_order,nx+bspline_order
    xgrid(i) = xmin + sign(1,i-1)*(real(i-1,dbl)/real(nx,dbl))**2 *(xmax-xmin)
  enddo
  !
  do i = 1, nxp
    xpgrid(i) = xpmin + real(i-1,dp)/real(nxp-1,dp) * (xpmax-xpmin)
  enddo
  !
  write(stdout,"(a,2i4)") "bspline order:  ", bspline_order
  write(stdout,"(a,2i4)") "       nbasis:  ", nbasis
  write(stdout,"(a,2i4)") "        xgrid:  [ndim =] [size =]", nx, size(xgrid)
  write(stdout,"(a,2i4)") "        xgrid:  [lbound =] [ubound =]", ubound(xgrid,1), lbound(xgrid,1)
  write(stdout,"(6f12.6)") xgrid(:)
  !
  write(stdout,"(/,a,2i4)") "    plot grid: ", nxp
  write(stdout,"(6f12.6)") xpgrid(1:6)
  write(stdout,"(6f12.6)") xpgrid(nxp-5:nxp)

  !
  ! bspline init & plot
  !
  call gauss_init(max(gauleg_order,14))
  !
  do ib = 1, nbasis
    !
    call bspline_init(lbound(xgrid,1),ubound(xgrid,1),xgrid,ib,gauleg_order,gauleg_order, &
                      bspline_order, bspline(ib) )
    !
    call bspline_eval(nxp,xpgrid,bplot(:,ib),bspline(ib),0)
  enddo

  !
  ! dump data
  !
  open(iun,file="DATA_bspline.dat")
  do ib = 1, nbasis
    do i = 1, nxp 
      write(iun,"(2f15.9)") xpgrid(i), bplot(i,ib)
    enddo
    write(iun,*)
  enddo
  close(iun)
  !
  open(iun,file="DATA_bars.dat")
  do i = 1, nx
    write(iun,"(2f15.9)") xgrid(i),  0.0d0
    write(iun,"(2f15.9)") xgrid(i), 10.0d0
    write(iun,"()")
  enddo 
  close(iun)
  !
  open(iun,file="DATA_dots.dat")
  do ib = 1, nbasis
    do i = 1, bspline(ib)%nx
      write(iun,"(2f15.9)") bspline(ib)%x(i), bspline(ib)%y(i)
    enddo
    write(iun,*)
  enddo
  close(iun)
  !
  open(iun,file="DATA_gauleg.dat")
  do ib = 1, nbasis
    write(iun,"(/,'# bspline ',i5)") ib
    do j = 1, bspline(ib)%nseg
    do i = 1, bspline(ib)%ngs
      write(iun,"(2f15.9)") bspline(ib)%xgs(i,j), 0.0d0
      write(iun,"(2f15.9)") bspline(ib)%xgs(i,j), 1.0d0
      write(iun,*)
    enddo
    enddo
    write(iun,*)
  enddo
  close(iun)

  !
  ! cleanup
  !
  deallocate(xgrid)
  deallocate(xpgrid)
  deallocate(bplot)
  do ib = 1, nbasis
    call bspline_deallocate( bspline(ib) )
  enddo
  deallocate(bspline)
  !
end program test

