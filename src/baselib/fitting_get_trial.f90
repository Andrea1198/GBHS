
!======================
subroutine fitting_get_trial(nx,xdata,ydata,grid_range,npoles,ipoles,ierr)
  !======================
  use constants
  implicit none
  !
  ! returns the indexes npoles points, 
  ! distributed as such to be denser where
  ! the input function (ydata) has peaks
  !
  ! ydata in input is not assumed to be positive
  !
  ! points are equally spaced along the curve y=f(x),
  ! L = int_a^b  sqrt( 1 + f'^2 ) dx
  !
  ! Actually, we can also use a different power, like 1.0 instead of 0.5
  ! ie. using: (1+f'^2)^power
  !
  ! on output: ierr=0  all good,  
  !            ierr/=0 problems found
  !
  integer, intent(in)  :: nx, npoles
  real(dp),intent(in)  :: xdata(nx),ydata(nx)
  integer, intent(in)  :: grid_range(2)
  integer, intent(out) :: ipoles(npoles)
  integer, intent(out) :: ierr
  !
  real(dp), allocatable :: dy(:)
  real(dp), allocatable :: lgrid(:)
  real(dp):: power
  integer :: i,j
  logical :: lfound
  real(dp):: dx,llen,lentot

  ierr=0
  if (npoles==0) return
  !
  power=1.0
  !
  ipoles(1:npoles)=0
  !
  allocate(dy(nx))
  allocate(lgrid(nx))
  !
  ! define the curve derivative (central differences)
  !
  do i = grid_range(1)+1,grid_range(2)-1
    dy(i)=(ydata(i+1)-ydata(i-1))/2.0_DP
  enddo
  dy(grid_range(1))=dy(grid_range(1)+1)
  dy(grid_range(2))=dy(grid_range(2)-1)
  !
  ! compute the curve partial lengths
  !
  lgrid(1:grid_range(1))=0.0
  do i = grid_range(1)+1,grid_range(2)
    !
    dx=xdata(i)-xdata(i-1)
    lgrid(i)=lgrid(i-1)+(1.0_DP + dy(i)**2)**power *dx
    !
  enddo
  ! not really needed
  lgrid(grid_range(2)+1:nx)=lgrid(grid_range(2))
  !
  lentot=lgrid(grid_range(2))
  llen=lentot/real(npoles+1,DP)
  !
  ipoles(1)=grid_range(1)-1
  do j = 1, npoles
    !
    ! get the point closes to j*llen
    lfound=.false.
    !
    search_loop:&
    do i = ipoles(j)+1, grid_range(2)-1
      if ( lgrid(i)-j*llen >= 0.0 ) then
        lfound=.true.
        ipoles(j)=i
        exit search_loop
      endif
    enddo search_loop
    !
    ! if needed, return with a problem
    if (.not.lfound) then
      ierr=j
      return
    endif
    !
  enddo

  deallocate(dy)
  deallocate(lgrid)
  return
  !
end subroutine fitting_get_trial
