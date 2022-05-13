
!======================
subroutine fitting_get_peaks(nx,xdata,ydata,npeaks,ipeaks)
  !======================
  use constants
  !
  ! returns the indexes of the largest npeaks
  ! (from input) of the give dataset
  ! Assumes ydata in input to be positive.
  !
  integer, intent(in)  :: nx, npeaks
  real(dp),intent(in)  :: xdata(nx),ydata(nx)
  integer, intent(out) :: ipeaks(npeaks)
  !
  real(dp), allocatable :: dy(:)
  integer,  allocatable :: ipeaks_tmp(:)
  integer :: i,ii,j,ind,npeaks_,imin
  real(dp):: val,vmin

  if (npeaks==0) return
  !
  ipeaks(1:npeaks)=0
  if (any(ydata<0.0)) return
  !
  allocate(dy(nx-1))
  allocate(ipeaks_tmp(nx-1))
  !
  ipeaks_tmp(:)=0
  ind=0
  !
  do i = 1, nx-1
    dy(i)=(ydata(i+1)-ydata(i))
  enddo
  !
  do i = 2,nx-1
    !
    ! find maxima
    !
    if ( dy(i)*dy(i-1)<0.0 .and. &
         ydata(i)>ydata(i-1) .and. ydata(i)>ydata(i+1) ) then
      ind=ind+1
      ipeaks_tmp(ind)=i
    endif
  enddo
  !
  npeaks_=0
  !
  do j = 1, ind
    i = ipeaks_tmp(j)
    val = ydata(i)
    !
    if ( npeaks_ < npeaks ) then
       npeaks_=npeaks_+1
       ipeaks(npeaks_) = i
    else
       imin=1
       vmin=ydata(ipeaks(imin))
       do ii = 2, npeaks
         if (ydata(ipeaks(ii))<vmin) then
           imin=ii
           vmin=ydata(ipeaks(ii))
         endif
       enddo
       !
       if (val > vmin ) then
          ipeaks(imin)=i
       endif
       !
    endif
  enddo

  deallocate(dy,ipeaks_tmp)
  return
  !
end subroutine fitting_get_peaks
