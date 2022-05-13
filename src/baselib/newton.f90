
!===============
module newton_module
  !===============
  !
  use kinds
  use constants
  implicit none
  !
  abstract interface
     !
     subroutine newton_func(ndim,x,F)
        use constants 
        integer,  intent(in) :: ndim
        real(dp), intent(in) :: x(ndim)
        real(dp), intent(out):: F(ndim)
     end subroutine newton_func
     !
     subroutine newton_jacob(ndim,x,J)
        use constants 
        integer,  intent(in) :: ndim
        real(dp), intent(in) :: x(ndim)
        real(dp), intent(out):: J(ndim,ndim)
     end subroutine newton_jacob
     !
  end interface

contains

!==========================
subroutine newton_solver(ndim,func,jacob,brhs,x0,alpha,niterx,thr,ierr,niter)
  !==========================
  !
  ! on exit:   
  !    ierr = 0    smooth exit
  !    ierr < 0    max number of iteration reached
  !    ierr > 0    problems inverting J
  !
  use kinds
  use constants 
  use util_module, only: mat_inv,mat_mul
  implicit none
  !
  integer,  intent(in)   :: ndim, niterx
  real(dp), intent(in)   :: brhs(ndim)
  real(dp), intent(inout):: x0(ndim)
  real(dp), intent(in)   :: alpha
  real(dp), intent(in)   :: thr
  procedure(newton_func) :: func
  procedure(newton_jacob):: jacob
  integer,  intent(out)  :: ierr
  integer, optional, intent(out) :: niter
  !
  integer  :: iter,niter_
  real(dp) :: F(ndim,1),J(ndim,ndim)
  real(dp) :: Jinv(ndim,ndim)
  real(dp) :: dx(ndim,1)
  real(dp) :: x_new(ndim), res
  
  ierr=0
  if (present(niter)) niter=0
  !
  iteration_loop:&
  do iter = 1, niterx
    !
    if (present(niter)) niter=iter
    !
    call func(ndim,x0,F(:,1))
    call jacob(ndim,x0,J)
    F(:,1)=F(:,1)-brhs(:)
    !
    call mat_inv(ndim,J,Jinv,ierr=ierr)
    if (ierr/=0) then
      ierr=abs(ierr)
      return
    endif
    !
    call mat_mul(dx,Jinv,'N',F,'N',ndim,1,ndim)
    x_new = x0 -dx(:,1)
    !
    res=dot_product(x_new-x0,x_new-x0)
    !
    if (res < thr) then
      exit iteration_loop
    else if (iter==niterx) then
      ierr=-1
    endif
    !
    ! mixing
    !
    x0 = alpha*x_new + (1.0d0-alpha)*x0
    !
  enddo iteration_loop
  !
  x0 = x_new
  !
end subroutine newton_solver
  
end module newton_module

