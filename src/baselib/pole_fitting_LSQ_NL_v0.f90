
!======================
subroutine pole_fitting_LSQ_NL(nx,xdata,ydata,npoles_fit,params,thr,niterx,ierr)
  !======================
  use constants
  use util_module,  only: mat_mul,mat_sv,mat_rank,mat_lsd
  !
  ! performs least-square fitting on input data.
  ! Fitting functions are poles (complex)
  ! ierr /=0 indicates a numerical problem
  !
  integer,    intent(in)    :: nx,npoles_fit
  real(dp),   intent(in)    :: xdata(nx)
  complex(dp),intent(in)    :: ydata(nx,1)
  real(dp),   intent(inout) :: params(*)
  real(dp),   intent(in)    :: thr
  integer(dp),intent(in)    :: niterx
  integer,    intent(out)   :: ierr
  !
  integer   :: ndim
  integer   :: i,j,ind,iter,niter
  logical   :: lconv
  real(dp)  :: res,res0,delta
  complex(dp) :: cpole
  complex(dp), allocatable :: ydata_fit(:,:)
  complex(dp), allocatable :: Xmat(:,:)
  complex(dp), allocatable :: Amat(:,:)
  complex(dp), allocatable :: Bmat(:,:)
  
  delta=0.01
  ierr=0
  ndim=1+2*npoles_fit
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(Xmat(nx,ndim))
  allocate(ydata_fit(nx,1))

  res0=0.0

  iter_loop:&
  do iter = 1, niterx

     !
     ! the first function is constant
     ! 
     Xmat(:,1)=1.0
     !
     do j = 2, ndim, 2
       !
       ind = 2+ (j/2-1)*3
       cpole = cmplx(params(ind),params(ind+1),dp)
       !
       do i = 1, nx
         Xmat(i,j+0) = 1.0d0/(xdata(i)-cpole)**2
         Xmat(i,j+1) = 1.0d0/(xdata(i)-cpole)
       enddo
       !
     enddo
     !
     call mat_mul(Bmat,Xmat,'C',ydata,'N',ndim,1,nx)
     call mat_mul(Amat,Xmat,'C',Xmat,'N', ndim,ndim,nx)
     !Amat=Amat/real(nx,dp)
     !Bmat=Bmat/real(nx,dp)

! XXX
WRITE(0,*) "rank A", mat_rank(ndim,ndim,Amat,1.0d-6)

     !
     ! solution
     !
     !call mat_sv(ndim,1,Amat,Bmat,ierr)
     call mat_lsd(ndim,ndim,1,Amat,Bmat,1.0d-8,IERR=ierr)
     if (ierr/=0) then
        call cleanup_loc()
        return
     endif

     !
     ! eval residuals
     !
     call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
     res = sqrt(real(dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1)),dp))/real(nx,dp)
!     res = sqrt(real(dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1)),dp))/real(ndim,dp)

     !
     ! reshape data
     ! assuming real data 
     ! (extension to complex quantities to be implemented)
     !
     params(1)=real(Bmat(1,1),dp)
     !
     do j = 2, ndim, 2
       ! 
       ind = 2+ (j/2-1)*3
       params(ind+2)=real(Bmat(j+1,1),dp)
       if (abs(real(Bmat(j+1,1),dp))< 1.0d-5) cycle
       !
       cpole=Bmat(j,1)/real(Bmat(j+1,1),dp)
       !
       if (real(cpole) <-delta) cpole=cmplx(-delta,aimag(cpole),dp)
       if (real(cpole) >delta)  cpole=cmplx(delta,aimag(cpole),dp)
       if (aimag(cpole)<-delta) cpole=cmplx(real(cpole,dp),-delta,dp)
       if (aimag(cpole)>delta)  cpole=cmplx(real(cpole,dp),delta,dp)
       !
       params(ind+0)=params(ind+0) +real(cpole)
       params(ind+1)=params(ind+1) +aimag(cpole)
       !
     enddo
     params(2+3*npoles_fit) = res
     write(0,"(a,i5, f21.9)") "iter res", iter, res

     !
     ! check convergence
     !
     if (iter>1 .and. abs(res-res0) < thr) then
        lconv=.true.
        exit iter_loop
     endif
     !
     res0=res
     !
  enddo iter_loop

  if (.not.lconv) ierr=1

  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     if (allocated(Amat)) deallocate(Amat)
     if (allocated(Bmat)) deallocate(Bmat)
     if (allocated(Xmat)) deallocate(Xmat)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
  end subroutine
  !
end subroutine pole_fitting_LSQ_NL
