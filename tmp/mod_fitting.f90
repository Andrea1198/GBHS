
!======================
module fitting_m
  !======================
  use kinds
  use constants
  use util_module,  only: mat_mul,mat_sv,mat_lsd
  use poly_m
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

  public :: pole_fitting_pywrap
  public :: pole_fitting_LSQ
  public :: pole_fitting_LSQ_cmplx
  public :: pole_fitting_LSQ_NL
  public :: pole_fitting_LSQ_NL_cmplx
!  public :: pole_fitting_pade
  public :: pole_fitting_pade_cmplx

contains

!======================
  function pole_fitting_LSQ_cmplx(nx,xdata,ydata_,npoles_fit,params,allow_offset) result(ierr)
  !======================
  implicit none
  !
  ! performs least-square fitting on input data.
  ! Fitting functions are poles (complex)
  ! ierr /=0 indicates a numerical problem
  !
  integer     :: nx,npoles_fit
  real(dp)    :: xdata(nx)
  complex(dp) :: ydata_(nx)
  complex(dp) :: params(*)
  integer     :: ierr
  logical, optional :: allow_offset
  !
  integer   :: ndim
  integer   :: i,j,ind
  logical  :: allow_offset_
  complex(dp) :: cpole,res
  complex(dp), allocatable :: ydata_fit(:,:)
  complex(dp), allocatable :: ydata(:,:)
  complex(dp), allocatable :: Xmat(:,:)
  complex(dp), allocatable :: Amat(:,:)
  complex(dp), allocatable :: Bmat(:,:)

  ierr=0
  ndim=npoles_fit+1
  !
  allow_offset_=.true.
  if (present(allow_offset)) allow_offset_=allow_offset
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(Xmat(nx,ndim))
  allocate(ydata_fit(nx,1))
  allocate(ydata(nx,1))
  ydata(:,1)=ydata_(:)

  if (allow_offset_) then
    Xmat(:,1)=1.0
  else
    Xmat(:,1)=0.0
    Xmat(nx,1)=1.0
  endif
  !
  !
  do j = 2, ndim
    !
    ind = 2+(j-2)*3
    cpole = cmplx(real(params(ind),dp),real(params(ind+1),dp),dp)
    !
    do i = 1, nx
      Xmat(i,j) = params(ind+2)/(xdata(i)-cpole)
    enddo
    !
  enddo
  !
  call mat_mul(Bmat,Xmat,'C',ydata,'N',ndim,1,nx)
  call mat_mul(Amat,Xmat,'C',Xmat,'N', ndim,ndim,nx)

  !
  ! solution
  !
  call mat_sv(ndim,1,Amat,Bmat,ierr)
  if (ierr/=0) then
     call cleanup_loc()
     return
  endif

  !
  ! eval residuals
  !
  call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
  res = dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1))/real(nx,dp)

  !
  ! reshape data
  ! assuming real data
  ! (extension to complex quantities to be implemented)
  !
  params(1)=Bmat(1,1)
  do j = 2, ndim
    ind = 2+(j-2)*3
    params(ind+2)=Bmat(j,1)*params(ind+2)
  enddo
  params(2+3*npoles_fit) = res

  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     if (allocated(Amat))  deallocate(Amat)
     if (allocated(Bmat))  deallocate(Bmat)
     if (allocated(Xmat))  deallocate(Xmat)
     if (allocated(ydata)) deallocate(ydata)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
  end subroutine
  !
end function pole_fitting_LSQ_cmplx

!======================
  function pole_fitting_LSQ(nx,xdata,ydata_,npoles_fit,params,allow_offset) result(ierr)
  !======================
  implicit none
  !
  ! performs least-square fitting on input data.
  ! Fitting functions are real lorentzians
  ! ierr /=0 indicates a numerical problem
  !
  integer  :: nx,npoles_fit
  real(dp) :: xdata(nx)
  real(dp) :: ydata_(nx)
  real(dp) :: params(*)
  logical, optional :: allow_offset
  integer  :: ierr
  !
  integer  :: ndim
  integer  :: i,j,ind
  logical  :: allow_offset_
  real(dp) :: res,dd,xx
  real(dp), allocatable :: ydata_fit(:,:)
  real(dp), allocatable :: ydata(:,:)
  real(dp), allocatable :: Xmat(:,:)
  real(dp), allocatable :: Amat(:,:)
  real(dp), allocatable :: Bmat(:,:)

  ierr=0
  ndim=npoles_fit+1
  !
  allow_offset_=.true.
  if (present(allow_offset)) allow_offset_=allow_offset
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(Xmat(nx,ndim))
  allocate(ydata_fit(nx,1))
  allocate(ydata(nx,1))
  ydata(:,1)=ydata_(:)
  !
  !
  ! the first function is constant
  if (allow_offset_) then
    Xmat(:,1)=1.0
  else
    Xmat(:,1)=0.0
    Xmat(nx,1)=1.0
  endif
  !
  do j = 2, ndim
    !
    ind = 2+(j-2)*3
    xx  = params(ind)
    dd  = params(ind+1)
    !
    do i = 1, nx
      Xmat(i,j)    = dd/((xdata(i)-xx)**2 + dd**2)
    enddo
    !
  enddo
  !
  call mat_mul(Bmat,Xmat,'T',ydata,'N',ndim,1,nx)
  call mat_mul(Amat,Xmat,'T',Xmat,'N', ndim,ndim,nx)

  !
  ! solution
  !
  call mat_sv(ndim,1,Amat,Bmat,ierr)
  if (ierr/=0) then
     call cleanup_loc()
     return
  endif

  !
  ! eval residuals
  !
  call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
  res = dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1))/real(nx,dp)

  !
  ! reshape data
  ! assuming real data
  ! (extension to complex quantities to be implemented)
  !
  params(1)=0.0d0
  if (allow_offset_) params(1)=Bmat(1,1)
  !
  do j = 2, ndim
    ind = 2+(j-2)*3
    params(ind+2)=Bmat(j,1)
  enddo
  params(2+3*npoles_fit) = res
  !

  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     if (allocated(Amat))  deallocate(Amat)
     if (allocated(Bmat))  deallocate(Bmat)
     if (allocated(Xmat))  deallocate(Xmat)
     if (allocated(ydata)) deallocate(ydata)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
  end subroutine
  !
end function pole_fitting_LSQ

!======================
  function pole_fitting_LSQ_NL_cmplx(nx,xdata,ydata_,npoles_fit,params,thr,niterx) result(ierr)
  !======================
  implicit none
  !
  ! performs least-square fitting on input data.
  ! Fitting functions are poles (complex)
  ! ierr /=0 indicates a numerical problem
  !
  integer     :: nx,npoles_fit
  real(dp)    :: xdata(nx)
  complex(dp) :: ydata_(nx)
  real(dp)    :: params(*)
  real(dp)    :: thr
  integer     :: niterx
  integer     :: ierr
  !
  integer   :: ndim
  integer   :: i,j,ind,iter
  logical   :: lconv
  real(dp)  :: res,res0,trust
  complex(dp) :: cpole
  complex(dp), allocatable :: ydata(:,:)
  complex(dp), allocatable :: ydata_fit(:,:)
  complex(dp), allocatable :: Xmat(:,:)
  complex(dp), allocatable :: Xmat_d1(:,:)
  complex(dp), allocatable :: Xmat_d2(:,:)
  complex(dp), allocatable :: beta(:,:), alpha(:,:)
  complex(dp), allocatable :: ctmp(:,:)
  complex(dp), allocatable :: vtmp(:,:)
  complex(dp), allocatable :: Amat(:,:)
  complex(dp), allocatable :: Bmat(:,:)

  trust=0.2
  !
  ierr=0
  ndim=1+npoles_fit
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(Xmat(nx,ndim))
  allocate(Xmat_d1(nx,ndim))
  allocate(Xmat_d2(nx,ndim))
  allocate(beta(ndim,1))
  allocate(alpha(ndim,ndim))
  allocate(ctmp(nx,1))
  allocate(vtmp(1,nx))
  allocate(ydata_fit(nx,1))
  allocate(ydata(nx,1))
  ydata(:,1)=ydata_(:)

  res0=0.0
  Xmat_d1=0.0
  Xmat_d2=0.0
  lconv=.false.

  iter_loop:&
  do iter = 1, niterx

     !
     ! the first function is constant
     !
     Xmat(:,1)=1.0
     !
     do j = 2, ndim
       !
       ind = 2+(j-2)*3
       cpole = cmplx(params(ind),params(ind+1),dp)
       !
       do i = 1, nx
         Xmat(i,j)    = 1.0d0/(xdata(i)-cpole)
         Xmat_d1(i,j) = 1.0d0*params(ind+2)/(xdata(i)-cpole)**2
         Xmat_d2(i,j) = 2.0d0*params(ind+2)/(xdata(i)-cpole)**3
       enddo
       !
     enddo
     !
     call mat_mul(Bmat,Xmat,'C',ydata,'N',ndim,1,nx)
     call mat_mul(Amat,Xmat,'C',Xmat,'N', ndim,ndim,nx)

     !
     ! LinSys solution
     !
     call mat_sv(ndim,1,Amat,Bmat,ierr)
     !call mat_lsd(ndim,ndim,1,Amat,Bmat,1.0d-8,IERR=ierr)
     !
     if (ierr/=0) then
        call cleanup_loc()
        return
     endif

     !
     ! eval residuals
     !
     call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
     res = sqrt(real(dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1)),dp)/real(nx,dp))

     !
     ! eval beta
     !
     ctmp(1:nx,1)=ydata(1:nx,1)-ydata_fit(1:nx,1)
     call mat_mul(beta,Xmat_d1,'C',ctmp,'N',ndim,1,nx)
     !
     beta=-beta

     !
     ! eval alpha
     !
     call mat_mul(alpha,Xmat_d1,'C',Xmat_d1,'N',ndim,ndim,nx)
     !
     ctmp(1:nx,1)=ydata(1:nx,1)-ydata_fit(1:nx,1)
     call mat_mul(vtmp,ctmp,'C',Xmat_d2,'N',1,ndim,nx)
     !
     do j = 2, ndim
        alpha(j,j)=alpha(j,j) -vtmp(1,j)
     enddo
     !
     call mat_sv(ndim,1,alpha,beta,IERR=ierr)
     !call mat_lsd(ndim,ndim,1,alpha,beta,1.0d-8,IERR=ierr)
     if (ierr/=0) then
        call cleanup_loc()
        return
     endif

     !
     ! reshape data
     ! assuming real data
     ! (extension to complex quantities to be implemented)
     !
     params(1)=real(Bmat(1,1),dp)
     !
     do j = 2, ndim
       !
       ind = 2+ (j-2)*3
       params(ind+2)=real(Bmat(j,1),dp)
       !
       params(ind+0)=params(ind+0)+real(beta(1,j),dp)
       params(ind+1)=params(ind+1)+aimag(beta(1,j))
       !
     enddo
     !
     params(2+3*npoles_fit) = res
     !write(stdout,"(a,i5, f21.9)") "iter res", iter, res

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

  if (.not.lconv) ierr=-1

  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     if (allocated(Amat))  deallocate(Amat)
     if (allocated(Bmat))  deallocate(Bmat)
     if (allocated(Xmat))  deallocate(Xmat)
     if (allocated(beta))  deallocate(beta)
     if (allocated(alpha)) deallocate(alpha)
     if (allocated(ctmp))  deallocate(ctmp)
     if (allocated(vtmp))  deallocate(vtmp)
     if (allocated(Xmat_d1))   deallocate(Xmat_d1)
     if (allocated(Xmat_d2))   deallocate(Xmat_d2)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
     if (allocated(ydata)) deallocate(ydata)
  end subroutine
  !
end function pole_fitting_LSQ_NL_cmplx

!======================
  function pole_fitting_LSQ_NL(nx,xdata,ydata_,npoles_fit,params,thr,niterx,allow_offset) result(ierr)
  !======================
  implicit none
  !
  ! performs least-square fitting on input data.
  ! Fitting functions are real Lorentzians
  ! ierr /=0 indicates a numerical problem
  !
  integer  :: nx,npoles_fit
  real(dp) :: xdata(nx)
  real(dp) :: ydata_(nx)
  real(dp) :: params(*)
  real(dp) :: thr
  integer  :: niterx
  integer  :: ierr
  logical,  optional :: allow_offset
  !
  character(24) :: subname='pole_fitting_LSQ_NL_real'
  character(256):: minimization_type
  integer   :: ndim,ndim_d1
  integer   :: i,j,jj,ind,iter
  logical   :: lconv
  logical   :: compute_alpha, compute_alpha_diag
  real(dp)  :: res,res0,res1,res_scal,trust,cost,xmx0
  real(dp)  :: xx,dd,rr
  logical   :: allow_offset_
  real(dp), allocatable :: ydata(:,:)
  real(dp), allocatable :: ydata_fit(:,:)
  real(dp), allocatable :: Xmat(:,:)
  real(dp), allocatable :: Xmat_d1(:,:)
  real(dp), allocatable :: Xmat_d2(:,:,:)
  real(dp), allocatable :: beta(:,:)
  real(dp), allocatable :: alpha(:,:), alpha_diag(:)
  real(dp), allocatable :: ctmp(:,:)
  real(dp), allocatable :: Amat(:,:)
  real(dp), allocatable :: Bmat(:,:)
  !
  real(dp), external :: ddot

  trust=0.005
  cost=0.075
  !
  allow_offset_=.true.
  if (present(allow_offset)) allow_offset_ = allow_offset
  !
  !minimization_type="steepest_descent"
  !minimization_type="full_hessian"
  minimization_type="diag_hessian"
  if (nx < npoles_fit ) call errore(subname, 'too many DOF, solutions are not uniques',10)
  !
  compute_alpha=.false.
  compute_alpha_diag=.false.
  if (trim(minimization_type)=="full_hessian") compute_alpha=.true.
  if (trim(minimization_type)=="diag_hessian") compute_alpha_diag=.true.
  !
  !
  ierr=0
  ndim=npoles_fit+1
  ndim_d1=2*npoles_fit
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(Xmat(nx,ndim))
  allocate(Xmat_d1(nx,ndim_d1))
  allocate(Xmat_d2(nx,ndim_d1,2))
  allocate(beta(ndim_d1,1))
  allocate(ctmp(nx,1))
  allocate(ydata_fit(nx,1))
  allocate(ydata(nx,1))
  ydata(:,1)=ydata_(:)
  !
  ! more
  !
  if (compute_alpha) then
    !
    allocate(alpha(ndim_d1,ndim_d1))
    !
  elseif (compute_alpha_diag) then
    !
    allocate(alpha_diag(ndim_d1))
    !
  endif

  res_scal=maxval(abs(ydata))
  res0=0.0
  Xmat_d1=0.0
  Xmat_d2=0.0
  lconv=.false.

  iter_loop:&
  do iter = 1, niterx

     !
     !===========
     ! LSQ
     !===========
     !
     ! the first function is constant
     if (allow_offset_) then
        Xmat(:,1)=1.0
     else
        Xmat(:,1)=0.0
        Xmat(nx,1)=1.0
     endif
     !
     do j = 2, ndim
       !
       ind = 2+(j-2)*3
       xx  = params(ind)
       dd  = params(ind+1)
       !
       do i = 1, nx
         Xmat(i,j)    = dd/((xdata(i)-xx)**2 + dd**2)
       enddo
       !
     enddo

     !
     ! build and solve the Lin Sys
     !
     call mat_mul(Bmat,Xmat,'T',ydata,'N',ndim,1,nx)
     call mat_mul(Amat,Xmat,'T',Xmat,'N', ndim,ndim,nx)
     !
     call mat_sv(ndim,1,Amat,Bmat,ierr)
     !call mat_lsd(ndim,ndim,1,Amat,Bmat,1.0d-8,IERR=ierr)
     !
     if (ierr/=0) then
        call cleanup_loc()
        return
     endif

     !
     ! update data I
     !
     params(1)=Bmat(1,1)
     params(1)=0.0
     Xmat_d2=0.0
     !
     do j = 2, ndim
       !
       ind = 2+(j-2)*3
       xx  = params(ind)
       dd  = params(ind+1)
       rr  = Bmat(j,1)
       params(ind+2)=rr
       !
       jj = 2*(j-2)
       !
       do i = 1, nx
         xmx0=xdata(i)-xx
         !
         Xmat_d1(i,jj+1) = rr*2.0d0*dd*(xmx0)/((xmx0)**2 + dd**2)**2
         Xmat_d1(i,jj+2) = rr*((xmx0)**2 -dd**2)/((xmx0)**2 + dd**2)**2
         !
         Xmat_d2(i,jj+1,1) = rr*2.0d0*dd*( 3.0d0*(xmx0)**2 -dd**2)/( xmx0**2+dd**2)**3
         Xmat_d2(i,jj+1,2) = rr*2.0d0*xmx0* ((xmx0)**2 -3.0d0*dd**2)/((xmx0)**2 + dd**2)**3
         Xmat_d2(i,jj+2,2) = Xmat_d2(i,jj+1,2)
         Xmat_d2(i,jj+2,1) = -rr*2.0d0*dd/((xmx0)**2 + dd**2)**2 &
                                -rr*4.0d0*dd*((xmx0)**2 -dd**2)/((xmx0)**2 + dd**2)**3
       enddo
       !
     enddo
     !
     ! eval residuals
     !
     call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
     res1 = ddot(nx,ydata(:,1)-ydata_fit(:,1),1,ydata(:,1)-ydata_fit(:,1),1)/(real(nx,dp)*res_scal**2)

     !
     !====================
     ! non-linear fitting
     !====================
     !
     ! eval beta
     !
     ctmp(1:nx,1)=ydata(1:nx,1)-ydata_fit(1:nx,1)
     call mat_mul(beta,Xmat_d1,'T',ctmp,'N',ndim_d1,1,nx)
     !
     beta=-beta

     !
     ! eval alpha
     !
     if (compute_alpha) then
       !
       CALL errore(subname,"Xmat_d2 data structure changed",10)
       !
       call mat_mul(alpha,Xmat_d1,'T',Xmat_d1,'N',ndim_d1,ndim_d1,nx)
       !
       ctmp(1:nx,1)=ydata(1:nx,1)-ydata_fit(1:nx,1)
       !
       do j = 1, ndim_d1, 2
          alpha(j,j)    =alpha(j,j)     -ddot(nx,ctmp(:,1),1,Xmat_d2(:,j,1),1)
          alpha(j+1,j+1)=alpha(j+1,j+1) -ddot(nx,ctmp(:,1),1,Xmat_d2(:,j+1,1),1)
          alpha(j,j+1)  =alpha(j,j+1)   -ddot(nx,ctmp(:,1),1,Xmat_d2(:,j,2),1)
          alpha(j+1,j)  =alpha(j,j+1)
       enddo
       !
     elseif (compute_alpha_diag) then
       !
       ctmp(1:nx,1)=ydata(1:nx,1)-ydata_fit(1:nx,1)
       !
       do j = 1, ndim_d1
          alpha_diag(j) = ddot(nx, Xmat_d1(:,j), 1, Xmat_d1(:,j), 1)
          alpha_diag(j) = alpha_diag(j) -ddot(nx,ctmp(:,1),1,Xmat_d2(:,j,1),1)
       enddo
       !
     endif
     !
     ! minimization strategy
     !
     select case(trim(minimization_type))
     case("stepest_descent")
       !
       ! rescale beta
       !
       beta = cost * beta
       !
     case("full_hessian")
       !
       !call mat_sv(ndim,1,alpha,beta,ierr)
       call mat_lsd(ndim,ndim,1,alpha,beta,1.0d-6,IERR=ierr)
       !
       if (ierr/=0) then
          call cleanup_loc()
          return
       endif
       !
       beta = 0.5*beta
       !
     case("diag_hessian")
       !
       do j = 1, ndim_d1
          !
          if ( alpha_diag(j) <= 1.0d-5 ) then
             beta(j,1) = cost * beta(j,1)
          else
             beta(j,1) = beta(j,1)/alpha_diag(j)
          endif
          !
       enddo
       !
     case default
       call errore(subname,"invalid minimization_type = "//trim(minimization_type),10)
     end select
     !
     ! trust radius
     !
     do j = 1, ndim_d1
       if (beta(j,1)<-trust) beta(j,1)=-trust
       if (beta(j,1) >trust) beta(j,1)=trust
     enddo
     !
     ! update data II
     !
     do j = 2, ndim
       !
       jj= 2*(j-2)
       ind = 2+ (j-2)*3
       !
       params(ind+0)=params(ind+0)-beta(jj+1,1)
       params(ind+1)=params(ind+1)-beta(jj+2,1)
       !
       !
       xx  = params(ind)
       dd  = params(ind+1)
       rr  = params(ind+2)
       !
       do i = 1, nx
         Xmat(i,j)    = dd/((xdata(i)-xx)**2 + dd**2)
       enddo
       !
     enddo
     !
     ! eval residuals
     !
     call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
     res = ddot(nx,ydata(:,1)-ydata_fit(:,1),1,ydata(:,1)-ydata_fit(:,1),1)/(real(nx,dp)*res_scal**2)
     !
     params(2+3*npoles_fit) = res
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
  !
  if (.not.allow_offset_) params(1)=0.0d0
  !
  if (.not.lconv) ierr=-1
  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     if (allocated(Amat))  deallocate(Amat)
     if (allocated(Bmat))  deallocate(Bmat)
     if (allocated(Xmat))  deallocate(Xmat)
     if (allocated(beta))  deallocate(beta)
     if (allocated(alpha)) deallocate(alpha)
     if (allocated(alpha_diag)) deallocate(alpha_diag)
     if (allocated(ctmp))  deallocate(ctmp)
     if (allocated(Xmat_d1))   deallocate(Xmat_d1)
     if (allocated(Xmat_d2))   deallocate(Xmat_d2)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
     if (allocated(ydata)) deallocate(ydata)
  end subroutine
  !
end function pole_fitting_LSQ_NL

#ifdef __TODO
!======================
  function pole_fitting_pade(nx,xdata,ydata_,npoles_fit,params,allow_offset) result(ierr)
  !======================
  implicit none
  !
  ! performs pade-approximant fitting on input data.
  ! Fitting functions are real lorentzians
  ! ierr /=0 indicates a numerical problem
  !
  integer  :: nx,npoles_fit
  real(dp) :: xdata(nx)
  real(dp) :: ydata_(nx)
  real(dp) :: params(*)
  logical, optional :: allow_offset
  integer  :: ierr
  !
  integer  :: ndim
  integer  :: i,j,ind
  logical  :: allow_offset_
  real(dp) :: res,dd,xx,rr
  real(dp), allocatable :: ydata_fit(:,:)
  real(dp), allocatable :: ydata(:,:)
  real(dp), allocatable :: Xmat(:,:)
  real(dp), allocatable :: Amat(:,:)
  real(dp), allocatable :: Bmat(:,:)

  ierr=0
  ndim=npoles_fit+1
  !
  allow_offset_=.true.
  if (present(allow_offset)) allow_offset_=allow_offset
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(Xmat(nx,ndim))
  allocate(ydata_fit(nx,1))
  allocate(ydata(nx,1))
  ydata(:,1)=ydata_(:)
  !
  !
  ! the first function is constant
  if (allow_offset_) then
    Xmat(:,1)=1.0
  else
    Xmat(:,1)=0.0
    Xmat(nx,1)=1.0
  endif
  !
  do j = 2, ndim
    !
    ind = 2+(j-2)*3
    xx  = params(ind)
    dd  = params(ind+1)
    !
    do i = 1, nx
      Xmat(i,j)    = dd/((xdata(i)-xx)**2 + dd**2)
    enddo
    !
  enddo
  !
  call mat_mul(Bmat,Xmat,'T',ydata,'N',ndim,1,nx)
  call mat_mul(Amat,Xmat,'T',Xmat,'N', ndim,ndim,nx)

  !
  ! solution
  !
  call mat_sv(ndim,1,Amat,Bmat,ierr)
  if (ierr/=0) then
     call cleanup_loc()
     return
  endif

  !
  ! eval residuals
  !
  call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
  res = dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1))/real(nx,dp)

  !
  ! reshape data
  ! assuming real data
  ! (extension to complex quantities to be implemented)
  !
  params(1)=0.0d0
  if (allow_offset_) params(1)=Bmat(1,1)
  !
  do j = 2, ndim
    ind = 2+(j-2)*3
    params(ind+2)=Bmat(j,1)
  enddo
  params(2+3*npoles_fit) = res
  !

  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     if (allocated(Amat))  deallocate(Amat)
     if (allocated(Bmat))  deallocate(Bmat)
     if (allocated(Xmat))  deallocate(Xmat)
     if (allocated(ydata)) deallocate(ydata)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
  end subroutine
  !
end function pole_fitting_pade
#endif

!======================
  function pole_fitting_pade_cmplx(nx,xdata,ydata_,npoles_fit,params,allow_offset,even) result(ierr)
  !======================
  implicit none
  !
  ! performs least-square fitting on input data.
  ! Fitting functions are poles (complex)
  ! ierr /=0 indicates a numerical problem
  !
  integer     :: nx,npoles_fit
  real(dp)    :: xdata(nx)
  complex(dp) :: ydata_(nx)
  complex(dp) :: params(*)
  integer     :: ierr
  logical, optional :: allow_offset
  logical, optional :: even
  !
  integer   :: ndim,countt
  integer   :: i,j,k,ind
  logical  :: allow_offset_
  logical  :: even_
  complex(dp) :: res,ctmp(4)
  complex(dp), allocatable :: ydata_fit(:,:)
  complex(dp), allocatable :: ydata(:,:)
  complex(dp), allocatable :: Amat(:,:)
  complex(dp), allocatable :: Bmat(:,:)
  complex(dp), allocatable :: p_roots(:), q_roots(:)
  !
  type(poly_t) :: p_pol, q_pol

  ierr=0
  ndim=2*npoles_fit
  !
  allow_offset_=.true.
  if (present(allow_offset)) allow_offset_=allow_offset
  !
  even_=.false.
  if (present(even)) even_=even
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(ydata_fit(nx,1))
  allocate(ydata(nx,1))
  !
  if(.not. even_) then
    allocate(q_roots(npoles_fit))
    allocate(p_roots(npoles_fit))
  else
    allocate(q_roots(2*npoles_fit))
    allocate(p_roots(2*npoles_fit))
  endif
  !
  ydata(:,1)=ydata_(:)

  do j = 1, npoles_fit
  do k = 1, npoles_fit
    !
    ctmp=0.0
    if(.not. even_) then
      do i = 1, nx
        ctmp(1)=ctmp(1)+xdata(i)**(j+k-2)
        ctmp(2)=ctmp(2)+xdata(i)**(j+k-2)*ydata(i,1)
        ctmp(3)=ctmp(3)+xdata(i)**(j+k-2)*conjg(ydata(i,1))
        ctmp(4)=ctmp(4)+xdata(i)**(j+k-2)*ydata(i,1)*conjg(ydata(i,1))
      enddo
    else
      do i = 1, nx
        ctmp(1)=ctmp(1)+xdata(i)**(2*(j+k-2))
        ctmp(2)=ctmp(2)+xdata(i)**(2*(j+k-2))*ydata(i,1)
        ctmp(3)=ctmp(3)+xdata(i)**(2*(j+k-2))*conjg(ydata(i,1))
        ctmp(4)=ctmp(4)+xdata(i)**(2*(j+k-2))*ydata(i,1)*conjg(ydata(i,1))
      enddo
    endif
    ctmp=ctmp/real(nx,dp)
    !
    Amat(j,k) = ctmp(1)
    Amat(j,k+npoles_fit) = -ctmp(2)
    Amat(j+npoles_fit,k) = ctmp(3)
    Amat(j+npoles_fit,k+npoles_fit) = -ctmp(4)
    !
  enddo
  enddo
  !
  do j = 1, npoles_fit
    !
    ctmp=0.0
    if(.not. even_) then
      do i = 1, nx
        ctmp(1)=ctmp(1)+xdata(i)**(j+npoles_fit-1)*ydata(i,1)
        ctmp(2)=ctmp(2)+xdata(i)**(j+npoles_fit-1)*ydata(i,1)*conjg(ydata(i,1))
      enddo
    else
      do i = 1, nx
        ctmp(1)=ctmp(1)+xdata(i)**(2*(j+npoles_fit-1))*ydata(i,1)
        ctmp(2)=ctmp(2)+xdata(i)**(2*(j+npoles_fit-1))*ydata(i,1)*conjg(ydata(i,1))
      enddo
    endif
    ctmp=ctmp/real(nx,dp)
    !
    Bmat(j,1)=ctmp(1)
    Bmat(j+npoles_fit,1)=ctmp(2)
    !
  enddo

  !
  ! solution
  !
  call mat_sv(ndim,1,Amat,Bmat,ierr)
  if (ierr/=0) then
     call cleanup_loc()
     return
  endif

  !
  ! eval residuals
  !
  if(.not. even_) then
    call poly_reset(p_pol,npoles_fit-1)
    call poly_reset(q_pol,npoles_fit)
    !
    p_pol%c(0:npoles_fit-1)=Bmat(1:npoles_fit,1)
    q_pol%c(0:npoles_fit-1)=Bmat(1+npoles_fit:2*npoles_fit,1)
    q_pol%c(npoles_fit)=1.0_dp
    !
  else
    call poly_reset(p_pol,2*(npoles_fit-1))
    call poly_reset(q_pol,2*(npoles_fit))
    !
    do i = 1,npoles_fit-1
      !
      p_pol%c(2*i-2) = Bmat(i,1)
      p_pol%c(2*i-1) = 0.0_dp
      q_pol%c(2*i-2) = Bmat(npoles_fit+i,1)
      q_pol%c(2*i-1) = 0.0_dp
      !
    enddo
    !
    p_pol%c(2*npoles_fit-2) = Bmat(npoles_fit,1)
    q_pol%c(2*npoles_fit-2) = Bmat(2*npoles_fit,1)
    q_pol%c(2*npoles_fit)   = 1.0_dp 
    q_pol%c(2*npoles_fit-1)   = 0.0_dp 
      !
  endif
  !
  ydata_fit(:,1) = poly_eval(xdata,nx,p_pol) / poly_eval(xdata,nx,q_pol)
  res = dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1))/real(nx,dp)

  !
  !
  ! reshape data
  !
  call poly_roots(q_roots,q_pol)
  call poly_roots(p_roots,p_pol)
  !
  ! XXX
  write(0,"(/,'P_ROOTS')")
  write(0,"(2f21.14)") p_roots
  write(0,"(/,'Q_ROOTS')")
  write(0,"(2f21.14)") q_roots
  ! XXX
  !
  params(1)=0.0d0
  if (allow_offset_) params(1)=0.0  ! to be fixed
  !
  if(.not. even_) then
    do j = 1, npoles_fit
      !
      ind = 2+(j-1)*3
      params(ind+0)=real(q_roots(j),dp)
      params(ind+1)=aimag(q_roots(j))
      !
      ctmp(1)=poly_eval(q_roots(j),p_pol)
      do i = 1, npoles_fit
        if (i==j) cycle
        ctmp(1)=ctmp(1)/(q_roots(j)-q_roots(i))
      enddo
      params(ind+2)=ctmp(1)
      !
    enddo
  else
    countt=0
    do j = 1, 2*npoles_fit
      !
      if(real(q_roots(j),dp)>=0) then
        countt = countt + 1
      else
        cycle
      endif
      !
      ind = 2+(countt-1)*3
      params(ind+0)=real(q_roots(j),dp)
      params(ind+1)=aimag(q_roots(j))
      !
      ctmp(1)=poly_eval(q_roots(j),p_pol)
      do i = 1, 2*npoles_fit
        if (i==j) cycle
        ctmp(1)=ctmp(1)/(q_roots(j)-q_roots(i))
      enddo
      params(ind+2)=ctmp(1)
      !
      if(countt==npoles_fit) exit
      !
    enddo
  endif
  !
  params(2+3*npoles_fit) = res

  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     if (allocated(Amat))  deallocate(Amat)
     if (allocated(Bmat))  deallocate(Bmat)
     if (allocated(ydata)) deallocate(ydata)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
     if (allocated(q_roots)) deallocate(q_roots)
     if (allocated(p_roots)) deallocate(p_roots)
     call poly_reset(p_pol)
     call poly_reset(q_pol)
  end subroutine
  !
end function pole_fitting_pade_cmplx

end module fitting_m

