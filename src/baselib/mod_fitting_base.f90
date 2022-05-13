!======================
module fitting_base_m
  !======================
  use kinds
  use constants
  use util_module,       only: mat_mul,mat_sv,mat_lsd
  use poly_m
  use func4fit_m
  use nonneg_leastsq 
  implicit none
  private

  logical :: fitting_verbosity=.false.



  public :: fitting_verbosity

  public :: dlsei_wrap 
  public :: nnls_wrap 
  public :: pole_fitting_LSQ
  public :: pole_fitting_LSQ_cmplx
  public :: pole_fitting_LSQ_cmplx_lsd
  public :: pole_fitting_LSQ_NL
  public :: pole_fitting_LSQ_NL_cmplx
  public :: pole_fitting_LSQ_NL_cmplx_c
  public :: pole_fitting_pade_cmplx
  !
  public :: pole_interp_pade_cmplx
  public :: pole_interp_pade_cmplx_dn
  public :: pole_interp_pade_cmplx_noLS
  public :: pole_rational_interp_cmplx
  !
  public :: rational_interp_A
  public :: rational_interp_B

contains

  subroutine dlsei_wrap(mAmat,nAmat,mCmat,Amat,Bmat,Cmat,res,ierr)
    implicit none
    real(dp)            :: Amat(:,:),Bmat(:,:),Cmat(:,:)
    integer, intent(in) :: mAmat,nAmat,mCmat
    real(dp),  intent(out):: res
    integer, intent(out):: ierr
    !
    real(dp), allocatable :: w(:,:),ws(:),x(:)
    real(dp), allocatable :: optionV(:) 
    integer,  allocatable :: ip(:)
    integer               :: MDW,N,ME,MA,MG,K,LIP,LW,i
    real(dp)              :: RNORME,RNORML
    !
    ! allocations

    N = nAmat
    ME = mCmat
    MA = mAmat
    MG = nAmat 
    K = max(MA+MG,N)
    MDW=ME+MA+MG
    LW = 2*(ME+N)+K+(MG+2)*(N+7)
    LIP = MG+2*N+2
    allocate(w(MDW,N+1))
    allocate(optionV(1))
    allocate(x(N))
    allocate(ws(LW))
    allocate(ip(LIP))
    optionV(1) = 1
    ip(1)=LW
    ip(2)=LIP
    w = 0._dp
    !
    ! set up work matrix
    !
    ! constraints (equalities)
    do i =1,ME
      w(i,1:N+1)=Cmat(i,1:N+1)
    enddo
    ! equations (approximated equalities) 
    do i =1,MA
      w(ME+i,1:N)=Amat(i,1:N)
      w(ME+i,N+1)=Bmat(i,1)
    enddo
    ! inequalities G=identity, H=0  (H is a vector)
    do i =1,MG
      w(ME+MA+i,i)=1._dp
      w(ME+MA+i,N+1)=0._dp
    enddo
    !
    CALL DLSEI (w, MDW, ME, MA, MG, N, optionV, x, RNORME,RNORML, ierr, ws, ip)
    !
    ! prune negative amps which are noise
    do i =1,N
      if(x(i)>0) then
        Bmat(i,1)=x(i)
      else
        Bmat(i,1)=0._dp
      endif
    enddo
    !
    res = RNORME
    !
    deallocate(w,optionV,x,ws,ip)
    !
    end subroutine dlsei_wrap

    subroutine nnls_wrap(mAmat,nAmat,Amat,Bmat,res,ierr)
      implicit none
      integer,   intent(in) :: mAmat,nAmat
      real(dp)              :: Amat(:,:),Bmat(:,:)
      real(dp),  intent(out):: res
      integer,   intent(out):: ierr
      !
      real(dp), allocatable :: Bvec(:)
      integer,  allocatable :: itmp(:)
      real(dp), allocatable :: work(:)
      !
      ! note that we suppose Bmat(:nAmat,1) to have at least nAmat lines (the system cannot be underdetermined)
      !
      allocate(Bvec(nAmat))
      allocate(itmp(nAmat))
      allocate(work(nAmat))
      !
      Bvec(:) = Bmat(:,1)
      call nnls(Amat, mAmat, nAmat, Bvec, Bmat(:nAmat,1), res, work, itmp, ierr)
      !
      deallocate(Bvec,itmp,work)
      !
    end subroutine nnls_wrap 


!======================
  function pole_fitting_LSQ_cmplx_lsd(nx,xdata,ydata_,npoles_fit,params) result(ierr)
  !======================
  implicit none
  !
  ! performs a modified (avoiding Xmat*Xmat) cmplx least-square fitting on input data.
  ! Fitting functions are poles (complex)
  ! ierr /=0 indicates a numerical problem
  !
  integer     :: nx,npoles_fit
  complex(dp) :: xdata(nx)
  complex(dp) :: ydata_(nx)
  complex(dp) :: params(*)
  integer     :: ierr
  !
  integer   :: ndim
  integer   :: i,j,ind
  complex(dp) :: cpole,res
  complex(dp), allocatable :: ydata_fit(:,:)
  complex(dp), allocatable :: ydata(:,:)
  complex(dp), allocatable :: Bmat(:,:)
  complex(dp), allocatable :: Xmat(:,:)

  ndim=npoles_fit
  !
  ierr=0
  !
  allocate(Xmat(nx,ndim))
  allocate(ydata_fit(nx,1))
  allocate(ydata(nx,1))
  allocate(Bmat(nx,1))
  ydata(:,1)=ydata_(:)
  Bmat = ydata

  !
  do j = 1, npoles_fit 
    !
    ind = 1+(j-1)*3
    cpole = cmplx(real(params(ind+1),dp),real(params(ind+2),dp),dp)
    !
    do i = 1, nx
      Xmat(i,j) = params(ind+3)/(xdata(i)-cpole)
    enddo
    !
  enddo
  !
  ! solution
  ! for rcond use epsilon(1._dp)*max(nx,ndim), default in numpy
  !
  call mat_lsd(nx, ndim, 1, Xmat, Bmat, epsilon(1._dp)*max(nx,ndim), IERR=ierr)
  !
  if (ierr/=0) then
     call cleanup_loc()
     return
  endif
  !
  !
  ! eval residuals
  !
  call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
  res = dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1)) ! /real(nx,dp)

  params(1)=0.0d0
  !
  do j = 1, npoles_fit
    ind = 1+(j-1)*3
    params(ind+3)=Bmat(j,1)*params(ind+3)
  enddo
  params(2+3*npoles_fit) = res

  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     if (allocated(Bmat))  deallocate(Bmat)
     if (allocated(Xmat))  deallocate(Xmat)
     if (allocated(ydata)) deallocate(ydata)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
  end subroutine
  !
end function pole_fitting_LSQ_cmplx_lsd


!======================
  function pole_fitting_LSQ_cmplx(nx,xdata,ydata_,npoles_fit,params,allow_offset,sum_amp) result(ierr)
  !======================
  implicit none
  !
  ! performs least-square fitting on input data.
  ! Fitting functions are poles (complex)
  ! ierr /=0 indicates a numerical problem
  !
  integer     :: nx,npoles_fit
  complex(dp) :: xdata(nx)
  complex(dp) :: ydata_(nx)
  complex(dp) :: params(*)
  integer     :: ierr
  logical,     optional :: allow_offset
  complex(dp), optional :: sum_amp 
  !
  integer   :: ndim
  integer   :: i,j,ind
  logical  :: allow_offset_
  logical  :: fix_sum_amp 
  complex(dp) :: cpole,res
  complex(dp), allocatable :: ydata_fit(:,:)
  complex(dp), allocatable :: ydata(:,:)
  complex(dp), allocatable :: Xmat(:,:)
  complex(dp), allocatable :: Amat(:,:)
  complex(dp), allocatable :: Bmat(:,:)
  complex(dp), allocatable :: Amat_(:,:)
  complex(dp), allocatable :: Bmat_(:,:)
  real(dp),    allocatable :: RAmat_(:,:)
  real(dp),    allocatable :: RBmat_(:,:)

  allow_offset_=.true.
  if (present(allow_offset)) allow_offset_=allow_offset
  fix_sum_amp=.false.
  if (present(sum_amp)) then
    if(huge(real(sum_amp,dp))/=real(sum_amp,dp)) fix_sum_amp = .true.
  endif
  !
  ierr=0
  !
  if (allow_offset_) then
    ndim=npoles_fit+1
  else
    ndim=npoles_fit
  endif
  !
  if (fix_sum_amp) then
    allocate(Amat_(ndim+1,ndim+1))
    allocate(Bmat_(ndim+1,1))
    allocate(RAmat_(ndim+1,ndim+1))
    allocate(RBmat_(ndim+1,1))
  else
    allocate(Amat_(ndim,ndim))
    allocate(Bmat_(ndim,1))
    allocate(RAmat_(ndim,ndim))
    allocate(RBmat_(ndim,1))
  endif
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(Xmat(nx,ndim))
  allocate(ydata_fit(nx,1))
  allocate(ydata(nx,1))
  ydata(:,1)=ydata_(:)

  !
  do j = 1, npoles_fit 
    !
    ind = 1+(j-1)*3
    cpole = cmplx(real(params(ind+1),dp),real(params(ind+2),dp),dp)
    !
    do i = 1, nx
      Xmat(i,j) = params(ind+3)/(xdata(i)-cpole)
    enddo
    !
  enddo
  if (allow_offset_) Xmat(:,ndim)=1.0
  !
  call mat_mul(Bmat,Xmat,'C',ydata,'N',ndim,1,nx)
  call mat_mul(Amat,Xmat,'C',Xmat,'N', ndim,ndim,nx)
  !
  Amat_(1:ndim,1:ndim) = Amat(1:ndim,1:ndim) 
  Bmat_(1:ndim,1) = Bmat(1:ndim,1) 
  !
  if (fix_sum_amp) then
    Amat_(ndim+1,:) = 1.0
    Bmat_(ndim+1,1) = sum_amp 
    Amat_(:,ndim+1) = 1.0
    Amat_(ndim+1,ndim+1) = 0.0 
  endif
  !
  RBmat_ = real(Bmat_,dp)
  RAmat_ = real(Amat_,dp)

  !
  ! solution
  !
  ! uncomment for real amplitudes and comment the next lines
  !
  if (fix_sum_amp) then
    !call mat_sv(ndim+1,1,RAmat_,RBmat_,ierr)
    call mat_sv(ndim+1,1,Amat_,Bmat_,ierr)
  else
    !call mat_sv(ndim,1,RAmat_,RBmat_,ierr)
    call mat_sv(ndim,1,Amat_,Bmat_,ierr)
  endif
  !
  if (ierr/=0) then
     call cleanup_loc()
     return
  endif
  !
  !
  ! uncomment for real amplitudes
  !
  Bmat(1:ndim,1) = Bmat_(1:ndim,1) 
  !Bmat(1:ndim,1) = RBmat_(1:ndim,1) 

  !
  ! eval residuals
  !
  call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
  res = dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1))/real(nx,dp)

  params(1)=0.0d0
  if (allow_offset_) params(1)=Bmat(ndim,1)
  !
  do j = 1, npoles_fit
    ind = 1+(j-1)*3
    params(ind+3)=Bmat(j,1)*params(ind+3)
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
     if (allocated(Amat_))  deallocate(Amat_)
     if (allocated(Bmat_))  deallocate(Bmat_)
     if (allocated(RAmat_)) deallocate(RAmat_)
     if (allocated(RBmat_)) deallocate(RBmat_)
     if (allocated(Xmat))  deallocate(Xmat)
     if (allocated(ydata)) deallocate(ydata)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
  end subroutine
  !
end function pole_fitting_LSQ_cmplx

!======================
  function pole_fitting_LSQ(nx,xdata,ydata_,npoles_fit,params,allow_offset,degree,nCnstr,xyCnstr) result(ierr)
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
  integer, optional :: degree 
  integer, optional :: nCnstr 
  real(dp),optional :: xyCnstr(:,:)
  integer  :: ierr
  !
  integer  :: ndim
  integer  :: i,j,ind,m
  logical  :: allow_offset_
  integer  :: degree_ 
  integer  :: nCnstr_ 
  real(dp) :: res,dd,xx
  real(dp), allocatable :: ydata_fit(:,:)
  real(dp), allocatable :: ydata(:,:)
  real(dp), allocatable :: Xmat(:,:)
  real(dp), allocatable :: Amat(:,:)
  real(dp), allocatable :: Bmat(:,:)
  real(dp), allocatable :: Cmat(:,:)
  real(dp), allocatable :: Bvec(:)
  integer,  allocatable :: itmp(:)
  real(dp), allocatable :: work(:)

  allow_offset_=.true.
  if (present(allow_offset)) allow_offset_=allow_offset
  degree_=1
  if (present(degree)) degree_=degree
  nCnstr_=0
  if (present(nCnstr)) then
    nCnstr_=nCnstr
    if (.not.present(xyCnstr)) then
      ierr=1
      return
    endif
  endif
  !
  ierr=0
  !
  if (allow_offset_) then
    ndim=npoles_fit+1
  else
    ndim=npoles_fit
  endif
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(Cmat(nCnstr_,ndim+1))
  allocate(Bvec(ndim))
  allocate(itmp(ndim))
  allocate(work(ndim))
  allocate(Xmat(nx,ndim))
  allocate(ydata_fit(nx,1))
  allocate(ydata(nx,1))
  ydata(:,1)=ydata_(:)
  !
  do j = 1, npoles_fit 
    !
    ind = 1+(j-1)*3
    xx  = params(ind+1)
    dd  = params(ind+2)
    !
    do i = 1, nx
      Xmat(i,j)    = dd**(2*degree_-1)/((xdata(i)-xx)**(2*degree_) + dd**(2*degree_))
    enddo
    ! constraints
    do m=1,nCnstr_ 
      ! special constrain for amplitude sum
      ! to be fixed more for a better usage
      if (xyCnstr(m,1)>xdata(nx)) then
        Cmat(m,j)=1._dp
        Cmat(m,ndim+1) =integrate_n_momentum(0,nx,abs(ydata(:,1)),xdata,INFINT=.true.)/PI
      else
        Cmat(m,j) = dd**(2*degree_-1)/((xyCnstr(m,1)-xx)**(2*degree_) + dd**(2*degree_))
        Cmat(m,ndim+1) = xyCnstr(m,2)
      endif
    enddo
    !
  enddo
  !
  if (allow_offset_) then
    ! to be checked
    Xmat(:,ndim)=1.0
  endif
  !
  call mat_mul(Bmat,Xmat,'T',ydata,'N',ndim,1,nx)
  call mat_mul(Amat,Xmat,'T',Xmat,'N', ndim,ndim,nx)

  !
  ! solution
  !
  !call mat_sv(ndim,1,Amat,Bmat,ierr)
  if (nCnstr_==0) then
    call nnls_wrap(ndim,ndim,Amat,Bmat,res,ierr)
    ierr = ierr - 1
  else
    call dlsei_wrap(ndim,ndim,nCnstr_,Amat,Bmat,Cmat,res,ierr) 
  endif
  if (ierr/=0) then
     call cleanup_loc()
     return
  endif

  !
  ! eval residuals
  !
  call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
  res = dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1))/real(nx,dp)

  params(1)=0.0d0
  if (allow_offset_) params(1)=Bmat(ndim,1)
  !
  do j = 1, npoles_fit
    ind = 1+(j-1)*3
    params(ind+3)=Bmat(j,1)
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
     if (allocated(Cmat))  deallocate(Cmat)
     if (allocated(Xmat))  deallocate(Xmat)
     if (allocated(ydata)) deallocate(ydata)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
     if (allocated(Bvec))  deallocate(Bvec)
     if (allocated(itmp))  deallocate(itmp)
     if (allocated(work))  deallocate(work)
  end subroutine
  !
end function pole_fitting_LSQ

!======================
  function pole_fitting_LSQ_NL_cmplx(nx,xdata,ydata_,npoles_fit,cparams,thr,niterx,allow_offset) result(ierr)
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
  complex(dp) :: cparams(*)
  real(dp)    :: thr
  integer     :: niterx
  integer     :: ierr
  logical, optional :: allow_offset
  !
  integer   :: ndim
  integer   :: i,j,ind,iter
  logical   :: lconv,allow_offset_
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
  allow_offset_=.true.
  if (present(allow_offset)) allow_offset_=allow_offset
  ndim = npoles_fit + 1
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
       cpole = cparams(ind)+CI*cparams(ind+1)
       !
       do i = 1, nx
         Xmat(i,j)    = 1.0d0/(xdata(i)-cpole)
         Xmat_d1(i,j) = 1.0d0*cparams(ind+2)/(xdata(i)-cpole)**2
         Xmat_d2(i,j) = 2.0d0*cparams(ind+2)/(xdata(i)-cpole)**3
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
     !
     if (ierr/=0) then
        call cleanup_loc()
        return
     endif

     !
     ! reshape data
     ! assuming real data
     ! (extension to complex quantities to be implemented)
     !
     cparams(1)=Bmat(1,1)
     !
     do j = 2, ndim
       !
       ind = 2+ (j-2)*3
       cparams(ind+2)=Bmat(j,1)
       !
       cparams(ind+0)=cparams(ind+0)+real(beta(j,1),dp)
       cparams(ind+1)=cparams(ind+1)+aimag(beta(j,1))
       !
     enddo
     !
     cparams(2+3*npoles_fit) = res
     !write(6,"(a,i5, f21.9)") "iter res", iter, res

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
  function pole_fitting_LSQ_NL_cmplx_c(nx,xdata,ydata_,npoles_fit,cparams,thr,niterx,allow_offset) result(ierr)
  !======================
  implicit none
  !
  ! performs least-square fitting on input data.
  ! Fitting functions are real Lorentzians
  ! ierr /=0 indicates a numerical problem
  !
  integer     :: nx,npoles_fit
  complex(dp) :: xdata(nx)
  complex(dp) :: ydata_(nx)
  complex(dp) :: cparams(2+3*npoles_fit)
  real(dp)    :: thr
  integer     :: niterx
  integer     :: ierr
  logical,  optional :: allow_offset
  !
  character(25) :: subname='pole_fitting_LSQ_NL_cmplx'
  character(256):: minimization_type
  integer       :: ndim,ndim_d1
  integer       :: i,j,ind,iter, niter_check_trust
  real(dp)      :: res0,res1,res2,res_scal
  real(dp)      :: trust,cost
  real(dp)      :: xx,dd
  complex(dp)   :: rr, cpole
  logical       :: allow_offset_
  logical       :: lconv
  logical       :: compute_alpha, compute_alpha_diag
  logical       :: do_reduce_trust, do_enlarge_trust, do_skip_LSQ
  complex(dp), allocatable :: ydata(:,:)
  complex(dp), allocatable :: ydata_fit(:,:)
  complex(dp), allocatable :: Xmat(:,:)
  complex(dp), allocatable :: Xmat_d1(:,:)
  complex(dp), allocatable :: Xmat_d2(:,:)
  complex(dp), allocatable :: beta(:,:)
  complex(dp), allocatable :: alpha(:,:), alpha_diag(:)
  complex(dp), allocatable :: ctmp(:,:)
  complex(dp), allocatable :: Amat(:,:)
  complex(dp), allocatable :: Bmat(:,:)
  complex(dp), allocatable :: cparams0(:)

  trust=0.02
  cost=0.05
  niter_check_trust=25
  !
  do_reduce_trust=.false.
  do_enlarge_trust=.false.
  do_skip_LSQ=.false.
  !
  !minimization_type="steepest_descent"
  minimization_type="diag_hessian"
  !minimization_type="full_hessian"
  !
  if (nx < npoles_fit ) call errore(subname, 'too many DOF, solutions are not uniques',10)
  !
  compute_alpha=.false.
  compute_alpha_diag=.false.
  if (trim(minimization_type)=="full_hessian") compute_alpha=.true.
  if (trim(minimization_type)=="diag_hessian") compute_alpha_diag=.true.
  !
  allow_offset_=.false.
  if (present(allow_offset)) allow_offset_=allow_offset
  !
  !
  ierr=0
  !
  if (allow_offset_) then
    ndim=npoles_fit+1
  else
    ndim=npoles_fit
  endif
  !
  ndim_d1=npoles_fit
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(Xmat(nx,ndim))
  allocate(Xmat_d1(nx,ndim_d1))
  allocate(Xmat_d2(nx,ndim_d1))
  allocate(beta(ndim_d1,1))
  allocate(ctmp(nx,1))
  allocate(ydata_fit(nx,1))
  allocate(ydata(nx,1))
  allocate(cparams0(2+3*npoles_fit))
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
  !
  res_scal=maxval(abs(ydata))
  res0=0.0
  lconv=.false.
  !
  iter_loop:&
  do iter = 1, niterx

     ! 
     ! reduce/enlarge trust radius if needed
     ! 
     if (mod(iter,niter_check_trust)==0) do_enlarge_trust=.true.
     !
     if (do_enlarge_trust.and..not.do_reduce_trust) then
       trust=trust*2.0_dp
       do_enlarge_trust=.false.
     endif
     !
     if (do_reduce_trust) then
       trust=trust/2.0_dp
       do_reduce_trust=.false.
     endif

     !
     !===========
     ! LSQ
     !===========
     !
     whether_skip_LSQ:&
     if (.not. do_skip_LSQ) then
       !
       do j = 1, npoles_fit 
         !
         ind = 1+(j-1)*3
         cpole = real(cparams(ind+1), dp) + CI * real(cparams(ind+2), dp)
         !
         do i = 1, nx
           Xmat(i,j)    = 1.0d0/(xdata(i)-cpole)
         enddo
         !
       enddo
       if (allow_offset_) Xmat(:,ndim)=1.0

       !
       ! build and solve the Lin Sys
       !
       call mat_mul(Bmat,Xmat,'C',ydata,'N',ndim,1,nx)
       call mat_mul(Amat,Xmat,'C',Xmat,'N', ndim,ndim,nx)
       !
       call mat_sv(ndim,1,Amat,Bmat,ierr)
       !
       if (ierr/=0) then
          call cleanup_loc()
          return
       endif

       !
       ! update data I
       !
       cparams(1)=0.0d0
       if (allow_offset_) cparams(1)=Bmat(ndim,1)
       !
       do j = 1, npoles_fit
         !
         ind = 1+(j-1)*3
         cparams(ind+3)=Bmat(j,1)
         !
       enddo
       !
       ! eval residuals
       !
       call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
       res1 = real(dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1)),dp)/(real(nx,dp)*res_scal**2)
       !
       if (fitting_verbosity) &
          write(6,"(a,i5, f21.9)") "iter res1", iter, res1
       !
     endif whether_skip_LSQ
     !
     do_skip_LSQ=.false.

     !
     !====================
     ! non-linear fitting
     !====================
     !
     ! evaluate derivatives, Xmat_d1, Xmat_d2
     !
     Xmat_d2=0.0
     Xmat_d1=0.0
     !
     do j = 1, npoles_fit
       !
       ind = 1+(j-1)*3
       !
       xx  = real(cparams(ind+1), dp)
       dd  = real(cparams(ind+2), dp)
       rr  = cparams(ind+3)
       !
       cpole = xx + CI*dd
       !
       do i = 1, nx
         !
         Xmat_d1(i,j) =       rr/(xdata(i) -cpole)**2
         Xmat_d2(i,j) = 2.0d0*rr/(xdata(i) -cpole)**3
         !                       
       enddo
       !
     enddo
     
     !
     ! eval    beta = d Func / d pole_j
     !
     ctmp(1:nx,1)=ydata(1:nx,1)-ydata_fit(1:nx,1)
     call mat_mul(beta,Xmat_d1,'C',ctmp,'N',ndim_d1,1,nx)
     !
     beta=-beta

     !
     ! eval alpha
     !
     if (compute_alpha) then
       !
       CALL errore(subname,"Xmat_d2 data structure changed",10)
       !
       call mat_mul(alpha,Xmat_d1,'C',Xmat_d1,'N',ndim_d1,ndim_d1,nx)
       !
       ctmp(1:nx,1)=ydata(1:nx,1)-ydata_fit(1:nx,1)
       !
       do j = 1, ndim_d1
          alpha(j,j) = alpha(j,j) -dot_product(Xmat_d2(:,j),ctmp(:,1))
       enddo
       !
     elseif (compute_alpha_diag) then
       !
       ctmp(1:nx,1)=ydata(1:nx,1)-ydata_fit(1:nx,1)
       !
       do j = 1, ndim_d1
          alpha_diag(j) = -dot_product(Xmat_d2(:,j),ctmp(:,1))
       enddo
       !
     endif
     !
     ! minimization strategy
     !
     select case(trim(minimization_type))
     case("steepest_descent")
       !
       ! rescale beta
       !
       beta = -cost * conjg(beta)
       !
     case("full_hessian")
       !
       !call mat_sv(ndim,1,alpha,beta,ierr)
       call mat_lsd(ndim_d1,ndim_d1,1,alpha,beta,1.0d-6,IERR=ierr)
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
          if ( abs(alpha_diag(j)) <= 1.0d-4 ) then
             beta(j,1) = -cost * conjg(beta(j,1))
          else
             beta(j,1) = -beta(j,1)/alpha_diag(j)
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
       if ( abs(beta(j,1)) > trust )  beta(j,1) = trust * beta(j,1)/abs(beta(j,1))
     enddo
     !
     ! update data II
     !
     cparams0=cparams
     !
     do j = 1, npoles_fit
       !
       ind = 1+ (j-1)*3
       !
       cparams(ind+1)=cparams(ind+1) + real(beta(j,1),dp)
       cparams(ind+2)=cparams(ind+2) +aimag(beta(j,1))
       !
       xx  = real( cparams(ind+1), dp)
       dd  = real( cparams(ind+2), dp)
       rr  = cparams(ind+3)
       !
       cpole = xx + CI*dd
       !
       do i = 1, nx
         Xmat(i,j)    = 1.0d0/(xdata(i)-cpole)
       enddo
       !
     enddo
     if (allow_offset_) Xmat(:,ndim)=1.0

     !
     ! eval residuals
     !
     call mat_mul(ydata_fit,Xmat,'N',Bmat,'N',nx,1,ndim)
     res2 = real(dot_product(ydata(:,1)-ydata_fit(:,1),ydata(:,1)-ydata_fit(:,1)),dp)/(real(nx,dp)*res_scal**2)
     !
     if (res2>res1) then
       ! reject the step
       cparams=cparams0
       ! reduce trust
       do_reduce_trust=.true.
       ! no need to perform the next LSQ
       do_skip_LSQ=.true.
     endif
     !
     if (fitting_verbosity) &
        write(6,"(a,i5, 3f21.9)") "iter res2", iter, res2, trust
     !
     cparams(2+3*npoles_fit) = res2
     !
     ! check convergence
     !
     if (iter>1 .and. abs(res2-res0) < thr) then
        lconv=.true.
        exit iter_loop
     endif
     !
     res0=res2
     !
  enddo iter_loop
  !
  if (.not.allow_offset_) cparams(1)=0.0d0
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
     if (allocated(ctmp))  deallocate(ctmp)
     if (allocated(alpha)) deallocate(alpha)
     if (allocated(alpha_diag)) deallocate(alpha_diag)
     if (allocated(Xmat_d1))    deallocate(Xmat_d1)
     if (allocated(Xmat_d2))    deallocate(Xmat_d2)
     if (allocated(cparams0))   deallocate(cparams0)
     if (allocated(ydata_fit))  deallocate(ydata_fit)
     if (allocated(ydata)) deallocate(ydata)
  end subroutine
  !
end function pole_fitting_LSQ_NL_cmplx_c


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

  call errore(subname,"implementation with allow shift bug done, needs debugging before usage",10)

  trust=0.005
  cost=0.075
  !
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
  allow_offset_=.true.
  if (present(allow_offset)) allow_offset_=allow_offset
  !
  ierr=0
  !
  if (allow_offset_) then
    ndim=npoles_fit+1
  else
    ndim=npoles_fit
  endif
  !
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
     do j = 1, npoles_fit 
       !
       ind = 1+(j-1)*3
       xx  = params(ind+1)
       dd  = params(ind+2)
       !
       do i = 1, nx
         Xmat(i,j)    = dd/((xdata(i)-xx)**2 + dd**2)
       enddo
       !
     enddo
     if (allow_offset_) Xmat(:,ndim)=1.0

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
     params(1)=0.0d0
     if (allow_offset_) params(1)=Bmat(1,ndim)
     !
     Xmat_d2=0.0
     !
     do j = 1, npoles_fit
       !
       ind = 1+(j-1)*3
       xx  = params(ind)
       dd  = params(ind+1)
       rr  = Bmat(j,1)
       params(ind+3)=rr
       !
       jj = 2*j-1
       !
       do i = 1, nx
         xmx0=xdata(i)-xx
         !
         Xmat_d1(i,jj) = rr*2.0d0*dd*(xmx0)/((xmx0)**2 + dd**2)**2
         Xmat_d1(i,jj+1) = rr*((xmx0)**2 -dd**2)/((xmx0)**2 + dd**2)**2
         !
         Xmat_d2(i,jj,1) = rr*2.0d0*dd*( 3.0d0*(xmx0)**2 -dd**2)/( xmx0**2+dd**2)**3
         Xmat_d2(i,jj,2) = rr*2.0d0*xmx0* ((xmx0)**2 -3.0d0*dd**2)/((xmx0)**2 + dd**2)**3
         Xmat_d2(i,jj+1,2) = Xmat_d2(i,j,2)
         Xmat_d2(i,jj+1,1) = -rr*2.0d0*dd/((xmx0)**2 + dd**2)**2 &
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
     do j = 1, npoles_fit
       !
       jj= 2*(j-1)
       ind = 1+ (j-1)*3
       !
       params(ind+1)=params(ind+1)-beta(jj+1,1)
       params(ind+2)=params(ind+2)-beta(jj+2,1)
       !
       xx  = params(ind+1)
       dd  = params(ind+2)
       rr  = params(ind+3)
       !
       do i = 1, nx
         Xmat(i,j)    = dd/((xdata(i)-xx)**2 + dd**2)
       enddo
       !
     enddo
     if (allow_offset_) Xmat(:,ndim)=1.0
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
  function pole_fitting_pade_cmplx(nx,xdata,ydata_,npoles_fit,cparams,allow_offset,even) result(ierr)
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
  complex(dp) :: cparams(*)
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
  complex(dp), allocatable :: roots(:)
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
    allocate(roots(npoles_fit))
  else
    allocate(roots(2*npoles_fit))
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
  call poly_roots(roots,q_pol)
  !
  cparams(1)=0.0d0
  if (allow_offset_) cparams(1)=0.0  ! to be fixed
  !
  if(.not. even_) then
    do j = 1, npoles_fit
      !
      ind = 2+(j-1)*3
      cparams(ind+0)=real(roots(j),dp)
      cparams(ind+1)=aimag(roots(j))
      !
      ctmp(1)=poly_eval(roots(j),p_pol)
      do i = 1, npoles_fit
        if (i==j) cycle
        ctmp(1)=ctmp(1)/(roots(j)-roots(i))
      enddo
      cparams(ind+2)=ctmp(1)
      !
    enddo
  else
    countt=0
    do j = 1, 2*npoles_fit
      !
      if(real(roots(j),dp)>=0) then
        countt = countt + 1
      else
        cycle
      endif
      !
      ind = 2+(countt-1)*3
      cparams(ind+0)=real(roots(j),dp)
      cparams(ind+1)=aimag(roots(j))
      !
      ctmp(1)=poly_eval(roots(j),p_pol)
      do i = 1, 2*npoles_fit
        if (i==j) cycle
        ctmp(1)=ctmp(1)/(roots(j)-roots(i))
      enddo
      cparams(ind+2)=ctmp(1)
      !
      if(countt==npoles_fit) exit
      !
    enddo
  endif
  !
  cparams(2+3*npoles_fit) = res

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
     if (allocated(roots)) deallocate(roots)
     call poly_reset(p_pol)
     call poly_reset(q_pol)
  end subroutine
  !
end function pole_fitting_pade_cmplx

!======================
  function pole_interp_pade_cmplx_dn(nx,xdata,ydata,npoles_fit,dn,cparams,allow_offset) result(ierr)
  !======================
  implicit none
  !
  ! performs interpolation on input data using pade' approx
  ! (equivalent to a sum over poles)
  ! ierr /=0 indicates a numerical problem
  !
  integer     :: nx,npoles_fit,dn
  real(dp)    :: xdata(nx)
  complex(dp) :: ydata(nx)
  complex(dp) :: cparams(*)
  integer     :: ierr
  logical, optional :: allow_offset
  !
  integer   :: ndim
  integer   :: i,j,ind
  logical  :: allow_offset_
  complex(dp) :: res,ctmp(1)
  complex(dp), allocatable :: ydata_fit(:)
  complex(dp), allocatable :: Amat(:,:)
  complex(dp), allocatable :: Bmat(:,:)
  complex(dp), allocatable :: roots(:)
  !
  integer :: p_ndeg, q_ndeg
  type(poly_t) :: p_pol, q_pol

  ierr=0
  !
  allow_offset_=.false.
  if (present(allow_offset)) allow_offset_=allow_offset
  !
  p_ndeg=npoles_fit-dn
  if (allow_offset_) p_ndeg=npoles_fit
  q_ndeg=npoles_fit
  !
  ! the highest order coefficient of q is equal to 1
  ! to get rid of a fictitious degree of freedom
  ! due to the cancellation of a constant coeff in the
  ! p/q ratio
  !
  ndim=p_ndeg+q_ndeg+1
  !
  if (nx/=ndim) then
    ierr=-1
    return
  endif
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(ydata_fit(nx))
  !
  allocate(roots(q_ndeg))
  !
  ! build the linear system
  !
  do j = 1, p_ndeg+1
    do i = 1, nx
       Amat(i,j) = xdata(i)**(j-1)
    enddo
  enddo
  do j = 1, q_ndeg
    do i = 1, nx
       Amat(i,j+p_ndeg+1) = -ydata(i)*xdata(i)**(j-1)
    enddo
  enddo
  do i = 1, nx
    Bmat(i,1)=ydata(i)*xdata(i)**(q_ndeg)
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
  call poly_reset(p_pol,p_ndeg)
  call poly_reset(q_pol,q_ndeg)
  !
  p_pol%c(0:p_ndeg)  =Bmat(1:p_ndeg+1,1)
  q_pol%c(0:q_ndeg-1)=Bmat(p_ndeg+2:p_ndeg+q_ndeg+1,1)
  q_pol%c(q_ndeg)=1.0_dp
  !
  ydata_fit(:) = poly_eval(xdata,nx,p_pol) / poly_eval(xdata,nx,q_pol)
  res = dot_product(ydata(:)-ydata_fit(:),ydata(:)-ydata_fit(:))/real(nx,dp)

  !
  ! reshape data
  !
  call poly_roots(roots,q_pol)
  !
  cparams(1)=0.0d0
  if (allow_offset_) cparams(1)=0.0  ! to be fixed
  !
  do j = 1, npoles_fit
    !
    ind = 2+(j-1)*3
    cparams(ind+0)=real(roots(j),dp)
    cparams(ind+1)=aimag(roots(j))
    !
    ctmp(1)=poly_eval(roots(j),p_pol)
    do i = 1, npoles_fit
      if (i==j) cycle
      ctmp(1)=ctmp(1)/(roots(j)-roots(i))
    enddo
    cparams(ind+2)=ctmp(1)
    !
  enddo
  !
  cparams(2+3*npoles_fit) = res

  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     if (allocated(Amat))  deallocate(Amat)
     if (allocated(Bmat))  deallocate(Bmat)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
     if (allocated(roots)) deallocate(roots)
     call poly_reset(p_pol)
     call poly_reset(q_pol)
  end subroutine
  !
end function pole_interp_pade_cmplx_dn


!======================
  function pole_interp_pade_cmplx(nx,xdata,ydata,npoles_fit,cparams,allow_offset) result(ierr)
  !======================
  implicit none
  !
  ! performs interpolation on input data using pade' approx
  ! (equivalent to a sum over poles)
  ! ierr /=0 indicates a numerical problem
  !
  integer     :: nx,npoles_fit
  real(dp)    :: xdata(nx)
  complex(dp) :: ydata(nx)
  complex(dp) :: cparams(*)
  integer     :: ierr
  logical, optional :: allow_offset
  !
  integer   :: ndim
  integer   :: i,j,ind
  logical  :: allow_offset_
  complex(dp) :: res,ctmp(1)
  complex(dp), allocatable :: ydata_fit(:)
  complex(dp), allocatable :: Amat(:,:)
  complex(dp), allocatable :: Bmat(:,:)
  complex(dp), allocatable :: roots(:)
  !
  integer :: p_ndeg, q_ndeg
  type(poly_t) :: p_pol, q_pol

  ierr=0
  !
  allow_offset_=.false.
  if (present(allow_offset)) allow_offset_=allow_offset
  !
  p_ndeg=npoles_fit-1
  if (allow_offset_) p_ndeg=npoles_fit
  q_ndeg=npoles_fit
  !
  ! the highest order coefficient of q is equal to 1
  ! to get rid of a fictitious degree of freedom
  ! due to the cancellation of a constant coeff in the
  ! p/q ratio
  !
  ndim=p_ndeg+q_ndeg+1
  !
  if (nx/=ndim) then
    ierr=-1
    return
  endif
  !
  allocate(Amat(ndim,ndim))
  allocate(Bmat(ndim,1))
  allocate(ydata_fit(nx))
  !
  allocate(roots(q_ndeg))
  !
  ! build the linear system
  !
  do j = 1, p_ndeg+1
    do i = 1, nx
       Amat(i,j) = xdata(i)**(j-1)
    enddo
  enddo
  do j = 1, q_ndeg
    do i = 1, nx
       Amat(i,j+p_ndeg+1) = -ydata(i)*xdata(i)**(j-1)
    enddo
  enddo
  do i = 1, nx
    Bmat(i,1)=ydata(i)*xdata(i)**(q_ndeg)
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
  call poly_reset(p_pol,p_ndeg)
  call poly_reset(q_pol,q_ndeg)
  !
  p_pol%c(0:p_ndeg)  =Bmat(1:p_ndeg+1,1)
  q_pol%c(0:q_ndeg-1)=Bmat(p_ndeg+2:p_ndeg+q_ndeg+1,1)
  q_pol%c(q_ndeg)=1.0_dp
  !
  ydata_fit(:) = poly_eval(xdata,nx,p_pol) / poly_eval(xdata,nx,q_pol)
  res = dot_product(ydata(:)-ydata_fit(:),ydata(:)-ydata_fit(:))/real(nx,dp)

  !
  ! reshape data
  !
  call poly_roots(roots,q_pol)
  !
  cparams(1)=0.0d0
  if (allow_offset_) cparams(1)=0.0  ! to be fixed
  !
  do j = 1, npoles_fit
    !
    ind = 2+(j-1)*3
    cparams(ind+0)=real(roots(j),dp)
    cparams(ind+1)=aimag(roots(j))
    !
    ctmp(1)=poly_eval(roots(j),p_pol)
    do i = 1, npoles_fit
      if (i==j) cycle
      ctmp(1)=ctmp(1)/(roots(j)-roots(i))
    enddo
    cparams(ind+2)=ctmp(1)
    !
  enddo
  !
  cparams(2+3*npoles_fit) = res

  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     if (allocated(Amat))  deallocate(Amat)
     if (allocated(Bmat))  deallocate(Bmat)
     if (allocated(ydata_fit)) deallocate(ydata_fit)
     if (allocated(roots)) deallocate(roots)
     call poly_reset(p_pol)
     call poly_reset(q_pol)
  end subroutine
  !
end function pole_interp_pade_cmplx

!======================
  function pole_interp_pade_cmplx_noLS(nx,xdata,ydata,npoles_fit,cparams,allow_offset) result(ierr)
  !======================
  implicit none
  !
  ! performs interpolation on input data using pade' approx
  ! it exploits Thiele's iterative algorithm to get a continued fraction
  ! interpolant and then an iterative procedure to compute the
  ! pade approximant. 
  !
  ! The procedure is explained (among many references) in PRB 54, R8285-R8288.
  !
  ! Then it uses polynomial algebra over the polynomial forms A and B of R=A/B to get coefficients
  ! of A and B. After, the SOP poles (roots of B) are found with the companion matrix,
  ! and residues with their definition (usual Heaviside cover-up formula).
  ! (equivalent to a sum over poles)
  !
  ! ierr /=0 indicates a numerical problem
  !
  integer     :: nx,npoles_fit
  complex(dp) :: xdata(nx)
  complex(dp) :: ydata(nx)
  complex(dp) :: cparams(*)
  integer     :: ierr
  logical, optional :: allow_offset
  !
  character(25) :: subname="pole_interp_pade_cmplx_LS"
  integer   :: ndim
  integer   :: i,j,ind
  logical  :: allow_offset_
  complex(dp) :: res,ctmp
  complex(dp), allocatable :: ydata_fit(:)
  complex(dp), allocatable :: aa(:), gg(:,:)
  complex(dp), allocatable :: roots(:)
  !
  integer :: p_ndeg, q_ndeg
  type(poly_t) :: p_pol, q_pol
  type(poly_t) :: p1,p2
  type(poly_t), allocatable :: A_pol(:),B_pol(:)

  ierr=0
  !
  allow_offset_=.false.
  if (present(allow_offset)) allow_offset_=allow_offset
  if (allow_offset_) call errore(subname,"offset not implemented 1",10)
  !
  p_ndeg=npoles_fit-1
  if (allow_offset_) p_ndeg=npoles_fit
  q_ndeg=npoles_fit
  !
  ! the highest order coefficient of q is equal to 1
  ! to get rid of a fictitious degree of freedom
  ! due to the cancellation of a constant coeff in the
  ! p/q ratio
  !
  ndim=p_ndeg+q_ndeg+1
  !
  if (nx/=ndim) then
    ierr=-1
    return
  endif
  !
  allocate(A_pol(0:ndim))
  allocate(B_pol(0:ndim))
  allocate(aa(ndim))
  allocate(gg(ndim,ndim))
  allocate(ydata_fit(nx))
  !
  allocate(roots(q_ndeg))
  !
  ! build gg and aa
  !
  gg=0.0
  gg(:,1) = ydata(:)
  aa(1)=gg(1,1)
  !
  do j = 2, ndim
    !
    do i = j, nx
      gg(i,j) = ( gg(j-1,j-1)-gg(i,j-1) ) / ( (xdata(i)-xdata(j-1)) * gg(i,j-1) )
    enddo
    aa(j) = gg(j,j)
    !
  enddo

  !
  ! evaluate the A_pol and B_pol coeff
  !
  call poly_reset(p1,1)
  !
  call poly_reset(A_pol(0),0)
  call poly_reset(A_pol(1),0)
  call poly_reset(B_pol(0),0)
  call poly_reset(B_pol(1),0)
  !
  A_pol(0)%c=0.0
  A_pol(1)%c=aa(1)
  B_pol(0)%c=1.0
  B_pol(1)%c=1.0
  !
  call poly_reset(p2)
  !
  do j = 2, ndim
    p1%c(0:1) = aa(j)*(/-xdata(j-1),cmplx(1.0_dp,0.0_dp,dp)/)
    call poly_mult(p2,  p1, A_pol(j-2))
    call poly_sum(A_pol(j),   A_pol(j-1), p2 )
    call poly_reset(p2)
  enddo
  !
  do j = 2, ndim
    p1%c(0:1) = aa(j)*(/-xdata(j-1),cmplx(1.0_dp,0.0_dp,dp)/)
    call poly_mult(p2,  p1, B_pol(j-2))
    call poly_sum(B_pol(j), B_pol(j-1), p2 )
    call poly_reset(p2)
  enddo

  !
  ! eval residuals
  !
  call poly_reset(p_pol,p_ndeg)
  call poly_reset(q_pol,q_ndeg)
  !
  p_pol=A_pol(ndim)
  q_pol=B_pol(ndim)
  p_pol%c=p_pol%c/q_pol%c(npoles_fit)
  q_pol%c=q_pol%c/q_pol%c(npoles_fit)
  !
  ydata_fit(:) = poly_eval(xdata,nx,p_pol) / poly_eval(xdata,nx,q_pol)
  res = dot_product(ydata(:)-ydata_fit(:),ydata(:)-ydata_fit(:))/real(nx,dp)

  !
  ! reshape data
  !
  call poly_roots(roots,q_pol)
  !
  cparams(1)=0.0d0
  if (allow_offset_) cparams(1)=0.0  ! to be fixed
  !
  do j = 1, npoles_fit
    !
    ind = 2+(j-1)*3
    cparams(ind+0)=real(roots(j),dp)
    cparams(ind+1)=aimag(roots(j))
    !
    ctmp=poly_eval(roots(j),p_pol)
    do i = 1, npoles_fit
      if (i==j) cycle
      ctmp=ctmp/(roots(j)-roots(i))
    enddo
    cparams(ind+2)=ctmp
    !
  enddo
  !
  cparams(2+3*npoles_fit) = res

  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     do i = 1, size(A_pol)
        call poly_reset(A_pol(i-1))
     enddo
     if (allocated(A_pol)) deallocate(A_pol)
     do i = 1, size(B_pol)
        call poly_reset(B_pol(i-1))
     enddo
     if (allocated(B_pol)) deallocate(B_pol)
     call poly_reset(p1)
     call poly_reset(p2)
     !
     if (allocated(ydata_fit)) deallocate(ydata_fit)
     if (allocated(roots))     deallocate(roots)
     if (allocated(aa))        deallocate(aa)
     if (allocated(gg))        deallocate(gg)
     call poly_reset(p_pol)
     call poly_reset(q_pol)
  end subroutine
  !
end function pole_interp_pade_cmplx_noLS

!======================
  function pole_rational_interp_cmplx(nx,xdata,ydata,npoles_fit,dn,cparams,allow_offset) result(ierr)
  !======================
  implicit none
  !
  ! perform rational interpolation on input data using the A_i or the B_i
  ! algorithm (see below)
  ! (equivalent to a sum over poles)
  ! ierr /=0 indicates a numerical problem
  !
  integer     :: nx,npoles_fit,dn
  real(dp)    :: xdata(nx)
  complex(dp) :: ydata(nx)
  complex(dp) :: cparams(*)
  integer     :: ierr
  logical, optional :: allow_offset
  !
  character(26) :: subname="pole_rational_interp_cmplx"
  integer   :: ndim,ntr,npol
  integer   :: i,j,ind
  logical  :: allow_offset_
  complex(dp) :: res,ctmp(1)
  complex(dp), allocatable :: ydata_fit(:)
  complex(dp), allocatable :: roots(:)
  !
  integer :: p_ndeg, q_ndeg
  type(poly_t) :: p_pol, q_pol
  type(poly_t), allocatable :: A_pol(:),B_pol(:)

  ierr=0
  !
  allow_offset_=.false.
  if (present(allow_offset)) allow_offset_=allow_offset
  if (allow_offset_) call errore(subname,"offset not implemented 1",10)
  !
  p_ndeg=npoles_fit+dn
  ! if (allow_offset_) p_ndeg=npoles_fit
  q_ndeg=npoles_fit
  !
  ndim=p_ndeg+q_ndeg+1
  !
  if (nx/=ndim) then
    ierr=-1
    return
  endif
  !
  if (dn<0 .and. npoles_fit<dn) then
    ierr=-2
    return
  endif
  !
  npol = nx*(nx+1)/2
  !
  allocate(A_pol(npol),B_pol(npol))
  allocate(ydata_fit(nx))
  allocate(roots(q_ndeg))
  !
  if (dn>=0) then
    !
    ntr = dn+1
    ierr = rational_interp_A(ntr,nx,xdata,ydata,A_pol,B_pol)
    !
  else
    !
    ntr = abs(dn)
    ierr = rational_interp_B(ntr,nx,xdata,ydata,A_pol,B_pol)
    ! XXX
    if (ierr/=0) then
      write(6,*) 'ntr=', ntr
      write(6,*) 'ierr=', ierr
      !stop('rational_interp_B error')
    endif
    ! XXX
    !
  endif

  !
  ! eval residues
  !
  call poly_reset(p_pol,p_ndeg)
  call poly_reset(q_pol,q_ndeg)
  !
  p_pol=A_pol(npol)
  q_pol=B_pol(npol)
  p_pol%c=p_pol%c/q_pol%c(q_ndeg)
  q_pol%c=q_pol%c/q_pol%c(q_ndeg)
  !
  ydata_fit(:) = poly_eval(xdata,nx,p_pol) / poly_eval(xdata,nx,q_pol)
  res = dot_product(ydata(:)-ydata_fit(:),ydata(:)-ydata_fit(:))/real(nx,dp)

  !
  ! reshape data
  !
  call poly_roots(roots,q_pol)
  !
  cparams(1)=0.0d0
  if (allow_offset_) cparams(1)=0.0  ! to be fixed
  !
  do j = 1, npoles_fit
    !
    ind = 2+(j-1)*3
    cparams(ind+0)=real(roots(j),dp)
    cparams(ind+1)=aimag(roots(j))
    !
    ctmp(1)=poly_eval(roots(j),p_pol)
    do i = 1, npoles_fit
      if (i==j) cycle
      ctmp(1)=ctmp(1)/(roots(j)-roots(i))
    enddo
    cparams(ind+2)=ctmp(1)
    !
  enddo
  !
  cparams(2+3*npoles_fit) = res

  !
  ! cleanup
  !
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     do i = 1, size(A_pol)
        call poly_reset(A_pol(i))
     enddo
     if (allocated(A_pol)) deallocate(A_pol)
     do i = 1, size(B_pol)
        call poly_reset(B_pol(i))
     enddo
     if (allocated(B_pol)) deallocate(B_pol)
     !
     if (allocated(ydata_fit)) deallocate(ydata_fit)
     if (allocated(roots)) deallocate(roots)
     call poly_reset(p_pol)
     call poly_reset(q_pol)
  end subroutine
  !
end function pole_rational_interp_cmplx 

!======================
  function rational_interp_A(ntr,nx,xdata,ydata,p_pol,q_pol) result(ierr)
  !======================
  implicit none
  !
  ! performs rational interpolation on input data using the A_i algorithm in
  ! F.M. Larkin, The Computer Journal 10, Issue 2, pp. 178-187
  ! ierr /=0 indicates a numerical problem
  !
  integer      :: ntr,nx
  real(dp)     :: xdata(nx)
  complex(dp)  :: ydata(nx)
  type(poly_t) :: p_pol(*), q_pol(*)
  integer      :: ierr
  !
  !character(17) :: subname="rational_interp_A"
  !
  integer   :: p_ndeg, q_ndeg, r
  integer   :: j, k, ir, ind, ind_km1, ind_jp1km1
  complex(dp)  :: alpha, beta

  ierr = 0
  !
  if (ntr>nx-1) then
    ierr=-1
    return
  endif

  !
  ! prepare polynomials
  !
  do j=1, nx
    !
    ind = j
    !
    call poly_reset(p_pol(ind),0)
    call poly_reset(q_pol(ind),0)
    !
    p_pol(ind)%c(0) = ydata(j)
    q_pol(ind)%c(0) = 1.0d0
    !
  enddo

  !
  ! execute ntr times the triangular rule
  !
  p_ndeg = 0
  q_ndeg = 0
  !
  do k=1, ntr
    !
    p_ndeg = p_ndeg+1
    !
    do j=1, nx-k
      !
      ind        = k*nx - k*(k-1)/2 + j
      ind_km1    = ind - (nx-(k-1))
      ind_jp1km1 = ind_km1 + 1
      !
      alpha = q_pol(ind_km1)%c(0)    / (q_pol(ind_km1)%c(0)+q_pol(ind_jp1km1)%c(0))
      beta  = q_pol(ind_jp1km1)%c(0) / (q_pol(ind_km1)%c(0)+q_pol(ind_jp1km1)%c(0))
      !
      call poly_reset(p_pol(ind),p_ndeg)
      call poly_reset(q_pol(ind),0)
      !
      p_pol(ind)%c(0) = alpha*(-xdata(j)*p_pol(ind_jp1km1)%c(0)) + &
                        beta *(xdata(j+k)*p_pol(ind_km1)%c(0))
      do ir=1, p_ndeg
        !
        ! p_pol(j,k-1)%c(p_ndeg) = 0 in principle, but is not allocated
        !
        if (ir==p_ndeg) then
          p_pol(ind)%c(ir) = alpha*(p_pol(ind_jp1km1)%c(ir-1)) + &
                            beta *(-p_pol(ind_km1)%c(ir-1))
        else
          p_pol(ind)%c(ir) = alpha*(p_pol(ind_jp1km1)%c(ir-1)-xdata(j)*p_pol(ind_jp1km1)%c(ir)) + &
                            beta *(xdata(j+k)*p_pol(ind_km1)%c(ir)-p_pol(ind_km1)%c(ir-1))
        endif
        !
      enddo
      !
      q_pol(ind)%c(0) = alpha*(-xdata(j)*q_pol(ind_jp1km1)%c(0)) + &
                        beta *(xdata(j+k)*q_pol(ind_km1)%c(0))
      !
    enddo
    !
  enddo

  !
  ! execute the rhombus rule for the remaining steps
  !
  do k=ntr+1, nx-1
    !
    if (mod(k+ntr,2) /= 0) then 
      !
      r = (k+ntr-1)/2
      q_ndeg = q_ndeg+1
      !
    else
      !
      r = (k-ntr)/2
      p_ndeg = p_ndeg+1
      !
    endif
    !
    do j=1, nx-k
      !
      ind        = k*nx - k*(k-1)/2 + j
      ind_km1    = ind  - (nx-(k-1))
      ind_jp1km1 = ind_km1 + 1
      !
      if (mod(k+ntr,2) /= 0) then 
        !
        alpha = p_pol(ind_km1)%c(r)    / (p_pol(ind_km1)%c(r)+p_pol(ind_jp1km1)%c(r))
        beta  = p_pol(ind_jp1km1)%c(r) / (p_pol(ind_km1)%c(r)+p_pol(ind_jp1km1)%c(r))
        !
      else
        !
        alpha = q_pol(ind_km1)%c(r)    / (q_pol(ind_km1)%c(r)+q_pol(ind_jp1km1)%c(r))
        beta  = q_pol(ind_jp1km1)%c(r) / (q_pol(ind_km1)%c(r)+q_pol(ind_jp1km1)%c(r))
        !
      endif
      !
      call poly_reset(p_pol(ind),p_ndeg)
      call poly_reset(q_pol(ind),q_ndeg)
      !
      p_pol(ind)%c(0) = alpha*(-xdata(j)*p_pol(ind_jp1km1)%c(0)) + &
                        beta *(xdata(j+k)*p_pol(ind_km1)%c(0))
      do ir=1, p_ndeg
        !
        if (p_pol(ind_km1)%order<p_ndeg .and. ir==p_ndeg) then
          p_pol(ind)%c(ir) = alpha*(p_pol(ind_jp1km1)%c(ir-1)) + &
                            beta *(-p_pol(ind_km1)%c(ir-1))
        else
          p_pol(ind)%c(ir) = alpha*(p_pol(ind_jp1km1)%c(ir-1)-xdata(j)*p_pol(ind_jp1km1)%c(ir)) + &
                            beta *(xdata(j+k)*p_pol(ind_km1)%c(ir)-p_pol(ind_km1)%c(ir-1))
        endif
        !
      enddo
      !
      q_pol(ind)%c(0) = alpha*(-xdata(j)*q_pol(ind_jp1km1)%c(0)) + &
                        beta *(xdata(j+k)*q_pol(ind_km1)%c(0))
      do ir=1, q_ndeg
        !
        if (q_pol(ind_km1)%order<q_ndeg .and. ir==q_ndeg) then
          q_pol(ind)%c(ir) = alpha*(q_pol(ind_jp1km1)%c(ir-1)) + &
                            beta *(-q_pol(ind_km1)%c(ir-1))
        else
          q_pol(ind)%c(ir) = alpha*(q_pol(ind_jp1km1)%c(ir-1)-xdata(j)*q_pol(ind_jp1km1)%c(ir)) + &
                            beta *(xdata(j+k)*q_pol(ind_km1)%c(ir)-q_pol(ind_km1)%c(ir-1))
        endif
        !
      enddo
      !
    enddo
    !
  enddo

  end function rational_interp_A

!======================
  function rational_interp_B(ntr,nx,xdata,ydata,p_pol,q_pol) result(ierr)
  !======================
  implicit none
  !
  ! performs rational interpolation on input data using the B_i algorithm in
  ! F.M. Larkin, The Computer Journal 10, Issue 2, pp. 178-187
  ! ierr /=0 indicates a numerical problem
  !
  integer      :: ntr,nx
  real(dp)     :: xdata(nx)
  complex(dp)  :: ydata(nx)
  type(poly_t) :: p_pol(*), q_pol(*)
  integer      :: ierr
  !
  !character(17) :: subname="rational_interp_A"
  !
  integer      :: p_ndeg, q_ndeg, r
  integer      :: j, k, ir, ind, ind_km1, ind_jp1km1
  complex(dp)  :: alpha, beta
  
  ierr=0
  !
  if (ntr>nx-1) then
    ierr=-1
    return
  endif

  !
  ! prepare polynomials
  !
  do j=1, nx
    !
    ind = j
    !
    call poly_reset(p_pol(ind),0)
    call poly_reset(q_pol(ind),0)
    !
    p_pol(ind)%c(0) = ydata(j)
    q_pol(ind)%c(0) = 1.0d0
    !
  enddo

  !
  ! execute ntr times the triangular rule
  !
  p_ndeg = 0
  q_ndeg = 0
  !
  do k=1, ntr
    !
    q_ndeg = q_ndeg+1
    !
    do j=1, nx-k
      !
      ind        = k*nx - k*(k-1)/2 + j
      ind_km1    = ind - (nx-(k-1))
      ind_jp1km1 = ind_km1 + 1
      !
      alpha = p_pol(ind_km1)%c(0)    / (p_pol(ind_km1)%c(0)+p_pol(ind_jp1km1)%c(0))
      beta  = p_pol(ind_jp1km1)%c(0) / (p_pol(ind_km1)%c(0)+p_pol(ind_jp1km1)%c(0))
      !
      call poly_reset(p_pol(ind),0)
      call poly_reset(q_pol(ind),q_ndeg)
      !
      p_pol(ind)%c(0) = alpha*(-xdata(j)*p_pol(ind_jp1km1)%c(0)) + &
                        beta *(xdata(j+k)*p_pol(ind_km1)%c(0))
      !
      q_pol(ind)%c(0) = alpha*(-xdata(j)*q_pol(ind_jp1km1)%c(0)) + &
                        beta *(xdata(j+k)*q_pol(ind_km1)%c(0))
      do ir=1, q_ndeg
        !
        ! q_pol(j,k-1)%c(q_ndeg) = 0 in principle, but is not allocated
        !
        if (ir==q_ndeg) then
          q_pol(ind)%c(ir) = alpha*(q_pol(ind_jp1km1)%c(ir-1)) + &
                            beta *(-q_pol(ind_km1)%c(ir-1))
        else
          q_pol(ind)%c(ir) = alpha*(q_pol(ind_jp1km1)%c(ir-1)-xdata(j)*q_pol(ind_jp1km1)%c(ir)) + &
                            beta *(xdata(j+k)*q_pol(ind_km1)%c(ir)-q_pol(ind_km1)%c(ir-1))
        endif
        !
      enddo
      !

      !
    enddo
    !
  enddo

  !
  ! execute the rhombus rule for the remaining steps
  !
  do k=ntr+1, nx-1
    !
    if (mod(k+ntr,2) == 0) then 
      !
      r = (k-ntr)/2
      q_ndeg = q_ndeg+1
      !
    else
      !
      r = (k+ntr-1)/2
      p_ndeg = p_ndeg+1
      !
    endif
    !
    do j=1, nx-k
      !
      ind        = k*nx - k*(k-1)/2 + j
      ind_km1    = ind  - (nx-(k-1))
      ind_jp1km1 = ind_km1 + 1
      !
      if (mod(k+ntr,2) == 0) then 
        !
        alpha = p_pol(ind_km1)%c(r)    / (p_pol(ind_km1)%c(r)+p_pol(ind_jp1km1)%c(r))
        beta  = p_pol(ind_jp1km1)%c(r) / (p_pol(ind_km1)%c(r)+p_pol(ind_jp1km1)%c(r))
        !
      else
        !
        alpha = q_pol(ind_km1)%c(r)    / (q_pol(ind_km1)%c(r)+q_pol(ind_jp1km1)%c(r))
        beta  = q_pol(ind_jp1km1)%c(r) / (q_pol(ind_km1)%c(r)+q_pol(ind_jp1km1)%c(r))
        !
      endif
      !
      call poly_reset(p_pol(ind),p_ndeg)
      call poly_reset(q_pol(ind),q_ndeg)
      !
      p_pol(ind)%c(0) = alpha*(-xdata(j)*p_pol(ind_jp1km1)%c(0)) + &
                        beta *(xdata(j+k)*p_pol(ind_km1)%c(0))
      do ir=1, p_ndeg
        !
        if (p_pol(ind_km1)%order<p_ndeg .and. ir==p_ndeg) then
          p_pol(ind)%c(ir) = alpha*(p_pol(ind_jp1km1)%c(ir-1)) + &
                            beta *(-p_pol(ind_km1)%c(ir-1))
        else
          p_pol(ind)%c(ir) = alpha*(p_pol(ind_jp1km1)%c(ir-1)-xdata(j)*p_pol(ind_jp1km1)%c(ir)) + &
                            beta *(xdata(j+k)*p_pol(ind_km1)%c(ir)-p_pol(ind_km1)%c(ir-1))
        endif
        !
      enddo
      !
      q_pol(ind)%c(0) = alpha*(-xdata(j)*q_pol(ind_jp1km1)%c(0)) + &
                        beta *(xdata(j+k)*q_pol(ind_km1)%c(0))
      do ir=1, q_ndeg
        !
        if (q_pol(ind_km1)%order<q_ndeg .and. ir==q_ndeg) then
          q_pol(ind)%c(ir) = alpha*(q_pol(ind_jp1km1)%c(ir-1)) + &
                            beta *(-q_pol(ind_km1)%c(ir-1))
        else
          q_pol(ind)%c(ir) = alpha*(q_pol(ind_jp1km1)%c(ir-1)-xdata(j)*q_pol(ind_jp1km1)%c(ir)) + &
                            beta *(xdata(j+k)*q_pol(ind_km1)%c(ir)-q_pol(ind_km1)%c(ir-1))
        endif
        !
      enddo
      !
    enddo
    !
  enddo

  end function rational_interp_B


end module fitting_base_m

