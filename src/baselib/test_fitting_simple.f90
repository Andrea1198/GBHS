
!==============
program test
!==============
use kinds
use constants
use fitting_m
implicit none
  !
  integer,     parameter :: stdout=6
  !
  integer :: ierr,narg
  integer :: i,j,ip
  integer :: nx, nx_aux, npoles,npoles_fit
  real(dp):: x, yinf,xmin,xmax,res
  complex(dp)   :: pole
  character(64) :: str
  real(dp),     allocatable :: xgrid_aux(:)
  real(dp),     allocatable :: xgrid(:), func_im(:)
  real(dp),     allocatable :: params(:)
  complex(dp),  allocatable :: params_cmplx(:)
  integer,      allocatable :: ipeaks(:)
  complex(dp),  allocatable :: poles(:), func(:), func_fit(:)
  complex(dp),  allocatable :: func_aux(:)
  real(dp),     allocatable :: pole_res(:)
  logical                       :: even
  !
  character(64) :: fitting_type
  character(64) :: subname="fitting_test"
  !
  integer, external :: python_initialize
  integer, external :: python_finalize

  !
  ! init python
  !
  ierr = python_initialize()
  !
  ! get cmd line input
  !
  narg=command_argument_count()
  if (narg<1) then
    write(stdout,"(' USAGE: ./test_fitting.x  npoles_fit [[fitting_type],[even]]')")
    write(stdout,"('     fitting_type = least-squares             | LSQ | ')")
    write(stdout,"('                  | least-squares-cmplx       | LSQ-cmplx | ')")
#ifdef __PYTHON
    write(stdout,"('                  | non-linear-py             | NL-PY | ')")
#endif
    write(stdout,"('                  | non-linear-fortran        | NL-F  | ')")
    write(stdout,"('                  | non-linear-fortran-cmplx  | NL-cmplx | ')")
    write(stdout,"('                  | pade                      | PA | ')")
    write(stdout,"('                  | pade-cmplx                | PA-cmplx | ')")
    write(stdout,"('                  | pade-interp-cmplx          | PA-INT-cmplx | ')")
    write(stdout,"('             even = logical                           ')")
    stop
  endif
  !
  call get_command_argument(1,str)
  read(str,*) npoles_fit
  !
  fitting_type="pade-cmplx"
#ifdef __PYTHON
  fitting_type="non-linear-py"
#endif
  if (narg>=2) call get_command_argument(2,fitting_type)
  !
  even=.false.
  !
  if (narg>=3) then
    call get_command_argument(3,str)
    read(str,*) even
  endif
  !  
  xmin = -10.0
  xmax =  10.0
  nx   = 1000
  yinf = 0.0
  !
  npoles = 4
  allocate(poles(npoles))
  allocate(pole_res(npoles))
  !
  poles(1) = -1.5 + CI * 0.10
  poles(2) = -1.0 + CI * 0.05
  poles(3) =  2.0 - CI * 0.05
  poles(4) =  2.8 - CI * 0.15
  !
  pole_res(1) = 0.5
  pole_res(2) = 1.5
  pole_res(3) = 1.7
  pole_res(4) = 0.3
  !
  if(even) then
    poles(1) = -1.5 + CI * 0.10
    poles(2) = -1.0 + CI * 0.05
    poles(3) =  1.0 - CI * 0.05
    poles(4) =  1.5 - CI * 0.10
    !
    pole_res(1) = -0.5
    pole_res(2) = -1.5
    pole_res(3) = 1.5
    pole_res(4) = 0.5
  endif
  
  !
  allocate(xgrid(nx), func(nx), func_im(nx))
  !
  do i = 1, nx
    xgrid(i) = xmin + real(i-1,dp)*(xmax-xmin)/real(nx-1,dp)
    x=xgrid(i)
    !
    func(i)=yinf
    do ip = 1, npoles
      func(i) = func(i) + pole_res(ip)/(x-poles(ip))
    enddo
    func_im(i) = aimag(func(i))
    !
  enddo

  !
  ! get pole position
  !
  allocate(ipeaks(npoles_fit))
  call fitting_get_peaks(nx,xgrid,abs(aimag(func)),npoles_fit,ipeaks)
  !
  if(fitting_type=="pade-cmplx") ipeaks=1
  
  !
  ! reporting
  !
  write(stdout,"('Fitting type: ',a)") trim(fitting_type)
  write(stdout,"('Fitting even: ',l4)") even
  write(stdout,"('Npoles ',i5)")npoles
  do i = 1, npoles
    write(stdout,"('  p[',i2,'] = ',3f15.9)") i, poles(i), pole_res(i) 
  enddo
  !
  write(stdout,"(/,'ipeaks ',i5)")npoles_fit
  do i = 1, npoles_fit
    if (ipeaks(i)==0) ipeaks(i)=1 
    write(stdout,"('ipk[',i2,'] = ',i5)") i, ipeaks(i) 
  enddo
  !
  write(stdout,"(/,'Trials ',i5)")npoles_fit
  do i = 1, npoles_fit
    write(stdout,"('  p[',i2,'] = ',3f15.9)") i, xgrid(ipeaks(i)),1.0d0/aimag(func(ipeaks(i)))
  enddo

  ! 
  ! fitting 
  ! 
  allocate(func_fit(nx),params((3*npoles_fit)+2))
  allocate(params_cmplx((3*npoles_fit)+2))
  !
  ! initial guess
  !
  params(:)=0.0
  params(1)=0.0
  !
  do i = 1, npoles_fit
    j=1+(i-1)*3
    params(j+1)=xgrid(ipeaks(i))
    params(j+2)=-0.05*sign(1.0d0,xgrid(ipeaks(i)))
    params(j+3)=1.0
  enddo
  write(stdout,"(/,'Fortran: p0 ',i5)") size(params)
  write(stdout,"(6f15.9)") params(:)
  write(stdout,"()") 
  
  !
  ! call to fitting wrappers
  !
  select case (trim(fitting_type) )
  case ("NL-PY","non-linear-py")
     !
     ierr = pole_fitting_pywrap(nx,xgrid,func_im,npoles_fit,params)
     !
  case ("LSQ","least-squares")
     !
     ierr = pole_fitting_LSQ(nx,xgrid,func_im,npoles_fit,params)
     !
  case ("LSQ-cmplx","least-squares-cmplx")
     !
     params_cmplx=params
     ierr = pole_fitting_LSQ_cmplx(nx,cmplx(xgrid,0._dp,dp),func,npoles_fit,params_cmplx)
     params=real(params_cmplx,dp)
     !
  case ("NL-F","non-linear")
     !
     ierr = pole_fitting_LSQ_NL(nx,xgrid,func_im,npoles_fit,params,thr=1.0d-8,niterx=100)
     !
  case ("NL-F-cmplx","NL-cmplx","non-linear-cmplx")
     !
     params_cmplx=params
     ierr = pole_fitting_LSQ_NL_cmplx_c(nx,cmplx(xgrid,0.0_dp,dp),func,npoles_fit,params_cmplx,thr=1.0d-8,niterx=100)
     params=real(params_cmplx,dp)
     !
!  case ("PA","pade")
!     !
!     ierr = pole_fitting_pade(nx,xgrid,func_im,npoles_fit,params)
     !
  case ("PA-cmplx","pade-cmplx")
     !
     params_cmplx=params
     ierr = pole_fitting_pade_cmplx(nx,xgrid,func,npoles_fit,params_cmplx,EVEN=even)
     params=real(params_cmplx,dp)
     !
  case ("PA-INT-cmplx","pade-interp-cmplx")
     !
     if (even) call errore("test_fitting_simple","EVEN=.true. not implemented",1)
     !
     params_cmplx=params
     nx_aux=npoles_fit*2
     allocate(xgrid_aux(nx_aux),func_aux(nx_aux))
     !
     do i = 1, nx_aux
       j = 1+(i-1)*int(nx/nx_aux)
       xgrid_aux(i) = xgrid(j)
       func_aux(i) = func(j)
     enddo
     !
     ierr = pole_interp_pade_cmplx_noLS(nx_aux,cmplx(xgrid_aux,0._dp,dp),func_aux,npoles_fit,params_cmplx)
     !
     deallocate(xgrid_aux,func_aux)
     params=real(params_cmplx,dp)
     !
  case default
     call errore(subname,"invalid fitting_type: "//trim(fitting_type),1)
  end select
  !
  if (ierr/=0) then
     call errore(subname,"fitting procedure failed",10)
  endif

  !
  ! raw result report
  !
  write(stdout,"('Fortran FIT: ')")
  write(stdout,"(3x,'   yinf = ',3f15.9)") params(1)
  do i = 2, size(params)-1, 3
     write(stdout,"(3x,'  p[',i3,'] = ',3f15.9)") (i-1)/3+1, params(i:i+2)
  enddo
  write(stdout,"(3x,'    ERR = ',3f15.9)") params(size(params))
  write(stdout,"()")

  !
  ! build the fitted function
  !
  do i = 1, nx
    !
    x=xgrid(i)
    func_fit(i)=params(1)
    !
    do ip = 1, npoles_fit
      pole=params(1+3*(ip-1)+1)+ CI*params(1+3*(ip-1)+2)
      res=params(1+3*(ip-1)+3)
      !
      func_fit(i) = func_fit(i) + res/(x-pole)
    enddo
    !
  enddo
  !
  ! dump data
  !
  open(10, file="data_fit.dat")
    write(10,"('# xgrid        func         func_fit  ')")
    do i = 1, nx
      write(10,"(5f15.9)") xgrid(i), func(i), func_fit(i)
    enddo
  close(10)

  !
  ! cleanup
  !
  deallocate(xgrid,func,func_im)
  deallocate(func_fit, params)
  deallocate(params_cmplx)
  deallocate(ipeaks)
  !
  ! shutdown
  !
  ierr = python_finalize()
  !
end program test

