
!==============
program test
!==============
use kinds
use constants
use fitting_m
implicit none
  !
  integer, parameter :: stdout=6
  !
  integer       :: ierr,narg
  integer       :: i,j,ip
  integer       :: nx, niterx
  integer       :: npoles, npoles_fit, nsampl_fit
  integer       :: grid_range(2)
  real(dp)      :: x, yinf, xmin,xmax
  real(dp)      :: xmin_sampl,xmax_sampl
  real(dp)      :: gamma_scale
  complex(dp)   :: pole,res
  character(64) :: str
  real(dp),     allocatable :: xgrid(:)
  complex(dp),  allocatable :: params_cmplx(:)
  complex(dp),  allocatable :: poles(:), func(:)
  complex(dp),  allocatable :: func_fit(:)
  integer,      allocatable :: isample_fit(:)
  complex(dp),  allocatable :: sample_fit(:)
  integer,      allocatable :: itrial_fit(:)
  complex(dp),  allocatable :: trial_fit(:)
  complex(dp),  allocatable :: func_aux(:)
  real(dp),     allocatable :: pole_res(:)
  !
  integer       :: npopulation
  integer       :: niter_gen
  real(dp)      :: crossover_frac,mutation_frac,select_frac
  real(dp)      :: rand_radius_re,rand_radius_im
  !
  complex(DP), external :: gaussian_pole
  !
  character(64) :: fitting_type
  character(64) :: subname="fitting_test"

  !
  ! get cmd line input
  !
  narg=command_argument_count()
  if (narg<1) then
    write(stdout,"(' USAGE: ./test_fitting_gau.x  npoles_fit [[fitting_type]')")
    write(stdout,"('     fitting_type = ')")
    write(stdout,"('                  | non-linear-cmplx  | NL-cmplx   | NL  ')")
    write(stdout,"('                  | pade-cmplx        | PA-cmplx   | pade | PA ')")
    write(stdout,"('                  | genetic-cmplx     | gen-cmplx  | genetic ')")
    stop
  endif
  !
  call get_command_argument(1,str)
  read(str,*) npoles_fit
  !
  fitting_type="non-linear-cmplx"
  if (narg>=2) call get_command_argument(2,fitting_type)
  !
  niterx = 500
  !
  xmin = -10.0
  xmax =  10.0
  nx   = 1000
  yinf = 0.0
  xmin_sampl = -5.0
  xmax_sampl =  5.0
  !
  npoles = 4
  allocate(poles(npoles))
  allocate(pole_res(npoles))
  !
  poles(1) = -1.5 + CI * 0.20
  poles(2) = -1.0 + CI * 0.10
  poles(3) =  2.0 - CI * 0.10
  poles(4) =  2.8 - CI * 0.30
  !
  pole_res(1) = 0.5
  pole_res(2) = 1.5
  pole_res(3) = 1.7
  pole_res(4) = 0.3
  !
  allocate(xgrid(nx), func(nx))
  !
  do i = 1, nx
    !
    xgrid(i) = xmin + real(i-1,dp)*(xmax-xmin)/real(nx-1,dp)
    x=xgrid(i)
    !
    func(i)=yinf
    do ip = 1, npoles
      func(i) = func(i) + pole_res(ip)* gaussian_pole(x,real(poles(ip),DP),aimag(poles(ip)))
    enddo
    !
  enddo

  !
  ! get trial poles and sampling points
  !
  nsampl_fit=2*npoles_fit
  !
  allocate(isample_fit(nsampl_fit))
  allocate(sample_fit(nsampl_fit))
  allocate(itrial_fit(npoles_fit))
  allocate(trial_fit(npoles_fit))
  !
  grid_range(1)=nx/4+1
  grid_range(2)=3*nx/4
  !
  call fitting_get_trial(nx,xgrid,abs(func),grid_range,npoles_fit,itrial_fit,ierr)
  if (ierr/=0) call errore("test_gau","error while searching for sampling points",1)
  !
  if(fitting_type=="pade-cmplx") isample_fit=1
  !
  ! use the distance among fitting points to guess the peak broadening
  !
  gamma_scale=0.1_DP
  !
  do i = 2, npoles_fit-1
    trial_fit(i)=cmplx( xgrid(itrial_fit(i)), &
                        gamma_scale*(xgrid(itrial_fit(i+1))-xgrid(itrial_fit(i-1)))/2.0_DP, kind=DP)
  enddo
  !
  ! handle extrema 
  ! 
  i=1
  trial_fit(i)=cmplx( xgrid(itrial_fit(i)), aimag(trial_fit(i+1)), kind=DP)
  i=npoles_fit
  trial_fit(i)=cmplx( xgrid(itrial_fit(i)), aimag(trial_fit(i-1)), kind=DP)
  !
  do i = 1, npoles_fit
    if ( real(trial_fit(i))>0 .and. aimag(trial_fit(i))>0 ) trial_fit(i)=conjg(trial_fit(i))
    if ( real(trial_fit(i))<0 .and. aimag(trial_fit(i))<0 ) trial_fit(i)=conjg(trial_fit(i))
  enddo

  !
  ! sampling
  !
  do j = 1, nsampl_fit
    sample_fit(j) = xmin_sampl + (j-1)*(xmax_sampl-xmin_sampl)/real(nsampl_fit-1,DP)
    !
    call locate(xgrid,nx,real(sample_fit(j),dp),isample_fit(j))
    sample_fit(j)=xgrid(isample_fit(j))
  enddo
  ! 
  allocate(func_aux(nsampl_fit))
  !
  do i = 1, nsampl_fit
    func_aux(i) = func(isample_fit(i))
  enddo
  
  !
  ! reporting
  !
  write(stdout,"('Fitting type: ',a)") trim(fitting_type)
  write(stdout,"('Npoles ',i5)")npoles
  do i = 1, npoles
    write(stdout,"('  p[',i2,'] = ',4f15.9)") i, poles(i), pole_res(i) 
  enddo
  !
  write(stdout,"(/,'Trials ',i5)")npoles_fit
  do i = 1, npoles_fit
    write(stdout,"('  p[',i2,'] = ',3f15.9)") i, trial_fit(i)
  enddo

  ! 
  ! fitting 
  ! 
  allocate(func_fit(nx))
  allocate(params_cmplx((3*npoles_fit)+2))
  !
  ! initial guess
  !
  params_cmplx(:)=0.0
  params_cmplx(1)=0.0
  !
  do i = 1, npoles_fit
    j=1+(i-1)*3
    params_cmplx(j+1)=real(trial_fit(i),DP)
    params_cmplx(j+2)=aimag(trial_fit(i))
    params_cmplx(j+3)=1.0
  enddo
  write(stdout,"(/,'Fortran: p0 ',i5)") size(params_cmplx)
  write(stdout,"(6f15.9)") params_cmplx(:)
  write(stdout,"()") 
  
  !
  ! call to fitting wrappers
  !
  ierr=0
  select case (trim(fitting_type) )
     !
  case ("NL-cmplx","non-linear-cmplx", "NL")
     !
     ierr = pole_fitting_LSQ_NL_cmplx_c(nsampl_fit,sample_fit,func_aux,npoles_fit,params_cmplx,&
                                        thr=1.0d-10,niterx=niterx,allow_offset=.false.)
     !
  case ("genetic-cmplx", "gen-cmplx", "genetic")
     !
     npopulation=400
     niter_gen=200
     crossover_frac=0.9
     mutation_frac=0.7
     select_frac=0.20
     rand_radius_re=1.5
     rand_radius_im=0.2
     !
     ierr = pole_fitting_genetic_cmplx(npopulation,niter_gen,crossover_frac,mutation_frac,select_frac,&
                                       rand_radius_re,rand_radius_im,nsampl_fit,sample_fit,func_aux,npoles_fit,params_cmplx,&
                                       thr=1.0d-10,niterx=niterx,allow_offset=.false.,efermi=0.0_DP)
     !
     npopulation=4000
     niter_gen=20
     crossover_frac=0.9
     mutation_frac=0.7
     select_frac=0.20
     rand_radius_re=0.4
     rand_radius_im=0.1
     !
     ierr = pole_fitting_genetic_cmplx(npopulation,niter_gen,crossover_frac,mutation_frac,select_frac,&
                                       rand_radius_re,rand_radius_im,nsampl_fit,sample_fit,func_aux,npoles_fit,params_cmplx,&
                                       thr=1.0d-10,niterx=niterx,allow_offset=.false.,efermi=0.0_DP)
     !
  case ("PA-cmplx","pade-cmplx","pade", "PA")
     !
     ierr = pole_interp_pade_cmplx_noLS(nsampl_fit,sample_fit,func_aux,npoles_fit,params_cmplx)
     !
  case default
     call errore(subname,"invalid fitting_type: "//trim(fitting_type),1)
  end select
  !
  if (ierr>0) then
     call errore(subname,"fitting procedure failed",10)
  elseif (ierr< 0) then
     call warning(subname,"fitting procedure: max number of iterations reached")
  endif

  !
  ! raw result report
  !
  write(stdout,"('Fortran FIT: ')")
  write(stdout,"(3x,'   yinf = ',3f15.9)") params_cmplx(1)
  do i = 2, size(params_cmplx)-1, 3
     write(stdout,"(3x,'  p[',i3,'] = ',2f15.9,3x,2f15.9)") (i-1)/3+1, real(params_cmplx(i:i+1)), params_cmplx(i+2)
  enddo
  write(stdout,"(3x,'    ERR = ',3f15.9)") real( params_cmplx(size(params_cmplx)) )
  write(stdout,"()")

  !
  ! build the fitted function
  !
  do i = 1, nx
    !
    x=xgrid(i)
    func_fit(i)=params_cmplx(1)
    !
    do ip = 1, npoles_fit
      pole=params_cmplx(1+3*(ip-1)+1)+ CI*params_cmplx(1+3*(ip-1)+2)
      res=params_cmplx(1+3*(ip-1)+3)
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
  open(10, file="data_sampling.dat")
    write(10,"('# xpoint    func     ')")
    do i = 1, nsampl_fit
      write(10,"(5f15.9)") sample_fit(i), func_aux(i)
    enddo
  close(10)

  !
  ! cleanup
  !
  deallocate(xgrid,func)
  deallocate(func_fit)
  deallocate(func_aux)
  deallocate(params_cmplx)
  deallocate(isample_fit)
  deallocate(sample_fit)
  deallocate(itrial_fit)
  deallocate(trial_fit)
  !
end program test

