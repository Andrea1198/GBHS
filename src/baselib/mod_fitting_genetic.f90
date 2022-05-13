!======================
module fitting_genetic_m
  !======================
  use kinds
  use constants
  use fitting_base_m,  only: pole_fitting_LSQ_NL_cmplx_c
  implicit none
  private

  public :: pole_fitting_genetic_cmplx

contains

!======================
  function pole_fitting_genetic_cmplx(npopulation,niter_gen,crossover_frac,mutation_frac,select_frac,radius_re,radius_im,&
                                      nx,xdata,ydata_,npoles_fit,cparams,thr,niterx,allow_offset,efermi) result(ierr)
  !======================
  !
  ! implements a poorman genetic algorithm to perform a global optimization,
  ! based on the pole_fitting_LSQ_NL_cmplx routine for convex optimization
  !
  use random_numbers_module, only : randy
  implicit none
  !
  integer, intent(in) :: npopulation
  integer, intent(in) :: niter_gen
  real(dp),intent(in) :: crossover_frac
  real(dp),intent(in) :: mutation_frac
  real(dp),intent(in) :: select_frac
  real(dp),intent(in) :: radius_re,radius_im
  !
  integer     :: nx,npoles_fit
  complex(dp) :: xdata(nx)
  complex(dp) :: ydata_(nx)
  complex(dp) :: cparams(2+3*npoles_fit)
  real(dp)    :: thr
  integer     :: niterx
  integer     :: ierr
  logical,  optional :: allow_offset
  real(dp), optional :: efermi
  !
  ! local vars
  !
  !character(26) :: subname='pole_fitting_genetic_cmplx'
  !
  integer     :: nparams
  integer     :: ngen_in_memory
  real(dp)    :: rand_radius_re, rand_radius_im
  complex(dp), allocatable :: cparams_gen(:,:,:)
  real(dp),    allocatable :: fitness_best(:)
  real(dp),    allocatable :: fitness_gen(:)
  !
  real(dp)    :: efermi_, arg
  logical     :: allow_offset_
  logical     :: skip_mutations
  integer     :: ipole, ipole_s, ipole_e
  integer     :: ip, ip1, ip2, isel, iter
  !
  real(dp), parameter :: THR_ZERO=1.0d-10
  real(dp), parameter :: HIGH_VAL=10.0


  !
  ! settings
  !
  rand_radius_re=radius_re
  rand_radius_im=radius_im
  nparams = 2+3*npoles_fit
  ngen_in_memory=2
  allocate(cparams_gen(nparams,npopulation,ngen_in_memory) )
  allocate(fitness_gen(npopulation))
  allocate(fitness_best(0:niter_gen))
  !
  efermi_=0.0
  if (present(efermi)) efermi_=efermi
  allow_offset_=.true.
  if (present(allow_offset)) allow_offset_=allow_offset

  !
  ! initialization
  !
  do ip = 1, npopulation
    !
    cparams_gen(:,ip,1)=cparams(:)
    !
    if (ip > 1 ) then
      call genetic_randomize(rand_radius_re,rand_radius_im,efermi,[1,npoles_fit],cparams_gen(:,ip,1))
    endif
    !
    ierr = pole_fitting_LSQ_NL_cmplx_c(nx,xdata,ydata_,npoles_fit,cparams_gen(:,ip,1),thr,niterx,allow_offset_)
    if (ierr > 0 ) cparams_gen(nparams,ip,1) = HIGH_VAL
    !
  enddo
  !
  call genetic_reorder(nparams,npopulation,cparams_gen(:,:,1),fitness_gen)
  fitness_best(0) = fitness_gen(1)

!write(*,*) "ITERATION 0"
!write(*,"(6f15.9)") fitness_gen(1:6)
!write(*,*)

  !
  ! main loop
  !
  ierr=-1
  !
  generation_loop:&
  do iter = 1, niter_gen

    !
    ! selection step (only the most fit stay)
    ! [1-isel] are kept as is
    !
    isel=nint( select_frac*real(npopulation) )
    !
    cparams_gen(:,1:isel,2)=cparams_gen(:,1:isel,1)

    !
    ! crossover (generate children)
    !
    do ip = isel+1, npopulation
      ! pick two random numbers in the range [1,crossover_frac*isel]
      ip1=nint(crossover_frac*real(isel,dp)*randy())
      ip2=nint(crossover_frac*real(isel,dp)*randy())
      if (ip1==0) ip1=1
      if (ip2==0) ip2=1
      !
      call genetic_crossover(nparams,cparams_gen(:,ip1,2),cparams_gen(:,ip2,2),cparams_gen(:,ip,2))
      !
      ierr = pole_fitting_LSQ_NL_cmplx_c(nx,xdata,ydata_,npoles_fit,cparams_gen(:,ip,2),thr,niterx,allow_offset_)
      if (ierr > 0 ) cparams_gen(nparams,ip,2) = HIGH_VAL
      !
    enddo
    !
    ! check early convergence
    !
    fitness_gen(:) = real( cparams_gen(nparams,:,2), dp )
    skip_mutations=.false.
    if (any( fitness_gen(:)< THR_ZERO )) skip_mutations=.true.

    !
    ! mutations
    !
    if (.not. skip_mutations ) then
      !
      do ip = 4, npopulation
        !
        ! check whether we have to include mutations
        arg=randy()
        if ( arg > mutation_frac ) cycle
        !
        ! identify pole for mutation
        ipole=nint(npoles_fit*randy())
        if (ipole==0) ipole=1
        ipole_s=max(1,ipole-2)
        ipole_e=min(npoles_fit,ipole+2)
        !
        call genetic_randomize(rand_radius_re,rand_radius_im,efermi,[ipole_s,ipole_e],cparams_gen(:,ip,2))
        !
        ierr = pole_fitting_LSQ_NL_cmplx_c(nx,xdata,ydata_,npoles_fit,cparams_gen(:,ip,2),thr,niterx,allow_offset_)
        if (ierr > 0 ) cparams_gen(nparams,ip,2) = HIGH_VAL
        !
      enddo
      !
    endif

    !
    ! finalize generation
    !
    call genetic_reorder(nparams,npopulation,cparams_gen(:,:,2),fitness_gen)
    fitness_best(iter) = fitness_gen(1)
    !
    cparams_gen(nparams,:,1)=cparams_gen(nparams,:,2)

!write(*,*) "ITERATION ", iter, " LAST "
!write(*,"(6f15.9)") fitness_gen(1:6)
!write(*,*)

    !
    ! check convergence
    !
    if ( fitness_gen(1) < THR_ZERO &
         !.or. abs(fitness_best(iter)-fitness_best(iter-1))<thr &
         ) then
      ierr=0
      exit generation_loop
    else
      ierr=-1
    endif
    !
  enddo generation_loop
  !
!write(*,*) "FITNESS_BEST "
!write(*,"(6f15.9)") fitness_best(1:min(iter,niter_gen))
!write(*,*)
  cparams=cparams_gen(:,1,1)
  call cleanup_loc()
  return

contains

  subroutine cleanup_loc()
     if (allocated(cparams_gen)) deallocate(cparams_gen)
     if (allocated(fitness_gen)) deallocate(fitness_gen)
     if (allocated(fitness_best)) deallocate(fitness_best)
  end subroutine
  !
  !=========
  subroutine genetic_randomize(rand_radius_re,rand_radius_im,efermi,pole_range,cparams_gen)
    !=========
    implicit none
    real(dp)  :: rand_radius_re,rand_radius_im,efermi
    integer   :: pole_range(2)
    complex(dp) :: cparams_gen(nparams)
    !
    integer :: i,ind
    real(dp):: rand_re,rand_im, xx,dd
    !
    do i = pole_range(1),pole_range(2)
      rand_re=rand_radius_re*(randy()-0.5_dp)*2.0_dp
      rand_im=rand_radius_im*(randy()-0.5_dp)*2.0_dp
      !
      ind=1+(i-1)*3
      cparams_gen(ind+1)=cparams_gen(ind+1)+rand_re
      cparams_gen(ind+2)=cparams_gen(ind+2)+rand_im
      xx=real(cparams_gen(ind+1), dp)
      dd=real(cparams_gen(ind+2), dp)
      if (xx<efermi .and. dd<0) cparams_gen(ind+2)=-cparams_gen(ind+2)
      if (xx>efermi .and. dd>0) cparams_gen(ind+2)=-cparams_gen(ind+2)
    enddo
    !
  end subroutine genetic_randomize
  !
  !=========
  subroutine genetic_reorder(nparams,npopulation,cparams_gen,fitness_gen)
    !=========
    implicit none
    integer     :: nparams, npopulation
    complex(dp) :: cparams_gen(nparams,npopulation)
    real(dp)    :: fitness_gen(npopulation)
    !
    integer     :: ip
    integer     :: imap(npopulation)
    complex(dp), allocatable :: cparams_aux(:,:)
    !
    fitness_gen(:)=real(cparams_gen(nparams,:), dp)
    !
    ! reordering of hpsort is in ascending order
    imap=0
    call hpsort(npopulation,fitness_gen,imap)
    !
    allocate(cparams_aux(nparams,npopulation))
    !
    do ip = 1, npopulation
      cparams_aux(:,ip) = cparams_gen(:,imap(ip))
    enddo
    cparams_gen=cparams_aux
    !
    deallocate(cparams_aux)
    !
  end subroutine genetic_reorder
  !
  !=========
  subroutine genetic_crossover(nparams,cparams_gen1,cparams_gen2,cparams_child)
    !=========
    implicit none
    integer     :: nparams
    complex(dp) :: cparams_gen1(nparams), cparams_gen2(nparams), cparams_child(nparams)
    !
    integer :: npoles, ip, ind, iarg
    !
    !npoles=(nparams-2)/3
    !do ip = 1, npoles
    !  ind=1+(ip-1)*3
    !  iarg=1+nint(randy())
    !  !
    !  if (iarg==1) then
    !     cparams_child(ind+1:ind+3)=cparams_gen1(ind+1:ind+3)
    !  else
    !     cparams_child(ind+1:ind+3)=cparams_gen2(ind+1:ind+3)
    !  endif
    !enddo
    cparams_child=(cparams_gen1+cparams_gen2)/2.0_dp
    !
  end subroutine genetic_crossover

end function pole_fitting_genetic_cmplx

end module fitting_genetic_m

