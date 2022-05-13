!
! module containing basic routines to
! define and use a bspline basis set of any order
!
!==============================
   module bspline_module
   !==============================
   use kinds
   implicit none
   private
   save
   !
   type bspline_type
       integer   :: k
       integer   :: nx
       integer   :: nseg
       integer   :: ngs
       integer   :: ngl
       real(dbl), allocatable :: x(:)
       real(dbl), allocatable :: dx(:)
       real(dbl), allocatable :: y(:)
       real(dbl), allocatable :: dy(:)
       real(dbl), allocatable :: xgs(:,:)
       real(dbl), allocatable :: wgs(:,:)
       real(dbl), allocatable :: ygs(:,:)
       real(dbl), allocatable :: dygs(:,:)
       real(dbl), allocatable :: xgl(:,:)
       real(dbl), allocatable :: wgl(:,:)
       real(dbl), allocatable :: ygl(:,:)
       real(dbl), allocatable :: dygl(:,:)
       real(dbl), allocatable :: pol(:,:)
       logical   :: alloc = .false.
   end type bspline_type
   !
   public :: bspline_type
   public :: bspline_allocate
   public :: bspline_deallocate
   public :: bspline_init
   public :: bspline_eval
   public :: bspline_define
   ! 
contains 
   !
   subroutine bspline_allocate( k, ngs, ngl, obj )
   implicit none
      integer, intent(in) :: k, ngs, ngl
      type(bspline_type)  :: obj
      !
      obj%k=k
      obj%nx=k+1
      obj%nseg=k
      obj%ngs=ngs
      obj%ngl=ngl
      allocate( obj%x(obj%nx) )
      allocate( obj%dx(obj%nseg) )
      allocate( obj%y(obj%nx) )
      allocate( obj%dy(obj%nx) )
      allocate( obj%xgs(obj%ngs,obj%nseg) )
      allocate( obj%wgs(obj%ngs,obj%nseg) )
      allocate( obj%ygs(obj%ngs,obj%nseg) )
      allocate( obj%dygs(obj%ngs,obj%nseg) )
      allocate( obj%xgl(obj%ngl,obj%nseg) )
      allocate( obj%wgl(obj%ngl,obj%nseg) )
      allocate( obj%ygl(obj%ngl,obj%nseg) )
      allocate( obj%dygl(obj%ngl,obj%nseg) )
      allocate( obj%pol(0:obj%k-1,obj%nseg) )
      obj%alloc=.true.
      !
   end subroutine bspline_allocate   
   !
   subroutine bspline_deallocate( obj )
   implicit none
      type(bspline_type) :: obj
      !
      if ( .not. obj%alloc ) return
      !
      obj%k=0
      obj%nx=0
      obj%nseg=0
      obj%ngs=0
      obj%ngl=0
      if (allocated(obj%x))    deallocate( obj%x )
      if (allocated(obj%dx))   deallocate( obj%dx )
      if (allocated(obj%y))    deallocate( obj%y )
      if (allocated(obj%dy))   deallocate( obj%dy )
      if (allocated(obj%xgs))  deallocate( obj%xgs )
      if (allocated(obj%wgs))  deallocate( obj%wgs )
      if (allocated(obj%ygs))  deallocate( obj%ygs )
      if (allocated(obj%dygs)) deallocate( obj%dygs )
      if (allocated(obj%xgl))  deallocate( obj%xgl )
      if (allocated(obj%wgl))  deallocate( obj%wgl )
      if (allocated(obj%ygl))  deallocate( obj%ygl )
      if (allocated(obj%dygl)) deallocate( obj%dygl )
      if (allocated(obj%pol))  deallocate( obj%pol )
      obj%alloc=.false.
      !
   end subroutine bspline_deallocate   
   !  
!====================================================
   subroutine bspline_init( istart, ndim, grid, ileft, ngs, ngl, k, obj )
   !====================================================
   !
   ! grid starts at negative indeces to give enough knot points
   ! for the definition of the first bspline through the use of
   ! bsplvd. The physical first knot is xgrid(1)
   !
   use gauss_module, only : xgauss, wgauss, ngaussx
   implicit none
      !
      integer,   intent(in) :: ndim, k, istart
      real(dbl), intent(in) :: grid(istart:ndim)
      integer,   intent(in) :: ngs, ngl, ileft
      type(bspline_type)    :: obj
      !
      character(256) :: subname="bspline_init"
      real(dbl) :: norm
      integer   :: is, i

      !
      ! checks
      !
      if ( ngs > ngl )     call errore(subname,"ngs larger than ngl",10)
      if ( ngs > ngaussx ) call errore(subname,"ngs too large",ngs)
      if ( ngl > ngaussx ) call errore(subname,"ngl too large",ngs)

      !
      ! let's consider that each bspline of deg k requires
      ! Nseg=k, ie k+1 points. In general we have
      !
      ! ndim = nspline + k
      !
      if ( ileft <= 0 ) call errore(subname,"invalid nspline",10)

      !
      ! alloc
      !
      call bspline_deallocate( obj )
      call bspline_allocate( k, ngs, ngl, obj )
     
      !
      ! grids, and bspline evaluations
      ! - x, knots
      ! - dx
      ! - xgs, xgl, gauss integration absissas
      ! - wgs, wgl, gauss weights
      !
      obj%x(1) = grid(ileft)
      !
      do is = 1, obj%nseg
          !
          obj%x(is+1) = grid(ileft+is)
          obj%dx(is)  = grid(ileft+is) -grid(ileft+is-1)
          !
          obj%xgs(1:ngs,is) = obj%dx(is)/2.0d0 * xgauss(1:ngs,ngs) &
                                  + (grid(ileft+is)+grid(ileft+is-1))/2.0d0
          obj%wgs(1:ngs,is) = obj%dx(is)/2.0d0 * wgauss(1:ngs,ngs)
          !
          obj%xgl(1:ngl,is) = obj%dx(is)/2.0d0 * xgauss(1:ngl,ngl) &
                                  + (grid(ileft+is)+grid(ileft+is-1))/2.0d0
          obj%wgl(1:ngl,is) = obj%dx(is)/2.0d0 * wgauss(1:ngl,ngl)
          !
      enddo
      !
      ! define the polinomyal coefficients
      !
      call bspline_define( istart, ndim, grid, ileft, obj )

      !
      ! eval the spline values
      !
      call bspline_eval( obj%nx, obj%x, obj%y,  obj, 0 )
      call bspline_eval( obj%nx, obj%x, obj%dy, obj, 1 )
      !
      norm=0.0d0
      !
      do is = 1, obj%nseg
          !
          call bspline_eval( ngs, obj%xgs(:,is), obj%ygs(:,is),  obj, 0 )
          call bspline_eval( ngs, obj%xgs(:,is), obj%dygs(:,is), obj, 1 )
          !
          call bspline_eval( ngl, obj%xgl(:,is), obj%ygl(:,is),  obj, 0 )
          call bspline_eval( ngl, obj%xgl(:,is), obj%dygl(:,is), obj, 1 )
          !
          do i = 1, ngl
              norm=norm +obj%wgl(i,is)*obj%ygl(i,is)**2
          enddo
          !
      enddo

      !
      ! normalize
      !
      obj%y   = obj%y/sqrt(norm)
      obj%dy  = obj%dy/sqrt(norm)
      obj%ygs = obj%ygs/sqrt(norm)
      obj%dygs= obj%dygs/sqrt(norm)
      obj%ygl = obj%ygl/sqrt(norm)
      obj%dygl= obj%dygl/sqrt(norm)
      obj%pol = obj%pol/sqrt(norm)
      !
      return
      !
   end subroutine bspline_init

!====================================================
   subroutine bspline_eval( ndim, x, y, obj, order )
   !====================================================
   !
   ! evaluate the values of a bspline function (obj)
   ! on the list of points x. order 0 stands for the
   ! bspline functions itself, order > 1 stands for
   ! its derivatives
   ! 
   !
   implicit none
      !
      integer,   intent(in)  :: ndim
      real(dbl), intent(in)  :: x(ndim)
      real(dbl), intent(out) :: y(ndim)
      type(bspline_type)  :: obj
      integer,   intent(in)  :: order
      !
      integer   :: k, i, ik, is, iss
      real(dbl) :: x0
      character(12) :: subname="bspline_eval"
      !
      ! this is the target bspline order
      k=obj%k
      if ( order < 0 ) call errore(subname,"order < 0", -order)
      !
      x0 = obj%x(1)
      !
      do i = 1, ndim
          !
          y(i) = 0.0d0
          !
          if ( x(i) <= obj%x(1) .or. x(i) >= obj%x(obj%nx) ) cycle
          !
          is = 1
          do iss = 1, obj%nseg
              if ( x(i) <= obj%x(iss+1) .and. x(i) >= obj%x(iss) ) then
                  is=iss
                  exit
              endif
          enddo
          !
          if ( order == 0 ) then
              !
              do ik = 0, k-1
                  y(i) = y(i) + obj%pol(ik,is) * (x(i)-x0)**ik
              enddo
              !
          elseif ( order == 1 ) then
              !
              do ik = 1, k-1
                  y(i) = y(i) + real(ik) * obj%pol(ik,is) * (x(i)-x0)**(ik-1)
              enddo
              !
          elseif ( order == 2 ) then
              !
              do ik = 2, k-1
                  y(i) = y(i) + real(ik*(ik-1)) * obj%pol(ik,is) * (x(i)-x0)**(ik-2)
              enddo
              !
          elseif ( order == 3 ) then
              !
              do ik = 3, k-1
                  y(i) = y(i) + real(ik*(ik-1)*(ik-2)) * obj%pol(ik,is) * (x(i)-x0)**(ik-3)
              enddo
              !
          else
              call errore(subname,"derivative not implemented",order)
          endif
          !
      enddo
      !
   end subroutine bspline_eval

!=========================================================
   subroutine bspline_define( istart, ndim, grid, ileft, obj )
   !=========================================================
   implicit none
      !
      integer,   intent(in) :: ndim, istart, ileft
      real(dbl), intent(in) :: grid(istart:ndim)
      type(bspline_type)    :: obj
      !
      character(14) :: subname="bspline_define"
      !
      integer   :: is, ik, ib, k, ind, i0
      real(dbl) :: dxk0, dxkp
      real(dbl), allocatable :: coeff0(:,:)
      real(dbl), allocatable :: coeff(:,:)
    
      if ( .not. obj%alloc ) call errore(subname,"obj not alloc", 10) 
    
      !
      ! local workspace
      !
      k=obj%k
      allocate( coeff0(0:k-1,k), coeff(0:k-1,k) )
    
      segment_loop:&
      do is = 1, obj%nseg
          !
          coeff0(:,:) = 0.0d0
          coeff(:,:)  = 0.0d0
          !
          coeff0(0,1) = 1.0d0
          !
          i0=ileft+is-1
          !
          do ik = 2, k
              !
              dxk0=grid(i0)   -grid(i0-ik+1) 
              dxkp=grid(i0+1) -grid(i0-ik+2) 
              !
              coeff(1:k-1,1) = coeff(1:k-1,1) -1.0d0/dxkp * coeff0(0:k-2,1) 
              coeff(0:k-2,1) = coeff(0:k-2,1) +(grid(i0+1)-grid(ileft))/dxkp * coeff0(0:k-2,1) 
              !
              do ib = 1, ik-1
                  !
                  dxk0=grid(i0+ib)   -grid(i0-ik+ib+1) 
                  dxkp=grid(i0+ib+1) -grid(i0-ik+ib+2) 
                  !
                  coeff(1:k-1,ib+1) = coeff(1:k-1,ib+1) +1.0d0/dxk0 * coeff0(0:k-2,ib) 
                  coeff(1:k-1,ib+1) = coeff(1:k-1,ib+1) -1.0d0/dxkp * coeff0(0:k-2,ib+1) 
                  !
                  coeff(0:k-2,ib+1) = coeff(0:k-2,ib+1) &
                                      -(grid(i0-ik+ib+1)-grid(ileft))/dxk0   * coeff0(0:k-2,ib) 
                  coeff(0:k-2,ib+1) = coeff(0:k-2,ib+1) &
                                      +(grid(i0+ib+1)-grid(ileft))/dxkp* coeff0(0:k-2,ib+1) 
                  !
              enddo
              !
              coeff0 = coeff
              coeff  = 0.0d0
              !
          enddo
          !
          ind = k-is+1
          obj%pol(0:k-1,is) = coeff0(0:k-1,ind)
          !
      enddo segment_loop
      !
      ! cleanup
      !
      deallocate( coeff0, coeff )
      !
      return
      !
    end subroutine bspline_define

end module bspline_module

