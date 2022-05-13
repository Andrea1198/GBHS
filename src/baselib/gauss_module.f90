!
! Gaussian integration module
!
!==============================
   module gauss_module
   !==============================
   use kinds
   use util_module
   implicit none
   private
   save
   !
   integer  :: nuniform=14
   integer  :: ngaussx
   real(dbl), allocatable :: xgauss(:,:)
   real(dbl), allocatable :: wgauss(:,:)
   logical  :: alloc
   !
   public :: ngaussx
   public :: xgauss
   public :: wgauss
   public :: alloc
   !
   public :: gauss_allocate
   public :: gauss_deallocate
   public :: gauss_init
   public :: gauss_abscissas
   public :: legendre_pol
   public :: legendre_pol_coeff
   public :: legendre_pol_roots
   !
contains 
   !
   subroutine gauss_allocate( ngaussx_ )
   implicit none
      integer, intent(in) :: ngaussx_
      !
      ngaussx=ngaussx_
      allocate( xgauss(ngaussx,ngaussx) )
      allocate( wgauss(ngaussx,ngaussx) )
      alloc=.true.
      !
   end subroutine gauss_allocate   
   !
   subroutine gauss_deallocate( )
   implicit none
      !
      if ( .not. alloc ) return
      !
      ngaussx=0
      deallocate( xgauss )
      deallocate( wgauss )
      alloc=.false.
      !
   end subroutine gauss_deallocate   
   !  
   subroutine gauss_init(ngaussx_)
   implicit none
      integer, intent(in) :: ngaussx_
      real(dbl) :: x, p, w
      integer   :: n, i
      !
      if ( alloc ) call gauss_deallocate()
      call gauss_allocate( ngaussx_ )
      !
      do n = 1, ngaussx_
          !
          !call colpnt( n, xgauss(:,n) )
          call gauss_abscissas( n, xgauss(:,n) )
          !
          do i = 1, min(n,nuniform)
              x=xgauss(i,n)
              p=legendre_pol(x,n+1)
              w=2.0d0*(1.0d0-x**2)/((n+1)*p)**2
              wgauss(i,n)=w
          enddo
          !
          do i = min(n,nuniform)+1, n
              x=xgauss(i,n)
              w=2.0d0/real(n,dbl)
              wgauss(i,n)=w
          enddo
          !
      enddo
      !
   end subroutine gauss_init
   !
   function legendre_pol( x, n)
   implicit none
      real(dbl) :: legendre_pol
      real(dbl) :: x
      integer   :: n
      !
      real(dbl) :: p, pm1, pm2
      integer   :: i
      ! 
      if ( n == 0 ) then 
          legendre_pol=1.0d0
          return
      elseif ( n == 1) then
          legendre_pol=x
          return
      endif
      !
      i=2
      pm1=x
      pm2=1.0d0
      !
      do while (i <= n)   
          p = (2.0d0*i-1.0d0)*x*pm1 -(1.0d0*i-1.0d0)*pm2    
          p = p/real(i,dbl)
          if ( i==n ) exit
          pm2=pm1
          pm1=p
          i=i+1
      enddo
      !
      legendre_pol=p
      !
   end function 
   !
   subroutine legendre_pol_coeff(coeff, n)
   implicit none
      integer,   intent(in)  :: n
      real(dbl), intent(out) :: coeff(0:n)
      !
      integer :: i
      real(dbl) :: c_m1(0:n)
      real(dbl) :: c_m2(0:n)
      !
      if (n==0) then
        coeff(0) = 1.0d0
        return
      elseif ( n== 1 ) then
        coeff(0) = 0.0
        coeff(1) = 1.0
        return
      endif
      !
      i=2
      c_m1=0.0
      c_m1(1) = 1.0
      c_m2=0.0
      c_m2(0) = 1.0
      !
      do while (i <= n)
         coeff(0:n) = -real(i-1,dbl)/real(i,dbl)*c_m2(0:n)
         coeff(1:n) = coeff(1:n) + real(2*i-1,dbl)/real(i,dbl)*c_m1(0:n-1) 
         !
         c_m2=c_m1
         c_m1=coeff
         i=i+1
      enddo
      !
   end subroutine
   !
   subroutine legendre_pol_roots(n, xv)
   implicit none
      !
      ! given the coefficients of Legendre polynomials,
      ! here we use the method of the companion matrix to
      ! get the roots
      !
      ! given a polynomial in the form:
      ! P(x) = c0 + c1*x + ...  + x^n
      !
      ! the companion matrix C(P), has P as characteristic
      ! polynomial, such that the eigevalues of C are the root of P 
      !
      !        (0   0   0   0  ..  -c0 )
      ! C(P) = (1   0   0   0  ..  -c1 )
      !        (0   1   0   0  ..  -c2 )
      !        (0   0   1   0  ..  -c3 )
      !        (0   0   0   .. ..  ..  )
      !        (0   0   0   0   1  -cn-1 )
      !
      integer,   intent(in)  :: n
      real(dbl), intent(out) :: xv(n)
      !
      integer :: i
      integer,      allocatable :: ind(:)
      real(dbl),    allocatable :: coeff(:)
      complex(dbl), allocatable :: zmat(:,:), z(:,:), w(:)
      !
      if (n==0) then 
        return
      elseif (n==1) then
        xv=1.0d0
        return
      endif
      !
      allocate(coeff(0:n))
      allocate(zmat(n,n))
      allocate(z(n,n),w(n))
      allocate(ind(n))
      !
      call legendre_pol_coeff(coeff,n)
      !
      zmat=0.0
      zmat(1,n) = -coeff(0)/coeff(n)
      do i = 1, n-1
        zmat(i+1,i) = 1.0d0
        zmat(i+1,n)= -coeff(i)/coeff(n)
      enddo
      !
      call zmat_diag(z,w,zmat,n,"r")
      xv=real(w,dbl)
      !
      ! sorting & symmetrization
      call hpsort(n,xv,ind)
      !
      do i = 1, n/2
        xv(n-i+1) = -xv(i)
      enddo
      !
      deallocate(coeff)
      deallocate(zmat)
      deallocate(z,w)
      deallocate(ind)
      !      
   end subroutine legendre_pol_roots
   !
   subroutine gauss_abscissas( n, xv )
   implicit none
      !
      ! Numerical values from 
      ! http://pomax.github.io/bezierinfo/legendre-gauss.html
      !
      integer,   intent(in)  :: n
      real(dbl), intent(out) :: xv(n)
      !
      integer :: i
      !
      select case ( n )
      case ( 1 )
         xv(1) =  0.0000000000000000d0
      case ( 2 )
         xv(1) = -0.5773502691896257d0
         xv(2) =  0.5773502691896257d0
      case ( 3 )
         xv(1) = -0.7745966692414834d0
         xv(2) =  0.0000000000000000d0
         xv(3) =  0.7745966692414834d0
      case ( 4 )
         xv(1) = -0.8611363115940526d0
         xv(2) = -0.3399810435848563d0
         xv(3) =  0.3399810435848563d0
         xv(4) =  0.8611363115940526d0
      case ( 5 )
         xv(1) = -0.9061798459386640d0
         xv(2) = -0.5384693101056831d0
         xv(3) =  0.0000000000000000d0
         xv(4) =  0.5384693101056831d0
         xv(5) =  0.9061798459386640d0
      case ( 6 )
         xv(1) = -0.9324695142031521d0
         xv(2) = -0.6612093864662645d0
         xv(3) = -0.2386191860831969d0
         xv(4) =  0.2386191860831969d0
         xv(5) =  0.6612093864662645d0
         xv(6) =  0.9324695142031521d0
      case ( 7 )
         xv(1) = -0.9491079123427585d0
         xv(2) = -0.7415311855993945d0
         xv(3) = -0.4058451513773972d0
         xv(4) =  0.0000000000000000d0
         xv(5) =  0.4058451513773972d0
         xv(6) =  0.7415311855993945d0
         xv(7) =  0.9491079123427585d0  
      case ( 8 )
         xv(1) = -0.9602898564975363d0
         xv(2) = -0.7966664774136267d0
         xv(3) = -0.5255324099163290d0
         xv(4) = -0.1834346424956498d0
         xv(5) =  0.1834346424956498d0
         xv(6) =  0.5255324099163290d0
         xv(7) =  0.7966664774136267d0
         xv(8) =  0.9602898564975363d0
      case ( 9 )
         xv(1) = -0.9681602395076261d0
         xv(2) = -0.8360311073266358d0
         xv(3) = -0.6133714327005904d0
         xv(4) = -0.3242534234038089d0
         xv(5) =  0.0000000000000000d0
         xv(6) =  0.3242534234038089d0
         xv(7) =  0.6133714327005904d0
         xv(8) =  0.8360311073266358d0
         xv(9) =  0.9681602395076261d0
      case ( 10 )
         xv(1) = -0.9739065285171717d0
         xv(2) = -0.8650633666889845d0
         xv(3) = -0.6794095682990244d0
         xv(4) = -0.4333953941292472d0
         xv(5) = -0.1488743389816312d0
         xv(6) =  0.1488743389816312d0
         xv(7) =  0.4333953941292472d0
         xv(8) =  0.6794095682990244d0
         xv(9) =  0.8650633666889845d0
         xv(10)=  0.9739065285171717d0
      case ( 11 )
         xv(1) = -0.9782286581460570d0
         xv(2) = -0.8870625997680953d0
         xv(3) = -0.7301520055740494d0
         xv(4) = -0.5190961292068118d0
         xv(5) = -0.2695431559523450d0
         xv(6) =  0.0000000000000000d0
         xv(7) =  0.2695431559523450d0
         xv(8) =  0.5190961292068118d0
         xv(9) =  0.7301520055740494d0
         xv(10)=  0.8870625997680953d0
         xv(11)=  0.9782286581460570d0
      case ( 12 )
         xv(1) = -0.9815606342467192d0
         xv(2) = -0.9041172563704749d0
         xv(3) = -0.7699026741943047d0
         xv(4) = -0.5873179542866175d0
         xv(5) = -0.3678314989981802d0
         xv(6) = -0.1252334085114689d0
         xv(7) =  0.1252334085114689d0
         xv(8) =  0.3678314989981802d0
         xv(9) =  0.5873179542866175d0
         xv(10)=  0.7699026741943047d0
         xv(11)=  0.9041172563704749d0
         xv(12)=  0.9815606342467192d0
      case ( 13 )
         xv(1) = -0.9841830547185881d0
         xv(2) = -0.9175983992229779d0
         xv(3) = -0.8015780907333099d0
         xv(4) = -0.6423493394403402d0
         xv(5) = -0.4484927510364469d0
         xv(6) = -0.2304583159551348d0
         xv(7) =  0.0000000000000000d0
         xv(8) =  0.2304583159551348d0
         xv(9) =  0.4484927510364469d0
         xv(10)=  0.6423493394403402d0
         xv(11)=  0.8015780907333099d0
         xv(12)=  0.9175983992229779d0
         xv(13)=  0.9841830547185881d0
      case ( 14 )
         xv(1) = -0.9862838086968123d0
         xv(2) = -0.9284348836635735d0
         xv(3) = -0.8272013150697650d0
         xv(4) = -0.6872929048116855d0
         xv(5) = -0.5152486363581541d0
         xv(6) = -0.3191123689278897d0
         xv(7) = -0.1080549487073437d0
         xv(8) =  0.1080549487073437d0
         xv(9) =  0.3191123689278897d0
         xv(10)=  0.5152486363581541d0
         xv(11)=  0.6872929048116855d0
         xv(12)=  0.8272013150697650d0
         xv(13)=  0.9284348836635735d0
         xv(14)=  0.9862838086968123d0
      case default
         !
         call legendre_pol_roots(n,xv)
         !
      end select
      !
   end subroutine gauss_abscissas
  
end module gauss_module

