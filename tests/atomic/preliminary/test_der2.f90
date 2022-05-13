
PROGRAM test
IMPLICIT NONE

  INTEGER, PARAMETER :: ndim=100
  INTEGER, PARAMETER :: dbl=KIND(1.0d0)

  REAL(dbl) :: xmin
  REAL(dbl) :: xmax
  REAL(dbl) :: zmesh
  REAL(dbl) :: dx
  REAL(dbl) :: rmax

  REAL(dbl), ALLOCATABLE :: rgrid(:)
  REAL(dbl) :: r0, rp, rpp, rm, x
  REAL(dbl) :: der1, der2
  
  INTEGER :: i, nr

!---------------------------------

  xmin   =  -7.0
  dx     =   0.01
  rmax   = 100.0
  zmesh  =   2.0
  
!--------------------------------

  xmax = LOG( rmax * zmesh )
  nr   = (xmax-xmin)/dx + 1
  nr   = 2*(nr/2)+1

  ALLOCATE( rgrid(nr) )

  DO i = 1, nr
      x = xmin + REAL(i-1,dbl)*dx
      rgrid(i)   = EXP(x) / zmesh
  ENDDO
  
  !
  ! evaluate the derivatives
  !
  OPEN ( 2, NAME="test_der2.dat" )
  DO i = 2, nr-1

      !   
      ! f(r)  = r * log(r) -r
      ! f'(r) = log(r)
      ! f''(r) = 1/r
      !   
      rm  = rgrid(i-1)
      r0  = rgrid(i)
      rp  = rgrid(i+1)
      rpp = rgrid(i+2)
      !   
      der1 = 8.0 / dx**2 / ( rm + 2.0*r0 + rp ) * ( & 
              ( rp * log(rp)-rp  -r0*log(r0)+r0 ) / ( r0 +rp ) + & 
              ( rm * log(rm)-rm  -r0*log(r0)+r0 ) / ( rm + r0 )  )
      !
!      der2 = 1.0/( r0 * dx )**2 * ( & 
!              ( rp * log(rp)-rp ) * ( 1.0 - dx/2.0 ) + &
!              ( rm * log(rm)-rm ) * ( 1.0 + dx/2.0 ) - &
!              2.0 * ( r0 * log(r0)-r0 ) )

!      der2 = 4.0 / dx**2 / ( r0 + rp ) * ( & 
!              ( rpp * log(rpp)-rpp ) / ( rp +rpp ) + & 
!              ( r0  * log(r0) -r0 )  / ( r0 +rp )  - & 
!              ( rp  * log(rp) -rp )  * ( 1.0d0/( r0 +rp ) + 1.0d0/(rp+rpp) ) )
     
      !   
      WRITE(2,"(10f25.9)") r0, 1.0_dbl/r0, der1, der2
      ! 
  ENDDO
  !
  CLOSE(2)

  DEALLOCATE( rgrid )

END PROGRAM test


