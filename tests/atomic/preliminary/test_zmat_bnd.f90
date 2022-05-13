
PROGRAM test
!
USE util_module
IMPLICIT NONE

  INTEGER, PARAMETER :: ndim=6
  INTEGER, PARAMETER :: dbl=KIND(1.0d0)

  INTEGER :: kl, ku

  COMPLEX(dbl) :: z(ndim, ndim)
  COMPLEX(dbl) :: zb(ndim, ndim)

  INTEGER :: i,j

!---------------------------------

  kl=2
  ku=1

  DO j=1,ndim
  DO i=1,ndim
      z(i,j) = REAL(j + 10*i)
      IF ( j-i > ku ) z(i,j) = 0 
      IF ( i-j > kl ) z(i,j) = 0 
  ENDDO
  ENDDO
  
!--------------------------------

  WRITE(6, "(    2x,'Z before packing')" )
  CALL write_matrix(ndim, z)

  CALL zmat_bnd_pack( zb, z, ndim, kl, ku)

  WRITE(6, "(2/, 2x,'Zb packed')" )
  CALL write_matrix(ndim, zb)

  CALL zmat_bnd_unpack( z, zb, ndim, kl, ku)

  WRITE(6, "(2/, 2x,'Z after unpacking')" )
  CALL write_matrix(ndim, z)

!--------------------------------


CONTAINS

SUBROUTINE write_matrix(n,zmat)
IMPLICIT NONE

   INTEGER :: n
   COMPLEX(dbl) :: zmat(n,n)
   !
   INTEGER :: i, j

   WRITE(6,*)
   DO i=1,ndim
      WRITE(6,"(100f7.2)") REAL(zmat(i,:))
   ENDDO

END SUBROUTINE

END PROGRAM test

