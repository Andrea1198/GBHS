!
! Copyright (C) 2004 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! <INFO>
!*********************************************
  MODULE util_module
  !*********************************************
  !
  USE kinds
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI
  !
  IMPLICIT NONE
  PRIVATE

! General purpose utilities
!
! routines in this module:
! SUBROUTINE  zmat_pack( zp, z, n)
! SUBROUTINE  zmat_unpack( z, zp, n)
! SUBROUTINE   mat_bnd_pack( zb, z, n, kl, ku)
! SUBROUTINE   mat_bnd_unpack( z, zb, n, kl, ku)
! SUBROUTINE   mat_bnd_getdims( z, n, kl, ku)
! SUBROUTINE   mat_herm( z, n)
! SUBROUTINE   mat_antiherm( z, n)
! SUBROUTINE   mat_svd( m, n, a, s, u, vt)
! SUBROUTINE   mat_sv ( n, nrhs, a, b [,ierr])
! SUBROUTINE   mat_lsd( m, n, nrhs, a, b, rcond [, sv, rank, ierr])
! SUBROUTINE   mat_mul( c, a, opa, b, opb, m, n, k)
! SUBROUTINE   mat_vec_mul( c, a, opa, b, m, n)
! SUBROUTINE   mat_hdiag( z, w, a, n[, uplo][, il, iu])
! SUBROUTINE   mat_inv( n, a, z [, kl, ku, ldab] [,det_a] [,ierr] )
! SUBROUTINE   mat_hsqrt( n, a, [,ldab, ka] z [,ierr] )
! SUBROUTINE  zmat_diag( z, w, a, n, side)
! SUBROUTINE   mat_chol( a, tmat, n, uplo, ierr)
! COMPLEX/REAL FUNCTION   mat_dotp( m, n, a, b)
! COMPLEX/REAL FUNCTION   mat_hdotp( m, a, b)
! COMPLEX/REAL FUNCTION   vec_dotp( m, a, b)
! LOGICAL FUNCTION  zmat_is_unitary( m, n, z [,side] [,toll])
! LOGICAL FUNCTION   mat_is_herm( n, z [,toll])
! LOGICAL FUNCTION   mat_is_symm( n, z [,toll])
! LOGICAL FUNCTION   mat_is_diag_dom( n, z [,toll])
! INTEGER FUNCTION   mat_rank( m, n, a, toll)
! 
! </INFO>
!

!
! banded matrix tools
INTERFACE mat_bnd_getdims
   MODULE PROCEDURE zmat_bnd_getdims
   MODULE PROCEDURE dmat_bnd_getdims
END INTERFACE
INTERFACE mat_bnd_pack
   MODULE PROCEDURE zmat_bnd_pack
   MODULE PROCEDURE dmat_bnd_pack
   MODULE PROCEDURE dmat_hbnd_pack
END INTERFACE
INTERFACE mat_bnd_unpack
   MODULE PROCEDURE zmat_bnd_unpack
   MODULE PROCEDURE dmat_bnd_unpack
   MODULE PROCEDURE dmat_hbnd_unpack
END INTERFACE
!
! hermitian, antohermitian parts
INTERFACE mat_herm
   MODULE PROCEDURE zmat_herm
   MODULE PROCEDURE dmat_herm
END INTERFACE
INTERFACE mat_antiherm
   MODULE PROCEDURE zmat_antiherm
END INTERFACE
!
INTERFACE mat_is_herm
   MODULE PROCEDURE zmat_is_herm
   MODULE PROCEDURE dmat_is_herm
END INTERFACE
!
INTERFACE mat_is_symm
   MODULE PROCEDURE zmat_is_symm
   MODULE PROCEDURE dmat_is_symm
END INTERFACE
!
INTERFACE mat_is_diagdom
   MODULE PROCEDURE zmat_is_diagdom
   MODULE PROCEDURE dmat_is_diagdom
END INTERFACE
!
! matrix multiplication
INTERFACE mat_mul
   MODULE PROCEDURE zmat_mul
   MODULE PROCEDURE dmat_mul
END INTERFACE
!
! matrix-vector multiplication
INTERFACE mat_vec_mul
    ! MODULE PROCEDURE dmat_vec_mul 
    MODULE PROCEDURE dmat_vec_mul
    MODULE PROCEDURE zmat_vec_mul
END INTERFACE
!
! singular value decomposition
INTERFACE mat_svd
   MODULE PROCEDURE zmat_svd
   MODULE PROCEDURE dmat_svd
END INTERFACE
!
! schur decomposition
INTERFACE mat_schur
   MODULE PROCEDURE zmat_schur
   !MODULE PROCEDURE dmat_schur
END INTERFACE
!
! simple linear system solver
INTERFACE mat_sv
   MODULE PROCEDURE zmat_sv
   MODULE PROCEDURE zmat_sv_1
   MODULE PROCEDURE dmat_sv
   MODULE PROCEDURE dmat_sv_1
END INTERFACE
!
INTERFACE mat_lsd
   MODULE PROCEDURE dmat_lsd
   MODULE PROCEDURE dmat_lsd_1
   MODULE PROCEDURE zmat_lsd
   MODULE PROCEDURE zmat_lsd_1
END INTERFACE
!
! rank calculation
INTERFACE mat_rank
   MODULE PROCEDURE zmat_rank
   MODULE PROCEDURE dmat_rank
END INTERFACE
!
! matrix diagonalization
INTERFACE mat_hdiag
   MODULE PROCEDURE zmat_hdiag
   MODULE PROCEDURE zmat_hdiag_pack
   MODULE PROCEDURE zmat_hdiag_gen
   MODULE PROCEDURE dmat_hdiag
   MODULE PROCEDURE dmat_hdiag_gen
   MODULE PROCEDURE dmat_bnd_hdiag
   MODULE PROCEDURE dmat_bnd_hdiag_ext
   MODULE PROCEDURE dmat_bnd_hdiag_gen_ext
END INTERFACE
!
INTERFACE zmat_diag
   MODULE PROCEDURE zmat_diag_oneside
   MODULE PROCEDURE zmat_diag_twoside
END INTERFACE
!
! matrix inversion
INTERFACE mat_inv
   MODULE PROCEDURE zmat_inv
   MODULE PROCEDURE zmat_bnd_inv
   MODULE PROCEDURE dmat_inv
   MODULE PROCEDURE dmat_bnd_inv
END INTERFACE
!
! matrix cholesky decomposition
INTERFACE mat_chol
   MODULE PROCEDURE zmat_chol
END INTERFACE
!
! matrix sqrt
INTERFACE mat_hsqrt
   MODULE PROCEDURE zmat_hsqrt
   MODULE PROCEDURE dmat_bnd_hsqrt
   MODULE PROCEDURE dmat_hsqrt
END INTERFACE
!
! matrix dot product
INTERFACE  mat_dotp
   MODULE PROCEDURE zmat_ge_dotp
   MODULE PROCEDURE dmat_ge_dotp
END INTERFACE
!
! matrix herm dot product
INTERFACE  mat_hdotp
   MODULE PROCEDURE zmat_hp_dotp
   MODULE PROCEDURE zmat_he_dotp
   MODULE PROCEDURE dmat_sy_dotp
END INTERFACE
!
! vec dot product
INTERFACE  vec_dotp
   MODULE PROCEDURE dvec_dotp
   MODULE PROCEDURE zvec_dotp
END INTERFACE



PUBLIC :: zmat_pack, zmat_unpack
PUBLIC ::  mat_bnd_pack,  mat_bnd_unpack
PUBLIC ::  mat_bnd_getdims
PUBLIC ::  mat_herm
PUBLIC ::  mat_antiherm
PUBLIC ::  mat_svd
PUBLIC ::  mat_sv
PUBLIC ::  mat_lsd
PUBLIC ::  mat_schur
PUBLIC ::  mat_mul
PUBLIC ::  mat_vec_mul
PUBLIC ::  mat_hdiag
PUBLIC ::  mat_inv
PUBLIC ::  mat_hsqrt
PUBLIC :: zmat_diag
PUBLIC ::  mat_chol
PUBLIC :: zmat_is_unitary
PUBLIC ::  mat_is_symm
PUBLIC ::  mat_is_diagdom
PUBLIC ::  mat_is_herm
PUBLIC ::  mat_rank
PUBLIC ::  mat_hdotp
PUBLIC ::  mat_dotp
PUBLIC ::  vec_dotp

CONTAINS

!
! Subroutines
!

!**********************************************************
   SUBROUTINE zmat_pack( zp, z, n )
   !**********************************************************
    IMPLICIT NONE
    COMPLEX(dbl), INTENT(OUT) :: zp(:)
    COMPLEX(dbl), INTENT(IN)  :: z(:,:)
    INTEGER,      INTENT(IN)  :: n
    !
    INTEGER :: i, j, ind

    ind = 1
    DO j = 1, n
      DO i = 1, n
        zp(ind) = z(i,j)
        ind = ind + 1
      END DO
    END DO

    RETURN
  END SUBROUTINE


!**********************************************************
   SUBROUTINE zmat_unpack( z, zp, n )
   !**********************************************************
   IMPLICIT NONE
   COMPLEX(dbl), INTENT(IN)  :: zp(:)
   COMPLEX(dbl), INTENT(OUT) :: z(:,:)
   INTEGER,      INTENT(IN)  :: n
   INTEGER :: i, j, ind

   ind = 1
   DO j = 1, n
     DO i = 1, n
       z(i,j) = zp(ind)
       ind = ind + 1
     END DO
   END DO

   RETURN
END SUBROUTINE


!**********************************************************
   SUBROUTINE zmat_bnd_pack( zb, z, n, kl, ku )
   !**********************************************************
    IMPLICIT NONE
    COMPLEX(dbl), INTENT(OUT)  :: zb(:,:)
    COMPLEX(dbl), INTENT(IN)   :: z(:,:)
    INTEGER,      INTENT(IN)   :: n
    INTEGER,      INTENT(IN)   :: kl, ku
    !
    ! kl :   # of subdiagonals
    ! ku :   # of supradiagonals
    ! ldzb must be >= 2*kl + ku +1
    !
    INTEGER :: i, j, ind

    IF ( SIZE(zb,1) < 2 *kl + ku + 1 ) CALL errore('zmat_bnd_pack','invalid ldzb',10) 
    !
    zb(:,:) = CZERO
    !
    DO j = 1, n
    DO i = MAX(1,j-ku),MIN(n,j+kl)
        ind = kl+ku+1+i-j
        zb(ind,j) = z(i,j)
    ENDDO
    ENDDO

    RETURN
  END SUBROUTINE zmat_bnd_pack

    
!**********************************************************
   SUBROUTINE zmat_bnd_unpack( z, zb, n, kl, ku )
   !**********************************************************
    IMPLICIT NONE
    COMPLEX(dbl), INTENT(OUT) :: z(:,:)
    COMPLEX(dbl), INTENT(IN)  :: zb(:,:)
    INTEGER,      INTENT(IN)  :: n
    INTEGER,      INTENT(IN)  :: kl, ku
    !
    ! kl :   # of subdiagonals
    ! ku :   # of supradiagonals
    ! ldzb must be >= 2*kl + ku +1
    !
    INTEGER :: i, j, ind

    IF ( SIZE(zb,1) < 2 *kl + ku + 1 ) CALL errore('zmat_bnd_unpack','invalid ldzb',10) 
    !
    z(:,:) = CZERO
    !
    DO j = 1, n
    DO i = MAX(1,j-ku),MIN(n,j+kl)
        ind = kl+ku+1+i-j
        z(i,j) = zb(ind,j)
    ENDDO
    ENDDO

    RETURN
  END SUBROUTINE zmat_bnd_unpack


!**********************************************************
   SUBROUTINE dmat_bnd_pack( zb, z, n, kl, ku )
   !**********************************************************
    IMPLICIT NONE
    REAL(dbl),    INTENT(OUT)  :: zb(:,:)
    REAL(dbl),    INTENT(IN)   :: z(:,:)
    INTEGER,      INTENT(IN)   :: n
    INTEGER,      INTENT(IN)   :: kl, ku
    !
    ! kl :   # of subdiagonals
    ! ku :   # of supradiagonals
    ! ldzb must be >= 2*kl + ku +1
    !
    INTEGER :: i, j, ind

    IF ( SIZE(zb,1) < 2 *kl + ku + 1 ) CALL errore('dmat_bnd_pack','invalid ldzb',10) 
    !
    zb(:,:) = ZERO
    !
    DO j = 1, n
    DO i = MAX(1,j-ku),MIN(n,j+kl)
        ind = kl+ku+1+i-j
        zb(ind,j) = z(i,j)
    ENDDO
    ENDDO

    RETURN
  END SUBROUTINE dmat_bnd_pack

    
!**********************************************************
   SUBROUTINE dmat_bnd_unpack( z, zb, n, kl, ku )
   !**********************************************************
    IMPLICIT NONE
    REAL(dbl),    INTENT(OUT) :: z(:,:)
    REAL(dbl),    INTENT(IN)  :: zb(:,:)
    INTEGER,      INTENT(IN)  :: n
    INTEGER,      INTENT(IN)  :: kl, ku
    !
    ! kl :   # of subdiagonals
    ! ku :   # of supradiagonals
    ! ldzb must be >= 2*kl + ku +1
    !
    INTEGER :: i, j, ind

    IF ( SIZE(zb,1) < 2 *kl + ku + 1 ) CALL errore('dmat_bnd_unpack','invalid ldzb',10) 
    !
    z(:,:) = ZERO
    !
    DO j = 1, n
    DO i = MAX(1,j-ku),MIN(n,j+kl)
        ind = kl+ku+1+i-j
        z(i,j) = zb(ind,j)
    ENDDO
    ENDDO

    RETURN
  END SUBROUTINE dmat_bnd_unpack


!**********************************************************
   SUBROUTINE dmat_hbnd_pack( zb, z, n, kd )
   !**********************************************************
    IMPLICIT NONE
    REAL(dbl),    INTENT(OUT)  :: zb(:,:)
    REAL(dbl),    INTENT(IN)   :: z(:,:)
    INTEGER,      INTENT(IN)   :: n
    INTEGER,      INTENT(IN)   :: kd
    !
    ! kd :   # of supradiagonals
    !
    INTEGER :: i, j, ind

    IF ( SIZE(zb,1) < kd+1 ) CALL errore('dmat_hbnd_pack','invalid ldzb',10) 
    !
    zb(:,:) = ZERO
    !
    DO j = 1, n
    DO i = MAX(1,j-kd),j
        ind = kd+1+i-j
        zb(ind,j) = z(i,j)
    ENDDO
    ENDDO

    RETURN
  END SUBROUTINE dmat_hbnd_pack

    
!**********************************************************
   SUBROUTINE dmat_hbnd_unpack( z, zb, n, kd )
   !**********************************************************
    IMPLICIT NONE
    REAL(dbl),    INTENT(OUT) :: z(:,:)
    REAL(dbl),    INTENT(IN)  :: zb(:,:)
    INTEGER,      INTENT(IN)  :: n
    INTEGER,      INTENT(IN)  :: kd
    !
    ! kd :   # of supradiagonals
    !
    INTEGER :: i, j, ind

    IF ( SIZE(zb,1) < kd + 1 ) CALL errore('dmat_hbnd_unpack','invalid ldzb',10) 
    !
    z(:,:) = ZERO
    !
    DO j = 1, n
    DO i = MAX(1,j-kd),j
        ind = kd+1+i-j
        z(i,j) = zb(ind,j)
    ENDDO
    ENDDO

    RETURN
  END SUBROUTINE dmat_hbnd_unpack


!**********************************************************
   SUBROUTINE zmat_bnd_getdims( a, thr, n, kl, ku )
   !**********************************************************
    IMPLICIT NONE
    COMPLEX(dbl), INTENT(IN)   :: a(:,:)
    REAL(dbl),    INTENT(IN)   :: thr
    INTEGER,      INTENT(IN)   :: n
    INTEGER,      INTENT(OUT)  :: kl, ku
    !
    ! kl :   # of subdiagonals
    ! ku :   # of supradiagonals
    !
    INTEGER :: i, j
    !
    kl = 0
    ku = 0
    !
    DO j = 1, n
    DO i = 1, j-1
        IF ( REAL( a(i,j)*CONJG(a(i,j) ) ) > thr .AND. &
             j-i > kl ) kl = j-i  
    ENDDO
    ENDDO

    DO j = 1, n
    DO i = j+1, n
        IF ( REAL( a(i,j)*CONJG(a(i,j) ) ) > thr .AND. &
             i-j > ku ) ku = i-j
    ENDDO
    ENDDO

    RETURN
  END SUBROUTINE zmat_bnd_getdims


!**********************************************************
   SUBROUTINE dmat_bnd_getdims( a, thr, n, kl, ku )
   !**********************************************************
    IMPLICIT NONE
    REAL(dbl),    INTENT(IN)   :: a(:,:)
    REAL(dbl),    INTENT(IN)   :: thr
    INTEGER,      INTENT(IN)   :: n
    INTEGER,      INTENT(OUT)  :: kl, ku
    !
    ! kl :   # of subdiagonals
    ! ku :   # of supradiagonals
    !
    INTEGER :: i, j
    !
    kl = 0
    ku = 0
    !
    DO j = 1, n
    DO i = 1, j-1
        IF ( a(i,j)*a(i,j)  > thr .AND. &
             j-i > kl ) kl = j-i  
    ENDDO
    ENDDO

    DO j = 1, n
    DO i = j+1, n
        IF ( a(i,j)*a(i,j) > thr .AND. &
             i-j > ku ) ku = i-j
    ENDDO
    ENDDO

    RETURN
  END SUBROUTINE dmat_bnd_getdims


!**********************************************************
   SUBROUTINE zmat_herm( z, n )
   !**********************************************************
   !
   ! overwrite the input matrix with its hermitean part
   !
   IMPLICIT NONE
   COMPLEX(dbl), INTENT(INOUT) :: z(:,:)
   INTEGER,      INTENT(IN)    :: n
   INTEGER :: i, j

   DO j = 1, n
   DO i = j, n
       z(i,j) = 0.5_dbl * ( z(i,j) + CONJG( z(j,i) ) )
       z(j,i) = CONJG( z(i,j) )
   ENDDO
   ENDDO

   RETURN
END SUBROUTINE zmat_herm


!**********************************************************
   SUBROUTINE dmat_herm( z, n )
   !**********************************************************
   !
   ! overwrite the input matrix with its hermitean part
   !
   IMPLICIT NONE
   REAL(dbl),    INTENT(INOUT) :: z(:,:)
   INTEGER,      INTENT(IN)    :: n
   INTEGER :: i, j

   DO j = 1, n
   DO i = j, n
       z(i,j) = 0.5_dbl * ( z(i,j) + z(j,i) )
       z(j,i) = z(i,j)
   ENDDO
   ENDDO

   RETURN
END SUBROUTINE dmat_herm


!**********************************************************
   SUBROUTINE zmat_antiherm( z, n )
   !**********************************************************
   !
   ! overwrite the input matrix with its antihermitean part
   !
   IMPLICIT NONE
   COMPLEX(dbl), INTENT(INOUT) :: z(:,:)
   INTEGER,      INTENT(IN)    :: n
   INTEGER :: i, j

   DO j = 1, n
   DO i = j, n
       z(i,j) = -CI * 0.5_dbl * ( z(i,j) - CONJG( z(j,i) ) )
       z(j,i) = CONJG( z(i,j) )
   ENDDO
   ENDDO

   RETURN
END SUBROUTINE


!**********************************************************
   SUBROUTINE dmat_svd(m, n, a, s, u, vt)
   !**********************************************************
   !
   !  computes the singular value decomposition (SVD) of a REAL(DP)
   !  M-by-N matrix A. The SVD is written
   !
   !       A = U * SIGMA * transpose(V)
   !
   !  where SIGMA is an M-by-N matrix which is zero except for its
   !  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
   !  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
   !  are the singular values of A; they are real and non-negative, and
   !  are returned in descending order. 
   !  Note that the routine returns V**T, not V.
   !
   IMPLICIT NONE
   INTEGER,   INTENT(IN)  :: m, n
   REAL(dbl), INTENT(IN)  :: a(:,:)
   REAL(dbl), INTENT(OUT) :: s(:)
   REAL(dbl), INTENT(OUT) :: u(:,:), vt(:,:)

   INTEGER   :: ierr, info, lwork
   !REAL(dbl) :: raux 
   REAL(dbl), ALLOCATABLE :: atmp(:,:), work(:)

   IF ( m <= 0 .OR. n<=0 ) CALL errore('dmat_svd','Invalid DIMs',1)
   IF ( m > SIZE(a,1) .OR. m > SIZE(u,1) .OR. m > SIZE(u,2) ) &
           CALL errore('dmat_svd','m too large',m)
   IF ( n > SIZE(a,2) .OR. n > SIZE(vt,1) .OR. n > SIZE(vt,2) ) &
           CALL errore('dmat_svd','n too large',n)
   IF ( SIZE(s) < MIN(m,n) ) CALL errore('dmat_svd','s dimension too small',1)

   !
   ! allocate local variables and workspace
   !
   ALLOCATE( atmp(m,n), STAT=ierr )
   IF (ierr/=0)  CALL errore('dmat_svd','allocating atmp',ABS(ierr))
   !
   ! save A (which is intent IN)
   atmp(:,:) = a(1:m, 1:n)

   !
   ! get lwork
   !
   lwork = MAX( 3*MIN(m,n) + MAX(m,n), 5*MIN(m,n) )
   !
!   lwork = -1
!   !
!   CALL DGESVD('A','A', m, n, atmp, m, s, u, SIZE(u,1), vt, SIZE(vt,1), &
!                raux, lwork, info)
!   !
!   lwork = NINT( raux )
!   !
   ALLOCATE( work(lwork), STAT=ierr )
   IF (ierr/=0)  CALL errore('dmat_svd','allocating work',ABS(ierr))


   ! use magma if possible
   CALL DGESVD('A','A', m, n, atmp, m, s, u, SIZE(u,1), vt, SIZE(vt,1), &
                work, lwork, info)

   IF ( info < 0 ) CALL errore('dmat_svd', 'DGESVD: info illegal value', -info )
   IF ( info > 0 ) CALL errore('dmat_svd', 'DGESVD: DBESQR not converged', info )
    
   DEALLOCATE( atmp, work, STAT=ierr)
   IF(ierr/=0) CALL errore('dmat_svd','deallocating atpm, work',ABS(ierr))

   RETURN
END SUBROUTINE dmat_svd


!**********************************************************
   SUBROUTINE zmat_svd(m, n, a, s, u, vt)
   !**********************************************************
   !
   !  computes the singular value decomposition (SVD) of a complex
   !  M-by-N matrix A. The SVD is written
   !
   !       A = U * SIGMA * conjugate-transpose(V)
   !
   !  where SIGMA is an M-by-N matrix which is zero except for its
   !  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
   !  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
   !  are the singular values of A; they are real and non-negative, and
   !  are returned in descending order. 
   !  Note that the routine returns V**H, not V.
   !
   IMPLICIT NONE
   INTEGER, INTENT(IN)       :: m, n
   COMPLEX(dbl), INTENT(IN)  :: a(:,:)
   REAL(dbl),    INTENT(OUT) :: s(:)
   COMPLEX(dbl), INTENT(OUT) :: u(:,:), vt(:,:)

   INTEGER    :: ierr, info, lwork
   !REAL(dbl)  :: raux
   REAL(dbl),    ALLOCATABLE :: rwork(:)
   COMPLEX(dbl), ALLOCATABLE :: atmp(:,:), work(:)

   IF ( m <= 0 .OR. n<=0 ) CALL errore('zmat_svd','Invalid DIMs',1)
   IF ( m > SIZE(a,1) .OR. m > SIZE(u,1) .OR. m > SIZE(u,2) ) &
           CALL errore('zmat_svd','m too large',m)
   IF ( n > SIZE(a,2) .OR. n > SIZE(vt,1) .OR. n > SIZE(vt,2) ) &
           CALL errore('zmat_svd','n too large',n)
   IF ( SIZE(s) < MIN(m,n) ) CALL errore('zmat_svd','s dimension too small',1)

   !
   ! allocate local variables and workspace
   !
   ALLOCATE( atmp(m,n), STAT=ierr )
   IF (ierr/=0)  CALL errore('zmat_svd','allocating atmp',ABS(ierr))
   !
   ALLOCATE( rwork(5 * MIN(m,n) ), STAT=ierr )
   IF (ierr/=0)  CALL errore('zmat_svd','allocating rwork',ABS(ierr))
   !
   ! save A (which is intent IN)
   atmp(:,:) = a(1:m, 1:n)
   
   !
   ! determine lwork
   !
   lwork = 2 * MIN(m,n) + MAX(m,n)
!   !
!   lwork = -1
!   !
!   CALL ZGESVD('A','A', m, n, atmp, m, s, u, SIZE(u,1), vt, SIZE(vt,1), &
!                raux, lwork, rwork, info)
!   !
!   lwork = NINT( raux )
!   ! 
   ALLOCATE( work(lwork), STAT=ierr )
   IF (ierr/=0)  CALL errore('zmat_svd','allocating work',ABS(ierr))


   ! use magma if possible
   CALL ZGESVD('A','A', m, n, atmp, m, s, u, SIZE(u,1), vt, SIZE(vt,1), &
                work, lwork, rwork, info)

   IF ( info < 0 ) CALL errore('zmat_svd', 'ZGESVD: info illegal value', -info )
   IF ( info > 0 ) CALL errore('zmat_svd', 'ZGESVD: ZBESQR not converged', info )
    
   DEALLOCATE( atmp, work, rwork, STAT=ierr)
   IF(ierr/=0) CALL errore('zmat_svd','deallocating atmp, work, rwork',ABS(ierr))

   RETURN
END SUBROUTINE zmat_svd


!**********************************************************
   SUBROUTINE dmat_sv(n, nrhs, a, b, ierr)
   !**********************************************************
   !
   !  Computes the solution of the real system of linear equations
   !     A * X = B,
   !  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
   !
   !  The LU decomposition with partial pivoting and row interchanges is
   !  used to factor A as
   !     A = P * L * U,
   !  where P is a permutation matrix, L is unit lower triangular, and U is
   !  upper triangular.  The factored form of A is then used to solve the
   !  system of equations A * X = B.
   !
   IMPLICIT NONE
   INTEGER,   INTENT(IN)       :: n, nrhs
   REAL(dbl), INTENT(IN)       :: a(:,:)
   REAL(dbl), INTENT(INOUT)    :: b(:,:)
   INTEGER, OPTIONAL, INTENT(out) :: ierr

   INTEGER :: ierr_, info
   INTEGER,      ALLOCATABLE :: ipiv(:)
   REAL(dbl), ALLOCATABLE    :: atmp(:,:)

   IF ( PRESENT(ierr) ) ierr=0
   !
   IF ( n > SIZE(a,1) .OR. n > SIZE(a,2) ) CALL errore('dmat_sv','matrix A too small',1)
   IF ( n > SIZE(b,1) ) CALL errore('dmat_sv','matrix B too small (I)',1)
   IF ( nrhs > SIZE(b,2) ) CALL errore('dmat_sv','matrix B too small (II)',1)

   ALLOCATE( atmp(n,n), ipiv(n), STAT=ierr_ )
     IF (ierr_/=0) CALL errore('dmat_sv','allocating atmp, ipiv',ABS(ierr_))
   !
   ! make a local copy of a
   atmp(:,:) = a(1:n,1:n) 

   ! use magma if possible
   CALL DGESV( n, nrhs, atmp, n, ipiv, b, SIZE(b,1), info)

   IF ( PRESENT(ierr) ) THEN
        IF (info/=0) ierr= info
   ELSE
        IF ( info < 0 ) CALL errore('dmat_sv', 'DGESV: info illegal value', -info )
        IF ( info > 0 ) CALL errore('dmat_sv', 'DGESV: singular matrix', info )
   ENDIF
    
   DEALLOCATE( atmp, ipiv, STAT=ierr_)
      IF(ierr_/=0) CALL errore('dmat_sv','deallocating atmp, ipiv',ABS(ierr_))

   RETURN
END SUBROUTINE dmat_sv


!**********************************************************
   SUBROUTINE dmat_sv_1(n, nrhs, a, b, ierr)
   !**********************************************************
   !
   ! Interface to dmat_sv when nrhs = 1
   !
   IMPLICIT NONE
   INTEGER,         INTENT(IN)    :: n, nrhs
   REAL(dbl),       INTENT(IN)    :: a(:,:)
   REAL(dbl),       INTENT(INOUT) :: b(:)
   INTEGER, OPTIONAL, INTENT(out) :: ierr
   !
   INTEGER  :: ierr_, info
   REAL(dbl),   ALLOCATABLE :: bl(:,:)

   IF ( PRESENT(ierr) ) ierr=0
   !
   IF ( nrhs /= 1) CALL errore('dmat_sv_1','more than 1 rhs ? ',ABS(nrhs)+1)
   IF ( n > SIZE(b,1) ) CALL errore('dmat_sv_1','vector b too small',2)
   !
   ALLOCATE( bl(n,1), STAT=ierr_ )
   IF (ierr_/=0) CALL errore('dmat_sv_1','allocating bl',ABS(ierr_))

   bl(:,1) = b(1:n)
   CALL dmat_sv( n, 1, a, bl, IERR=info)
   !
   IF ( PRESENT(ierr) ) THEN 
       ierr=info
   ELSE
       IF ( info /=0 ) CALL errore('dmat_sv_1','info/=0',ABS(info))
   ENDIF
   !
   b(1:n) = bl(1:n,1)

   DEALLOCATE( bl, STAT=ierr_)
   IF (ierr_/=0) CALL errore('dmat_sv_1','deallocating bl',ABS(ierr_))

   RETURN
END SUBROUTINE dmat_sv_1


!**********************************************************
   SUBROUTINE dmat_lsd(m, n, nrhs, a, b, rcond, sv, rank, ierr)
   !**********************************************************
   !
   !  Computes the minimum-norm solution of the real system of linear equations
   !     A * X = B,
   !  where A is an M-by-N matrix and X and B are N-by-NRHS matrices.
   !
   !  The solutions found are:   min 2-norm| B-AX |
   !
   !  The SVD method is used, see the manual of DGELSD
   !  Singular values and rank can be provided in output
   !
   !  Only square matrices are treated
   !
   IMPLICIT NONE
   INTEGER,   INTENT(IN)       :: m, n, nrhs
   REAL(dbl), INTENT(IN)       :: rcond
   REAL(dbl), INTENT(IN)       :: a(:,:)
   REAL(dbl), INTENT(INOUT)    :: b(:,:)
   REAL(dbl), OPTIONAL, INTENT(OUT) :: sv(:)
   INTEGER,   OPTIONAL, INTENT(OUT) :: rank
   INTEGER,   OPTIONAL, INTENT(OUT) :: ierr
   !
   CHARACTER(8)              :: subname="dmat_lsd"
   INTEGER                   :: ierr_, rank_, info
   INTEGER,      ALLOCATABLE :: iwork(:)
   REAL(dbl),    ALLOCATABLE :: atmp(:,:), work(:), sv_(:)
   !
   REAL(dbl)         :: workl
   INTEGER           :: smlsiz, minmn, nlvl, lwork, liwork
   INTEGER, EXTERNAL :: ILAENV



   IF ( PRESENT(ierr) ) ierr=0
   !
   IF ( m > SIZE(a,1) .OR. n > SIZE(a,2) ) CALL errore(subname,'matrix A too small',1)
   IF ( n > SIZE(b,1) ) CALL errore(subname,'matrix B too small (I)',1)
   IF ( nrhs > SIZE(b,2) ) CALL errore(subname,'matrix B too small (II)',1)

   !
   ! dimensions
   !
   minmn  = MIN(m,n)
   smlsiz = ILAENV( 9, 'DGELSD', ' ', 0, 0, 0, 0 )
   nlvl   = MAX( INT( LOG( DBLE( minmn ) / DBLE( smlsiz+1 ) ) /  &
            LOG( 2.0d0 ) ) + 1, 0 )
   !
   liwork = MAX(1, 3 * minmn * nlvl + 11 * minmn)
   !

   !
   ! workspace
   !
   ALLOCATE( atmp(m,n), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating atmp',ABS(ierr_))
   ALLOCATE( sv_(MIN(m,n)), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating sv_',ABS(ierr_))
   ALLOCATE( iwork(liwork), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating iwork',ABS(ierr_))
   !
   ! get lwork
   !
   info   =  0
   lwork  = -1
   !
   ! use magma if possible
   ! get proper dimensions here
   CALL DGELSD( m, n, nrhs, atmp, m, b, SIZE(b,1), sv_, rcond, rank_, &
                workl, lwork, iwork, info)
   !
   lwork  = INT( workl )
   !
   ALLOCATE( work(lwork), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating work',ABS(ierr_))

   !
   ! make a local copy of a
   atmp(:,:) = a(1:m,1:n) 

   !
   ! use magma if possible
   CALL DGELSD( m, n, nrhs, atmp, m, b, SIZE(b,1), sv_, rcond, rank_, &
                work, lwork, iwork, info)

   IF ( PRESENT(ierr) ) THEN
        IF (info/=0) ierr= info
   ELSE
        IF ( info < 0 ) CALL errore(subname, 'DGELSD: info illegal value', -info )
        IF ( info > 0 ) CALL errore(subname, 'DGELSD: algorithm failure', info )
   ENDIF
   !
   IF ( PRESENT(sv) )     sv(1:minmn) = sv_(1:minmn)
   IF ( PRESENT(rank) )   rank        = rank_
    

   DEALLOCATE( atmp, STAT=ierr_)
   IF(ierr_/=0) CALL errore(subname,'deallocating atmp',ABS(ierr_))
   DEALLOCATE( sv_, STAT=ierr_)
   IF(ierr_/=0) CALL errore(subname,'deallocating atmp',ABS(ierr_))
   DEALLOCATE( work, iwork, STAT=ierr_)
   IF(ierr_/=0) CALL errore(subname,'deallocating work, iwork',ABS(ierr_))

   RETURN
END SUBROUTINE dmat_lsd


!**********************************************************
   SUBROUTINE dmat_lsd_1(m, n, nrhs, a, b, rcond, sv, rank, ierr)
   !**********************************************************
   !
   ! Interface to dmat_lsd when nrhs = 1
   !
   IMPLICIT NONE
   INTEGER,   INTENT(IN)       :: m, n, nrhs
   REAL(dbl), INTENT(IN)       :: rcond
   REAL(dbl), INTENT(IN)       :: a(:,:)
   REAL(dbl), INTENT(INOUT)    :: b(:)
   REAL(dbl), OPTIONAL, INTENT(OUT) :: sv(:)
   INTEGER,   OPTIONAL, INTENT(OUT) :: rank
   INTEGER,   OPTIONAL, INTENT(OUT) :: ierr
   !
   CHARACTER(10)             :: subname="dmat_lsd_1"
   INTEGER                   :: ierr_, rank_, info
   REAL(dbl),    ALLOCATABLE :: bl(:,:)
   REAL(dbl),    ALLOCATABLE :: sv_(:)

   IF ( PRESENT(ierr) ) ierr=0
   !
   IF ( nrhs /= 1) CALL errore(subname,'more than 1 rhs ? ',ABS(nrhs)+1)
   IF ( n > SIZE(b,1) ) CALL errore(subname,'vector b too small',2)
   !
   ALLOCATE( bl(n,1), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating bl',ABS(ierr_))
   ALLOCATE( sv_(MIN(m,n)), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating sv_',ABS(ierr_))

   bl(:,1) = b(1:n)
   CALL dmat_lsd( m, n, 1, a, bl, rcond, SV=sv_, RANK=rank_, IERR=info)
   !
   IF ( PRESENT(ierr) ) THEN 
       ierr=info
   ELSE
       IF (info /=0) CALL errore(subname,'info/=0',ABS(info))
   ENDIF
   IF ( PRESENT(sv) )    sv(1:MIN(m,n))  = sv_(:)
   IF ( PRESENT(rank) )  rank            = rank_
   !
   b(1:n) = bl(1:n,1)

   DEALLOCATE( bl, STAT=ierr_)
   IF (ierr_/=0) CALL errore(subname,'deallocating bl',ABS(ierr_))
   DEALLOCATE( sv_, STAT=ierr_)
   IF (ierr_/=0) CALL errore(subname,'deallocating sv_',ABS(ierr_))

   RETURN
END SUBROUTINE dmat_lsd_1


!**********************************************************
   SUBROUTINE zmat_lsd(m, n, nrhs, a, b, rcond, sv, rank, ierr)
   !**********************************************************
   !
   !  Computes the minimum-norm solution of the real system of linear equations
   !     A * X = B,
   !  where A is an M-by-N matrix and X and B are N-by-NRHS matrices.
   !
   !  The solutions found are:   min 2-norm| B-AX |
   !
   !  The SVD method is used, see the manual of ZGELSD
   !  Singular values and rank can be provided in output
   !
   !  Only square matrices are treated
   !
   IMPLICIT NONE
   INTEGER,      INTENT(IN)       :: m, n, nrhs
   REAL(dbl),    INTENT(IN)       :: rcond
   COMPLEX(dbl), INTENT(IN)       :: a(:,:)
   COMPLEX(dbl), INTENT(INOUT)    :: b(:,:)
   REAL(dbl),    OPTIONAL, INTENT(OUT) :: sv(:)
   INTEGER,      OPTIONAL, INTENT(OUT) :: rank
   INTEGER,      OPTIONAL, INTENT(OUT) :: ierr
   !
   CHARACTER(8)              :: subname="zmat_lsd"
   INTEGER                   :: ierr_, rank_, info
   INTEGER,      ALLOCATABLE :: iwork(:)
   REAL(dbl),    ALLOCATABLE :: rwork(:), sv_(:)
   COMPLEX(dbl), ALLOCATABLE :: atmp(:,:), work(:)
   !
   REAL(dbl)         :: rworkl
   COMPLEX(dbl)      :: workl
   INTEGER           :: smlsiz, minmn, nlvl, lwork, lrwork, liwork
   INTEGER, EXTERNAL :: ILAENV



   IF ( PRESENT(ierr) ) ierr=0
   !
   IF ( m > SIZE(a,1) .OR. n > SIZE(a,2) ) CALL errore(subname,'matrix A too small',1)
   IF ( n > SIZE(b,1) ) CALL errore(subname,'matrix B too small (I)',1)
   IF ( nrhs > SIZE(b,2) ) CALL errore(subname,'matrix B too small (II)',1)

   !
   ! dimensions
   !
   minmn  = MIN(m,n)
   smlsiz = ILAENV( 9, 'ZGELSD', ' ', 0, 0, 0, 0 )
   nlvl   = MAX( INT( LOG( DBLE( minmn ) / DBLE( smlsiz+1 ) ) /  &
            LOG( 2.0d0 ) ) + 1, 0 )
   !
   !mnthr  = ILAENV( 6, 'ZGELSD', ' ', m, n, nrhs, -1 )
   liwork = MAX(1, 3 * minmn * nlvl + 11 * minmn)
   !

   !
   ! workspace
   !
   ALLOCATE( atmp(m,n), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating atmp',ABS(ierr_))
   ALLOCATE( sv_(MIN(m,n)), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating sv_',ABS(ierr_))
   ALLOCATE( iwork(liwork), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating iwork',ABS(ierr_))
   !
   ! get lwork
   !
   info   =  0
   lwork  = -1
   !
   ! use magma if possible
   ! get proper dimensions here
   CALL ZGELSD( m, n, nrhs, atmp, m, b, SIZE(b,1), sv_, rcond, rank_, &
                workl, lwork, rworkl, iwork, info)
   !
   lwork  = INT( workl )
   lrwork = INT( rworkl )
   !
   ALLOCATE( work(lwork), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating work',ABS(ierr_))
   ALLOCATE( rwork(lrwork), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating rwork',ABS(ierr_))

   !
   ! make a local copy of a
   atmp(:,:) = a(1:m,1:n) 

   !
   ! use magma if possible
   CALL ZGELSD( m, n, nrhs, atmp, m, b, SIZE(b,1), sv_, rcond, rank_, &
                work, lwork, rwork, iwork, info)

   IF ( PRESENT(ierr) ) THEN
        IF (info/=0) ierr= info
   ELSE
        IF ( info < 0 ) CALL errore(subname, 'ZGELSD: info illegal value', -info )
        IF ( info > 0 ) CALL errore(subname, 'ZGELSD: algorithm failure', info )
   ENDIF
   !
   IF ( PRESENT(sv) )     sv(1:minmn) = sv_(1:minmn)
   IF ( PRESENT(rank) )   rank        = rank_
    

   DEALLOCATE( atmp, STAT=ierr_)
   IF(ierr_/=0) CALL errore(subname,'deallocating atmp',ABS(ierr_))
   DEALLOCATE( sv_, STAT=ierr_)
   IF(ierr_/=0) CALL errore(subname,'deallocating atmp',ABS(ierr_))
   DEALLOCATE( work, rwork, iwork, STAT=ierr_)
   IF(ierr_/=0) CALL errore(subname,'deallocating work, rwork, iwork',ABS(ierr_))

   RETURN
END SUBROUTINE zmat_lsd


!**********************************************************
   SUBROUTINE zmat_lsd_1(m, n, nrhs, a, b, rcond, sv, rank, ierr)
   !**********************************************************
   !
   ! Interface to zmat_lsd when nrhs = 1
   !
   IMPLICIT NONE
   INTEGER,   INTENT(IN)       :: m, n, nrhs
   REAL(dbl), INTENT(IN)       :: rcond
   COMPLEX(dbl), INTENT(IN)    :: a(:,:)
   COMPLEX(dbl), INTENT(INOUT) :: b(:)
   REAL(dbl), OPTIONAL, INTENT(OUT) :: sv(:)
   INTEGER,   OPTIONAL, INTENT(OUT) :: rank
   INTEGER,   OPTIONAL, INTENT(OUT) :: ierr
   !
   CHARACTER(10)             :: subname="dmat_lsd_1"
   INTEGER                   :: ierr_, rank_, info
   COMPLEX(dbl), ALLOCATABLE :: bl(:,:)
   REAL(dbl),    ALLOCATABLE :: sv_(:)

   IF ( PRESENT(ierr) ) ierr=0
   !
   IF ( nrhs /= 1) CALL errore(subname,'more than 1 rhs ? ',ABS(nrhs)+1)
   IF ( n > SIZE(b,1) ) CALL errore(subname,'vector b too small',2)
   !
   ALLOCATE( bl(n,1), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating bl',ABS(ierr_))
   ALLOCATE( sv_(MIN(m,n)), STAT=ierr_ )
   IF (ierr_/=0) CALL errore(subname,'allocating sv_',ABS(ierr_))

   bl(:,1) = b(1:n)
   CALL zmat_lsd( m, n, 1, a, bl, rcond, SV=sv_, RANK=rank_, IERR=info)
   !
   IF ( PRESENT(ierr) ) THEN 
       ierr=info
   ELSE
       IF (info /=0) CALL errore(subname,'info/=0',ABS(info))
   ENDIF
   IF ( PRESENT(sv) )    sv(1:MIN(m,n))  = sv_(:)
   IF ( PRESENT(rank) )  rank            = rank_
   !
   b(1:n) = bl(1:n,1)

   DEALLOCATE( bl, STAT=ierr_)
   IF (ierr_/=0) CALL errore(subname,'deallocating bl',ABS(ierr_))
   DEALLOCATE( sv_, STAT=ierr_)
   IF (ierr_/=0) CALL errore(subname,'deallocating sv_',ABS(ierr_))

   RETURN
END SUBROUTINE zmat_lsd_1


!**********************************************************
   SUBROUTINE zmat_sv(n, nrhs, a, b, ierr)
   !**********************************************************
   !
   !  Computes the solution of the complex system of linear equations
   !     A * X = B,
   !  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
   !
   !  The LU decomposition with partial pivoting and row interchanges is
   !  used to factor A as
   !     A = P * L * U,
   !  where P is a permutation matrix, L is unit lower triangular, and U is
   !  upper triangular.  The factored form of A is then used to solve the
   !  system of equations A * X = B.
   !
   IMPLICIT NONE
   INTEGER, INTENT(IN)         :: n, nrhs
   COMPLEX(dbl), INTENT(IN)    :: a(:,:)
   COMPLEX(dbl), INTENT(INOUT) :: b(:,:)
   INTEGER, OPTIONAL, INTENT(out) :: ierr

   INTEGER :: ierr_, info
   INTEGER,      ALLOCATABLE :: ipiv(:)
   COMPLEX(dbl), ALLOCATABLE :: atmp(:,:)

   IF ( PRESENT(ierr) ) ierr=0
   !
   IF ( n > SIZE(a,1) .OR. n > SIZE(a,2) ) CALL errore('zmat_sv','matrix A too small',1)
   IF ( n > SIZE(b,1) ) CALL errore('zmat_sv','matrix B too small (I)',1)
   IF ( nrhs > SIZE(b,2) ) CALL errore('zmat_sv','matrix B too small (II)',1)

   ALLOCATE( atmp(n,n), ipiv(n), STAT=ierr_ )
     IF (ierr_/=0) CALL errore('zmat_sv','allocating atmp, ipiv',ABS(ierr_))
   !
   ! make a local copy of a
   atmp(:,:) = a(1:n,1:n) 

   ! use magma if possible
   CALL ZGESV( n, nrhs, atmp, n, ipiv, b, SIZE(b,1), info)

   IF ( PRESENT(ierr) ) THEN
        IF (info/=0) ierr= info
   ELSE
        IF ( info < 0 ) CALL errore('zmat_sv', 'ZGESV: info illegal value', -info )
        IF ( info > 0 ) CALL errore('zmat_sv', 'ZGESV: singular matrix', info )
   ENDIF

   DEALLOCATE( atmp, ipiv, STAT=ierr_)
      IF(ierr_/=0) CALL errore('zmat_sv','deallocating atmp, ipiv',ABS(ierr_))

   RETURN
END SUBROUTINE zmat_sv


!**********************************************************
   SUBROUTINE zmat_sv_1(n, nrhs, a, b, ierr)
   !**********************************************************
   !
   ! Interface to zmat_sv when nrhs = 1
   !
   IMPLICIT NONE
   INTEGER, INTENT(IN)         :: n, nrhs
   COMPLEX(dbl), INTENT(IN)    :: a(:,:)
   COMPLEX(dbl), INTENT(INOUT) :: b(:)
   INTEGER, OPTIONAL, INTENT(out) :: ierr
   !
   INTEGER  :: ierr_, info
   COMPLEX(dbl),   ALLOCATABLE :: bl(:,:)

   IF ( PRESENT(ierr) ) ierr = 0
   !
   IF ( nrhs /= 1) CALL errore('zmat_sv_1','more than 1 rhs ?',ABS(nrhs)+1)
   IF ( n > SIZE(b,1) ) CALL errore('zmat_sv_1','vector b too small',2)
   !
   ALLOCATE( bl(n,1), STAT=ierr_ )
   IF (ierr_/=0) CALL errore('zmat_sv_1','allocating bl',ABS(ierr_))

   bl(:,1) = b(1:n)
   CALL zmat_sv( n, 1, a, bl, IERR=info)
   !
   IF ( PRESENT(ierr) ) THEN 
       ierr=info
   ELSE
       IF ( info /=0 ) CALL errore('zmat_sv_1','info/=0',ABS(info))
   ENDIF
   !
   b(1:n) = bl(1:n,1)

   DEALLOCATE( bl, STAT=ierr_)
   IF (ierr_/=0) CALL errore('zmat_sv_1','deallocating bl',ABS(ierr_))

   RETURN
END SUBROUTINE zmat_sv_1


!**********************************************************
   SUBROUTINE zmat_mul( c, a, opa, b, opb, m, n, k )
   !**********************************************************
   IMPLICIT NONE
   COMPLEX(dbl), INTENT(IN)  :: a(:,:)
   COMPLEX(dbl), INTENT(IN)  :: b(:,:)
   COMPLEX(dbl), INTENT(OUT) :: c(:,:)
   CHARACTER,    INTENT(IN)  :: opa, opb
   INTEGER,      INTENT(IN)  :: m, n, k
   !
   INTEGER :: i, j, l, ierr
   COMPLEX(dbl), ALLOCATABLE :: c_(:,:)
   !
   ! According to BLAS convention:
   ! C = opa(A) * opb(B)      op* = 'N' normal, 'C' complx conjg (i.e. herm conjg)
   !                          C is m*n,   opa(A) is m*k, opb(B) = k*n
   !
   IF ( m <= 0 .OR. n<=0 .OR. k<=0) CALL errore('zmat_mul','Invalid DIM',1)
   IF( opb /= 'N' .AND. opb /= 'C' ) &
     CALL errore('zmat_mul','argument value not allowed', 5 )

   ALLOCATE( c_( SIZE(c,1),SIZE(c,2) ), STAT=ierr )
   IF ( ierr/=0 ) CALL errore('zmat_mul','allocating c_',ABS(ierr))

   !
   ! this filter value has to be checked. Here we try to use BLAS
   ! as much as we can.
   !
   IF( k <= 3 ) THEN
       IF( ( opb == 'N' ) .AND. ( opa == 'N' ) ) THEN
           !
           IF ( m > SIZE(c,1) .OR. m > SIZE(a,1) ) CALL  errore('zmat_mul','Invalid C,A',m)
           IF ( n > SIZE(c,2) .OR. n > SIZE(b,2) ) CALL  errore('zmat_mul','Invalid C,B',n)
           IF ( k > SIZE(a,2) .OR. k > SIZE(b,1) ) CALL  errore('zmat_mul','Invalid A,B',k)
           !
           DO j = 1, n
           DO i = 1, m
               c_(i,j) = CZERO
               DO l = 1, k
                   c_(i,j) = c_(i,j) + a(i,l) * b(l,j)
               ENDDO
           ENDDO
           ENDDO
       ELSE IF( ( opa == 'N' ) .AND. ( opb == 'C' ) ) THEN
           !
           IF ( m > SIZE(c,1) .OR. m > SIZE(a,1) ) CALL  errore('zmat_mul','Invalid C,A',m)
           IF ( n > SIZE(c,2) .OR. n > SIZE(b,1) ) CALL  errore('zmat_mul','Invalid C,B',n)
           IF ( k > SIZE(a,2) .OR. k > SIZE(b,2) ) CALL  errore('zmat_mul','Invalid A,B',k)
           !
           DO j = 1, n
           DO i = 1, m
               c_(i,j) = CZERO
               DO l = 1, k
                   c_(i,j) = c_(i,j) + a(i,l) * CONJG( b(j,l) )
               ENDDO
           ENDDO
           ENDDO
       ELSE IF( ( opa == 'C' ) .AND. ( opb == 'N' ) ) THEN
           !
           IF ( m > SIZE(c,1) .OR. m > SIZE(a,2) ) CALL  errore('zmat_mul','Invalid C,A',m)
           IF ( n > SIZE(c,2) .OR. n > SIZE(b,2) ) CALL  errore('zmat_mul','Invalid C,B',n)
           IF ( k > SIZE(a,1) .OR. k > SIZE(b,1) ) CALL  errore('zmat_mul','Invalid A,B',k)
           !
           DO j = 1, n
           DO i = 1, m
               c_(i,j) = CZERO
               DO l = 1, k
                   c_(i,j) = c_(i,j) + CONJG( a(l,i) )* b(l,j) 
               ENDDO
           ENDDO
           ENDDO
       ELSE IF( ( opb == 'C' ) .AND. ( opa == 'C' ) ) THEN
           !
           IF ( m > SIZE(c,1) .OR. m > SIZE(a,2) ) CALL  errore('zmat_mul','Invalid C,A',m)
           IF ( n > SIZE(c,2) .OR. n > SIZE(b,1) ) CALL  errore('zmat_mul','Invalid C,B',n)
           IF ( k > SIZE(a,1) .OR. k > SIZE(b,2) ) CALL  errore('zmat_mul','Invalid A,B',k)
           !
           DO j = 1, n
           DO i = 1, m
               c_(i,j) = CZERO
               DO l = 1, k
                   c_(i,j) = c_(i,j) + CONJG( a(l,i) )* CONJG( b(j,l) )
               ENDDO
           ENDDO
           ENDDO
       ENDIF
   ELSE

       CALL ZGEMM( opa, opb, m, n, k, CONE, a, SIZE(a,1), &
                   b, SIZE(b,1), CZERO, c_, SIZE(c,1) )
   ENDIF
   !
   c(:,:) = c_(:,:)
   !
   DEALLOCATE( c_, STAT=ierr)
   IF ( ierr/=0 ) CALL errore('zmat_mul','deallocating c_', ABS(ierr))
   !
END SUBROUTINE zmat_mul


!**********************************************************
   SUBROUTINE dmat_mul( c, a, opa, b, opb, m, n, k )
   !**********************************************************
   IMPLICIT NONE
   REAL(dbl), INTENT(IN)  :: a(:,:)
   REAL(dbl), INTENT(IN)  :: b(:,:)
   REAL(dbl), INTENT(OUT) :: c(:,:)
   CHARACTER, INTENT(IN)  :: opa, opb
   INTEGER,   INTENT(IN)  :: m, n, k
   !
   INTEGER :: i, j, l, ierr
   REAL(dbl), ALLOCATABLE :: c_(:,:)
   !
   ! According to BLAS convention:
   ! C = opa(A) * opb(B)      op* = 'N' normal, 'T' transpose 
   !                          C is m*n,   opa(A) is m*k, opb(B) = k*n
   !
   IF ( m <= 0 .OR. n<=0 .OR. k<=0) CALL errore('dmat_mul','Invalid DIM',1)
   IF( opb /= 'N' .AND. opb /= 'T' ) &
     CALL errore('dmat_mul','argument value not allowed ', 5 )

   ALLOCATE( c_( SIZE(c,1),SIZE(c,2) ), STAT=ierr )
   IF ( ierr/=0 ) CALL errore('dmat_mul','allocating c_',ABS(ierr))

   !
   ! this filter value has to be checked. Here we try to use BLAS
   ! as much as we can.
   !
   IF( k <= 3 ) THEN
       IF( ( opb == 'N' ) .AND. ( opa == 'N' ) ) THEN
           !
           IF ( m > SIZE(c,1) .OR. m > SIZE(a,1) ) CALL  errore('dmat_mul','Invalid C,A',m)
           IF ( n > SIZE(c,2) .OR. n > SIZE(b,2) ) CALL  errore('dmat_mul','Invalid C,B',n)
           IF ( k > SIZE(a,2) .OR. k > SIZE(b,1) ) CALL  errore('dmat_mul','Invalid A,B',k)
           !
           DO j = 1, n
           DO i = 1, m
               c_(i,j) = ZERO
               DO l = 1, k
                   c_(i,j) = c_(i,j) + a(i,l) * b(l,j)
               ENDDO
           ENDDO
           ENDDO
       ELSE IF( ( opa == 'N' ) .AND. ( opb == 'T' ) ) THEN
           !
           IF ( m > SIZE(c,1) .OR. m > SIZE(a,1) ) CALL  errore('dmat_mul','Invalid C,A',m)
           IF ( n > SIZE(c,2) .OR. n > SIZE(b,1) ) CALL  errore('dmat_mul','Invalid C,B',n)
           IF ( k > SIZE(a,2) .OR. k > SIZE(b,2) ) CALL  errore('dmat_mul','Invalid A,B',k)
           !
           DO j = 1, n
           DO i = 1, m
               c_(i,j) = ZERO
               DO l = 1, k
                   c_(i,j) = c_(i,j) + a(i,l) * b(j,l)
               ENDDO
           ENDDO
           ENDDO
       ELSE IF( ( opa == 'T' ) .AND. ( opb == 'N' ) ) THEN
           !
           IF ( m > SIZE(c,1) .OR. m > SIZE(a,2) ) CALL  errore('dmat_mul','Invalid C,A',m)
           IF ( n > SIZE(c,2) .OR. n > SIZE(b,2) ) CALL  errore('dmat_mul','Invalid C,B',n)
           IF ( k > SIZE(a,1) .OR. k > SIZE(b,1) ) CALL  errore('dmat_mul','Invalid A,B',k)
           !
           DO j = 1, n
           DO i = 1, m
               c_(i,j) = ZERO
               DO l = 1, k
                   c_(i,j) = c_(i,j) + a(l,i) * b(l,j) 
               ENDDO
           ENDDO
           ENDDO
       ELSE IF( ( opb == 'T' ) .AND. ( opa == 'T' ) ) THEN
           !
           IF ( m > SIZE(c,1) .OR. m > SIZE(a,2) ) CALL  errore('dmat_mul','Invalid C,A',m)
           IF ( n > SIZE(c,2) .OR. n > SIZE(b,1) ) CALL  errore('dmat_mul','Invalid C,B',n)
           IF ( k > SIZE(a,1) .OR. k > SIZE(b,2) ) CALL  errore('dmat_mul','Invalid A,B',k)
           !
           DO j = 1, n
           DO i = 1, m
               c_(i,j) = ZERO
               DO l = 1, k
                   c_(i,j) = c_(i,j) + a(l,i) * b(j,l) 
               ENDDO
           ENDDO
           ENDDO
       ENDIF
   ELSE

       CALL DGEMM( opa, opb, m, n, k, ONE, a, SIZE(a,1), &
                   b, SIZE(b,1), ZERO, c_, SIZE(c_,1) )
   ENDIF
   !
   c(:,:) = c_(:,:)
   !
   DEALLOCATE( c_, STAT=ierr)
   IF ( ierr/=0 ) CALL errore('dmat_mul','deallocating c_', ABS(ierr))
   !
END SUBROUTINE dmat_mul


!**********************************************************
   SUBROUTINE dmat_vec_mul( c, a, opa, b, m, n )
   !**********************************************************
   IMPLICIT NONE
   REAL(dbl), INTENT(IN)  :: a(:,:)
   REAL(dbl), INTENT(IN)  :: b(:)
   REAL(dbl), INTENT(OUT) :: c(:)
   CHARACTER, INTENT(IN)  :: opa
   INTEGER,   INTENT(IN)  :: m, n
   !
   INTEGER :: i, j, l, ierr
   REAL(dbl), ALLOCATABLE :: c_(:)
   !
   ! According to BLAS convention:
   ! C = opa(A) * opb(B)      op* = 'N' normal, 'T' transpose 
   !                          C is m*n,   opa(A) is m*k, opb(B) = k*n
   !
   IF ( m <= 0 .OR. n<=0 ) CALL errore('dmat_vec_mul','Invalid DIM',1)
   IF( opa /= 'N' .AND. opa /= 'T' ) &
     CALL errore('dmat_vec_mul','argument value not allowed ', 5 )

   ALLOCATE( c_( SIZE(c) ), STAT=ierr )
   IF ( ierr/=0 ) CALL errore('dmat_vec_mul','allocating c_',ABS(ierr))

   CALL DGEMV( opa, m, n, ONE, a, SIZE(a,1), b, 1, ZERO, c_, 1 )
   !
   c = c_
   !
   DEALLOCATE( c_, STAT=ierr)
   IF ( ierr/=0 ) CALL errore('dmat_vec_mul','deallocating c_', ABS(ierr))
   !
END SUBROUTINE dmat_vec_mul


!**********************************************************
   SUBROUTINE zmat_vec_mul( c, a, opa, b, m, n )
   !**********************************************************
   IMPLICIT NONE
   COMPLEX(dbl), INTENT(IN)  :: a(:,:)
   COMPLEX(dbl), INTENT(IN)  :: b(:)
   COMPLEX(dbl), INTENT(OUT) :: c(:)
   CHARACTER, INTENT(IN)  :: opa
   INTEGER,   INTENT(IN)  :: m, n
   !
   INTEGER :: i, j, l, ierr
   COMPLEX(dbl), ALLOCATABLE :: c_(:)
   !
   ! According to BLAS convention:
   ! C = opa(A) * opb(B)      op* = 'N' normal, 'T' transpose 
   !                          C is m*n,   opa(A) is m*k, opb(B) = k*n
   !
   IF ( m <= 0 .OR. n<=0 ) CALL errore('dmat_vec_mul','Invalid DIM',1)
   IF( opa /= 'N' .AND. opa /= 'T' ) &
     CALL errore('dmat_vec_mul','argument value not allowed ', 5 )

   ALLOCATE( c_( SIZE(c) ), STAT=ierr )
   IF ( ierr/=0 ) CALL errore('dmat_vec_mul','allocating c_',ABS(ierr))

   CALL ZGEMV( opa, m, n, CONE, a, SIZE(a,1), b, 1, CZERO, c_, 1 )
   !
   c = c_
   !
   DEALLOCATE( c_, STAT=ierr)
   IF ( ierr/=0 ) CALL errore('dmat_vec_mul','deallocating c_', ABS(ierr))
   !
END SUBROUTINE zmat_vec_mul


!**********************************************************
   SUBROUTINE zmat_hdiag( z, w, a, n )
   !**********************************************************
   !
   ! utility to diagonalize complex hermitean matrices
   !
   IMPLICIT NONE
   COMPLEX(dbl), INTENT(IN)  :: a(:,:)
   COMPLEX(dbl), INTENT(OUT) :: z(:,:)
   REAL(dbl),    INTENT(OUT) :: w(:)
   INTEGER,      INTENT(in)  :: n

   INTEGER :: i, j, ierr, info
   COMPLEX(dbl), ALLOCATABLE :: ap(:)
   COMPLEX(dbl), ALLOCATABLE :: work(:)
   REAL(dbl), ALLOCATABLE :: rwork(:)
   INTEGER, ALLOCATABLE :: ifail(:)
   INTEGER, ALLOCATABLE :: iwork(:)

   ! get the dimension of the problem
   IF ( n <= 0 ) CALL errore('zmat_hdiag','Invalid N',ABS(n)+1)
   IF ( n > SIZE(a,1) .OR. n > SIZE(a,2) ) &
        CALL errore('zmat_hdiag','Invalid A dimensions',ABS(n)+1)
   IF ( n > SIZE(z,1) .OR. n > SIZE(z,2) ) &
        CALL errore('zmat_hdiag','Invalid Z dimensions',ABS(n)+1)
   
   ALLOCATE( ap(n*(n+1)/2), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating ap',ABS(ierr))
   ALLOCATE( work(2*n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating work',ABS(ierr))
   ALLOCATE( rwork(7*n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating rwork',ABS(ierr))
   ALLOCATE( ifail(n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating ifail',ABS(ierr))
   ALLOCATE( iwork(5*n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating iwork',ABS(ierr))

   DO j = 1, n
   DO i = 1, j
      ap(i + ( (j-1)*j)/2 ) = a(i,j)
   ENDDO
   ENDDO

   ! use magma if possible
   CALL ZHPEVX( 'v', 'a', 'u', n, ap, ZERO, ZERO, 0, 0, -ONE, i, w, &
                 z, SIZE(z,1), work, rwork, iwork, ifail, info )

   IF ( info < 0 ) CALL errore('zmat_hdiag', 'zhpevx: info illegal value', -info )
   IF ( info > 0 ) &
        CALL errore('zmat_hdiag', 'zhpevx: eigenvectors not converged', info )
    
   DEALLOCATE( ap, work, rwork, iwork, ifail, STAT=ierr)
      IF(ierr/=0) CALL errore('zmat_hdiag','deallocating ap...ifail',ABS(ierr))

   RETURN
END SUBROUTINE zmat_hdiag


!**********************************************************
   SUBROUTINE zmat_hdiag_pack( z, w, ap, n, uplo )
   !**********************************************************
   !
   ! utility to diagonalize complex hermitean matrices, in packed form
   !
   IMPLICIT NONE
   COMPLEX(dbl), INTENT(IN)  :: ap(:)
   COMPLEX(dbl), INTENT(OUT) :: z(:,:)
   REAL(dbl),    INTENT(OUT) :: w(:)
   INTEGER,      INTENT(IN)  :: n
   CHARACTER,    INTENT(IN)  :: uplo

   INTEGER :: i, ierr, info
   COMPLEX(dbl), ALLOCATABLE :: work(:)
   REAL(dbl), ALLOCATABLE :: rwork(:)
   INTEGER, ALLOCATABLE :: ifail(:)
   INTEGER, ALLOCATABLE :: iwork(:)

   ! get the dimension of the problem
   IF ( n <= 0 ) CALL errore('zmat_hdiag','Invalid N',ABS(n)+1)
   IF ( n*(n+1)/2 > SIZE(ap) .OR. n*(n+1)/2 > SIZE(ap) ) &
        CALL errore('zmat_hdiag','Invalid A dimensions',ABS(n)+1)
   IF ( n > SIZE(z,1) .OR. n > SIZE(z,2) ) &
        CALL errore('zmat_hdiag','Invalid Z dimensions',ABS(n)+1)
   
   ALLOCATE( work(2*n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating work',ABS(ierr))
   ALLOCATE( rwork(7*n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating rwork',ABS(ierr))
   ALLOCATE( ifail(n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating ifail',ABS(ierr))
   ALLOCATE( iwork(5*n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating iwork',ABS(ierr))

   ! use magma if possible
   CALL ZHPEVX( 'v', 'a', uplo, n, ap, ZERO, ZERO, 0, 0, -ONE, i, w, &
                 z, SIZE(z,1), work, rwork, iwork, ifail, info )

   IF ( info < 0 ) CALL errore('zmat_hdiag', 'zhpevx: info illegal value', -info )
   IF ( info > 0 ) &
        CALL errore('zmat_hdiag', 'zhpevx: eigenvectors not converged', info )
    
   DEALLOCATE( work, rwork, iwork, ifail, STAT=ierr)
      IF(ierr/=0) CALL errore('zmat_hdiag','deallocating work...ifail',ABS(ierr))

   RETURN
END SUBROUTINE zmat_hdiag_pack


!**********************************************************
   SUBROUTINE zmat_hdiag_gen( z, w, a, b, n )
   !**********************************************************
   !
   ! utility to solve the generalized eigenvalue problem with
   ! A,B hermitean matrices:
   !  A * z = w * B * z
   !
   IMPLICIT NONE
   COMPLEX(dbl), INTENT(IN)  :: a(:,:)
   COMPLEX(dbl), INTENT(IN)  :: b(:,:)
   COMPLEX(dbl), INTENT(OUT) :: z(:,:)
   REAL(dbl),    INTENT(OUT) :: w(:)
   INTEGER,      INTENT(in)  :: n

   INTEGER :: i, j, ierr, info
   COMPLEX(dbl), ALLOCATABLE :: ap(:), bp(:)
   COMPLEX(dbl), ALLOCATABLE :: work(:)
   REAL(dbl),    ALLOCATABLE :: rwork(:)
   INTEGER,      ALLOCATABLE :: ifail(:)
   INTEGER,      ALLOCATABLE :: iwork(:)

   ! get the dimension of the problem
   IF ( n <= 0 ) CALL errore('zmat_hdiag_gen','Invalid N',ABS(n)+1)
   IF ( n > SIZE(a,1) .OR. n > SIZE(a,2) ) &
        CALL errore('zmat_hdiag_gen','Invalid A dimensions',ABS(n)+1)
   IF ( n > SIZE(b,1) .OR. n > SIZE(b,2) ) &
        CALL errore('zmat_hdiag_gen','Invalid B dimensions',ABS(n)+1)
   IF ( n > SIZE(z,1) .OR. n > SIZE(z,2) ) &
        CALL errore('zmat_hdiag_gen','Invalid Z dimensions',ABS(n)+1)
   
   ALLOCATE( ap(n*(n+1)/2), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating ap',ABS(ierr))
   ALLOCATE( bp(n*(n+1)/2), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating bp',ABS(ierr))
   ALLOCATE( work(2*n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating work',ABS(ierr))
   ALLOCATE( rwork(7*n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating rwork',ABS(ierr))
   ALLOCATE( iwork(5*n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating iwork',ABS(ierr))
   ALLOCATE( ifail(n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_hdiag','allocating ifail',ABS(ierr))

   DO j = 1, n
   DO i = 1, j
      ap(i + ( (j-1)*j)/2 ) = a(i,j)
      bp(i + ( (j-1)*j)/2 ) = b(i,j)
   ENDDO
   ENDDO

   ! use magma if possible
   CALL ZHPGVX( 1, 'v', 'a', 'u', n, ap, bp, ZERO, ZERO, 0, 0, -ONE, i, w, &
                 z, SIZE(z,1), work, rwork, iwork, ifail, info )

   IF ( info < 0 ) CALL errore('zmat_hdiag', 'zhpevx: info illegal value', -info )
   IF ( info > 0 ) &
        CALL errore('zmat_hdiag', 'zhpevx: eigenvectors not converged', info )
    
   DEALLOCATE( ap, bp, work, rwork, iwork, ifail, STAT=ierr)
      IF(ierr/=0) CALL errore('zmat_hdiag','deallocating ap...ifail',ABS(ierr))

   RETURN
END SUBROUTINE zmat_hdiag_gen


!**********************************************************
   SUBROUTINE dmat_hdiag( z, w, a, n )
   !**********************************************************
   !
   ! utility to diagonalize real symmetric matrices
   !
   IMPLICIT NONE
   REAL(dbl),    INTENT(IN)  :: a(:,:)
   REAL(dbl),    INTENT(OUT) :: z(:,:)
   REAL(dbl),    INTENT(OUT) :: w(:)
   INTEGER,      INTENT(in)  :: n

   INTEGER :: i, j, ierr, info
   REAL(dbl), ALLOCATABLE :: ap(:)
   REAL(dbl), ALLOCATABLE :: work(:)
   INTEGER,   ALLOCATABLE :: ifail(:)
   INTEGER,   ALLOCATABLE :: iwork(:)

   ! get the dimension of the problem
   IF ( n <= 0 ) CALL errore('dmat_hdiag','Invalid N',ABS(n)+1)
   IF ( n > SIZE(a,1) .OR. n > SIZE(a,2) ) &
        CALL errore('dmat_hdiag','Invalid A dimensions',ABS(n)+1)
   IF ( n > SIZE(z,1) .OR. n > SIZE(z,2) ) &
        CALL errore('dmat_hdiag','Invalid Z dimensions',ABS(n)+1)
   
   ALLOCATE( ap(n*(n+1)/2), STAT=ierr )
   IF(ierr/=0) CALL errore('dmat_hdiag','allocating ap',ABS(ierr))
   ALLOCATE( work(8*n), STAT=ierr )
   IF(ierr/=0) CALL errore('dmat_hdiag','allocating work',ABS(ierr))
   ALLOCATE( ifail(n), STAT=ierr )
   IF(ierr/=0) CALL errore('dmat_hdiag','allocating ifail',ABS(ierr))
   ALLOCATE( iwork(5*n), STAT=ierr )
   IF(ierr/=0) CALL errore('dmat_hdiag','allocating iwork',ABS(ierr))

   DO j = 1, n
   DO i = 1, j
      ap(i + ( (j-1)*j)/2 ) = a(i,j)
   ENDDO
   ENDDO

   ! use magma if possible
   CALL DSPEVX( 'v', 'a', 'u', n, ap(1), ZERO, ZERO, 0, 0, -ONE, i, w(1), &
                 z(1,1), SIZE(z,1), work(1), iwork(1), ifail(1), info )

   IF ( info < 0 ) CALL errore('dmat_hdiag', 'zhpevx: info illegal value', -info )
   IF ( info > 0 ) &
        CALL errore('dmat_hdiag', 'zhpevx: eigenvectors not converged', info )
    
   DEALLOCATE( ap, work, iwork, ifail, STAT=ierr)
      IF(ierr/=0) CALL errore('dmat_hdiag','deallocating ap...ifail',ABS(ierr))

   RETURN
END SUBROUTINE dmat_hdiag


!**********************************************************
   SUBROUTINE dmat_hdiag_gen( z, w, a, b, n )
   !**********************************************************
   !
   ! utility to diagonalize real symmetric matrices
   ! in the general form A x = w B x
   !
   IMPLICIT NONE
   REAL(dbl),    INTENT(IN)  :: a(:,:), b(:,:)
   REAL(dbl),    INTENT(OUT) :: z(:,:)
   REAL(dbl),    INTENT(OUT) :: w(:)
   INTEGER,      INTENT(in)  :: n

   INTEGER :: i, j, ierr, info
   REAL(dbl), ALLOCATABLE :: ap(:), bp(:)
   REAL(dbl), ALLOCATABLE :: work(:)
   INTEGER,   ALLOCATABLE :: ifail(:)
   INTEGER,   ALLOCATABLE :: iwork(:)

   ! get the dimension of the problem
   IF ( n <= 0 ) CALL errore('dmat_hdiag','Invalid N',ABS(n)+1)
   IF ( n > SIZE(a,1) .OR. n > SIZE(a,2) ) &
        CALL errore('dmat_hdiag','Invalid A dimensions',ABS(n)+1)
   IF ( n > SIZE(z,1) .OR. n > SIZE(z,2) ) &
        CALL errore('dmat_hdiag','Invalid Z dimensions',ABS(n)+1)
   
   ALLOCATE( ap(n*(n+1)/2), STAT=ierr )
   IF(ierr/=0) CALL errore('dmat_hdiag','allocating ap',ABS(ierr))
   ALLOCATE( bp(n*(n+1)/2), STAT=ierr )
   IF(ierr/=0) CALL errore('dmat_hdiag','allocating bp',ABS(ierr))
   ALLOCATE( work(8*n), STAT=ierr )
   IF(ierr/=0) CALL errore('dmat_hdiag','allocating work',ABS(ierr))
   ALLOCATE( ifail(n), STAT=ierr )
   IF(ierr/=0) CALL errore('dmat_hdiag','allocating ifail',ABS(ierr))
   ALLOCATE( iwork(5*n), STAT=ierr )
   IF(ierr/=0) CALL errore('dmat_hdiag','allocating iwork',ABS(ierr))

   DO j = 1, n
   DO i = 1, j
      ap(i + ( (j-1)*j)/2 ) = a(i,j)
      bp(i + ( (j-1)*j)/2 ) = b(i,j)
   ENDDO
   ENDDO

   ! use magma if possible
   CALL DSPGVX( 1, 'v', 'a', 'u', n, ap, bp, ZERO, ZERO, 0, 0, -ONE, i, w, &
                 z, SIZE(z,1), work, iwork, ifail, info )

   IF ( info < 0 ) CALL errore('dmat_hdiag_gen', 'DSPGVX: info illegal value', -info )
   IF ( info > 0 ) CALL errore('dmat_hdiag_gen', 'DSPGVX: eigenvectors not converged', info )
    
   DEALLOCATE( ap, bp, work, iwork, ifail, STAT=ierr)
   IF(ierr/=0) CALL errore('dmat_hdiag_gen','deallocating ap...ifail',ABS(ierr))

   RETURN
END SUBROUTINE dmat_hdiag_gen


!**********************************************************
   SUBROUTINE dmat_bnd_hdiag( z, w, ab, ldab, n, kd )
   !**********************************************************
   !
   ! utility to diagonalize real symmetric matrices
   !
   IMPLICIT NONE
   REAL(dbl),    INTENT(IN)  :: ab(:,:)
   REAL(dbl),    INTENT(OUT) :: z(:,:)
   REAL(dbl),    INTENT(OUT) :: w(:)
   INTEGER,      INTENT(in)  :: ldab, n, kd

   CHARACTER(14) :: subname="dmat_bnd_hdiag"
   INTEGER :: i, j, m, ierr, info
   REAL(dbl), ALLOCATABLE :: work(:), qmat(:,:)
   INTEGER,   ALLOCATABLE :: ifail(:)
   INTEGER,   ALLOCATABLE :: iwork(:)

   ! get the dimension of the problem
   IF ( n <= 0 ) CALL errore(subname,'Invalid N',ABS(n)+1)
   !
   IF ( ldab < kd + 1)  CALL errore (subname, 'invalid ldab', 10)
   IF ( SIZE(z,1) < n)  CALL errore (subname, 'invalid ldz', 10)
   IF ( SIZE(z,2) < n)  CALL errore (subname, 'invalid ldz II', 11)
   !
   ALLOCATE( qmat(n,n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating qmat',ABS(ierr))
   ALLOCATE( work(7*n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating work',ABS(ierr))
   ALLOCATE( ifail(n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating ifail',ABS(ierr))
   ALLOCATE( iwork(5*n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating iwork',ABS(ierr))


   ! use magma if possible
   CALL DSBEVX( 'v', 'a', 'u', n, kd, ab, ldab, qmat, SIZE(qmat,1), ZERO, ZERO, 0, 0, ZERO, &
                 m, w, z, SIZE(z,1), work, iwork, ifail, info )

   IF ( info < 0 ) CALL errore(subname, 'DSBEVX: info illegal value', -info )
   IF ( info > 0 ) CALL errore(subname, 'DSBEVX: eigenvectors not converged', info )
    
   DEALLOCATE( qmat, work, iwork, ifail, STAT=ierr)
      IF(ierr/=0) CALL errore(subname,'deallocating qmat...ifail',ABS(ierr))

   RETURN
END SUBROUTINE dmat_bnd_hdiag


!**********************************************************
   SUBROUTINE dmat_bnd_hdiag_ext( z, w, ab, ldab, n, kd, il, iu )
   !**********************************************************
   !
   ! utility to diagonalize real symmetric matrices,
   ! expert algorithm
   !
   IMPLICIT NONE
   REAL(dbl),    INTENT(IN)  :: ab(:,:)
   REAL(dbl),    INTENT(OUT) :: z(:,:)
   REAL(dbl),    INTENT(OUT) :: w(:)
   INTEGER,      INTENT(in)  :: ldab, n, kd
   INTEGER,      INTENT(in)  :: il, iu

   CHARACTER(18) :: subname="dmat_bnd_hdiag_ext"
   INTEGER :: i, j, m, ierr, info
   REAL(dbl), ALLOCATABLE :: work(:), qmat(:,:)
   INTEGER,   ALLOCATABLE :: ifail(:)
   INTEGER,   ALLOCATABLE :: iwork(:)

   ! get the dimension of the problem
   IF ( n <= 0 ) CALL errore(subname,'Invalid N',ABS(n)+1)
   !
   IF ( ldab < kd + 1)  CALL errore (subname, 'invalid ldab', 10)
   IF ( SIZE(z,1) < n)  CALL errore (subname, 'invalid ldz I', 10)
   IF ( SIZE(z,2) < iu-il+1)  &
                        CALL errore (subname, 'invalid ldz II', 11)
   IF ( il <= 0 )       CALL errore (subname, 'invalid il', 12)
   IF ( iu <= 0 )       CALL errore (subname, 'invalid iu', 12)
   !
   ALLOCATE( qmat(n,n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating qmat',ABS(ierr))
   ALLOCATE( work(7*n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating work',ABS(ierr))
   ALLOCATE( ifail(n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating ifail',ABS(ierr))
   ALLOCATE( iwork(5*n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating iwork',ABS(ierr))


   ! use magma if possible

! XXX
!#if defined __OPENMP
!   CALL DSBEVX_omp( 'v', 'i', 'u', n, kd, ab, ldab, qmat, SIZE(qmat,1), 0.0d0, 0.0d0, il, iu, 0.0d0, &
!                    m, w, z, SIZE(z,1), work, iwork, ifail, info )
!#else
   CALL DSBEVX( 'v', 'i', 'u', n, kd, ab, ldab, qmat, SIZE(qmat,1), 0.0d0, 0.0d0, il, iu, 0.0d0, &
                 m, w, z, SIZE(z,1), work, iwork, ifail, info )
!#endif

   IF ( info < 0 ) CALL errore(subname, 'DSBEVX: info illegal value', -info )
   IF ( info > 0 ) CALL errore(subname, 'DSBEVX: eigenvectors not converged', info )
    
   DEALLOCATE( qmat, work, iwork, ifail, STAT=ierr)
      IF(ierr/=0) CALL errore(subname,'deallocating qmat...ifail',ABS(ierr))

   RETURN
END SUBROUTINE dmat_bnd_hdiag_ext


!**********************************************************
   SUBROUTINE dmat_bnd_hdiag_gen_ext( z, w, ab, ldab, bb, ldbb, n, ka, kb, il, iu )
   !**********************************************************
   !
   ! utility to diagonalize real symmetric matrices,
   ! expert algorithm, generalized eigevalue problem AZ = w BZ
   !
   ! According to the manual of DSBGVX:
   ! The eigenvectors are normalized so Z^T *B *Z = I.
   !
   IMPLICIT NONE
   REAL(dbl),    INTENT(IN)  :: ab(:,:)
   REAL(dbl),    INTENT(IN)  :: bb(:,:)
   REAL(dbl),    INTENT(OUT) :: z(:,:)
   REAL(dbl),    INTENT(OUT) :: w(:)
   INTEGER,      INTENT(in)  :: ldab, ldbb, n, ka, kb
   INTEGER,      INTENT(in)  :: il, iu

   CHARACTER(18) :: subname="dmat_bnd_hdiag_ext"
   INTEGER :: i, j, m, ierr, info
   REAL(dbl), ALLOCATABLE :: work(:), qmat(:,:)
   INTEGER,   ALLOCATABLE :: ifail(:)
   INTEGER,   ALLOCATABLE :: iwork(:)

   ! get the dimension of the problem
   IF ( n <= 0 ) CALL errore(subname,'Invalid N',ABS(n)+1)
   !
   IF ( ldab < ka + 1)  CALL errore (subname, 'invalid ldab', 10)
   IF ( ldbb < kb + 1)  CALL errore (subname, 'invalid ldbb', 10)
   IF ( SIZE(z,1) < n)  CALL errore (subname, 'invalid ldz I', 10)
   IF ( SIZE(z,2) < iu-il+1)  &
                        CALL errore (subname, 'invalid ldz II', 11)
   IF ( il <= 0 )       CALL errore (subname, 'invalid il', 12)
   IF ( iu <= 0 )       CALL errore (subname, 'invalid iu', 12)
   !
   ALLOCATE( qmat(n,n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating qmat',ABS(ierr))
   ALLOCATE( work(7*n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating work',ABS(ierr))
   ALLOCATE( ifail(n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating ifail',ABS(ierr))
   ALLOCATE( iwork(5*n), STAT=ierr )
   IF(ierr/=0) CALL errore(subname,'allocating iwork',ABS(ierr))

   CALL DSBGVX( 'v', 'i', 'u', n, ka, kb, ab, ldab, bb, ldbb, qmat, SIZE(qmat,1), &
                 0.0d0, 0.0d0, il, iu, 0.0d0, m, w, z, SIZE(z,1), work, iwork, ifail, info )

   IF ( info < 0 ) CALL errore(subname, 'DSBGVX: info illegal value', -info )
   IF ( info > 0 ) CALL errore(subname, 'DSBGVX: eigenvectors not converged', info )
    
   DEALLOCATE( qmat, work, iwork, ifail, STAT=ierr)
   IF(ierr/=0) CALL errore(subname,'deallocating qmat...ifail',ABS(ierr))

   RETURN
END SUBROUTINE dmat_bnd_hdiag_gen_ext


!**********************************************************
   SUBROUTINE zmat_diag_oneside( z, w, a, n, side )
   !**********************************************************
   !
   ! utility to diagonalize complex non-hermitean matrices
   !
   IMPLICIT NONE
   COMPLEX(dbl),        INTENT(IN)  :: a(:,:)
   COMPLEX(dbl),        INTENT(OUT) :: z(:,:)
   COMPLEX(dbl),        INTENT(OUT) :: w(:)
   INTEGER,             INTENT(in)  :: n
   CHARACTER,           INTENT(in)  :: side

   INTEGER   :: ierr, info, lwork
   CHARACTER :: jobvl, jobvr
   COMPLEX(dbl), ALLOCATABLE :: work(:), vl(:,:), vr(:,:)
   COMPLEX(dbl), ALLOCATABLE :: a_(:,:)
   REAL(dbl),    ALLOCATABLE :: rwork(:)

   ! get the dimension of the problem
   IF ( n <= 0 ) CALL errore('zmat_diag','Invalid N',ABS(n)+1)
   IF ( n > SIZE(a,1) .OR. n > SIZE(a,2) ) &
        CALL errore('zmat_diag','Invalid A dimensions',ABS(n)+1)
   IF ( n > SIZE(z,1) .OR. n > SIZE(z,2) ) &
        CALL errore('zmat_diag','Invalid Z dimensions',ABS(n)+1)

   SELECT CASE ( side )
   CASE ( 'L', 'l' )
       jobvl = 'V'
       jobvr = 'N'
   CASE ( 'R', 'r' )
       jobvl = 'N'
       jobvr = 'V'
   CASE DEFAULT
       CALL errore('zmat_diag','Invalid side',3)
   END SELECT

   lwork = 2 * n
   ALLOCATE( work(lwork), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_diag','allocating work',ABS(ierr))
   ALLOCATE( rwork(2*n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_diag','allocating rwork',ABS(ierr))
   ALLOCATE( vl(n,n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_diag','allocating vl',ABS(ierr))
   ALLOCATE( vr(n,n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_diag','allocating vr',ABS(ierr))

   ALLOCATE( a_(n,n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_diag','allocating a_',ABS(ierr))

   a_(1:n,1:n) = a(1:n,1:n)

   ! use magma if possible
   CALL ZGEEV( jobvl, jobvr, n, a_, SIZE(a_,1), w, vl, n, vr, n, work, &
               lwork, rwork, info )

   IF ( info < 0 ) CALL errore('zmat_diag', 'zgeev: info illegal value', -info )
   IF ( info > 0 ) CALL errore('zmat_diag', 'zgeev: eigenvectors not converged', info )

   SELECT CASE ( side )
   CASE ( 'L', 'l' )
       z(1:n,1:n) = vl
   CASE ( 'R', 'r' )
       z(1:n,1:n) = vr
   END SELECT

   DEALLOCATE( a_, work, rwork, vl, vr, STAT=ierr)
      IF(ierr/=0) CALL errore('zmat_diag','deallocating a_, work--vr',ABS(ierr))

   RETURN
END SUBROUTINE zmat_diag_oneside


!**********************************************************
   SUBROUTINE zmat_diag_twoside( zl, zr, w, a, n)
   !**********************************************************
   !
   ! utility to diagonalize complex non-hermitean matrices
   !
   IMPLICIT NONE
   COMPLEX(dbl),        INTENT(IN)  :: a(:,:)
   COMPLEX(dbl),        INTENT(OUT) :: zl(:,:)
   COMPLEX(dbl),        INTENT(OUT) :: zr(:,:)
   COMPLEX(dbl),        INTENT(OUT) :: w(:)
   INTEGER,             INTENT(in)  :: n

   INTEGER   :: ierr, info, lwork
   COMPLEX(dbl), ALLOCATABLE :: work(:)
   COMPLEX(dbl), ALLOCATABLE :: a_(:,:)
   REAL(dbl),    ALLOCATABLE :: rwork(:)

   ! get the dimension of the problem
   IF ( n <= 0 ) CALL errore('zmat_diag','Invalid N',ABS(n)+1)
   IF ( n > SIZE(a,1) .OR. n > SIZE(a,2) ) &
        CALL errore('zmat_diag','Invalid A dimensions',ABS(n)+1)
   IF ( n > SIZE(zl,1) .OR. n > SIZE(zl,2) .OR.  n > SIZE(zr,1) .OR. n > SIZE(zr,2) ) &
        CALL errore('zmat_diag','Invalid Z dimensions',ABS(n)+1)

   lwork = 2 * n
   ALLOCATE( work(lwork), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_diag','allocating work',ABS(ierr))
   ALLOCATE( rwork(2*n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_diag','allocating rwork',ABS(ierr))

   ALLOCATE( a_(n,n), STAT=ierr )
      IF(ierr/=0) CALL errore('zmat_diag','allocating a_',ABS(ierr))

   a_(1:n,1:n) = a(1:n,1:n)

   ! use magma if possible
   CALL ZGEEV( 'V', 'V', n, a_, SIZE(a_,1), w, zl, n, zr, n, work, &
               lwork, rwork, info )

   IF ( info < 0 ) CALL errore('zmat_diag', 'zgeev: info illegal value', -info )
   IF ( info > 0 ) CALL errore('zmat_diag', 'zgeev: eigenvectors not converged', info )

   DEALLOCATE( a_, work, rwork, STAT=ierr)
      IF(ierr/=0) CALL errore('zmat_diag','deallocating a_, work--vr',ABS(ierr))

   RETURN
END SUBROUTINE zmat_diag_twoside


!**********************************************************
   SUBROUTINE zmat_inv( n, a, z, det_a, ierr )
   !**********************************************************
   !
   ! compute Z = inv( A ), and, if required, 
   ! the determinant of A
   ! The current implementation is generalized from invmat.f90 
   ! of the Quantum-Espresso distribution
   !
   IMPLICIT NONE
   INTEGER,                INTENT(IN)  :: n
   COMPLEX(dbl),           INTENT(IN)  :: a(:,:)
   COMPLEX(dbl),           INTENT(OUT) :: z(:,:)
   COMPLEX(dbl), OPTIONAL, INTENT(OUT) :: det_a
   INTEGER,      OPTIONAL, INTENT(OUT) :: ierr
   !
   INTEGER           :: i, nb, info, ierr_, ldz, lwork, ipiv (n)
   INTEGER, EXTERNAL :: ILAENV
   ! info=0: inversion was successful
   ! ldz   : leading dimension (the same as n)
   ! ipiv  : work space for pivoting 
   COMPLEX(dbl), ALLOCATABLE :: work(:)
   !
   !
   IF ( PRESENT( ierr ) ) ierr=0
   !
   ldz = SIZE(z,1)
   z(1:n,1:n) = a(1:n,1:n)

   !
   ! perform matrix inversion according to LAPACK
   !
   ! First get the optimum LWORK
   ! check with the use of magma
   !
   nb = ILAENV( 1, 'ZGETRI', ' ', n, -1, -1, -1 )
   !
   lwork = n * nb
   !
   ALLOCATE( work( lwork ), STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore ('zmat_inv', 'allocating work', ABS (ierr_) )
   ! 
   !
   ! use magma if possible
   CALL ZGETRF (n, n, z, ldz, ipiv, info)
   !
   IF ( PRESENT(ierr) ) THEN
       ! 
       IF ( info/=0 ) THEN
           ierr = info
           RETURN
       ENDIF
       ! 
   ELSE
       IF ( info/=0 ) CALL errore ('zmat_inv', 'error in ZGETRF', ABS (info) )
   ENDIF

   !
   ! compute the determinan if required
   !
   IF ( PRESENT( det_a ) ) THEN
      !
      det_a = ONE
      DO i = 1, n
         det_a = det_a * z(i,i)
      ENDDO
      !
   ENDIF
   !
   ! use magma if possible
   CALL ZGETRI (n, z, ldz, ipiv, work, lwork, info)
   !
   IF ( PRESENT(ierr) ) THEN
       IF ( info/=0 ) ierr = info
   ELSE
       IF ( info/=0 ) CALL errore ('zmat_inv', 'error in ZGETRI', ABS (info) )
   ENDIF
   !
   ! 
   DEALLOCATE( work, STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore ('zmat_inv', 'deallocating work', ABS (ierr_) )
   !
END SUBROUTINE zmat_inv 


!**********************************************************
   SUBROUTINE zmat_bnd_inv( n, kl, ku, ab, ldab, z, ierr )
   !**********************************************************
   !
   ! compute Z = inv( AB ) for a banded matrix AB 
   ! AB should be entry in a band matrix format
   !
   IMPLICIT NONE
   INTEGER,                INTENT(IN)  :: n
   INTEGER,                INTENT(IN)  :: kl, ku, ldab
   COMPLEX(dbl),           INTENT(IN)  :: ab(ldab,n)
   COMPLEX(dbl),           INTENT(OUT) :: z(n,n)
   INTEGER,      OPTIONAL, INTENT(OUT) :: ierr
   !
   INTEGER           :: i, info, ldz, ipiv (n)
   !INTEGER           :: ierr_, lwork
   INTEGER, EXTERNAL :: ILAENV
   ! info=0: inversion was successful
   ! ldz   : leading dimension (the same as n)
   ! ipiv  : work space for pivoting 
   !
   !
   IF ( PRESENT( ierr ) ) ierr=0
   !
   ldz = n
   !
   ! perform matrix inversion according to LAPACK
   !
   IF ( ldab < 2 * kl + ku + 1) CALL errore ('zmat_bnd_inv', 'invalid ldab', 10 )
   !
   z = CZERO
   DO i = 1, n
       z(i,i) = CONE
   ENDDO
   !
   ! use magma if possible
   CALL ZGBSV (n, kl, ku, n, ab, ldab, ipiv, z, ldz, info)


   IF ( PRESENT(ierr) ) THEN
       ! 
       IF ( info/=0 ) THEN
           ierr = info
           RETURN
       ENDIF
       ! 
   ELSE
       IF ( info/=0 ) CALL errore ('zmat_bnd_inv', 'error in ZGBTRF', ABS (info) )
   ENDIF
   !
END SUBROUTINE zmat_bnd_inv 


!**********************************************************
   SUBROUTINE dmat_inv( n, a, z, det_a, ierr )
   !**********************************************************
   !
   ! compute Z = inv( A ), and, if required, 
   ! the determinant of A
   ! The current implementation is generalized from invmat.f90 
   ! of the Quantum-Espresso distribution
   !
   IMPLICIT NONE
   INTEGER,             INTENT(IN)    :: n
   REAL(dbl),           INTENT(IN)    :: a(:,:)
   REAL(dbl),           INTENT(OUT)   :: z(:,:)
   REAL(dbl), OPTIONAL, INTENT(OUT)   :: det_a
   INTEGER,   OPTIONAL, INTENT(OUT)   :: ierr
   !
   INTEGER           :: i, nb, info, ierr_, ldz, lwork, ipiv (n)
   INTEGER, EXTERNAL :: ILAENV
   ! info=0: inversion was successful
   ! ldz   : leading dimension (the same as n)
   ! ipiv  : work space for pivoting 
   REAL(dbl), ALLOCATABLE :: work(:)
   !
   !
   IF ( PRESENT( ierr ) ) ierr = 0
   !
   ldz = SIZE(z,1)
   z(1:n,1:n) = a(1:n,1:n)
   !
   ! perform matrix inversion according to LAPACK
   !
   ! First get the optimum LWORK
   ! check with magma
   !
   nb = ILAENV( 1, 'DGETRI', ' ', n, -1, -1, -1 )
   lwork = n * nb
   !
   !
   ALLOCATE( work( lwork ), STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore ('dmat_inv', 'allocating work', ABS(ierr_) )
   ! 
   !
   ! use magma if possible
   CALL DGETRF (n, n, z, ldz, ipiv, info)
   !
   IF ( PRESENT(ierr) ) THEN
       ! 
       IF ( info/=0 ) THEN
           ierr = info
           RETURN
       ENDIF
       ! 
   ELSE
       IF ( info/=0 ) CALL errore ('dmat_inv', 'error in DGETRF', ABS (info) )
   ENDIF

   !
   ! compute the determinan if required
   !
   IF ( PRESENT( det_a ) ) THEN
       !
       det_a = ONE
       DO i = 1, n
          det_a = det_a * z(i,i)
       ENDDO
       !
   ENDIF
   !
   ! use magma if possible
   CALL DGETRI (n, z, ldz, ipiv, work, lwork, info)
   !
   IF ( PRESENT(ierr) ) THEN
       IF ( info/=0 ) ierr = info
   ELSE
       IF ( info/=0 ) CALL errore ('dmat_inv', 'error in ZGETRI', ABS (info) )
   ENDIF
   !
   ! 
   DEALLOCATE( work, STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore ('dmat_inv', 'deallocating work', ABS (ierr_) )
   !
END SUBROUTINE dmat_inv 


!**********************************************************
   SUBROUTINE dmat_bnd_inv( n, kl, ku, ab, ldab, z, ierr )
   !**********************************************************
   !
   ! compute Z = inv( AB ) for a banded matrix AB 
   ! AB should be entry in a band matrix format
   !
   IMPLICIT NONE
   INTEGER,                INTENT(IN)  :: n
   INTEGER,                INTENT(IN)  :: kl, ku, ldab
   REAL(dbl),              INTENT(IN)  :: ab(ldab,n)
   REAL(dbl),              INTENT(OUT) :: z(n,n)
   INTEGER,      OPTIONAL, INTENT(OUT) :: ierr
   !
   INTEGER           :: i, info, ldz, ipiv (n)
   !INTEGER           :: ierr_, lwork
   INTEGER, EXTERNAL :: ILAENV
   ! info=0: inversion was successful
   ! ldz   : leading dimension (the same as n)
   ! ipiv  : work space for pivoting 
   !
   !
   IF ( PRESENT( ierr ) ) ierr=0
   !
   ldz = n
   !
   ! perform matrix inversion according to LAPACK
   !
   IF ( ldab < 2 * kl + ku + 1) CALL errore ('dmat_bnd_inv', 'invalid ldab', 10 )
   !
   z = ZERO
   DO i = 1, n
       z(i,i) = ONE
   ENDDO
   !
   ! use magma if possible
   CALL DGBSV (n, kl, ku, n, ab, ldab, ipiv, z, ldz, info)


   IF ( PRESENT(ierr) ) THEN
       ! 
       IF ( info/=0 ) THEN
           ierr = info
           RETURN
       ENDIF
       ! 
   ELSE
       IF ( info/=0 ) CALL errore ('dmat_bnd_inv', 'error in DGBTRF', ABS (info) )
   ENDIF
   !
END SUBROUTINE dmat_bnd_inv 


!**********************************************************
   SUBROUTINE zmat_hsqrt( n, a, z, ierr )
   !**********************************************************
   !
   ! compute Z = sqrt( A ) for a Hermitean matrix A
   !
   IMPLICIT NONE
   INTEGER,                INTENT(IN)  :: n
   COMPLEX(dbl),           INTENT(IN)  :: a(:,:)
   COMPLEX(dbl),           INTENT(OUT) :: z(:,:)
   INTEGER,      OPTIONAL, INTENT(OUT) :: ierr
   !
   COMPLEX(dbl), ALLOCATABLE :: work(:,:), work1(:,:)
   REAL(dbl),    ALLOCATABLE :: w(:)
   INTEGER  :: ierr_, i, j
   !
   IF ( PRESENT( ierr ) ) ierr=0
   !
   ALLOCATE( work(n,n), w(n), STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore("zmat_hsqrt","allocating work, w", ABS(ierr_))
   ALLOCATE( work1(n,n), STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore("zmat_hsqrt","allocating work1", ABS(ierr_))
   !
   ! diagonalize A
   CALL mat_hdiag( work, w, a, n )
   !
   DO i = 1, n
       IF ( w(i) <= 0.0d0 ) THEN
           ierr=1
           RETURN
       ELSE
           w(i) = SQRT(w(i))
       ENDIF
   ENDDO
   !
   DO j = 1, n
   DO i = 1, n
       !   
       work1(i,j) = work(i,j) * w(j)
       !   
   ENDDO
   ENDDO
   !   
   CALL mat_mul( z, work1, 'N', work, 'C', n, n, n)   
   !
   DEALLOCATE( work, work1, w, STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore("zmat_hsqrt","deallocating work -- w", ABS(ierr_))
   !
   RETURN
   !
END SUBROUTINE zmat_hsqrt


!**********************************************************
   SUBROUTINE dmat_hsqrt( n, a, z, ierr )
   !**********************************************************
   !
   ! compute Z = sqrt( A ) for a real hermitean (symmetric) matrix A
   !
   IMPLICIT NONE
   INTEGER,                INTENT(IN)  :: n
   REAL(dbl),              INTENT(IN)  :: a(:,:)
   REAL(dbl),              INTENT(OUT) :: z(:,:)
   INTEGER,      OPTIONAL, INTENT(OUT) :: ierr
   !
   REAL(dbl),  ALLOCATABLE :: work(:,:), work1(:,:)
   REAL(dbl),  ALLOCATABLE :: w(:)
   INTEGER  :: ierr_, i, j
   !
   IF ( PRESENT( ierr ) ) ierr=0
   !
   ALLOCATE( work(n,n), w(n), STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore("dmat_hsqrt","allocating work, w", ABS(ierr_))
   ALLOCATE( work1(n,n), STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore("dmat_hsqrt","allocating work1", ABS(ierr_))
   !
   ! diagonalize A
   CALL mat_hdiag( work, w, a, n )
   !
   DO i = 1, n
       IF ( w(i) <= 0.0d0 ) THEN
           ierr=1
           RETURN
       ELSE
           w(i) = SQRT(w(i))
       ENDIF
   ENDDO
   !
   DO j = 1, n
   DO i = 1, n
       !   
       work1(i,j) = work(i,j) * w(j)
       !   
   ENDDO
   ENDDO
   !   
   CALL mat_mul( z, work1, 'N', work, 'T', n, n, n)   
   !
   DEALLOCATE( work, work1, w, STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore("dmat_hsqrt","deallocating work -- w", ABS(ierr_))
   !
   RETURN
   !
END SUBROUTINE dmat_hsqrt


!**********************************************************
   !SUBROUTINE dmat_bnd_hsqrt( n, ab, ldab, ka, zb, ldzb, kz, ierr )
   SUBROUTINE dmat_bnd_hsqrt( n, ab, ldab, ka, z, ierr )
   !**********************************************************
   !
   ! compute Z = sqrt( A ) for a real BANDED hermitean (symmetric) matrix A
   !
   IMPLICIT NONE
   INTEGER,                INTENT(IN)  :: n, ldab, ka
   REAL(dbl),              INTENT(IN)  :: ab(:,:)
   REAL(dbl),              INTENT(OUT) :: z(:,:)
   !INTEGER,                INTENT(OUT) :: kz
   INTEGER,      OPTIONAL, INTENT(OUT) :: ierr
   !
   REAL(dbl),  ALLOCATABLE :: work(:,:), work1(:,:)
   REAL(dbl),  ALLOCATABLE :: w(:)
   INTEGER  :: ierr_, i, j
   !
   IF ( PRESENT( ierr ) ) ierr=0
   !
   ALLOCATE( work(n,n), w(n), STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore("dmat_bnd_hsqrt","allocating work, w", ABS(ierr_))
   !
   ! diagonalize A
   CALL dmat_bnd_hdiag( work, w, ab, ldab, n, ka )
   !
   ALLOCATE( work1(n,n), STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore("dmat_bnd_hsqrt","allocating work1", ABS(ierr_))
   !ALLOCATE( z(n,n), STAT=ierr_ )
   !IF ( ierr_/=0 ) CALL errore("dmat_bnd_hsqrt","allocating z", ABS(ierr_))
   !
   DO i = 1, n
       IF ( w(i) <= 0.0d0 ) THEN
           ierr=1
           RETURN
       ELSE
           w(i) = SQRT(w(i))
       ENDIF
   ENDDO
   !
   DO j = 1, n
   DO i = 1, n
       !   
       work1(i,j) = work(i,j) * w(j)
       !   
   ENDDO
   ENDDO
   !   
   CALL mat_mul( z, work1, 'N', work, 'T', n, n, n)   
   !
   DEALLOCATE( work, work1, w, STAT=ierr_ )
   IF ( ierr_/=0 ) CALL errore("dmat_hsqrt","deallocating work -- w", ABS(ierr_))
   !
   RETURN
   !
END SUBROUTINE dmat_bnd_hsqrt


!**********************************************************
   FUNCTION  zmat_is_herm( n, z, toll )
   !**********************************************************
   IMPLICIT NONE
   LOGICAL               :: zmat_is_herm
   INTEGER               :: n
   COMPLEX(dbl)          :: z(:,:)
   REAL(dbl), OPTIONAL   :: toll
   !
   ! n : actual dimension of Z
   ! check if a complex matrix is hermitean.
   !
   REAL(dbl)     :: toll_
   CHARACTER(12) :: subname='zmat_is_herm'
   REAL(dbl)     :: rtmp
   COMPLEX(dbl)  :: ztmp
   INTEGER       :: i, j
   

   toll_ = TINY(ZERO)  
   IF ( PRESENT(toll) ) toll_ = toll
   IF ( toll_ <= 0.0 ) CALL errore(subname,'Invalid TOLL',1)
  
   IF ( n > SIZE( z, 1) ) CALL errore(subname,'Invalid n I',n)
   IF ( n > SIZE( z, 2) ) CALL errore(subname,'Invalid n II',n)

   zmat_is_herm = .TRUE.
   !
   main_loop: &
   DO j = 1, n
   DO i = 1, j
       !
       ztmp = z(i,j) - CONJG( z(j,i) )
       rtmp = REAL( ztmp * CONJG( ztmp ), dbl ) 
       !
       IF ( rtmp > toll_ ) THEN
          zmat_is_herm = .FALSE. 
          EXIT main_loop
       ENDIF
       !
   ENDDO
   ENDDO main_loop
       
   RETURN
END FUNCTION zmat_is_herm


!**********************************************************
   FUNCTION  dmat_is_herm( n, A, toll )
   !**********************************************************
   IMPLICIT NONE
   LOGICAL               :: dmat_is_herm
   INTEGER               :: n
   REAL(dbl)             :: A(:,:)
   REAL(dbl), OPTIONAL   :: toll
   !
   ! n : actual dimension of Z
   ! check if a real matrix is hermitean.
   !
   REAL(dbl)     :: toll_
   CHARACTER(12) :: subname='dmat_is_herm'
   REAL(dbl)     :: rtmp
   INTEGER       :: i, j
   

   toll_ = TINY(ZERO)  
   IF ( PRESENT(toll) ) toll_ = toll
   IF ( toll_ <= 0.0 ) CALL errore(subname,'Invalid TOLL',1)
  
   IF ( n > SIZE( A, 1) ) CALL errore(subname,'Invalid n I',n)
   IF ( n > SIZE( A, 2) ) CALL errore(subname,'Invalid n II',n)

   dmat_is_herm = .TRUE.
   !
   main_loop: &
   DO j = 1, n
   DO i = 1, j
       !
       rtmp = A(i,j) - A(j,i)
       rtmp = rtmp**2
       !
       IF ( rtmp > toll_ ) THEN
          dmat_is_herm = .FALSE. 
          EXIT main_loop
       ENDIF
       !
   ENDDO
   ENDDO main_loop
       
   RETURN
END FUNCTION dmat_is_herm


!**********************************************************
   FUNCTION  zmat_is_symm( n, z, toll )
   !**********************************************************
   IMPLICIT NONE
   LOGICAL               :: zmat_is_symm
   INTEGER               :: n
   COMPLEX(dbl)          :: z(:,:)
   REAL(dbl), OPTIONAL   :: toll
   !
   ! n : actual dimension of Z
   ! check if a complex matrix is symmetric.
   !
   REAL(dbl)     :: toll_
   CHARACTER(12) :: subname='zmat_is_symm'
   REAL(dbl)     :: rtmp
   COMPLEX(dbl)  :: ztmp
   INTEGER       :: i, j
   

   toll_ = TINY(ZERO)  
   IF ( PRESENT(toll) ) toll_ = toll
   IF ( toll_ <= 0.0 ) CALL errore(subname,'Invalid TOLL',1)
  
   IF ( n > SIZE( z, 1) ) CALL errore(subname,'Invalid n I',n)
   IF ( n > SIZE( z, 2) ) CALL errore(subname,'Invalid n II',n)

   zmat_is_symm = .TRUE.
   !
   main_loop: &
   DO j = 1, n
   DO i = 1, j
       !
       ztmp = z(i,j) - z(j,i)
       rtmp = REAL( ztmp * CONJG( ztmp ), dbl ) 
       !
       IF ( rtmp > toll_ ) THEN
          zmat_is_symm = .FALSE. 
          EXIT main_loop
       ENDIF
       !
   ENDDO
   ENDDO main_loop
       
   RETURN
END FUNCTION zmat_is_symm


!**********************************************************
   FUNCTION  dmat_is_symm( n, A, toll )
   !**********************************************************
   IMPLICIT NONE
   LOGICAL               :: dmat_is_symm
   INTEGER               :: n
   REAL(dbl)             :: A(:,:)
   REAL(dbl), OPTIONAL   :: toll
   !
   ! n : actual dimension of Z
   ! check if a real matrix is hermitean.
   !
   REAL(dbl)     :: toll_
   CHARACTER(12) :: subname='dmat_is_symm'
   REAL(dbl)     :: rtmp
   INTEGER       :: i, j
   

   toll_ = TINY(ZERO)  
   IF ( PRESENT(toll) ) toll_ = toll
   IF ( toll_ <= 0.0 ) CALL errore(subname,'Invalid TOLL',1)
  
   IF ( n > SIZE( A, 1) ) CALL errore(subname,'Invalid n I',n)
   IF ( n > SIZE( A, 2) ) CALL errore(subname,'Invalid n II',n)

   dmat_is_symm = .TRUE.
   !
   main_loop: &
   DO j = 1, n
   DO i = 1, j
       !
       rtmp = A(i,j) - A(j,i)
       rtmp = rtmp**2
       !
       IF ( rtmp > toll_ ) THEN
          dmat_is_symm = .FALSE. 
          EXIT main_loop
       ENDIF
       !
   ENDDO
   ENDDO main_loop
       
   RETURN
END FUNCTION dmat_is_symm

!**********************************************************
   FUNCTION  zmat_is_diagdom( n, z, dom_type )
   !**********************************************************
   IMPLICIT NONE
   LOGICAL                :: zmat_is_diagdom
   INTEGER                :: n
   COMPLEX(dbl)           :: z(:,:)
   CHARACTER(*), OPTIONAL :: dom_type
   !
   ! n : actual dimension of Z
   ! check if a complex matrix is diagonally dominant.
   !
   CHARACTER(15) :: subname='zmat_is_diagdom'
   REAL(dbl), ALLOCATABLE :: z_abs(:,:)
   REAL(dbl)              :: rtmp
   INTEGER                :: i, ierr
   

   IF ( n > SIZE( z, 1) ) CALL errore(subname,'Invalid n I',n)
   IF ( n > SIZE( z, 2) ) CALL errore(subname,'Invalid n II',n)

   ALLOCATE(z_abs(n,n))
   IF (ierr/=0)  CALL errore(subname,'allocating z_abs',ABS(ierr))
   z_abs = abs(z)

   zmat_is_diagdom = .TRUE.
   !
   SELECT CASE (TRIM(dom_type))
   CASE("row")
     !
     DO i = 1, n
         !
         rtmp = sum(z_abs(1:i-1,i))+sum(z_abs(i+1:n,i))
         !
         IF ( rtmp > z_abs(i,i) ) THEN
            zmat_is_diagdom = .FALSE. 
            EXIT
         ENDIF
         !
     ENDDO
     !
   CASE("column")
     !
     DO i = 1, n
         !
         rtmp = sum(z_abs(i,1:i-1))+sum(z_abs(i,i+1:n))
         !
         IF ( rtmp > z_abs(i,i) ) THEN
            zmat_is_diagdom = .FALSE. 
            EXIT
         ENDIF
         !
     ENDDO
     !
   CASE DEFAULT
     !
     call errore(subname,"invalid dom_type",10)
     !
   END SELECT
       
   RETURN
END FUNCTION zmat_is_diagdom


!**********************************************************
   FUNCTION  dmat_is_diagdom( n, A, dom_type )
   !**********************************************************
   IMPLICIT NONE
   LOGICAL                :: dmat_is_diagdom
   INTEGER                :: n
   REAL(dbl)              :: A(:,:)
   CHARACTER(*), OPTIONAL :: dom_type
   !
   ! n : actual dimension of A
   ! check if a real matrix is diagonally dominant.
   !
   CHARACTER(15) :: subname='dmat_is_diagdom'
   !
   REAL(dbl), ALLOCATABLE :: A_abs(:,:)
   REAL(dbl)              :: rtmp
   INTEGER                :: i, ierr
   

   IF ( n > SIZE( A, 1) ) CALL errore(subname,'Invalid n I',n)
   IF ( n > SIZE( A, 2) ) CALL errore(subname,'Invalid n II',n)

   ALLOCATE(A_abs(n,n))
   IF (ierr/=0)  CALL errore(subname,'allocating A_abs',ABS(ierr))
   A_abs = abs(A)

   dmat_is_diagdom = .TRUE.
   !
   SELECT CASE (TRIM(dom_type))
   CASE("row")
     !
     DO i = 1, n
         !
         rtmp = sum(A_abs(1:i-1,i))+sum(A_abs(i+1:n,i))
         !
         IF ( rtmp > A_abs(i,i) ) THEN
            dmat_is_diagdom = .FALSE. 
            EXIT
         ENDIF
         !
     ENDDO
     !
   CASE("column")
     !
     DO i = 1, n
         !
         rtmp = sum(A_abs(i,1:i-1))+sum(A_abs(i,i+1:n))
         !
         IF ( rtmp > A_abs(i,i) ) THEN
            dmat_is_diagdom = .FALSE. 
            EXIT
         ENDIF
         !
     ENDDO
     !
   CASE DEFAULT
     !
     call errore(subname,"invalid dom_type",10)
     !
   END SELECT
       
   RETURN
END FUNCTION dmat_is_diagdom


!**********************************************************
   FUNCTION  zmat_is_unitary( m, n, z, side, toll )
   !**********************************************************
   IMPLICIT NONE
   LOGICAL                :: zmat_is_unitary
   INTEGER                :: m,n
   COMPLEX(dbl)           :: z(:,:)
   CHARACTER(*), OPTIONAL :: side 
   REAL(dbl),    OPTIONAL :: toll
   !
   ! m, n : actual dimensions of A
   ! check if a complex matrix is unitary.
   ! SIDE='left'  only   A^{\dag} * A = I
   ! SIDE='right' only   A * A^{\dag} = I
   ! SIDE='both'  both sides (DEFAULT)
   !
   REAL(dbl)     :: toll_
   CHARACTER(10) :: side_
   CHARACTER(15) :: subname='zmat_is_unitary'
   INTEGER       :: dim1,dim2
   INTEGER       :: i, j, ierr
   COMPLEX(dbl), ALLOCATABLE  :: result(:,:)
   
   zmat_is_unitary = .TRUE. 

   toll_ = TINY(ZERO)  
   side_ = 'both'
   IF ( PRESENT(side) ) side_ = TRIM(side)
   IF ( PRESENT(toll) ) toll_ = toll
   IF ( toll_ <= 0 ) CALL errore(subname,'Invalid TOLL',1)
  
   IF ( m > SIZE( z, 1) ) CALL errore(subname,'Invalid m',m)
   IF ( n > SIZE( z, 2) ) CALL errore(subname,'Invalid n',n)
   dim1 = m
   dim2 = n
   IF ( dim1 <= 0) CALL errore(subname,'Invalid dim1',ABS(dim1)+1)
   IF ( dim2 <= 0) CALL errore(subname,'Invalid dim2',ABS(dim2)+1)


   !
   ! check side LEFT
   !
   IF ( TRIM(side_) == 'both' .OR. TRIM(side_) == 'BOTH' .OR. &
        TRIM(side_) == 'left' .OR. TRIM(side_) == 'LEFT'  ) THEN 

       ALLOCATE( result(dim2,dim2), STAT=ierr )
          IF ( ierr /= 0 ) CALL errore(subname,'allocating result',ABS(ierr))
       ! 
       ! matrix mult
       CALL zmat_mul( result, z, 'C', z, 'N', dim2,dim2,dim1)

       DO j=1,dim2
       DO i=1,dim2
           IF ( i==j ) THEN
                IF ( ABS( result(i,j) -CONE ) > toll_ ) zmat_is_unitary = .FALSE.
           ELSE
                IF ( ABS( result(i,j) ) > toll_ ) zmat_is_unitary = .FALSE.
           ENDIF
       ENDDO
       ENDDO

       DEALLOCATE( result, STAT=ierr)
          IF ( ierr /= 0 ) CALL errore(subname,'deallocating result',ABS(ierr))
   ENDIF
       
   !
   ! check side RIGHT
   !
   IF ( TRIM(side_) == 'both' .OR. TRIM(side_) == 'BOTH' .OR. &
        TRIM(side_) == 'right'.OR. TRIM(side_) == 'RIGHT' ) THEN 

       ALLOCATE( result(dim1,dim1), STAT=ierr )
          IF ( ierr /= 0 ) CALL errore(subname,'allocating result',ABS(ierr))
       ! 
       ! matrix mult
       CALL zmat_mul( result, z, 'N', z, 'C', dim1,dim1,dim2)

       DO j=1,dim1
       DO i=1,dim1
           IF ( i==j ) THEN
                IF ( ABS( result(i,j) -CONE ) > toll_ ) zmat_is_unitary = .FALSE.
           ELSE
                IF ( ABS( result(i,j) ) > toll_ ) zmat_is_unitary = .FALSE.
           ENDIF
       ENDDO
       ENDDO

       DEALLOCATE( result, STAT=ierr)
          IF ( ierr /= 0 ) CALL errore(subname,'deallocating result',ABS(ierr))
   ENDIF
   RETURN
END FUNCTION zmat_is_unitary


!**********************************************************
   FUNCTION  zmat_rank( m, n, a, toll )
   !**********************************************************
   IMPLICIT NONE
   INTEGER          :: zmat_rank
   INTEGER          :: m,n
   COMPLEX(dbl)     :: a(:,:)
   REAL(dbl)        :: toll
   !
   INTEGER :: i,ierr 
   REAL(dbl),    ALLOCATABLE :: s(:)
   COMPLEX(dbl), ALLOCATABLE :: atmp(:,:), u(:,:), vt(:,:)

   IF ( m > SIZE(a,1) ) CALL errore('zmat_rank','Invalid m',ABS(m)+1)
   IF ( n > SIZE(a,2) ) CALL errore('zmat_rank','Invalid n',ABS(n)+1)

   ALLOCATE( atmp(m,n), u(m,m), vt(n,n), s(MIN(m,n)), STAT=ierr )
   IF ( ierr /=0 ) CALL errore('zmat_rank','allocating atmp--s',ABS(ierr))

   !
   ! local copy
   atmp(:,:) = a(1:m,1:n)
   !
   ! svd decomposition
   CALL mat_svd(m,n,atmp,s,u,vt)
   !
   ! TO BE SUBSTITUTED TO A CALL TO LOCATE
   zmat_rank = MIN(m,n)
   DO i=1, MIN(m,n)
       IF ( ABS(s(i)) < toll ) THEN
          zmat_rank = i-1
          EXIT
       ENDIF
   ENDDO
   !
   DEALLOCATE( atmp, u, vt, s, STAT=ierr )
   IF ( ierr /=0 ) CALL errore('zmat_rank','deallocating atmp--s',ABS(ierr))
   !
END FUNCTION zmat_rank


!**********************************************************
   FUNCTION  dmat_rank( m, n, a, toll )
   !**********************************************************
   IMPLICIT NONE
   INTEGER            :: dmat_rank
   INTEGER            :: m,n
   REAL(dbl)          :: a(:,:)
   REAL(dbl)          :: toll
   !
   INTEGER :: i,ierr 
   REAL(dbl), ALLOCATABLE :: s(:)
   REAL(dbl), ALLOCATABLE :: atmp(:,:), u(:,:), vt(:,:)


   IF ( m > SIZE(a,1) ) CALL errore('dmat_rank','Invalid m',ABS(m)+1)
   IF ( n > SIZE(a,2) ) CALL errore('dmat_rank','Invalid n',ABS(n)+1)

   ALLOCATE( atmp(m,n), u(m,m), vt(n,n), s(MIN(m,n)), STAT=ierr )
   IF ( ierr /=0 ) CALL errore('dmat_rank','allocating atmp--s',ABS(ierr))

   !
   ! local copy
   atmp(:,:) = a(1:m,1:n)

   !
   ! svd decomposition
   CALL mat_svd(m,n,atmp,s,u,vt)
   !
   ! TO BE SUBSTITUTED TO A CALL TO LOCATE
   dmat_rank = MIN(m,n)
   DO i=1, MIN(m,n)
       IF ( ABS(s(i)) < toll ) THEN
          dmat_rank = i-1
          EXIT
       ENDIF
   ENDDO
   !
   DEALLOCATE( atmp, u, vt, s, STAT=ierr )
   IF ( ierr /=0 ) CALL errore('dmat_rank','deallocating atmp--s',ABS(ierr))
   !
END FUNCTION dmat_rank


!**********************************************************
   FUNCTION  zmat_ge_dotp( m, n, a, b)
   !**********************************************************
   IMPLICIT NONE
   COMPLEX(dbl)       :: zmat_ge_dotp
   INTEGER            :: m, n
   COMPLEX(dbl)       :: a(:,:), b(:,:)
   !
   COMPLEX(dbl) :: dotp
   INTEGER      :: i, j

   IF ( m > SIZE(a,1) ) CALL errore('zmat_dotp','Invalid m I',ABS(m)+1)
   IF ( n > SIZE(a,2) ) CALL errore('zmat_dotp','Invalid n I',ABS(n)+1)
   IF ( m > SIZE(b,1) ) CALL errore('zmat_dotp','Invalid m II',ABS(m)+1)
   IF ( n > SIZE(b,2) ) CALL errore('zmat_dotp','Invalid n II',ABS(n)+1)

   dotp = CZERO
   !
   DO j = 1, n
   DO i = 1, m
      dotp = dotp + a(i,j) * CONJG ( b(i,j) )
   ENDDO
   ENDDO
   !
   zmat_ge_dotp = dotp
   RETURN
   !
END FUNCTION zmat_ge_dotp


!**********************************************************
   FUNCTION  dmat_ge_dotp( m, n, a, b)
   !**********************************************************
   IMPLICIT NONE
   REAL(dbl)       :: dmat_ge_dotp
   INTEGER         :: m, n
   REAL(dbl)       :: a(:,:), b(:,:)
   !
   REAL(dbl)    :: dotp
   INTEGER      :: i, j

   IF ( m > SIZE(a,1) ) CALL errore('dmat_ge_dotp','Invalid m I',ABS(m)+1)
   IF ( n > SIZE(a,2) ) CALL errore('dmat_ge_dotp','Invalid n I',ABS(n)+1)
   IF ( m > SIZE(b,1) ) CALL errore('dmat_ge_dotp','Invalid m II',ABS(m)+1)
   IF ( n > SIZE(b,2) ) CALL errore('dmat_ge_dotp','Invalid n II',ABS(n)+1)

   dotp = ZERO
   !
   DO j = 1, n
   DO i = 1, m
      dotp = dotp + a(i,j) * b(i,j)
   ENDDO
   ENDDO
   !
   dmat_ge_dotp = dotp
   RETURN
END FUNCTION dmat_ge_dotp


!**********************************************************
   FUNCTION  zmat_he_dotp( m, a, b)
   !**********************************************************
   IMPLICIT NONE
   COMPLEX(dbl)       :: zmat_he_dotp
   INTEGER            :: m
   COMPLEX(dbl)       :: a(:,:), b(:,:)
   !
   COMPLEX(dbl) :: dotp
   INTEGER      :: i, j

   IF ( m > SIZE(a,1) ) CALL errore('zmat_hdotp','Invalid m1',ABS(m)+1)
   IF ( m > SIZE(a,2) ) CALL errore('zmat_hdotp','Invalid m2',ABS(m)+1)
   IF ( m > SIZE(b,1) ) CALL errore('zmat_hdotp','Invalid m3',ABS(m)+1)
   IF ( m > SIZE(b,2) ) CALL errore('zmat_hdotp','Invalid m4',ABS(m)+1)

   dotp = CZERO
   !
   DO j = 1, m
   DO i = 1, j-1
      dotp = dotp + 2.0_dbl * CONJG( a(i,j) ) * b(i,j)
   ENDDO
   ENDDO
   !
   DO i = 1, m
      dotp = dotp + 1.0_dbl * CONJG( a(i,i) ) * b(i,i)
   ENDDO
   !
   zmat_he_dotp = dotp
   RETURN
   !
END FUNCTION zmat_he_dotp


!**********************************************************
   FUNCTION  dmat_sy_dotp( m, a, b)
   !**********************************************************
   IMPLICIT NONE
   REAL(dbl)          :: dmat_sy_dotp
   INTEGER            :: m
   REAL(dbl)          :: a(:,:), b(:,:)
   !
   REAL(dbl) :: dotp
   INTEGER   :: i, j

   IF ( m > SIZE(a,1) ) CALL errore('dmat_sy_dotp','Invalid m1',ABS(m)+1)
   IF ( m > SIZE(a,2) ) CALL errore('dmat_sy_dotp','Invalid m2',ABS(m)+1)
   IF ( m > SIZE(b,1) ) CALL errore('dmat_sy_dotp','Invalid m3',ABS(m)+1)
   IF ( m > SIZE(b,2) ) CALL errore('dmat_sy_dotp','Invalid m4',ABS(m)+1)

   dotp = ZERO
   !
   DO j = 1, m
   DO i = 1, j-1
      dotp = dotp + 2.0_dbl * a(i,j) * b(i,j)
   ENDDO
   ENDDO
   !
   DO i = 1, m
      dotp = dotp + 1.0_dbl * a(i,i) * b(i,i)
   ENDDO
   !
   dmat_sy_dotp = dotp
   RETURN
   !
END FUNCTION dmat_sy_dotp


!**********************************************************
   FUNCTION  zmat_hp_dotp( m, uplo, ap, bp)
   !**********************************************************
   !
   ! UPLO: 'U', 'u', provide the upper triangular parts of A, B,
   !       'L', 'l', provide the lower triangular parts
   !
   IMPLICIT NONE
   COMPLEX(dbl)       :: zmat_hp_dotp
   INTEGER            :: m
   CHARACTER          :: uplo
   COMPLEX(dbl)       :: ap(:), bp(:)
   !
   COMPLEX(dbl) :: dotp
   INTEGER      :: i, j, l

   IF ( m*(m+1)/2 > SIZE(ap) ) CALL errore('zmat_hdotp','Invalid m A',ABS(m)+1)
   IF ( m*(m+1)/2 > SIZE(ap) ) CALL errore('zmat_hdotp','Invalid m B',ABS(m)+1)

   dotp = CZERO
   l    = 0
   !
   SELECT CASE ( uplo)
   CASE ( 'U', 'u' )
       !
       DO j = 1, m
           !
           l = (j-1) * j / 2
           !
           DO i = 1, j-1
               l = l+1
               dotp = dotp + 2.0_dbl * CONJG( ap(l) ) * bp(l)
           ENDDO
           !
           l = l+1
           dotp = dotp + 1.0_dbl * CONJG ( ap(l) ) * bp(l)
           !
       ENDDO
       !
   CASE ( 'L', 'l' )
       !
       DO j = 1, m
           !
           l = (j-1) * ( 2*m -j +2) /2 +1
           dotp = dotp + 1.0_dbl * CONJG ( ap(l) ) * bp(l)
           !
           DO i = j+1, m
               l = l+1
               dotp = dotp + 2.0_dbl * CONJG( ap(l) ) * bp(l)
           ENDDO
           !
           !
       ENDDO
       !
   CASE DEFAULT
       CALL errore('zmat_hdotp','invalid uplo: '//uplo, 2 )
   END SELECT
    
   !
   zmat_hp_dotp = dotp
   RETURN
   !
END FUNCTION zmat_hp_dotp


!**********************************************************
   FUNCTION  dvec_dotp( m, va, vb)
   !**********************************************************
   IMPLICIT NONE
   REAL(dbl)       :: dvec_dotp
   INTEGER         :: m
   REAL(dbl)       :: va(:), vb(:)
   !
   REAL(dbl)    :: dotp
   INTEGER      :: i
   REAL(dbl), external :: DDOT

   IF ( m > SIZE(va) ) CALL errore('dvec_dotp','Invalid m I',ABS(m)+1)
   IF ( m > SIZE(vb) ) CALL errore('dvec_dotp','Invalid m II',ABS(m)+1)
   !
   dotp =  DDOT(m,va,1,vb,1)
   !
   !dotp = ZERO
   ! DO i = 1, m
   !     dotp = dotp + va(i) * vb(i)
   ! ENDDO
   !
   dvec_dotp = dotp
   RETURN
END FUNCTION dvec_dotp


!**********************************************************
   FUNCTION  zvec_dotp( m, va, vb)
   !**********************************************************
   IMPLICIT NONE
   COMPLEX(dbl)    :: zvec_dotp
   INTEGER         :: m
   COMPLEX(dbl)    :: va(:), vb(:)
   !
   COMPLEX(dbl) :: dotp
   INTEGER      :: i
   COMPLEX(dbl), external :: ZDOTC

   IF ( m > SIZE(va) ) CALL errore('zvec_dotp','Invalid m I',ABS(m)+1)
   IF ( m > SIZE(vb) ) CALL errore('zvec_dotp','Invalid m II',ABS(m)+1)

   dotp = CZERO
   !
   dotp =  ZDOTC(m,va,1,vb,1)
   ! DO i = 1, m
   !     dotp = dotp + CONJG(va(i)) * vb(i)
   ! ENDDO
   !
   zvec_dotp = dotp
   RETURN
END FUNCTION zvec_dotp


!**********************************************************
   SUBROUTINE zmat_chol( a, tmat, n, uplo, ierr )
   !**********************************************************
   !
   ! utility to perform the Cholesky decomposition of
   ! positive defined matrices
   !
   implicit none
   complex(dbl),  intent(IN)  :: a(:,:)
   complex(dbl),  intent(OUT) :: tmat(:,:)
   integer,       intent(IN)  :: n
   character,     intent(IN)  :: uplo
   integer,       intent(OUT) :: ierr
   !
   integer :: rank
   integer :: piv(n)
   real(dbl) :: w(n), work(2*n)
   complex(dbl) :: z(n,n)

   ierr=1
   if (size(a,1)< n .or. size(a,2)<n) return
   ierr=2
   if (size(tmat,1)< n .or. size(tmat,2)<n) return
   ierr=3
   if (uplo/="U" .and. uplo/="L") return
   !
   ierr=0
   tmat=0.0
   tmat(1:n,1:n)=a(1:n,1:n)
   !
   !call ZPOTRF(uplo,n,tmat,size(tmat,1),ierr)
   call ZPSTRF(uplo,n,tmat,size(tmat,1),piv,rank,-1.0_dbl,work,ierr)
   !
END SUBROUTINE zmat_chol


!**********************************************************
   SUBROUTINE zmat_schur( z, TT, w, a, n )
   !**********************************************************
   !
   ! Utility to decompose in Schur form a (generally non-hermitean) A matrix
   !
   ! A = Z * TT * Z^dag
   !
   IMPLICIT NONE
   COMPLEX(dbl),        INTENT(IN)  :: a(:,:)
   COMPLEX(dbl),        INTENT(OUT) :: z(:,:)
   COMPLEX(dbl),        INTENT(OUT) :: TT(:,:)
   COMPLEX(dbl),        INTENT(OUT) :: w(:)
   INTEGER,             INTENT(in)  :: n

   INTEGER   :: ierr, info, lwork, sdim
   COMPLEX(dbl), ALLOCATABLE :: work(:)
   COMPLEX(dbl), ALLOCATABLE :: a_(:,:)
   REAL(dbl),    ALLOCATABLE :: rwork(:)
   LOGICAL,      ALLOCATABLE :: bwork(:)

   ! get the dimension of the problem
   IF ( n <= 0 ) CALL errore('zmat_schur','Invalid N',ABS(n)+1)
   IF ( n > SIZE(a,1) .OR. n > SIZE(a,2) ) &
        CALL errore('zmat_schur','Invalid A dimensions',ABS(n)+1)
   IF ( n > SIZE(z,1) .OR. n > SIZE(z,2) ) &
        CALL errore('zmat_schur','Invalid Z dimensions',ABS(n)+1)

   lwork = 2 * n
   ALLOCATE( work(lwork), STAT=ierr )
   IF(ierr/=0) CALL errore('zmat_schur','allocating work',ABS(ierr))
   ALLOCATE( rwork(n), STAT=ierr )
   IF(ierr/=0) CALL errore('zmat_schur','allocating rwork',ABS(ierr))
   ALLOCATE( bwork(1), STAT=ierr )
   IF(ierr/=0) CALL errore('zmat_schur','allocating rwork',ABS(ierr))
   ALLOCATE( a_(n,n), STAT=ierr )
   IF(ierr/=0) CALL errore('zmat_schur','allocating a_',ABS(ierr))
   !
   a_(1:n,1:n) = a(1:n,1:n)
   !
   ! main driver
   !
   CALL ZGEES( 'V', 'N', .true., n, a_, SIZE(a_,1), sdim, w, z, SIZE(z,1), &
               work, lwork, rwork, bwork, info)
   IF ( info < 0 ) CALL errore('zmat_schur', 'zgees: info illegal value', abs(info) )
   !
   TT(1:n,1:n)=a_
   !
   DEALLOCATE(a_)
   DEALLOCATE(bwork)
   DEALLOCATE(rwork)
   DEALLOCATE(work)
   !
END SUBROUTINE zmat_schur

END MODULE util_module


