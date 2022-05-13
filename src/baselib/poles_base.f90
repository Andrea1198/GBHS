! 
! Copyright (C) 2012 Andrea Ferretti
! 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 
! <INFO>
!*********************************************
   MODULE poles_base_module
   !*********************************************
   !
   USE kinds,             ONLY : dbl
   USE constants,         ONLY : CI
   !
   IMPLICIT NONE
   PRIVATE

! This module contains the definition of POLE type, and POLES_LIST
! 
! routines in this module:
! SUBROUTINE  pole_allocate()
! </INFO>
!


   TYPE pole_type
      !
      REAL(dbl)                   :: rpole             ! real part of the pole position
      REAL(dbl)                   :: ipole             ! imaginary part of the pole
      COMPLEX(dbl)                :: cpole             ! cmplx pole (rpole + i ipole)
      LOGICAL                     :: is_occupied
      REAL(dbl)                   :: weight
      !
      LOGICAL                     :: is_zero
      LOGICAL                     :: full_opr
      INTEGER                     :: ndim              ! dim of the space
      COMPLEX(dbl), POINTER       :: R(:,:)            ! spectral func operator 
      !
      INTEGER                     :: nvec              ! # of vectors to represent R
      COMPLEX(dbl), POINTER       :: R_vec(:,:)        ! vectors used to represent R
      REAL(dbl),    POINTER       :: R_eig(:)          ! eigs corresponding to the above vects
      !
      TYPE( pole_type ), POINTER  :: next
      TYPE( pole_type ), POINTER  :: prev
      !
      LOGICAL                     :: alloc 
   END TYPE pole_type
      

!
! end of declarations
!
   INTERFACE ASSIGNMENT(=)
      MODULE PROCEDURE pole_assignment
   END INTERFACE

   PUBLIC ::  pole_type, ASSIGNMENT(=)
   PUBLIC ::  pole_init
   !
   PUBLIC ::  pole_allocate
   PUBLIC ::  pole_deallocate
   PUBLIC ::  pole_update
   !
   PUBLIC ::  pole_convolution
   PUBLIC ::  pole_merge
!   !
!   PUBLIC ::  pole_insert_list
!   PUBLIC ::  pole_remove_list


CONTAINS

!
! Subroutines
!   


!**********************************************************
   SUBROUTINE pole_init( pole )
   !**********************************************************
      IMPLICIT NONE
      TYPE (pole_type)  :: pole
      !
      pole%ndim = 0
      pole%weight = 0.0d0
      NULLIFY( pole%R ) 
      pole%nvec = 0
      NULLIFY( pole%R_vec ) 
      NULLIFY( pole%R_eig ) 
      !
      NULLIFY( pole%prev )
      NULLIFY( pole%next )
      !
      pole%alloc = .FALSE.
      !
  END SUBROUTINE pole_init


!**********************************************************
   SUBROUTINE pole_allocate( pole, ndim, is_zero, full_opr, nvec)
   !**********************************************************
      IMPLICIT NONE
      TYPE (pole_type),  INTENT(INOUT)  :: pole
      INTEGER,           INTENT(IN)     :: ndim
      LOGICAL, OPTIONAL, INTENT(IN)     :: is_zero
      LOGICAL, OPTIONAL, INTENT(IN)     :: full_opr
      INTEGER, OPTIONAL, INTENT(IN)     :: nvec
      !
      CHARACTER(13) :: subname = "pole_allocate"
      LOGICAL       :: is_zero_  = .FALSE.
      LOGICAL       :: full_opr_ = .TRUE.
      INTEGER       :: ierr
      
      IF ( pole%alloc ) CALL errore(subname,"pole already alloc", 10)
      !
      IF ( PRESENT( full_opr) ) full_opr_ = full_opr
      IF ( PRESENT( is_zero) )  is_zero_  = is_zero
      !
      IF ( is_zero_ .AND. ( PRESENT(full_opr) .OR. PRESENT(nvec) ) ) &
          CALL errore(subname,"not sure you want to do that",10)
      !
      pole%weight   = 0.0d0
      pole%full_opr = full_opr_
      pole%ndim     = ndim
      pole%is_zero  = is_zero_
      !
      !
      IF ( .NOT. is_zero_ ) THEN
          !
          IF ( full_opr_ ) THEN
              !
              ALLOCATE( pole%R(ndim,ndim), STAT=ierr )
              IF ( ierr/=0 ) CALL errore(subname,"allocating R",ABS(ierr))
              !
          ELSE
              IF ( .NOT. PRESENT( nvec ) ) CALL errore(subname,"missing nvec",10)
              !
              pole%nvec = nvec
              !
              ALLOCATE( pole%R_vec(ndim,nvec), STAT=ierr )
              IF ( ierr/=0 ) CALL errore(subname,"allocating R_vec",ABS(ierr))
              ALLOCATE( pole%R_eig(nvec), STAT=ierr )
              IF ( ierr/=0 ) CALL errore(subname,"allocating R_eig",ABS(ierr))
              !
          ENDIF
          !
      ENDIF
      !  
      pole%alloc = .TRUE.
      !
  END SUBROUTINE pole_allocate


!**********************************************************
   SUBROUTINE pole_deallocate( pole )
   !**********************************************************
      IMPLICIT NONE
      TYPE (pole_type),  INTENT(INOUT)  :: pole
      !
      CHARACTER(15) ::   subname = "pole_deallocate"
      INTEGER       :: ierr
      
      IF ( .NOT. pole%alloc ) RETURN
      !
      IF ( ASSOCIATED( pole%R ) ) THEN
          DEALLOCATE( pole%R, STAT=ierr)
          IF (ierr/=0) CALL errore(subname,"deallocating R",ABS(ierr))
      ENDIF
      IF ( ASSOCIATED( pole%R_vec ) ) THEN
          DEALLOCATE( pole%R_vec, STAT=ierr)
          IF (ierr/=0) CALL errore(subname,"deallocating R_vec",ABS(ierr))
      ENDIF
      IF ( ASSOCIATED( pole%R_eig ) ) THEN
          DEALLOCATE( pole%R_eig, STAT=ierr)
          IF (ierr/=0) CALL errore(subname,"deallocating R_eig",ABS(ierr))
      ENDIF
      !
      pole%alloc = .FALSE.
      CALL pole_init( pole )
      !
  END SUBROUTINE pole_deallocate


!**********************************************************
   SUBROUTINE pole_assignment( pole1, pole2 )
   !**********************************************************
      IMPLICIT NONE
      TYPE (pole_type),  INTENT(INOUT) :: pole1
      TYPE (pole_type),  INTENT(IN)    :: pole2
      !
      CHARACTER(15) :: subname="pole_assignment"
  
      IF ( .NOT. pole2%alloc ) CALL errore(subname,'pole2 not allocated',1)
      !
      IF ( pole1%alloc )       CALL pole_deallocate(pole1)
      CALL pole_allocate( pole1, pole2%ndim, pole2%is_zero, pole2%full_opr, pole2%nvec )
      !
      pole1%rpole = pole2%rpole
      pole1%ipole = pole2%ipole
      pole1%cpole = pole2%cpole
      pole1%is_occupied = pole2%is_occupied
      pole1%is_zero     = pole2%is_zero
      pole1%weight      = pole2%weight
      !
      IF ( .NOT. pole1%is_zero ) THEN
          !
          IF ( pole2%full_opr ) THEN 
              pole1%R = pole2%R
          ELSE
              pole1%R_vec = pole2%R_vec
              pole1%R_eig = pole2%R_eig
          ENDIF
          !
      ENDIF
      !
      pole1%next => pole2%next
      pole1%prev => pole2%prev
      !
  END SUBROUTINE pole_assignment


!**********************************************************
   SUBROUTINE pole_update( pole, cpole, is_zero, ndim, R, nvec, R_vec, R_eig)
   !**********************************************************
      IMPLICIT NONE
      TYPE (pole_type),     INTENT(INOUT) :: pole
      !
      COMPLEX(dbl),  OPTIONAL, INTENT(IN) :: cpole
      LOGICAL,       OPTIONAL, INTENT(IN) :: is_zero
      INTEGER,       OPTIONAL, INTENT(IN) :: ndim, nvec
      COMPLEX(dbl),  OPTIONAL, INTENT(IN) :: R(:,:)
      COMPLEX(dbl),  OPTIONAL, INTENT(IN) :: R_vec(:,:)
      REAL(dbl),     OPTIONAL, INTENT(IN) :: R_eig(:)
      !
      CHARACTER(11) :: subname="pole_update"
  
      IF ( .NOT. pole%alloc ) CALL errore(subname,'pole not allocated',1)
      !
      IF ( PRESENT( cpole ) ) THEN
          !
          pole%rpole = REAL( cpole, dbl)
          pole%ipole = AIMAG( cpole )
          pole%cpole = cpole
          !
          pole%is_occupied = .FALSE.
          IF ( pole%ipole >= 0.0_dbl ) pole%is_occupied = .TRUE.
          !
      ENDIF
      !
      IF ( PRESENT( is_zero ) ) THEN
          pole%is_zero = .TRUE.
          pole%weight  = 0.0d0
          IF ( ASSOCIATED( pole%R ) )     DEALLOCATE( pole%R )
          IF ( ASSOCIATED( pole%R_vec ) ) DEALLOCATE( pole%R_vec )
          IF ( ASSOCIATED( pole%R_eig ) ) DEALLOCATE( pole%R_eig )
      ENDIF
      !
      IF ( PRESENT(ndim) ) THEN
          IF ( pole%ndim /= ndim )   CALL errore(subname,"invalid ndim",10)
      ENDIF
      IF ( PRESENT(nvec) ) THEN
          IF ( pole%nvec /= nvec )   CALL errore(subname,"invalid nvec",10)
      ENDIF
      !
      IF ( PRESENT( R ) ) THEN
          !
          IF ( .NOT. PRESENT(ndim) ) CALL errore(subname,"missing ndim",10)
          !
          pole%R(1:ndim,1:ndim) = R(1:ndim,1:ndim)
          pole%weight = 1.0d0
          !
      ENDIF
      !
      IF ( PRESENT( R_vec ) .OR. PRESENT( R_eig ) ) THEN
          !
          IF ( .NOT. PRESENT(nvec) ) CALL errore(subname,"missing nvec",11)
          IF ( .NOT. PRESENT(ndim) ) CALL errore(subname,"missing ndim",11)
          !
          IF ( PRESENT(R_vec) ) pole%R_vec(1:ndim,1:nvec) = R_vec(1:ndim,1:nvec)
          IF ( PRESENT(R_eig) ) pole%R_eig(1:nvec)        = R_eig(1:nvec)
          !
          pole%weight = 1.0d0
          !
      ENDIF
      !          
      RETURN 
      !
  END SUBROUTINE pole_update


!**********************************************************
   SUBROUTINE pole_convolution( pole_c, pole1, pole2, ctype, side )
   !**********************************************************
   !
   ! compute the convolution of pole1 and pole2 according to:
   !
   ! vtype      
   ! "N":       pole_c(w) = \int dw'/(2 pi i)   pole1(w')   * pole2(w-w')
   ! "C":       pole_c(w) = \int dw'/(2 pi i)   pole1(w')   * pole2(w+w')
   !
   ! side:  "Nornal" as above
   !        "Opposite"   pole2 * pole1
   ! This is important only for the spatial degrees of freedom
   !
      IMPLICIT NONE
      TYPE (pole_type),     INTENT(INOUT) :: pole_c
      TYPE (pole_type),     INTENT(IN)    :: pole1
      TYPE (pole_type),     INTENT(IN)    :: pole2
      CHARACTER(1),         INTENT(IN)    :: ctype 
      CHARACTER(1),         INTENT(IN)    :: side
      !
      CHARACTER(16) :: subname="pole_convolution"
      LOGICAL       :: conv_is_zero = .FALSE.
      COMPLEX(dbl)  :: cpole
      REAL(dbl)     :: csign
      

      IF ( pole_c%alloc ) CALL pole_deallocate( pole_c )
      !
      ! few checks
      !
      IF ( pole1%ndim /= pole2%ndim ) CALL errore(subname,"ndim differs",10)

      !
      ! check whether we need to perform the convolution
      !
      IF ( pole1%is_zero )  conv_is_zero = .TRUE.
      IF ( pole2%is_zero )  conv_is_zero = .TRUE.
      !
      SELECT CASE ( side )
      CASE ( "N",'n' )
         IF ( pole1%ipole * pole2%ipole <= 0.0_dbl )  conv_is_zero = .TRUE.
      CASE ( "O",'o' )
         IF ( pole1%ipole * pole2%ipole >= 0.0_dbl )  conv_is_zero = .TRUE.
      CASE DEFAULT 
         CALL errore(subname,"invalid SIDE: "//side, 10 )
      END SELECT
      !
      IF ( conv_is_zero ) THEN
          !
          CALL pole_allocate( pole_c, pole1%ndim, IS_ZERO=.TRUE. )
          !
          RETURN
          !
      ENDIF

      ! 
      ! perform the real calculation
      ! 
      SELECT CASE ( ctype )
      CASE ( "N", 'n')
          !
          cpole = pole1%cpole + pole2%cpole
          !
      CASE ( "C", 'c')
          !
          cpole = pole2%cpole - pole1%cpole
          !
      CASE DEFAULT 
         CALL errore(subname,"invalid CTYPE: "//ctype, 10 )
      END SELECT


      !
      ! build R1 * R2
      !
      IF ( .NOT. pole1%full_opr ) CALL errore(subname,"no-full-opr  not yet implemented", 10)
      IF ( .NOT. pole2%full_opr ) CALL errore(subname,"no-full-opr  not yet implemented", 11)
      !
      CALL pole_allocate( pole_c, pole1%ndim, FULL_OPR=.TRUE. )
      !
      SELECT CASE ( side )
      CASE ( "N", "n" ) 
          !
          ! R = R1 * R2
          !
          CALL mat_mul( pole_c%R, pole1%R, 'N', pole2%R, 'N', pole_c%ndim, pole_c%ndim, pole_c%ndim) 
          !
      CASE ( "O", "o" )
          !
          ! R = R2 * R1
          !
          CALL mat_mul( pole_c%R, pole2%R, 'N', pole1%R, 'N', pole_c%ndim, pole_c%ndim, pole_c%ndim) 
          !
      END SELECT
      !
      csign = SIGN( pole1%ipole, 1.0_dbl )
      !
      pole_c%R = pole_c%R * csign
  

      RETURN 
      !
  END SUBROUTINE pole_convolution


!**********************************************************
   SUBROUTINE pole_merge( pole_m, pole1, pole2 )
   !**********************************************************
   !
   ! merge two poles into a single one:
   !
   ! R_m = R_1 + R_2
   ! r_m = weighed average (r_1, r_2 )
   ! i_m = such that covers both original poles
   ! 
   ! basically the position of the new pole is optimized to 
   ! mimick the original two poles
   !
      IMPLICIT NONE
      TYPE (pole_type),     INTENT(INOUT) :: pole_m
      TYPE (pole_type),     INTENT(IN)    :: pole1
      TYPE (pole_type),     INTENT(IN)    :: pole2
      !
      COMPLEX(dbl)  :: cpole
      CHARACTER(10) :: subname='pole_merge'
      INTEGER       :: ierr, ndim
      !REAL(dbl)     :: cost
      COMPLEX(dbl), ALLOCATABLE :: R(:,:)


      IF ( .NOT. pole1%alloc ) CALL errore(subname,'pole1 not alloc',10)
      IF ( .NOT. pole2%alloc ) CALL errore(subname,'pole2 not alloc',10)
      !
      IF ( pole1%ndim /= pole2%ndim ) CALL errore(subname,'invalid grids',10)
      IF ( .NOT. pole1%full_opr )     CALL errore(subname,'full opr only',10)
      IF ( .NOT. pole2%full_opr )     CALL errore(subname,'full opr only',10)
      !
      IF ( pole1%is_occupied .NEQV. pole2%is_occupied ) &
           CALL errore(subname,'mixing occupied and empty states',10)
      !
      IF ( pole_m%alloc ) CALL pole_deallocate( pole_m )
      !
      ! treat trivial cases
      !
      IF ( pole1%is_zero .AND. pole2%is_zero ) CALL pole_allocate( pole_m, pole1%ndim, IS_ZERO=.TRUE.)
      IF ( pole1%is_zero )  pole_m = pole2
      IF ( pole2%is_zero )  pole_m = pole1
      ! 
      IF ( pole1%is_zero .OR. pole2%is_zero ) RETURN

      !
      ! std case
      !
      ndim = pole1%ndim
      !
      CALL pole_allocate( pole_m, ndim )
      !
      ALLOCATE( R(ndim,ndim), STAT=ierr )
      IF ( ierr/=0 ) CALL errore(subname,'allocating R',ABS(ierr))
      !
      !cost = 1.0d0/( pole1%weight + pole2%weight )
      !
      R(:,:) =  pole1%weight * pole1%R(:,:) + pole2%weight * pole2%R(:,:) 
      !
      ! add a weight here XXX
      !
      cpole = (pole1%rpole + pole2%rpole)/2.0d0 + CI * &
              ( MAX( pole1%rpole +ABS(pole1%ipole), pole2%rpole +ABS(pole2%ipole) ) - &
                MIN( pole1%rpole -ABS(pole1%ipole), pole2%rpole -ABS(pole2%ipole) ) )/2.0d0
      !
      CALL pole_update( pole_m, CPOLE=cpole, IS_ZERO=.FALSE., R=R )
      !
      DEALLOCATE( R, STAT=ierr )
      IF ( ierr/=0 ) CALL errore(subname,'deallocating R',ABS(ierr))
      ! 
      RETURN 
      !
  END SUBROUTINE pole_merge

END MODULE poles_base_module

