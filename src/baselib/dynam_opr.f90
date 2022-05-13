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
   MODULE dynam_opr_module
   !*********************************************
   !
   USE kinds,             ONLY : dbl
   USE poles_base_module
   !
   IMPLICIT NONE
   PRIVATE

! This module contains the definition of DYNAM_OPR_TYPE
! 
! routines in this module:
! SUBROUTINE  
! </INFO>
!


   TYPE dynam_opr_type
      !
      INTEGER                     :: npoles
      INTEGER                     :: npolesx
      REAL(dbl),         POINTER  :: rpoles(:)
      !
      TYPE( pole_type ), POINTER  :: root
      TYPE( pole_type ), POINTER  :: last
      !
      LOGICAL                     :: alloc 
      !
   END TYPE dynam_opr_type
      
   INTEGER, PARAMETER :: NUM_MAX_POLES = 100

!
! end of declarations
!
   PUBLIC ::  dynam_opr_type
   PUBLIC ::  dynam_opr_init
   PUBLIC ::  dynam_opr_deallocate
   !
   PUBLIC ::  dynam_opr_pole_insert
!   PUBLIC ::  dynam_opr_pole_remove
!   PUBLIC ::  dynam_opr_sort


CONTAINS

!
! Subroutines
!   


!**********************************************************
   SUBROUTINE dynam_opr_init( opr, npolesx )
   !**********************************************************
      IMPLICIT NONE
      TYPE (dynam_opr_type)         :: opr
      INTEGER, OPTIONAL, INTENT(IN) :: npolesx
      !
      CHARACTER(14) :: subname="dynam_opr_init"
      INTEGER :: npolesx_
      INTEGER :: ierr
      !
      npolesx_ = NUM_MAX_POLES
      IF ( PRESENT( npolesx ) ) npolesx_ = npolesx
      !
      opr%npoles  = 0
      opr%npolesx = npolesx_
      !
      ALLOCATE( opr%rpoles(npolesx_), STAT=ierr )
      IF ( ierr/=0 ) CALL errore(subname,"allocating rpoles",ABS(ierr))
      !
      NULLIFY( opr%root )
      NULLIFY( opr%last )
      !
      opr%alloc = .TRUE.
      !
  END SUBROUTINE dynam_opr_init


!**********************************************************
   SUBROUTINE dynam_opr_deallocate( opr )
   !**********************************************************
      IMPLICIT NONE
      TYPE (dynam_opr_type),  INTENT(INOUT)  :: opr
      !
      CHARACTER(20) :: subname = "dynam_opr_deallocate"
      INTEGER       :: ierr
      
      IF ( .NOT. opr%alloc ) RETURN
      !
      IF ( ASSOCIATED( opr%rpoles ) ) THEN
          DEALLOCATE( opr%rpoles, STAT=ierr)
          IF (ierr/=0) CALL errore(subname,"deallocating rpoles",ABS(ierr))
      ENDIF
      !
      DO WHILE ( opr%npoles > 0 )
          CALL dynam_opr_pole_remove( opr, opr%last )
      ENDDO
      !
      opr%alloc = .FALSE.
      !
  END SUBROUTINE dynam_opr_deallocate


!**********************************************************
   SUBROUTINE dynam_opr_pole_insert( opr, pole )
   !**********************************************************
      IMPLICIT NONE
      TYPE (dynam_opr_type),  INTENT(INOUT)  :: opr
      TYPE (pole_type),       INTENT(IN)     :: pole
      !
      CHARACTER(21)  :: subname="dynam_opr_pole_insert"
      !
      IF ( .NOT. opr%alloc )   CALL errore(subname,"opr not alloc",10)
      IF ( .NOT. pole%alloc )  CALL errore(subname,"pole not alloc",10)
      !
      IF ( opr%npoles == opr%npolesx ) &
          CALL errore(subname,"max number of poles reached", opr%npolesx )
      !
      IF ( opr%npoles == 0 ) THEN
          !
          ALLOCATE( opr%root )    
          opr%last => opr%root
          !
      ELSE
          !
          ALLOCATE( opr%last%next )    
          opr%last => opr%last%next
          !
      ENDIF
      !
      opr%last = pole
      !
      opr%npoles = opr%npoles + 1
      !
      opr%rpoles( opr%npoles ) = pole%rpole
      !
      RETURN
      !
   END SUBROUTINE dynam_opr_pole_insert
      
END MODULE dynam_opr_module

