!
! Copyright (C) 2004 Andrea Ferretti
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License\'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! <INFO>
!*********************************************
   MODULE energy_grid_module
!*********************************************
   !
   USE io_module,          ONLY : ionode, ionode_id
   USE ptk_module,         ONLY : ptk_bcast
   USE constants,          ONLY : ZERO
   USE parameters,         ONLY : nstrx, toll => toll_egrid
   USE utilities_module,   ONLY : locate_nearest_index, vect_compare
   USE parser_module,      ONLY : change_case
   USE iotk_module
   IMPLICIT NONE
   PRIVATE

! This MODULE contains the definition of ENERGY_GRID type
! and some routines to manage it;
! 
! routines in this module:
! SUBROUTINE energy_grid_allocate(grid,grid_type,obj[,iEfermi][,delta])
! SUBROUTINE energy_grid_setup(emin,emax,nint .OR. de,obj[,efermi][,delta])
! SUBROUTINE energy_grid_update(obj[,efermi][,delta])
! SUBROUTINE energy_grid_assignment(obj1,obj2)
! FUNCTION   energy_grid_compare(obj1,obj2) 
! SUBROUTINE energy_grid_deallocate(obj)
! SUBROUTINE energy_grid_write(unit,name,obj[,fmt])
! SUBROUTINE energy_grid_read(unit,name,obj,found)
!
! </INFO>


! <INFO>
! TYPE definitions
!
   TYPE energy_grid
      INTEGER              :: nint            ! number of intervals
      INTEGER              :: iEfermi         ! Fermi energy index: IF NOT SET -> 0
      REAL                 :: Efermi          ! Fermi energy: IF NOT SET -> 0.0
      INTEGER              :: ieps            ! index for eps
      REAL                 :: eps             ! the grid point closest to zero 
      REAL                 :: Emin            ! energy minimum
      REAL                 :: Emax            ! energy maximum
      REAL                 :: de              ! energy step if the grid is homogeneus
      REAL                 :: delta           ! imaginary part of the energy
      REAL, POINTER        :: grid(:)         ! real energy grid
      CHARACTER(15)        :: grid_type       ! if 'UNIFORM'       -> uniform grid
                                              !    'NOT_UNIFORM'   -> general grid 
      LOGICAL              :: alloc=.FALSE.   ! allocation status
   END TYPE energy_grid
!
! </INFO>
! end TYPE definition
!   

   INTERFACE ASSIGNMENT(=)
      MODULE PROCEDURE energy_grid_assignment
   END INTERFACE
   INTERFACE OPERATOR(==)
      MODULE PROCEDURE energy_grid_compare
   END INTERFACE

   INTERFACE energy_grid_setup
      MODULE PROCEDURE energy_grid_setup_nint
      MODULE PROCEDURE energy_grid_setup_de
   END INTERFACE

   PUBLIC ::  energy_grid,             &
              ASSIGNMENT(=),           &
              OPERATOR(==),            &
              energy_grid_allocate,    &
              energy_grid_deallocate,  &
              energy_grid_write,       &
              energy_grid_read,        &
              energy_grid_setup,       &
              energy_grid_update


!
! end of declarations
!

CONTAINS 
   
!
! Subroutines
!

!**********************************************************
   SUBROUTINE energy_grid_allocate(grid,grid_type,obj,iEfermi,delta)
   !**********************************************************
      IMPLICIT NONE 
       
      REAL, DIMENSION(:), INTENT(in)   :: grid
      CHARACTER(LEN=*),   INTENT(in)   :: grid_type
      TYPE(energy_grid),  INTENT(out)  :: obj
      INTEGER, OPTIONAL,  INTENT(in)   :: iEfermi 
      REAL   , OPTIONAL,  INTENT(in)   :: delta
      INTEGER                          :: nint,ierr

      nint=SIZE(grid(:)) 
      ALLOCATE( obj%grid(nint), STAT=ierr )
      IF ( ierr /= 0 ) &
          CALL errore('energy_grid_allocate','allocation of obj%grid failed',1)
      IF ( PRESENT( iEfermi) .AND. ( iEfermi < 0 .OR. iEfermi > nint ) )  &
          CALL errore('energy_grid_allocate','iEfermi out range',1)
         

      obj%iEfermi=0   
      IF ( PRESENT(iEfermi) )  obj%iEfermi=iEfermi
      obj%Efermi=0.0
      IF ( obj%iEfermi /= 0 )  obj%Efermi=grid(obj%iEfermi)
      obj%delta=0.0   
      IF ( PRESENT(delta) )  obj%delta=delta
      obj%nint=nint 
      obj%grid(:)= grid(:)
      obj%Emin=grid(1)
      obj%Emax=grid(nint)

      !
      ! set grid_type related quantities
      !
      SELECT CASE ( grid_type )
      CASE ( 'UNIFORM' )
          obj%de=( obj%Emax - obj%Emin )/REAL(nint-1)        
      CASE ( 'NOT_UNIFORM' )
          obj%de=0.0
      CASE DEFAULT
          CALL errore('energy_grid_allocate','Invalid grid type',1)
      END SELECT
      obj%grid_type=grid_type 

      !
      ! if emin <= zero <= emax,
      ! setup the eps value as the grid point closest to zero (with its sign)
      ! otherwise eps is set to MOD( grid(1), de )
      !
      IF ( obj%emin <= 0.0 .AND. obj%emax >= 0.0 ) THEN
          !
          CALL locate_nearest_index( obj%grid, 0.0, obj%ieps, IERR=ierr)
          IF ( ierr/=0 ) CALL errore('energy_grid_allocate', 'searching for ZERO',ABS(ierr))
          !
          obj%eps =  obj%grid( obj%ieps ) 
          !
      ELSE
          !
          obj%ieps = -1
          !
          obj%eps  = 0.0 
          IF ( TRIM(obj%grid_type) == 'UNIFORM' ) &
               obj%eps = MOD( obj%grid(1), obj%de )
          !
      ENDIF

      !
      obj%alloc=.TRUE.
      !
   END SUBROUTINE energy_grid_allocate


!**********************************************************
   SUBROUTINE energy_grid_setup_nint(emin,emax,nint,obj,efermi,delta)
   !**********************************************************
      IMPLICIT NONE 
      REAL,              INTENT(in)     :: emin
      REAL,              INTENT(in)     :: emax
      INTEGER,           INTENT(in)     :: nint
      TYPE(energy_grid), INTENT(inout)  :: obj
      REAL,    OPTIONAL, INTENT(in)     :: efermi
      REAL,    OPTIONAL, INTENT(in)     :: delta
  
      REAL                              :: de
      REAL, ALLOCATABLE                 :: grid(:)
      INTEGER                           :: i, iefermi, ierr

      IF (  obj%alloc ) CALL energy_grid_deallocate(obj)
      IF ( emin > emax ) &
          CALL errore('energy_grid_setup','EMIN is greater than EMAX',1)
      IF ( nint <= 1 )   &
          CALL errore('energy_grid_setup','Invalid NINT on input',1)
 
      ALLOCATE( grid(nint) )
      de = ( emax - emin ) / REAL(nint -1)
      DO i=1,nint
         grid(i) = emin + REAL(i-1) * de
      ENDDO
 
      iefermi = 0 
      IF ( PRESENT(efermi) ) THEN 
          !
          CALL locate_nearest_index(grid, efermi, iefermi, TOLL=de*0.5, IERR=ierr)
          IF (ierr/=0 ) CALL errore('energy_grid_setup','setting efermi',ABS(ierr))
          !
      ENDIF

      IF ( PRESENT(delta) ) THEN
          CALL energy_grid_allocate(grid,"UNIFORM",obj,IEFERMI=iefermi,DELTA=delta)
      ELSE
          CALL energy_grid_allocate(grid,"UNIFORM",obj,IEFERMI=iefermi)
      ENDIF
      !
      DEALLOCATE( grid )
      !
   END SUBROUTINE energy_grid_setup_nint


!**********************************************************
   SUBROUTINE energy_grid_setup_de(emin,emax,de,obj,efermi,delta)
   !**********************************************************
      IMPLICIT NONE 
      REAL,              INTENT(in)     :: emin
      REAL,              INTENT(in)     :: emax
      REAL,              INTENT(in)     :: de
      TYPE(energy_grid), INTENT(inout)  :: obj
      REAL,    OPTIONAL, INTENT(in)     :: efermi
      REAL,    OPTIONAL, INTENT(in)     :: delta
  
      INTEGER                           :: ninte
      REAL, ALLOCATABLE                 :: grid(:)
      INTEGER                           :: i

      IF (  obj%alloc ) CALL energy_grid_deallocate(obj)
      IF ( emin > emax ) &
          CALL errore('energy_grid_setup','EMIN is greater than EMAX',1)
      IF ( de < 0.0 )   &
          CALL errore('energy_grid_setup','Invalid DE on input',1)
 
      ninte = NINT( (emax - emin)/de ) + 1

      IF ( ninte <= 1 ) CALL errore('energy_grid_setup','Invalid DE on input',1)

      IF ( PRESENT(efermi) .AND. PRESENT(delta) ) THEN
          CALL energy_grid_setup_nint(emin,emax,ninte,obj,EFERMI=efermi,DELTA=delta)
      ELSEIF ( PRESENT(efermi) ) THEN
          CALL energy_grid_setup_nint(emin,emax,ninte,obj,EFERMI=efermi)
      ELSEIF ( PRESENT(delta) )  THEN
          CALL energy_grid_setup_nint(emin,emax,ninte,obj,DELTA=delta)
      ELSE
          CALL energy_grid_setup_nint(emin,emax,ninte,obj)
      ENDIF
   END SUBROUTINE energy_grid_setup_de


!**********************************************************
   SUBROUTINE energy_grid_update(obj,iefermi,efermi,delta)
   !**********************************************************
      IMPLICIT NONE 
      TYPE(energy_grid), INTENT(inout)  :: obj
      INTEGER, OPTIONAL, INTENT(in)     :: iefermi
      REAL,    OPTIONAL, INTENT(in)     :: efermi
      REAL,    OPTIONAL, INTENT(in)     :: delta
      CHARACTER(18)   :: sub_name='energy_grid_update'
      REAL            :: mean_de
      INTEGER         :: i, ierr

      IF ( .NOT. obj%alloc )   &
          CALL errore(sub_name,'Energy_grid NOT allocated',1)
      IF ( PRESENT(iefermi) ) THEN
          obj%iefermi=iefermi
          obj%efermi = obj%grid(iefermi)
      ENDIF
      !
      IF ( PRESENT(iefermi) .AND. PRESENT(efermi) ) &
          CALL errore(sub_name,'IEFERMI and EFERMI cannot be both present',1)
      !
      IF ( PRESENT(efermi) ) THEN
          !
          mean_de = (obj%emax - obj%emin)/REAL(obj%nint -1)
          !
          CALL locate_nearest_index(obj%grid, efermi, i, TOLL=mean_de, IERR=ierr)
          IF ( ierr/=0 ) CALL errore(sub_name,'setting efermi',ABS(ierr))
          !
          obj%iefermi=i
          obj%efermi = obj%grid(i)
          !
      ENDIF
      IF ( PRESENT(delta) ) THEN
          obj%delta = delta
      ENDIF
   END SUBROUTINE energy_grid_update


!**********************************************************
   SUBROUTINE energy_grid_deallocate(obj)
   !**********************************************************
      IMPLICIT NONE 
      TYPE(energy_grid), INTENT(inout)      :: obj
      IF ( .NOT. obj%alloc )   &
          CALL errore('energy_grid_deallocate','Energy_grid NOT allocated',1)
      IF ( ASSOCIATED(obj%grid) ) DEALLOCATE(obj%grid)
      NULLIFY(obj%grid)
      obj%iEfermi=0
      obj%Efermi=0.0
      obj%ieps=0
      obj%eps=0.0
      obj%delta=0.0
      obj%nint=0
      obj%Emin=0.0
      obj%Emax=0.0
      obj%de=0.0
      obj%grid_type=""
      obj%alloc=.FALSE.
   END SUBROUTINE energy_grid_deallocate


!**********************************************************
   SUBROUTINE energy_grid_assignment(obj1,obj2)
   !**********************************************************
      IMPLICIT NONE 
      TYPE(energy_grid), INTENT(inout)      :: obj1
      TYPE(energy_grid), INTENT(in)         :: obj2

      IF ( .NOT. obj2%alloc )   &
          CALL errore('energy_grid_assignment','RHS Energy_grid NOT allocated',1)
      IF ( obj1%alloc )  CALL energy_grid_deallocate(obj1)
   
      IF ( .NOT. ASSOCIATED(obj2%grid) ) &
          CALL errore('energy_grid_assignment','Unexpected error on RHS Energy_grid',1)
      CALL energy_grid_allocate(obj2%grid,obj2%grid_type,obj1,IEFERMI=obj2%iefermi,  &
                                                              DELTA=obj2%delta       )
   END SUBROUTINE energy_grid_assignment


!**********************************************************
   LOGICAL FUNCTION energy_grid_compare(obj1,obj2)
   !**********************************************************
      IMPLICIT NONE 
      TYPE(energy_grid), INTENT(in)         :: obj1
      TYPE(energy_grid), INTENT(in)         :: obj2
      !
      ! If one of both obj's are not allocated, RETURN false
      !

      energy_grid_compare=.FALSE.
      !
      IF ( .NOT. obj1%alloc ) RETURN 
      IF ( .NOT. obj2%alloc ) RETURN
      IF ( obj1%nint /= obj2%nint ) RETURN
      IF ( obj1%grid_type /= obj2%grid_type )  RETURN
      IF ( obj1%iEfermi /= obj2%iEfermi )  RETURN
      IF ( obj1%delta /= obj2%delta )  RETURN
      CALL vect_compare(obj1%grid,obj2%grid,energy_grid_compare,TOLL=toll)
   END FUNCTION  energy_grid_compare     
      

!**********************************************************
   SUBROUTINE energy_grid_write(unit,name,obj,fmt)
   !**********************************************************
   !
   ! allowed fmt are:
   ! 'data' | 'summary'
   !
   IMPLICIT NONE 
      !
      INTEGER,                  INTENT(in) :: unit
      CHARACTER(*),             INTENT(in) :: name 
      CHARACTER(*), OPTIONAL,   INTENT(in) :: fmt
      TYPE(energy_grid),        INTENT(in) :: obj

      CHARACTER(15)                        :: fmt_, formatted
      INTEGER                              :: ierr
      CHARACTER(nstrx)                     :: str

      IF ( .NOT. obj%alloc ) &
          CALL errore('energy_grid_write','Energy_grid NOT allocated',1)

      IF ( ionode ) THEN
          !
          fmt_="data"
          IF ( PRESENT(fmt) ) THEN
               !
               ! set input value
               !
               fmt_ = TRIM(fmt)
               !
               ! check whether unit is formatted or not
               !
               INQUIRE(unit,FORMATTED=formatted)
               CALL change_case(formatted,'upper')
               !
               IF ( TRIM(formatted) /= "YES" .AND. TRIM(formatted) /= "UNKNOWN" ) &
                  CALL errore('energy_grid_write','FMT data not allowed for binary write', 2)
               !
          ENDIF
          !
          CALL change_case(fmt_,'lower')
      
          SELECT CASE ( TRIM(fmt_) )
          CASE ( "data" ) 
              str=" "
              CALL iotk_write_attr(str,"type","energy_grid")
              CALL iotk_write_begin(unit,name,ATTR=str)
              str=" "
              CALL iotk_write_attr(str,"nint",obj%nint)
              CALL iotk_write_attr(str,"iefermi",obj%iefermi)
              IF ( obj%iefermi /= 0) CALL iotk_write_attr(str,"efermi",obj%efermi)
              CALL iotk_write_attr(str,"de",obj%de)
              CALL iotk_write_attr(str,"emin",obj%emin)
              CALL iotk_write_attr(str,"emax",obj%emax)
              CALL iotk_write_attr(str,"delta",obj%delta)
              CALL iotk_write_attr(str,"grid_type",TRIM(obj%grid_type))
              CALL iotk_write_empty(unit,"DATA",ATTR=str)
    
              CALL iotk_write_dat(unit,"GRID",obj%grid,FMT="(6f15.9)",IERR=ierr)
                 IF ( ierr /= 0 ) &
                     CALL errore('energy_grid_write','Unable to write data',abs(ierr))
              CALL iotk_write_end(unit,name)
          CASE ( "summary" ) 
              IF ( TRIM(obj%grid_type) == "UNIFORM" ) THEN
                  WRITE(unit,"(7x,a6,' Grid [eV]:  Uniform grid  ( step dE = ',f10.4, &
                              &' )')") TRIM(name), obj%de
              ELSE
                  WRITE(unit,"(7x,a6,' Grid [eV]:  Variable-step grid' )") TRIM(name)
              ENDIF
              WRITE(unit,"(12x,'emin = ',f10.4,',',7x, 'Fermi energy = ',f10.4,',',/,  "// &
                         " 12x,'emax = ',f10.4,',',7x, 'EFermi index = ',i5,',',/,    " // &
                         " 12x,'nint = ',i5,   ',',12x,'       E_eps = ',f10.4,',',/,  "// &
                         " 37x,                        '       delta = ',f10.4, /)" )  &
                         obj%emin, obj%efermi, obj%emax, obj%iefermi, obj%nint, obj%eps, obj%delta
    
          CASE DEFAULT
              CALL errore('energy_grid_write','Invalid FMT = '//TRIM(fmt_), 2)
          END SELECT
          !
      ENDIF
      !
   END SUBROUTINE energy_grid_write


!**********************************************************
   SUBROUTINE energy_grid_read(unit,name,obj,found)
   !**********************************************************
      IMPLICIT NONE 
      INTEGER,                  INTENT(in) :: unit
      CHARACTER(*),             INTENT(in) :: name 
      TYPE(energy_grid),     INTENT(inout) :: obj
      LOGICAL,                 INTENT(out) :: found

      CHARACTER(17)                        :: sub_name='energy_grid_read'
      CHARACTER(nstrx)                     :: val_str, grid_type, str
      INTEGER                              :: nint, iefermi, ierrl
      REAL                                 :: delta
      REAL, ALLOCATABLE                    :: grid(:)

      IF ( obj%alloc ) CALL energy_grid_deallocate(obj)
      
      IF ( ionode ) THEN 
           CALL iotk_scan_begin(unit,name,ATTR=str,FOUND=found,IERR=ierrl)
           IF (ierrl > 0) CALL errore(sub_name,'Error in scanning for '//TRIM(name),ierrl)
      ENDIF
      !
      CALL ptk_bcast( found, ionode_id )
      IF ( .NOT. found ) RETURN
      
      IF ( ionode ) THEN
         !
         CALL iotk_scan_attr(str,"type",val_str,IERR=ierrl)  
         IF (ierrl /= 0) CALL errore(sub_name,'Wrong input format in attr TYPE',abs(ierrl))
         IF ( TRIM(val_str) /= 'energy_grid' ) CALL errore(sub_name,'Invalid TYPE on input',1)
         !
         CALL iotk_scan_empty(unit,"DATA",str,IERR=ierrl)
         IF (ierrl /= 0) CALL errore(sub_name,'reading block DATA',abs(ierrl))
         !
         CALL iotk_scan_attr(str,"nint",nint,IERR=ierrl)  
         IF (ierrl /= 0) CALL errore(sub_name,'reading attr NINT',abs(ierrl))
         !
         CALL iotk_scan_attr(str,"iefermi",iefermi,IERR=ierrl)  
         IF (ierrl /= 0) CALL errore(sub_name, 'reading attr IEFERMI',abs(ierrl))
         !
         CALL iotk_scan_attr(str,"delta",delta,IERR=ierrl)  
         IF (ierrl /= 0) CALL errore(sub_name,'reading attr DELTA',abs(ierrl))
         !
         CALL iotk_scan_attr(str,"grid_type",grid_type,IERR=ierrl)  
         IF (ierrl /= 0) CALL errore(sub_name, 'reading attr GRID_TYPE',abs(ierrl))
         !
      ENDIF
      !
      CALL ptk_bcast( nint,      ionode_id )
      CALL ptk_bcast( iefermi,   ionode_id )
      CALL ptk_bcast( delta,     ionode_id )
      CALL ptk_bcast( grid_type, ionode_id )
      !
      ALLOCATE( grid(nint), STAT=ierrl ) 
      IF ( ierrl/=0 ) CALL errore(sub_name, 'allocating grid',ABS(ierrl))
      !
      IF ( ionode ) THEN
         !
         grid(:) = ZERO
         !
         CALL iotk_scan_dat(unit,'GRID',grid(:),FOUND=found,IERR=ierrl)
         IF ( .NOT. found )  CALL errore(sub_name,'GRID dat NOT found on input',1)
         IF (ierrl > 0) CALL errore(sub_name,'Wrong input format in GRID',abs(ierrl))
         !
         CALL iotk_scan_end(unit,name,IERR=ierrl)   
         IF ( ierrl /= 0 ) CALL errore(sub_name,'Unable to end '//TRIM(name),1)
         !
      ENDIF
      !
      CALL ptk_bcast( grid, ionode_id )
      ! 
      CALL energy_grid_allocate(grid,TRIM(grid_type),obj,IEFERMI=iefermi,DELTA=delta)        
      !
      DEALLOCATE( grid, STAT=ierrl)
      IF ( ierrl/=0 ) CALL errore(sub_name, 'deallocating grid', ABS(ierrl))
      !
   END SUBROUTINE energy_grid_read

END MODULE energy_grid_module


