!
! Copyright (C) 2010   A. Ferretti
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License\'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!**********************************************************
   SUBROUTINE cleandir( dirname, fileprefix, filesuffix )
   !**********************************************************
   !
   ! remove files
   !
   !      dirname/<fileroot><num><filesuffix>
   !
   USE kinds
   USE files_module,        ONLY : file_exist, file_delete
   USE parser_base_module,  ONLY : int2char
   USE timing_module
   !
   IMPLICIT NONE

   !
   ! input variables
   !
   CHARACTER(*),   INTENT(IN) :: dirname
   CHARACTER(*),   INTENT(IN) :: fileprefix, filesuffix

   !
   ! local variables
   !
   INTEGER        :: nfail_max = 50
   INTEGER        :: nfail
   CHARACTER(9)   :: subname='cleandir'
   CHARACTER(256) :: filename
   LOGICAL        :: found
   INTEGER        :: i, ierr

!
!------------------------------
! main body
!------------------------------
!
   i = 0
   nfail = 0
   found = .TRUE.
   
   DO WHILE ( nfail <= nfail_max ) 
       i = i+1
       !
       filename = TRIM( dirname ) // "/" // TRIM(fileprefix) // &
                  TRIM(int2char(i)) // TRIM(filesuffix)
       !
       found = file_exist( filename )
       !
       IF ( .NOT. found ) THEN 
           nfail = nfail + 1
           IF ( nfail > nfail_max ) EXIT
       ELSE
           nfail = 0
           CALL file_delete( filename )
       ENDIF
       ! 
   ENDDO
   !
   RETURN
   ! 
END SUBROUTINE cleandir

