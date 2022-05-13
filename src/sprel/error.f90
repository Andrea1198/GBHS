!================================================
! SUBROUTINE FOR ERRORS IN ALLOCATION AND OTHERS
!================================================
! ANDREA (MR) P !
!================
SUBROUTINE error(subname, message, ierr)
        IMPLICIT NONE
        CHARACTER(64)     :: subname, message
        INTEGER           :: ierr, stdout=6
      
        WRITE(stdout, *) "Error while ", message, " in ", subname 
        WRITE(stdout, *) "Error code : ", ierr
END SUBROUTINE
        
