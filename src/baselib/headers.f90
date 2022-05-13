!
! Copyright (C) 2007 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License\'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!**********************************************************
   SUBROUTINE headers(iun, msg, htype)
   !**********************************************************
   !
   ! Print out the given header message msg
   !
   USE io_global_module, ONLY : ionode
   IMPLICIT NONE

   !
   ! input variables
   !
   INTEGER,      INTENT(IN) :: iun
   CHARACTER(*), INTENT(IN) :: msg
   CHARACTER(*), INTENT(IN) :: htype

   !
   ! local variables
   !
   INTEGER :: msglen
   CHARACTER(256) :: str
   CHARACTER(64)  :: subname='headers'

!
!------------------------------
! main body
!------------------------------
!

   IF ( ionode ) THEN
      !
      msglen = LEN_TRIM( msg )
      !
      SELECT CASE ( TRIM(htype) )
      CASE ( "main", "major" )
          WRITE( str, *) "(2x,'=  ',a,", 70-4-msglen, "x, '=')"
          WRITE( iun, "(/,2x,70('='))" )
          IF ( msglen > 65 ) CALL errore(subname,"msg too long",msglen)
          WRITE( iun, FMT=TRIM(str) ) TRIM(msg)
          WRITE( iun, "(2x,70('='),/)" )
      CASE ( "medium")
          WRITE( str, *) "(2x,'   ',a,", 60-4-msglen, "x, ' ')"
          IF ( msglen > 55 ) CALL errore(subname,"msg too long",msglen)
          WRITE( iun, FMT=TRIM(str) ) TRIM(msg)
          WRITE( iun, "(2x,50('='),/)" )
      CASE ( "minor")
          WRITE( str, *) "(2x,'   ',a,", 40-4-msglen, "x, ' ')"
          WRITE( iun, "(2x,40('-'))" )
          IF ( msglen > 35 ) CALL errore(subname,"msg too long",msglen)
          WRITE( iun, FMT=TRIM(str) ) TRIM(msg)
          WRITE( iun, "(2x,40('-'),/)" )
      CASE DEFAULT
          CALL errore(subname,"invalid htype",10)
      END SELECT
      !
   ENDIF
 
END SUBROUTINE headers

