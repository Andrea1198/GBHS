!
! Copyright (C) 2006 Andrea Ferretti
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE mp_tools
IMPLICIT NONE
PRIVATE

  PUBLIC :: divide_et_impera

CONTAINS

!*********************************************
   SUBROUTINE divide_et_impera(ns_global, ne_global, ns_local, ne_local, mpime, nproc, contig, indx)
   !*********************************************
   !
   ! given the global indexes of a loop, divide the 
   ! loop over the NPROC tasks, and determine the 
   ! extrema for the local process identified by MPIME
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN)     :: ns_global, ne_global
   INTEGER, INTENT(IN)     :: mpime, nproc
   INTEGER, INTENT(OUT)    :: ns_local, ne_local
   LOGICAL, OPTIONAL, INTENT(IN) :: contig
   INTEGER, OPTIONAL, INTENT(OUT):: indx(*)
   !
   ! local variables
   !
   CHARACTER(16)  :: subname='divide_et_impera'
   LOGICAL        :: contig_
   INTEGER, ALLOCATABLE :: dim_local(:)
   INTEGER        :: dim_global, dim_remind, i, ierr
   !
!
!----------------
! main body
!----------------
!
   IF ( nproc < 1 ) CALL errore(subname,'invalid nproc_',1)
   IF ( mpime < 0 ) CALL errore(subname,'invalid mpime',1)
   IF ( mpime >= nproc ) CALL errore(subname,'mpime too large',mpime)
   IF ( ne_global < ns_global ) CALL errore(subname,'invalid global indeces',1)
   !
   dim_global = ne_global - ns_global +1
   !
   contig_=.true.
   if (present(contig)) contig_=contig
   if (.not.contig_.and..not.present(indx)) CALL errore(subname,"non gontig indexes and indx not present",10)
   !
   ALLOCATE( dim_local( 0:nproc-1), STAT=ierr )
   IF (ierr/=0) CALL errore(subname,'allocating dim_local',ABS(ierr))
   !
   IF (contig_) THEN
     !
     dim_local(:)  = dim_global / nproc
     dim_remind = MOD( dim_global, nproc )
     !
     DO i=1, dim_remind
        dim_local(i-1) = dim_local(i-1) + 1
     ENDDO 
     !
     IF ( SUM(dim_local(:)) /= dim_global ) CALL errore(subname,'invalid sum rule',1)
     !
     ns_local  = ns_global + SUM( dim_local(0:mpime -1 ) )
     ne_local  = ns_local  + dim_local( mpime ) -1  
     !
     if (present(indx)) then
       !
       do i = ns_local,ne_local
         indx(i-ns_local+1)=i
       enddo
       !
     endif
     !
   else
     !
     ns_local = 1
     ne_local = 0
     do i = ns_global, ne_global
       if (mod(i-1,nproc)==mpime) then
          ne_local=ne_local+1 
          indx(ne_local)=i
       endif
     enddo
     !
   ENDIF
   !
   DEALLOCATE( dim_local, STAT=ierr )
   IF (ierr/=0) CALL errore(subname,'allocating dim_local',ABS(ierr))

END SUBROUTINE divide_et_impera 

END MODULE mp_tools
