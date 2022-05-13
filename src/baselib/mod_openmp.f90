!
!        Copyright (C) 2000-2018 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM,AF
! 
! This file is distributed under the terms of the GNU 
! General Public License. You can redistribute it and/or 
! modify it under the terms of the GNU General Public 
! License as published by the Free Software Foundation; 
! either version 2, or (at your option) any later version.
!
! This program is distributed in the hope that it will 
! be useful, but WITHOUT ANY WARRANTY; without even the 
! implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License 
! for more details.
!
! You should have received a copy of the GNU General Public 
! License along with this program; if not, write to the Free 
! Software Foundation, Inc., 59 Temple Place - Suite 330,Boston, 
! MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
!
#if defined __OPENMP  && ! defined _OPENMP
#   define _OPENMP
#endif
!
module openmp_m
 !
 implicit none
 !
 ! ... Threads
 !     the default value of 0 means that the actual number of threads is obtained
 !     from the environment var OMP_NUM_THREADS
 !     any other value given from input will overwrite the environment
 !
 integer     :: n_threads       = 0   ! Max number of threads possible
 integer     :: n_threads_now   = 0
 integer     :: n_threads_limit  = 0
 logical     :: master_thread   = .TRUE.
 !
#if defined _OPENMP
 !
!$omp threadprivate(master_thread)
!$omp threadprivate(n_threads_now)
 !
 integer, external :: omp_get_num_threads
 integer, external :: omp_get_thread_num
 integer, external :: omp_get_thread_limit
 integer, external :: omp_get_max_threads
 logical, external :: omp_get_nested
 !
 logical     :: omp_is_off    = .FALSE.
 !
#else
 logical     :: omp_is_off    = .TRUE.
#endif
 !
 contains
   !
   subroutine OPENMP_initialize( )
     !
#if !defined _OPENMP
     !
     n_threads     = 1
     n_threads_limit = 1
     n_threads_now = 1
     !
#endif
     !
#if defined _OPENMP
     !
     if (omp_is_off) then
       !
       n_threads=1
       n_threads_now=1
       n_threads_limit=1
       !
     else
       !
       n_threads=omp_get_max_threads()
       n_threads_now=n_threads
       n_threads_limit=omp_get_thread_limit()
       !
       ! set defaults against OMP_THREAD_LIMIT being not set
       !
       if(n_threads_limit<0.or.n_threads_limit>1024) n_threads_limit=0
       !
     endif
     !
     call omp_set_dynamic(.false.)
     call omp_set_nested(.false.)
#endif
     !
   end subroutine
   !
   subroutine OPENMP_set_threads(n_threads_in,use_nested)
     !
     integer, optional :: n_threads_in
     logical, optional :: use_nested
     !
     logical :: use_nested_=.false.
     !
     if (omp_is_off) then
       if (present(n_threads_in)) then
         n_threads_in=1
       else
         n_threads=1
       endif
       n_threads_now=1
#if defined _OPENMP
       call omp_set_num_threads( n_threads_now )
#endif
       !
     else
#if defined _OPENMP
       if (present(n_threads_in)) then
         if (n_threads_in==0) n_threads_in=n_threads
         call omp_set_num_threads( n_threads_in )
         n_threads_now=n_threads_in
       else
         call omp_set_num_threads( n_threads )
         n_threads_now=n_threads
       endif
       !
#  if defined _NESTING
       if (present(use_nested)) use_nested_=use_nested
       !
       if (use_nested_) then
          call omp_set_dynamic(.true.)
          call omp_set_nested(.true.)
          call omp_set_max_active_levels(2)
          !     
       endif
#  endif 
       !
#else
       if (present(n_threads_in)) then
          n_threads_in=1
       else
          n_threads=1
       endif
       n_threads_now=1
       !
#endif
     endif
     !
   end subroutine
   !
   subroutine OPENMP_update(master_)
     !
     logical, optional :: master_
     logical :: master__
     !
     if (omp_is_off) then
       master__=.TRUE.
       n_threads_now=1
     else
#if defined _OPENMP
       master__= ( omp_get_thread_num() == 0)
       n_threads_now= omp_get_num_threads()
#else
       master__=.TRUE.
       n_threads_now=1       
#endif
     endif
     if (present(master_)) master_=master__
     !
   end subroutine
   !
   function OPENMP_get_thread_num() result (num)
     !
     integer :: num
     !
#if defined _OPENMP
     num = omp_get_thread_num()
#else
     num = 0
#endif
     !
   end function OPENMP_get_thread_num
   !
end module openmp_m

