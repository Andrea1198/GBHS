!
! Copyright (C) 2018  AGWX Team
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License\'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
! Authors: AF, TC
!
! Some parts of the module below have been taken from QE
!
! Copyright (C) 2004-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!==============================
   module wrappers_m
   !==============================
   use kinds
   use iso_c_binding
   implicit none
   private
   !
   INTERFACE
     FUNCTION mkdir(dirname,mode) BIND(C,name="mkdir") RESULT(r)
       USE iso_c_binding
       CHARACTER(kind=c_char),INTENT(in)  :: dirname(*)
       INTEGER(c_int),VALUE  ,INTENT(in)  :: mode
       INTEGER(c_int)        :: r
     END FUNCTION
   END INTERFACE
   !
   public :: execute_command_line_wrapper
   public :: f_mkdir, f_mkdir_safe, f_rmdir
   !
contains

!============================================
subroutine execute_command_line_wrapper(msg,ierr)
  !============================================
  implicit none
  character(len=*) :: msg
  integer, optional :: ierr
  !
  integer :: ierr_
  !
  ! execute_command_line  not supported before ifort 15
  !
  ierr_=0
#if defined __INTEL_COMPILER
#   if  __INTEL_COMPILER >= 1500 
      call execute_command_line(msg,exitstat=ierr_)
#   else
      ierr_=0
      call system(msg)
#   endif
#else
    call execute_command_line(msg,exitstat=ierr_)
#endif
  if (present(ierr)) ierr=ierr_ 
  !
end subroutine execute_command_line_wrapper

function f_mkdir(dirname) RESULT(res)
  character(*),intent(in)  :: dirname
  integer :: res
  !
  res=0
  call execute_command_line_wrapper("mkdir -p "//trim(dirname)//" 2> /dev/null", res)
  !
end function

function f_rmdir(dirname) RESULT(res)
  character(*),intent(in)  :: dirname
  integer :: res
  !
  res=0
  call execute_command_line_wrapper("if [ -d "//trim(dirname)//" ]; then rm -rf "&
                                     //trim(dirname)// "; fi 2> /dev/null", res)
  !
end function



!  f_mkdir from QE (the one below), causes an ERROR if dirname already exists: use f_mkdir_safe instead
!
!  FUNCTION f_mkdir(dirname, mode) RESULT(r)
!    CHARACTER(*),INTENT(in)  :: dirname
!    INTEGER,INTENT(in) :: mode
!    INTEGER(c_int) :: r
!    INTEGER(c_int) :: c_mode
!    c_mode = INT(mode, kind=c_int)
!    r= mkdir(TRIM(dirname)//C_NULL_CHAR, c_mode)
!  END FUNCTION

  !
  ! safe mkdir from clib/c_mkdir.c that creates a directory, if necessary, 
  ! and checks permissions. It can be called in parallel.
  ! Returns: 0 = all ok
  !          1 = error
  !         -1 = the directory already existed and is properly writable
  FUNCTION f_mkdir_safe(dirname) RESULT(r)
    INTERFACE
    FUNCTION mkdir_safe(dirname) BIND(C,name="c_mkdir_safe") RESULT(r)
      USE iso_c_binding
      CHARACTER(kind=c_char),INTENT(in)  :: dirname(*)
      INTEGER(c_int)         :: r
    END FUNCTION mkdir_safe
    END INTERFACE
    CHARACTER(*),INTENT(in)  :: dirname
    INTEGER(c_int) :: r
    r= mkdir_safe(TRIM(dirname)//C_NULL_CHAR)
  END FUNCTION

end module wrappers_m
