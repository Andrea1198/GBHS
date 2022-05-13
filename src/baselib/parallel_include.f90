!
! Copyright (C) 2004 PWSCF-CP-FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
!   SISSA Code Interface -- Carlo Cavazzoni
!------------------------------------------------------------------------------C
      MODULE parallel_include

         USE kinds
         LOGICAL tparallel

#if defined __MPI
!
!     Include file for MPI
!
         INCLUDE 'mpif.h'
         DATA tparallel /.true./
#else
         integer, parameter :: MPI_COMM_WORLD=0
         DATA tparallel /.false./
#endif

      END MODULE parallel_include
