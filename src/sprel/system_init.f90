!=============================================
! HERE ALL SYSTEM CHARACTERISTICS ARE STORED
!=============================================
! ANDREA (MR) P
!=============================================

USE kinds
USE constants

MODULE system
        implicit none
        REAL(dbl), ALLOCATABLE  :: pos(:)
        REAL(dbl)               :: time, pr_time
        REAL(dbl), ALLOCATABLE  :: vel(:)
        REAL(dbl), ALLOCATABLE  :: acc(:)
        LOGICAL                 :: alloc
        !
        PUBLIC  :: pos
        PUBLIC  :: time
        PUBLIC  :: pr_time
CONTAINS
        SUBROUTINE system_allocate(pos_, vel_, acc_):
                !==========================
                ! init the system
                !==========================
                INTEGER :: ierr
                !
                ALLOCATE(pos(ndim), vel(ndim), acc(ndim), STAT=ierr)
                IF (ierr/=0) CALL error( rout_name, 'allocating pos, vel, acc', ABS(ierr) )
                !
                !
        END SUBROUTINE system_allocate
        !
        SUBROUTINE system_deallocate():
                IF (alloc) DEALLOCATE(pos, vel, acc)
                alloc = .false.
        END SUBROUTINE system_deallocate
        !
        SUBROUTINE system_init(pos_, vel_, acc_)
                !
                REAL(dbl), OPTIONAL, INTENT(IN) :: pos_, vel_, acc_
                !
                CALL system_allocate() 
                !
                IF (PRESENT(pos_)) pos = pos_
                IF (PRESENT(vel_)) vel = vel_
                IF (PRESENT(pos_)) acc = acc_
                !
        END SUBROUTINE system_init
