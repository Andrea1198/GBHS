!===================================
! MOCULE FOR VARIABLES OF THE MODEL
!===================================
! ANDREA (MR) P
!===================================
MODULE sprel_data
        USE kinds
        !        
        INTEGER                         :: ndim
        REAL(dbl)                       :: gamm
        REAL(dbl),      ALLOCATABLE     :: pos(:), vel(:), acc(:)
        REAL(dbl)                       :: norm_vel
        !
        PUBLIC  :: gamm
        PUBLIC  :: pos, vel, acc
        PUBLIC  :: ndim
        PUBLIC  :: norm_vel
        !
        PUBLIC  :: allocate_datas
        PUBLIC  :: deallocate_datas
        PUBLIC  :: init_gamma

CONTAINS

        SUBROUTINE allocate_datas(ndim_, pos_, vel_, acc_)
                USE kinds
                USE lin_alg
                IMPLICIT NONE
                !
                INTEGER,                INTENT(IN) :: ndim_
                REAL(dbl), OPTIONAL,    INTENT(IN) :: pos_(ndim_), vel_(ndim_), acc_(ndim_)
                
                !
                INTEGER         :: ierr
                REAL(dbl)       :: norm_vel
                CHARACTER(14)   :: subname='allocate_datas'
                !
                !
                ndim = ndim_
                !
                ALLOCATE(pos(ndim), vel(ndim), acc(ndim), stat=ierr)
                IF (ierr/=0) CALL error(subname, 'allocating pos, vel, acc', ABS(ierr))
                !
                IF (PRESENT(pos_)) pos = pos_
                IF (PRESENT(vel_)) vel = vel_
                IF (PRESENT(acc_)) acc = acc_
                !
                CALL norma2(vel, size(vel), norm_vel)
                CALL compute_gamma(gamm, norm_vel)
                !
        END SUBROUTINE allocate_datas

        SUBROUTINE deallocate_datas()
                IMPLICIT NONE
                !
                gamm = 0.0d0

        END SUBROUTINE deallocate_datas

        SUBROUTINE compute_gamma(v, gamm)
                USE kinds
                USE constants
                !
                IMPLICIT NONE
                !
                REAL(dbl),     INTENT(OUT) :: gamm
                REAL(dbl),      INTENT(IN) :: v
                !
                gamm = sqrt(1.0d0 - v / C_SI)
        END SUBROUTINE compute_gamma
END MODULE sprel_data
