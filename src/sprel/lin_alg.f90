MODULE lin_alg

        USE kinds
        USE constants

        IMPLICIT NONE
        PRIVATE

        INTERFACE mat_prod
                MODULE PROCEDURE square_mat_prod
                MODULE PROCEDURE general_mat_prod
        END INTERFACE

        PUBLIC  :: mat_prod
        PUBLIC  :: norma2
        PUBLIC  :: trans

CONTAINS
SUBROUTINE norma2(vec, ndim, norm)
        USE kinds
        !
        IMPLICIT NONE
        INTEGER,        INTENT(IN) :: ndim
        REAL(dbl),      INTENT(IN) :: vec(ndim)
        REAL(dbl),     INTENT(OUT) :: norm
        !
        INTEGER         :: i
        CHARACTER(5)    :: subname='norm2'        
        !
        IF (SIZE(vec) /= ndim) CALL error(subname, 'size and space dimension are different', 1)
        !
        norm = 0.0d0
        DO i = 1, ndim
                norm = norm + vec(i) * vec(i)
        ENDDO
        !
END SUBROUTINE norma2

SUBROUTINE trans(matrix, ndim1, ndim2, t_matrix)
        USE kinds
        !
        IMPLICIT NONE
        !
        INTEGER,        INTENT(IN) :: ndim1, ndim2
        REAL(dbl),      INTENT(IN) :: matrix(ndim1, ndim2)
        REAL(dbl),     INTENT(OUT) :: t_matrix(ndim2, ndim1)
        !
        INTEGER :: i,j
        !
        DO i = 1, ndim1
        DO j = 1, ndim2
                t_matrix(j,i) = matrix(i,j)
        ENDDO
        ENDDO
        !
END SUBROUTINE trans

SUBROUTINE square_mat_prod(mat1, mat2, ndim, mat_f)
        USE kinds
        !
        IMPLICIT NONE
        !
        INTEGER,                INTENT(IN) :: ndim
        REAL(dbl),              INTENT(IN) :: mat1(ndim, ndim), mat2(ndim, ndim)
        REAL(dbl),             INTENT(OUT) :: mat_f(ndim, ndim)
        !
        INTEGER :: i, j, k
        !
        DO i = 1, ndim
        DO j = 1, ndim
                mat_f(i,j) = 0.0d0
                DO k = 1, ndim
                        mat_f(i,j) = mat_f(i,j) + mat1(i,k) * mat2(k,j)
                ENDDO
        ENDDO
        ENDDO
END SUBROUTINE square_mat_prod

SUBROUTINE general_mat_prod(mat1, mat2, ndim1, ndim2, ndim3, ndim4, mat_f)
        USE kinds
        !
        IMPLICIT NONE
        !
        INTEGER,                INTENT(IN) :: ndim1, ndim2, ndim3, ndim4
        REAL(dbl),              INTENT(IN) :: mat1(ndim1, ndim2), mat2(ndim3, ndim4)
        REAL(dbl),             INTENT(OUT) :: mat_f(ndim1, ndim4)
        !
        CHARACTER(8) :: subname='mat_prod'
        INTEGER :: i, j, k
        !
        ! CHECK DIMENSIONS
        !
        IF (ndim2/=ndim3) CALL error(subname, "matrix dimensions don't match")
        !
        DO i = 1, ndim1
        DO j = 1, ndim4
                mat_f(i,j) = 0.0d0
                DO k = 1, ndim3
                        mat_f(i,j) = mat_f(i,j) + mat1(i,k) * mat2(k,j)
                ENDDO
        ENDDO
        ENDDO
        !
END SUBROUTINE general_mat_prod

END MODULE
