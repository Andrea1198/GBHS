
!**********************************************************
   PROGRAM test_omp
   !**********************************************************
   !
   ! test the integration and matrix multiplication made by epsilon_drv to
   ! obtain alpha. Expected result is (pi/4)*(A**3)
   !
   USE kinds
   USE random_numbers_module,   ONLY: rndm
   USE util_module,             ONLY: mat_mul
   USE openmp_m,                ONLY: n_threads,n_threads_now,n_threads_limit,&
                                      master_thread, OPENMP_update, &
                                      OPENMP_set_threads

   !
   IMPLICIT NONE
      !
      INTEGER, PARAMETER :: stdout=6
      INTEGER            :: ierr
      !
      INTEGER, PARAMETER :: ndim=500
      REAL(dbl), ALLOCATABLE :: A(:,:),B(:,:),C1(:,:),C2(:,:)
      !
      INTEGER :: i, j, k

      !
      ! startup
      !
      CALL startup("1.0","test_omp")
      CALL headers( stdout, "test OMP", "main" )

      !
      ! preliminary computations
      !
      ALLOCATE(A(ndim,ndim), B(ndim,ndim))
      ALLOCATE(C1(ndim,ndim), C2(ndim,ndim))
      !
      A = 0.0_dbl
      B = 0.0_dbl
      !
      do i=1,ndim
      do j=1,ndim
        A(i,j) = rndm()
        B(i,j) = rndm()
      enddo
      enddo

      !
      ! test initialization values
      !
      CALL headers( stdout, "initialization values", "minor" )
      !
      WRITE( stdout, "(2x,'n_threads      ',i5)") n_threads
      WRITE( stdout, "(2x,'n_threads_now  ',i5)") n_threads_now
      WRITE( stdout, "(2x,'n_threads_limit',i5)") n_threads_limit
      WRITE( stdout, "(/)")

      !
      ! test non-parallel region  
      !  
      CALL headers( stdout, "test OPENMP_set_threads(8)", "minor" )
      !
      call OPENMP_set_threads(8)
      !
      WRITE( stdout, "(2x,'n_threads      ',i5)") n_threads
      WRITE( stdout, "(2x,'n_threads_now  ',i5)") n_threads_now
      WRITE( stdout, "(2x,'n_threads_limit',i5)") n_threads_limit
      WRITE( stdout, "(/)")

      !
      ! test parallel region
      !
      ! CALL headers( stdout, "parallel region", "minor" )
      !
      !$omp parallel default(shared), private(i,j,k)
      !
      CALL OPENMP_update(master_thread)
      !
      !$omp do collapse(2)
      DO i = 1, ndim
      DO J = 1, ndim
        C1(i,j)=0.0
        DO k = 1, ndim
          C1(i,j) = C1(i,j) + A(i,k) * B(k,j)
        ENDDO
      ENDDO
      ENDDO
      !omp end do
      ! if (master_thread) WRITE( stdout, "(2x,'Outside loop:')")
      ! if (master_thread) WRITE( stdout, "(2x,'n_threads      ',i5)") n_threads
      ! if (master_thread) WRITE( stdout, "(2x,'n_threads_now  ',i5)") n_threads_now
      ! if (master_thread) WRITE( stdout, "(2x,'n_threads_limit',i5)") n_threads_limit
      ! if (master_thread) WRITE( stdout, "(/)")
      !$omp end parallel

      !
      ! check result
      !
      CALL mat_mul(C2,A,'N',B,'N',ndim,ndim,ndim)
      !
      DO j = 1, ndim
      DO i = 1, ndim
         if (ABS(C1(i,j)-C2(i,j)) >1.0d-8 ) &
           CALL errore("test_omp", "C1 is wrong", i+j)
      ENDDO
      ENDDO

      DEALLOCATE(A,B,C1,C2)

   END PROGRAM test_omp

