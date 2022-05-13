!======================================================
! SUBROUTINE FOR PROPER TIME COMPUTATION
!======================================================
! ANDREA (MR) P
!======================================================



SUBROUTINE compute_pr_time(t, tau)
        use kinds,              only : dbl
        use sprel_data,        only : gamm
        use constants
        !
        implicit none
        !
        real(dbl),   intent(in) :: t
        real(dbl),  intent(out) :: tau
        !
        tau = gamm * t
end subroutine compute_pr_time

