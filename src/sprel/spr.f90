PROGRAM special_rel
        USE kind
        USE constants
        USE system_init         ONLY : pos, vel, acc
        USE compute_pr_time
        
        IMPLICIT NONE
        !
        INTEGER, PARAMETER :: stdout = 6
        
        WRITE(stdout, *) pos, vel, acc
END PROGRAM special_rel
