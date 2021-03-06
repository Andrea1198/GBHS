!
! Copyright (C) 2006 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License\'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*********************************************
   MODULE smearing_base_module
!*********************************************
   USE kinds,         ONLY : dbl
   USE constants,     ONLY : ONE, TWO, PI, SQRTPI, SQRT2
   IMPLICIT NONE
   PRIVATE 
!
! Contains smearing definitions
! 
   
   PUBLIC :: smearing_func


CONTAINS

!
! subroutines
!

!**********************************************************
   FUNCTION smearing_func(x, smearing_type )
   !**********************************************************
   !
   ! evaluate the smearing function at x = ( eps - eps0 )/delta
   ! where eps, eps0 are energies and delta is the broadening parameter
   !
   ! the smearing function F is defined in such a way that
   ! 1/delta * F( x ) 
   ! 
   ! is normalized to one when integrated in energy 
   !
   ! allowed values for SMEARING_TYPE are:
   !  "lorentzian"
   !  "gaussian"
   !  "methfessel-paxton"
   !  "fermi-dirac"
   !  "marzari-vanderbilt"
   !
   IMPLICIT NONE
   !
   ! input variables
   !
   REAL(dbl)          :: smearing_func
   REAL(dbl)          :: x
   CHARACTER(*)       :: smearing_type
   !
   ! local variables
   !
   CHARACTER(13)      :: subname="smearing_func"
   REAL(dbl)          :: cost, d1, d2

!
!-------------------
! main body
!-------------------
!
       !
       ! dummy init
       smearing_func = ONE
       !
       SELECT CASE (TRIM(smearing_type))
       !
       CASE ( "lorentzian", "selfenergy0" )   
            !
            cost = ONE / PI
            !
            smearing_func = cost * ONE/( ONE + x**2 )  
    
       CASE ( "selfenergy", "selfenergy1" )   
            !
            d1 = 0.6_dbl
            !
            cost = d1**2 / PI
            !
            smearing_func = cost * ONE/( x**4 + (1.0_dbl -2.0_dbl*d1**2) * x**2 + d1**4 )
    
       CASE ( "selfenergy2" )   
            !
            d1 = 0.6_dbl
            d2 = 0.3_dbl
            !
            cost = d1**2 * d2**2 / PI
            !
            smearing_func = cost * ONE/( x**6 + (1.0_dbl -2.0_dbl*d1**2 -2.0d0*d2**2) * x**4 &
                                               - (2.0d0 * d2**2 -( d1**2+d2**2)**2 ) *x**2 + d2**4 ) 

       CASE ( "gaussian" )
            !
            cost = ONE / SQRTPI
            !
            smearing_func = cost * EXP( -x**2 )
            
       CASE ( "fermi-dirac", "fd" )
            !
            cost = ONE / TWO 
            !
            smearing_func = cost * ONE / ( ONE + COSH(x) )

       CASE ( "methfessel-paxton", "mp" )
            !
            ! Only N=1 is implemented at the moment
            !
            cost = ONE / SQRTPI
            !
            smearing_func = cost * EXP( -x**2 ) * ( 3.0_dbl/2.0_dbl - x**2 )

       CASE ( "marzari-vanderbilt", "mv" )
            !
            cost = ONE / SQRTPI 
            !
            smearing_func = cost * EXP( -(x- ONE/SQRT2 )**2 ) * ( TWO - SQRT2 * x )

       CASE DEFAULT
            CALL errore(subname, 'invalid smearing_type = '//TRIM(smearing_type),1)
       
       END SELECT

   END FUNCTION smearing_func

END MODULE smearing_base_module

