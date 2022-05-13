!
! Copyright (C) 2016 A. Ferretti
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*******************************
module math_common_m
   !*******************************
   !
   implicit none
   private
   !
   public :: binomial
   public :: factorial
   !
contains
!
!********************************************************
   integer function factorial(n)
   !********************************************************
   !
   implicit none
   integer :: n
   !
   integer :: i
   factorial=1
   if (n==0.or.n==1) return
   if (n<0) call errore("factorial","invalid n<0",10)
   !
   do i = 2,n
     factorial=factorial*i
   enddo
   return
   !
end function

!********************************************************
   integer function binomial(n,k)
   !********************************************************
   !
   implicit none
   integer :: n,k
   !
   integer :: i
   binomial=1
   if (n<0) call errore("binomial","invalid n<0",10)
   if (k<0) call errore("binomial","invalid k<0",10)
   if (n<k) call errore("binomial","invalid n<k",10)
   if (n==0.or.n==1) return
   if (k==0) return
   if (n==k) return
   !
   do i = k+1,n
     binomial=binomial*i
   enddo
   binomial=binomial/factorial(n-k)
   return
   !
end function

end module math_common_m
