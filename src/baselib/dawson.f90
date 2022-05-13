
!===================
function gaussian_pole(x,x0,delta)
  !===================
  use kinds,     only : DP
  use constants, only : PI, CI
  implicit none
  !
  complex(DP) :: gaussian_pole
  real(DP)    :: x,x0,delta
  real(DP), external :: dawson
  !
  real(DP)    :: arg
  !  
  arg=(x-x0)/delta
  gaussian_pole =  2.0_DP/delta*dawson(arg) &
                  +CI*sqrt(PI)/delta * exp(-arg**2)
  return
end function gaussian_pole
!
!===================
function dawson(x)
  !===================
  use kinds, only : DP
  implicit none
  integer  :: NMAX
  real(DP) :: dawson,x,H,A1,A2,A3
  parameter (NMAX=6,H=0.4,A1=2./3.,A2=0.4,A3=2./7.)
  !
  ! Returns Dawson’s integral F(x) = exp(−x2) int_0^x exp(t2)dt for any real x.
  !
  integer  :: i,init,n0
  real(DP) :: d1,d2,e1,e2,sum,x2,xp,xx,c(NMAX)
  save     :: init,c
  data init/0/   !Flag is 0 if we need to initialize, else 1. 

  if(init == 0)then
    init=1
    do i=1,NMAX
      c(i)=exp(-((2.*float(i)-1.)*H)**2) 
    enddo
  endif
  if(abs(x) < 0.2)then    !Use series expansion.
    x2=x**2
    dawson=x*(1.-A1*x2*(1.-A2*x2*(1.-A3*x2)))
  else   ! Use sampling theorem representation.
    xx=abs(x) 
    n0=2*nint(0.5*xx/H) 
    xp=xx-float(n0)*H 
    e1=exp(2.*xp*H) 
    e2=e1**2 
    d1=float(n0+1) 
    d2=d1-2.
    sum=0.
    do i=1,NMAX
      sum=sum+c(i)*(e1/d1+1./(d2*e1)) 
      d1=d1+2.
      d2=d2-2.
      e1=e2*e1
    enddo
    dawson=0.5641895835*sign(exp(-xp**2),x)*sum 
  endif
  return 
end function dawson

