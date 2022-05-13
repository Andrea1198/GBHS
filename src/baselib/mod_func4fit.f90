!======================
module func4fit_m
  !======================
  use kinds
  use constants
  implicit none
  private
  !
  interface lor_ndegree 
     module procedure lor_ndegree_single 
     module procedure lor_ndegree_grid 
     module procedure lor_ndegree_elemental
  end interface
  !
  public :: lor_ndegree
  public :: integrate_n_momentum 
contains
  !======================
  elemental function lor_ndegree_elemental(degree,x,amp,pole) result(r)
    integer,    intent(in) :: degree
    real(dp),   intent(in) :: x
    complex(dp),intent(in) :: amp,pole
    !
    real(dp)               :: r
    integer                :: i 
    !
    r=amp*aimag(pole)**(2*degree-1)/((x-real(pole,dp)**(2*degree)+aimag(pole))**(2*degree))
    !
  end function lor_ndegree_elemental

  !======================
  function lor_ndegree_single(degree,x,npoles,amp,pole) result(r)
    integer,    intent(in) :: degree,npoles
    real(dp),   intent(in) :: x
    complex(dp),intent(in) :: amp(:),pole(:)
    !
    real(dp)               :: r
    integer                :: i 
    !
    r = 0._dp
    do i = 1,npoles
      r=r+lor_ndegree_elemental(degree,x,amp(i),pole(i))
    enddo
    !
  end function lor_ndegree_single

  !======================
  function lor_ndegree_grid(nx,degree,x,npoles,amp,pole) result(r)
    integer,    intent(in) :: nx,degree,npoles
    real(dp),   intent(in) :: x(:)
    complex(dp),intent(in) :: amp(:),pole(:)
    !
    real(dp)               :: r(nx)
    integer                :: i 
    !
    r = 0._dp
    do i = 1,nx
      r(i) = lor_ndegree_single(degree,x(i),npoles,amp,pole)
    enddo
    !
  end function lor_ndegree_grid

  function integrate_n_momentum(n,nx,obj,xgrd,mueff,infint) result(res)
  !------------------------------------
    implicit none
    integer,           intent(in) :: n
    integer,           intent(in) :: nx
    real(dp),          intent(in) :: obj(:)
    real(dp),          intent(in) :: xgrd(:)
    real(dp),optional, intent(in) :: mueff
    logical ,optional, intent(in) :: infint
    !
    integer                       :: ix
    real(dp)                      :: res
    real(dp)                      :: mueff_ 
    real(dp)                      :: mq(2) 
    real(dp)                      :: m,q 
    real(dp)                      :: x_1,x_2,y_1,y_2 
    logical                       :: infint_
    real(dp)                      :: wgt_x(nx)
    !
    mueff_ = xgrd(nx)
    if(present(mueff)) mueff_ = mueff
    infint_ = .false.
    if(present(infint)) infint_ = infint
    !
    wgt_x(1)=0.5_dp*abs(xgrd(2)-xgrd(1))
    do ix=2,nx-1
      wgt_x(ix)=0.5_dp*abs(xgrd(ix-1)-xgrd(ix+1))
    enddo
    wgt_x(nx)=0.5_dp*abs(xgrd(nx)-xgrd(nx-1))
    !
    res = 0._dp
    do ix = 1, nx
      if(xgrd(ix)<=mueff_) then
        res = res + xgrd(ix)**n * obj(ix)* wgt_x(ix) 
      else
        exit
      endif
    enddo
    !
    !
    if(nx<2) return
    if(.not. infint_) return
    !
    ! tail modelling -infty
    !
    x_1 = abs(xgrd(1))
    x_2 = abs(xgrd(2))
    y_1 = abs(obj(1))
    y_2 = abs(obj(2))
    mq = xgrid_2point_linear_interp(log(x_1),log(x_2),log(y_1),log(y_2))
    !
    m = mq(1)
    q = mq(2)
    !
    if ( -m > n+1 ) then
      res=res-sign(1._dp,obj(1))*exp(q)*x_1**(n+m+1._dp)/(n+m+1._dp)
    else
      res=huge(res)
      res=res+res
    endif
    !
    if (xgrd(nx) /= mueff_ .or. xgrd(nx)<=0._dp) return
    !
    ! tail modelling +infty
    !
    x_1 = abs(xgrd(nx-1))
    x_2 = abs(xgrd(nx))
    y_1 = abs(obj(nx-1))
    y_2 = abs(obj(nx))
    mq = xgrid_2point_linear_interp(log(x_1),log(x_2),log(y_1),log(y_2))
    !
    m = mq(1)
    q = mq(2)
    !
    if ( -m > n+1 ) then
      res=res-sign(1._dp,obj(nx))*exp(q)*x_2**(n+m+1._dp)/(n+m+1._dp)
    else
      res=huge(res)
      res=res+res
    endif
    !
    return
    !
  end function integrate_n_momentum 

  function xgrid_2point_linear_interp(x_1,x_2,y_1,y_2) result(res)
  !------------------------------------
    implicit none
    real(dp),        intent(in) :: x_1,x_2
    real(dp),        intent(in) :: y_1,y_2
    real(dp)                    :: res(2)
    !
    !
    res(1) = (y_2 - y_1)/(x_2-x_1)
    res(2) = (x_2 * y_1 - y_2 * x_1)/(x_2-x_1)
    !
    return
    !
  end function xgrid_2point_linear_interp 
end module func4fit_m 
