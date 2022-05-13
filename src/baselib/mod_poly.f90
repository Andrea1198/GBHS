
!===============
module poly_m
  !===============
  use constants
  use util_module
  implicit none
  private

  type poly_t
    integer :: order
    complex(dp), allocatable :: c(:)
    logical :: alloc=.false.
  end type

  interface assignment(=)
    module procedure poly_assign
  end interface
  !
  interface poly_eval
    module procedure poly_eval_c0d
    module procedure poly_eval_c1d
    module procedure poly_eval_r0d
    module procedure poly_eval_r1d
  end interface

  public :: poly_t
  public :: assignment(=)
  public :: poly_reset
  public :: poly_assign
  public :: poly_sum
  public :: poly_mult
  public :: poly_adjust
  public :: poly_write
  public :: poly_roots
  public :: poly_eval

contains

  recursive subroutine poly_reset(obj,order)
     implicit none
     type(poly_t) :: obj
     integer, optional :: order
     !
     if (present(order)) then
        if (obj%alloc) call poly_reset(obj)
        allocate(obj%c(0:order))
        obj%order=order
        obj%alloc=.true.
     else
        if (allocated(obj%c)) deallocate(obj%c) 
        obj%order=0
        obj%alloc=.false.
     endif
     !
  end subroutine

  subroutine poly_assign(obj1,obj2)
     implicit none
     type(poly_t), intent(inout) :: obj1
     type(poly_t), intent(in)    :: obj2
     !
     if (.not.obj2%alloc) return
     call poly_reset(obj1,obj2%order)
     obj1%c(0:)=obj2%c(0:)
     ! 
  end subroutine poly_assign

  subroutine poly_sum(res,obj1,obj2,alpha1,alpha2)
     ! res = alpha1 * p1 + alpha2 * p2
     implicit none
     type(poly_t) :: res
     type(poly_t) :: obj1, obj2
     complex(dp), optional  :: alpha1, alpha2
     !
     complex(dp)  :: alpha1_, alpha2_

     if (.not.obj1%alloc.or..not.obj2%alloc) return
     call poly_reset(res,max(obj1%order,obj2%order))
     !
     alpha1_=1.0
     alpha2_=1.0
     if (present(alpha1)) alpha1_ = alpha1
     if (present(alpha2)) alpha2_ = alpha2
     !
     res%c=0.0
     res%c(0:obj1%order)=res%c(0:obj1%order)+alpha1_*obj1%c(0:obj1%order)
     res%c(0:obj2%order)=res%c(0:obj2%order)+alpha2_*obj2%c(0:obj2%order)
     return
     !
  end subroutine poly_sum

  subroutine poly_mult(res,obj1,obj2)
     ! res = obj1 * obj2
     implicit none
     type(poly_t) :: res
     type(poly_t) :: obj1, obj2
     integer :: i,j

     if (.not.obj1%alloc.or..not.obj2%alloc) return
     call poly_reset(res,obj1%order+obj2%order)
     !
     res%c=0.0
     do j = 0, obj2%order
     do i = 0, obj1%order
       res%c(i+j) = res%c(i+j) + obj1%c(i)*obj2%c(j)
     enddo
     enddo
     !
  end subroutine poly_mult

  subroutine poly_roots(vroots, obj)
     implicit none
     type(poly_t) :: obj
     complex(dp)  :: vroots(obj%order)
     !
     integer :: i,ndim
     complex(dp), allocatable :: zmat(:,:), cmat(:,:) 

     if (.not.obj%alloc) return
     !
     vroots=0.0
     ndim=obj%order
     if (abs(obj%c(ndim)) < 1.0d-6) return
     !
     ! build the companion matrix
     !
     allocate(cmat(ndim,ndim))
     allocate(zmat(ndim,ndim))
     !
     cmat=0.0
     do i = 1, ndim-1
        cmat(i+1,i) = 1.0
     enddo
     cmat(1:ndim,ndim) = -obj%c(0:ndim-1)/obj%c(ndim)
     !
     call zmat_diag(zmat,vroots,cmat,ndim,"r")
     deallocate(zmat,cmat)
     !
  end subroutine poly_roots

  subroutine poly_adjust(obj)
     implicit none
     type(poly_t) :: obj
     !
     type(poly_t) :: p_tmp
     integer :: i,ndim

     if (.not.obj%alloc) return
     ndim = obj%order
     do i = obj%order, 0, -1
       if (abs(obj%c(i)) > 0.0_dp) then
         ndim=i
         exit
       endif
     enddo
     !
     p_tmp=obj
     call poly_reset(obj)
     call poly_reset(obj,ndim)
     obj%c(0:ndim)=p_tmp%c(0:ndim)
     !
     call poly_reset(p_tmp)
     !
  end subroutine poly_adjust
  !
  function poly_eval_c0d(x,obj) result(res)
     implicit none
     complex(dp)  :: x
     type(poly_t) :: obj
     complex(dp)  :: res
     !
     integer :: i

     res=0.0
     if (.not.obj%alloc) return
     do i = obj%order, 0, -1
       res = res + obj%c(i)*x**i
     enddo
     !
  end function poly_eval_c0d
  !
  function poly_eval_r0d(x,obj) result(res)
     implicit none
     real(dp)     :: x
     type(poly_t) :: obj
     complex(dp)  :: res
     !
     integer :: i

     res=0.0
     if (.not.obj%alloc) return
     do i = obj%order, 0, -1
       res = res + obj%c(i)*x**i
     enddo
     !
  end function poly_eval_r0d
  !
  function poly_eval_c1d(x,n,obj) result(res)
     implicit none
     integer      :: n
     complex(dp)  :: x(n)
     type(poly_t) :: obj
     complex(dp)  :: res(n)
     !
     integer :: i

     res=0.0
     if (.not.obj%alloc) return
     do i = obj%order, 0, -1
       res(1:n) = res(1:n) + obj%c(i)*x(1:n)**i
     enddo
     !
  end function poly_eval_c1d
  !
  function poly_eval_r1d(x,n,obj) result(res)
     implicit none
     integer      :: n
     real(dp)     :: x(n)
     type(poly_t) :: obj
     complex(dp)  :: res(n)
     !
     integer :: i

     res=0.0
     if (.not.obj%alloc) return
     do i = obj%order, 0, -1
       res(1:n) = res(1:n) + obj%c(i)*x(1:n)**i
     enddo
     !
  end function poly_eval_r1d
  !
  subroutine poly_write(iun,obj,label)
     implicit none
     integer :: iun
     type(poly_t) :: obj
     character(*), optional :: label
     !
     character(128):: str
     integer :: i
     
     str="P(x) = "
     if (present(label)) str=trim(label)//" = "
     !
     do i = obj%order, 0, -1
        write(str,"(a,' (',2f15.9,') x^',i1)") trim(str), obj%c(i),i
        write(iun,"(2x,a)") trim(str)
        str= "     + "
     enddo
     !
  end subroutine poly_write

end module poly_m
