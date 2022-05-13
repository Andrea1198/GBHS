
!=================
program test_poly
  !=================
  use constants
  use poly_m
  implicit none

  type(poly_t) :: p1
  type(poly_t) :: p2
  type(poly_t) :: p3
  complex(dp), allocatable :: vzero(:)

  call poly_reset(p1,1)
  call poly_reset(p2,2)

  p1%c(0) = -2.0
  p1%c(1) = 1.0
  p2%c(0) = -1.0
  p2%c(1) = 0.0
  p2%c(2) = 1.0

  write(6,"(2x,a)") "p1 = x -2"
  call poly_write(6,p1)
  write(6,"(2x,a)") "p1 = x^2 -1"
  call poly_write(6,p2)

  write(6,"(/,2x,a)") "p3 = p1"
  p3=p1
  call poly_write(6,p3)

  call poly_sum(p3,  p1,p2)
  write(6,"(/, 2x,a)") "p3 = p1 + p2"
  call poly_write(6,p3)

  call poly_mult(p3, p1,p2)
  write(6,"(/, 2x,a)") "p3 = p1*p2"
  call poly_write(6,p3)

  allocate(vzero(p3%order))
  call poly_roots(vzero,p3)
  write(6,"(/, 2x,a)") "Roots p3"
  write(6,"( 10f15.9)") vzero
  deallocate(vzero)

  call poly_reset(p3)
  call poly_reset(p3,3)
  allocate(vzero(3))

  p3%c(0) = 0.0
  p3%c(1) = 0.0
  p3%c(2) = 0.0
  p3%c(3) = 1.0
  write(6,"(/, 2x,a)") "p3 = X^3"
  call poly_write(6,p3)
  call poly_roots(vzero,p3)
  write(6,"(/, 2x,a)") "Roots p3"
  write(6,"( 10f15.9)") vzero
  deallocate(vzero)

end program test_poly
