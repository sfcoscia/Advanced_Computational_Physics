program thiele
  use numtype
  use thiele_approx
  implicit none
  real(dp), external :: sinf
  real(dp) :: x, dx, y, dxx, xmax, xmin
  integer :: i, n


  n = 20
  dx = 0.5_dp

  do i = n, 1, -1
    zn(i) = i * dx !zn is already defined in thiele_approx
    fn(i) = sinf(zn(i))
    write(1,*) zn(i), fn(i)
  end do

  call thiele_coef(n, zn, fn, an)
  do i = n/2, -n/2, -1
    x = dx/2 + i * dx + tiny
    y = thiele_cf(x,n,zn,an)
    write(3,*) x, y
  end do

  xmax = n*dx 
  xmin = -n/2*dx
  dxx = (xmax-xmin)/200

  do i = 0, 200
    x = xmin + i*dxx
    write(4,*) x, sinf(x)
  end do


end program thiele

function sinf(x) result(f)
  use numtype
  implicit none
  real(dp) :: x, f

  f = sin(x)/x

end function sinf


