! Here we are going to apply the Runge Kutta Method (order 4) to solve 
! Planetary motion 
! which is dy/dx = f(x,y)
!
!2025-09-25

!create a module for our physics quantities
module setup_planet

    use numtype
    implicit none

    integer, parameter :: n_eq = 6
    !Here we set our physics parameters
    real(dp), parameter :: gravity = 6.6743e-11_dp, &
      mass_sun = 1.981e+30_dp, mass_earth = 6.1e+24_dp

end module setup_planet

!Here is where the program starts
program planet

    use numtype
    use setup_planet
    implicit none

    real(dp) :: t,dt,tmax
    real(dp), dimension(n_eq) :: y

    ! These our the initial parameters
    t = 0
    tmax = 60 * 60 * 24 * 365*10 !10 years of orbit
    dt = 60*60*24*7 ! 30 days
    
    y(1:3) = (/ 1.496e+11_dp, 0._dp, 0._dp /) ! x, y, z
    y(4:6) = (/ 0._dp, 29.783e+3_dp, 0._dp /) ! vx, vy, vz

    ! This is our main program loop
    do while(t<=tmax)

        write(1,*) y(1), y(2)
        write(2, *) t, dt
        call rk45step(t,dt,y)
    end do

end program planet


subroutine rk45step(x,h,y)

    use numtype
    use setup_planet
    implicit none
    real(dp), intent(inout) :: x !The x variable can change since we want
    !x_n in, x_n+1 out
    real(dp), intent(inout) :: h
    real(dp), intent(inout), dimension(n_eq) :: y !define y as a vector to have I and Q
    real(dp), dimension(n_eq) :: k1, k2, k3, k4, k5, k6, y1, y2 !k and dy also need to be vectors
    real(dp), parameter :: tiny = 1.e-20_dp, epsilon = 1.e-5_dp
    real(dp) :: rr, delta
    ! interface
    interface !matches up variables with the kv function
              !look up how it works later
        function kv(x, h, y) result(k)
            use numtype
            use setup_planet
            implicit none
            real(dp), intent(in) :: x 
            real(dp), intent(in) :: h
            real(dp), intent(in), dimension(n_eq) :: y
            real(dp), dimension(n_eq) :: k
        end function kv

    end interface

    !k coefficients
    k1 = kv(x,         h,  y)
    k2 = kv(x+h/4,     h,  y + k1/4)
    k3 = kv(x + 3*h/8,   h,  y + 3*k1/32 + 9*k2/32)
    k4 = kv(x + 12*h/13, h,  y + 1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197)
    k5 = kv(x + h,     h,  y + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104)
    k6 = kv(x + h/2,     h,  y - 8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40)

    ! update y and x
    y1 = y + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5
    y2 = y + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55
    
    rr = sqrt(dot_product(y2-y1,y2-y1))/h + tiny ! truncation error (see Wikipedia)

    if (rr < epsilon) then
        x = x + h
        y = y2
        delta = 0.92_dp*(epsilon/rr)**(0.2_dp) ! optimal exponent value
        h = delta * h
    else
        delta = 0.92_dp*(epsilon/rr)**(0.25_dp) ! optimal exponent value
        h = delta * h
    end if   


end subroutine rk45step

function kv(x, h, y) result(k)

    use numtype
    use setup_planet
    implicit none
    real(dp), intent(in) :: x 
    real(dp), intent(in) :: h
    real(dp), intent(in), dimension(n_eq) :: y
    real(dp), dimension(n_eq) :: f, k

   !physics
    real(dp) :: r
   ! r is the distance between the Earth and the Sun 
    r = sqrt(y(1)**2+y(2)**2+y(3)**2)
   
    f(1:3) = y(4:6) 
    f(4:6) = -gravity * mass_sun/r**3 * y(1:3) !y(1:3) is x,y,z
   !end physics
    k = h * f

 
end function kv

