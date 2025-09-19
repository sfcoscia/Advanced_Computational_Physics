! Here we are going to apply the Runge Kutta Method (order 4) to solve 
! the RC Circuit diff equation
! dQ/dt = (V-Q/C)/I = f(t,Q)
! which is dy/dx = f(x,y)
!
! and for I where dI/dt = - I/(RC)
! I(0) = emf/R
! 09/18/2025

!create a module for our physics quantities
module setup_rc

    use numtype
    implicit none

    integer, parameter :: n_eq = 2
    !Here we set our physics parameters
    real(dp), parameter :: emf = 10._dp, capacity = 1.5_dp, resistance = 2._dp

end module setup_rc

!Here is where the program starts
program rc

    use numtype
    use setup_rc
    implicit none

    real(dp) :: t,dt,tmax
    real(dp), dimension(n_eq) :: y

    ! These our the initial parameters
    t = 0
    tmax = 25._dp
    dt = 0.1_dp
    y(1) = 0._dp ! charge
    y(2) = emf/resistance ! current

    ! This is our main program loop
    do while(t<tmax)
        write(1,*) t, y(1)
        write(2, *) t, y(2)
        call rk4step(t,dt,y)

    end do

end program rc


subroutine rk4step(x,h,y)

    use numtype
    use setup_rc
    implicit none
    real(dp), intent(inout) :: x !The x variable can change since we want
    !x_n in, x_n+1 out
    real(dp), intent(in) :: h
    real(dp), intent(inout), dimension(n_eq) :: y !define y as a vector to have I and Q
    real(dp), dimension(n_eq) :: k1,k2,k3,k4, dy !k and dy also need to be vectors
    ! interface
    interface !matches up variables with the kv function
              !look up how it works later
        function kv(x, h, y) result(k)
            use numtype
            use setup_rc
            implicit none
            real(dp), intent(in) :: x 
            real(dp), intent(in) :: h
            real(dp), intent(in), dimension(n_eq) :: y
            real(dp), dimension(n_eq) :: k
        end function kv

    end interface

    !k coefficients
    k1 = kv(x, h, y)
    k2 = kv(x+h/2,h,y+k1/2)
    k3 = kv(x+h/2,h,y+k2/2)
    k4 = kv(x+h,h,y+k3)
    ! average k is dy
    dy = (k1+ 2*k2 + 2*k3 + k4)/6
    ! update y and x
    y = y + dy
    x = x + h

end subroutine rk4step

function kv(x, h, y) result(k)

    use numtype
    use setup_rc
    implicit none
    real(dp), intent(in) :: x 
    real(dp), intent(in) :: h
    real(dp), intent(in), dimension(n_eq) :: y
    real(dp), dimension(n_eq) :: f, k

    !physics
    !y = Q (charge)
    f(1) = (emf - y(1)/capacity)/resistance
    f(2) = -y(2)/(capacity * resistance)
    !end physics
    k = h*f

 
end function

