! Here we are going to apply the Runge Kutta Method (order 4) to solve 
! the RC Circuit diff equation
! dQ/dt = (V-Q/C)/I = f(t,Q)
! which is dy/dx = f(x,y)
!
! 09/18/2025

!create a module for our physics quantities
module setup_rc

    use numtype
    implicit none
    !Here we set our physics parameters
    real(dp), parameter :: emf = 10._dp, capacity = 1.5_dp, resistance = 2._dp

end module setup_rc

!Here is where the program starts
program rc

    use numtype
    use setup_rc
    implicit none

    real(dp) :: t,dt,tmax
    real(dp) :: charge

    ! These our the initial parameters
    t = 0
    tmax = 25._dp
    dt = 0.1_dp
    charge = 0._dp

    ! This is our main program loop
    do while(t<tmax)
        write(1,*) t, charge
        call rk4step(t,dt,charge)

    end do

end program rc


subroutine rk4step(x,h,y)

    use numtype
    implicit none
    real(dp), intent(inout) :: x !The x variable can change since we want
    !x_n in, x_n+1 out
    real(dp), intent(in) :: h
    real(dp), intent(inout) :: y
    real(dp) :: k1,k2,k3,k4, dy
    ! interface
    interface !matches up variables with the kv function
              !look up how it works later
        function kv(x, h, y) result(k)
            use numtype
            implicit none
            real(dp), intent(in) :: x 
            real(dp), intent(in) :: h
            real(dp), intent(in) :: y
            real(dp) :: k
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
    real(dp), intent(in) :: y
    real(dp) :: k1,k2,k3,k4,dy
    real(dp) :: f, k

    !physics
    !y = Q (charge)
    f = (emf - y/capacity)/resistance
    !end physics
    k = h*f

 
end function

