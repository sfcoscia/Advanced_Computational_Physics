! Here we are going to apply the Runge Kutta Method (order 4) to solve 
! projectile motion with drag diff equation
! which is dy/dx = f(x,y)
!
!2025-09-23

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
    tmax = 60 * 60 * 24 * 365 * 10 !10 years of orbit
    dt = 60*60*24 ! 1 day
    
    y(1:3) = (/ 1.496e+11_dp, 0._dp, 0._dp /) ! x, y, z
    y(4:6) = (/ 0._dp, 29.783e+3_dp._dp, 0._dp /) ! vx, vy, vz


    ! Initial kinetic energy:
    ke_i = mass/2 * v0**2

    ! This is our main program loop
    do while(y(3)>0)

        write(1,*) y(1), y(3)
        write(2, *) t, y(5)
        call rk4step(t,dt,y)

    end do

    ke_f = (y(2)**2 + y(4)**2 )/(2*mass) ! KE = p^2/2m
    print *, ke_f - ke_i, y(5)

end program planet


subroutine rk4step(x,h,y)

    use numtype
    use setup_planet
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
            use setup_planet
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
    use setup_planet
    implicit none
    real(dp), intent(in) :: x 
    real(dp), intent(in) :: h
    real(dp), intent(in), dimension(n_eq) :: y
    real(dp), dimension(n_eq) :: f, k

    !physics
    
    real(dp) :: vv(2), vrel(2), v !vrel is for the wind
    
    vv(1) = y(2)/mass; vv(2) = y(4)/mass ! v = p/m
    vrel = vv-wind                       ! relative velocity
    v = sqrt(vrel(1)**2 + vrel(2)**2)

    f(1) = vv(1)                          ! v_x
    f(2) = -drag * vrel(1) * v              ! F_x
    f(3) = vv(2)                          ! v_y
    f(4) = -mass*gravity - drag*vrel(2)*v   ! F_y
    f(5) = f(2)*vv(1) + f(4)*vv(2)           ! Power
    !end physics
    k = h * f

 
end function

