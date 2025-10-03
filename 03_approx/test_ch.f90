
program chebytest

    use chebyshev
    implicit none 
    real(dp), external :: func, dfunc
    real(dp) :: x, dx, ya, yb, eps, dz, zz
    integer :: nn, np, i, maxf 

    nn = 11
    ya = -pi/4
    yb = 9*pi/4

    call chebyex(func,nn,cheb,ya,yb)

    np = 100
    dx = (yb-ya)/np
    do i = 0, np
        x = ya + i*dx 
        write(1,*) x, 0._dp, func(x), cheby(x,cheb,nn,ya,yb)
    end do 

    call chebyzero(nn,cheb,ya,yb,z0,iz0) 
    print *,z0(1:iz0)
    eps = 1.e-15_dp
    dz = 0.01
    maxf = 10
    do i = 1, iz0 
        zz = z0(i)
        call root_polish(func,zz,dz,eps,maxf)
        print *,i, z0(i),func(z0(i)), zz,func(zz)
    end do 

    call chebyderiv(cheb,nn,chder,ya,yb) 
    call chebyderiv(chder,nn-1,chder2,ya,yb)

    np = 100
    dx = (yb-ya)/np
    do i = 0, np
        x = ya + i*dx 
        write(2,*) x, 0._dp, dfunc(x), cheby(x,chder,nn-1,ya,yb)
    end do    

    call chebyzero(nn-1,chder,ya,yb,z0,iz0) 
    do i = 1, iz0
        print *," f'(x)=0 ", i,z0(i), cheby( z0(i),chder2,nn-2,ya,yb)
    end do

end program chebytest


function func(x) result(f)
    use numtype
    implicit none
    real(dp) :: x, f 

    f = x*cos(x)

end function func



function dfunc(x) result(df)
    use numtype
    implicit none
    real(dp) :: x, df 

    df = cos(x) - x*sin(x)

end function dfunc

