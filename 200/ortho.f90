!A program to execute the horner scheme for a polynomial
program ortho_test

    use numtype
    implicit none

    real(dp) :: opol(0:100), x, y
    integer :: nm, n, i

    nm = 4
    
    do i = 0,100
        x = -2 + i * 4._dp/100
        call hermite_poly(nm,x, opol)
        write(7, *) x, opol(0:nm)

    end do

     do i = 0,100
        x = -1 + i * 2._dp/100
        call legendre_poly(nm,x, opol)
        write(8, *) x, opol(0:nm)

    end do
    
     do i = 0,100
        x = -1 + i * 2._dp/100
        write(18, *) x, Legendre_P(0,x), Legendre_P(1,x), Legendre_P(2,x), Legendre_P(3,x), Legendre_P(4,x)

    end do

     do i = 0,100
        x = 0+i*pi/100

        write(19, *) x, 15._dp/2*(7*(cos(x))**2-1)*(sin(x)**2), Legendre_A(4,2,cos(x)) 
    end do

    print *, fact(3), fact(4)

    contains 

        ! Another way to calculate the Legendre Polynomials for P_n
        recursive function Legendre_A(l,m,x) result(s0)   !Legendre_P
            use numtype
            implicit none
            integer, intent(in) :: l, m
            real(dp), intent(in) :: x
            real(dp) :: s0

            !logic statements!
            if( abs(x) > 1 ) then 
                stop ' |x|>1 '
            else if(abs(m) > l) then 
                s0 = 0
            else if (l<0) then
                !flip the sign of l and then put it in the Legendre function
                s0 = Legendre_A(-l-1,m,x)
            else if (m<0) then
                ! here we have to flip the sign of m as well so m --> -m
               s0 = (-1)**(-m) * fact(l-(-m))/fact(l+(-m)) * Legendre_A(l,(-m),x)
            else if (x < 0) then 
                s0 = (-1)**(l-m) * Legendre_A(l,m,-x)
            else if (l == 0 .and. m==0 ) then
                s0 = 1
            else if (l==m) then
                s0 = -(2*l-1)*sqrt(1-x**2)*Legendre_A(l-1,l-1,x)
            else 
                s0 = ((2*l-1) * x * Legendre_A(l-1,m,x)-(l-1+m)*Legendre_A(l-2,m,x))/(l-m)
            end if    !end the logic statement

        end function Legendre_A



        ! Another way to calculate the Legendre Polynomials for P_n
        recursive function Legendre_P(n,x) result(s0)   !Legendre_P
            use numtype
            implicit none
            integer, intent(in) :: n
            real(dp), intent(in) :: x
            real(dp) :: s0

            !logic statements!
            if(n < 0) then 
                s0 = 0
            else if (n==0) then
                s0 = 1._dp
            else 
                s0 = ((2*n-1) * x * Legendre_P(n-1,x) - (n-1) * Legendre_P(n-2,x))/n
            end if    !end the logic statement

    
        end function Legendre_P

        subroutine legendre_poly(nm, x, lpoly) !Legendre polynomials P_n+1(x)
            use numtype
            implicit none
            integer, intent(in) :: nm
            real(dp), intent(in) :: x
            real(dp), dimension(0:nm) :: lpoly
            integer :: n

            lpoly(0) = 1
            lpoly(1) = x
            do n = 1, nm-1
                lpoly(n+1) = ((2*n+1) * x * lpoly(n) - n * lpoly(n-1))/(n+1)
            end do

        end subroutine legendre_poly

        subroutine hermite_poly(nm, x, hpoly) !Hermite polynomials H_n(x)
            use numtype
            implicit none
            integer, intent(in) :: nm
            real(dp), intent(in) :: x
            real(dp), dimension(0:nm) :: hpoly
            integer :: n

            hpoly(0) = 1
            hpoly(1) = 2*x
            do n = 1, nm-1
                hpoly(n+1) = 2 * x * hpoly(n) - 2 * n * hpoly(n-1)
            end do

        end subroutine hermite_poly


        !calculate n factorial recursively
        recursive function fact(n) result(s0)   !n!
            use numtype
            implicit none
            integer, intent(in) :: n
            real(dp) :: s0

            !logic statements!
            if(n < 0) then 
                stop ' something went wrong '
            else if (n==0) then
                s0 = 1._dp
            else 
                s0 = n*fact(n-1)
            end if    !end the logic statement

    
        end function fact

end program ortho_test

