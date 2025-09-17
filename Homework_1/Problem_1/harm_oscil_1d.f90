!A program to execute the horner scheme for a polynomial
program harmonic_oscillator_1d

    use numtype
    implicit none

    real(dp) :: opol(0:100), x, y
    integer :: nm, n, i

    nm = 5
    
    do i = 0,100
        x = -4 + i * 8._dp/100
        !First calculate the Hermite polynomials
        call hermite_poly(nm,x, opol)
        !Now multiply the polynomials by the relevant coefficients
        call wavefunction(nm,x,opol)
        write(7, *) x, opol(0:nm)
    end do

    print *, 'computing 1D Harmonic Oscillator functions'

    contains 


        subroutine hermite_poly(nm, x, hpoly) !Hermite polynomials H_n(x)
            use numtype
            implicit none
            integer, intent(in) :: nm
            real(dp), intent(in) :: x
            real(dp), dimension(0:nm) :: hpoly
            integer :: n

            !First few Hermite polynomails
            hpoly(0) = 1
            hpoly(1) = 2*x
            ! Hermite polynomials to arbitrary degree n
            do n = 1, nm-1
                hpoly(n+1) = 2 * x * hpoly(n) - 2 * n * hpoly(n-1)
            end do

        end subroutine hermite_poly

        subroutine wavefunction(nm, x, wfunc) !Wavefunctions for each n
            use numtype
            implicit none
            integer, intent(in) :: nm
            real(dp), intent(in) :: x
            real(dp), dimension(0:nm) :: wfunc
            integer :: n
            real(dp) :: coeff, beta
            beta = 1
            do n = 0, nm
              !Here we calculate the coefficients of the wavefunctions and multiply them
              !by the input hermite polynomials
              coeff = sqrt(beta/(sqrt(pi)*2**n*fact(n)))*exp(-beta**2*x**2/2)
              wfunc(n) = coeff*wfunc(n)
            end do

        end subroutine wavefunction


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
                !Calculate the factorial
                s0 = n*fact(n-1)
            end if    !end the logic statement

    
        end function fact

end program harmonic_oscillator_1d

