!A program to execute the horner scheme for a polynomial
program harmonic_oscillator_1d

    use numtype
    implicit none

    real(dp) :: opol(0:100), x, y
    integer :: nm, n, i

    nm = 5
    
    do i = 0,100
        x = -2 + i * 4._dp/100
        call wavefunction(nm,x, opol)
        write(7, *) x, opol(0:nm)*sqrt(1/(sqrt(pi)*2**nm*fact(nm)))*exp(-x**2/2)
    end do

    print *, 'computing 1D Harmonic Oscillator functions'

    contains 

        subroutine wavefunction(nm, x, wfunc) !Hermite polynomials H_n(x)
            use numtype
            implicit none
            integer, intent(in) :: nm
            real(dp), intent(in) :: x
            real(dp), dimension(0:nm) :: wfunc
            integer :: n
            real(dp) :: coeff
             
            wfunc(0) = 1*sqrt(1/(2*sqrt(pi)))*exp(-x**2)/2
            wfunc(1) = 2*x*sqrt(1/(sqrt(pi)*2**2*fact(2)))*exp(-x**2/2)
            do n = 1, nm-1
              coeff = sqrt(1/(sqrt(pi)*2**n*fact(n)))*exp(-x**2/2)
              wfunc(n+1) = coeff*(2 * x * wfunc(n) - 2 * n * wfunc(n-1))
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
                s0 = n*fact(n-1)
            end if    !end the logic statement

    
        end function fact

end program harmonic_oscillator_1d

