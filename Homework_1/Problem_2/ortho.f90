!A program to execute the horner scheme for a polynomial
program harmonic_oscillator_nd

    use numtype
    implicit none

    real(dp) :: opol(0:100), x, y
    integer :: nm, n, i

    nm = 5
    
      
     do i = 0,100
        x = -2 + i * 1._dp/10
        write(18, *) x, Laguerre_P(3,0._dp,x) 
    end do

     print *, 'computing 1D Harmonic Oscillator functions'

    contains 

      
        recursive function Laguerre_P(n,a,x) result(s0)   !Legendre_P
            use numtype
            implicit none
            integer, intent(in) :: n
            real(dp), intent(in) :: x, a
            real(dp) :: s0

            !logic statements!
            if(n < 0) then 
                s0 = 0
            else if (n==0) then
                s0 = 1._dp
            else if (n==0) then
                s0 = 1 + a - x
            else 
                s0 = ((2*n-1+a-x) * Laguerre_P(n-1,a,x) - (n-1+a) * Laguerre_P(n-2,a,x))/n
            end if    !end the logic statement

    
        end function Laguerre_P

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

end program harmonic_oscillator_nd

