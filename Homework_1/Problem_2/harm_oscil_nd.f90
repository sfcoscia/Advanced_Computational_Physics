!A program to calculate the N-dimensional radial wavefunction for the harmonic
!oscillator
program harmonic_oscillator_nd

    use numtype
    implicit none

    real(dp) :: wfunc(0:100), x, y
    integer :: nm, n, i
    real(dp) :: v

    nm = 5
    v = 1
    
     ! Test of the Laguerre polynomial function 
     ! do i = 0,100
     !   x = -2 + i * 1._dp/10
     !   write(17, *) x, Laguerre_P(3,0._dp,x) 
     ! end do

    do i = 0,200
        x = 0 + i * 4._dp/100
        ! Run for l = 0,1 , D = 2,3
        call wavefunction(nm,0,v,x,2,wfunc)
        write(02, *) x, wfunc(0:nm) 
        call wavefunction(nm,0,v,x,3,wfunc)
        write(03, *) x, wfunc(0:nm) 
        call wavefunction(nm,1,v,x,2,wfunc)
        write(12, *) x, wfunc(0:nm) 
        call wavefunction(nm,1,v,x,3,wfunc)
        write(13, *) x, wfunc(0:nm) 
    end do


     print *, 'Computing N-D Harmonic Oscillator functions ...'

    contains 
        ! This is a function to calculate the wavefunction from the Laguerre polynomials
        subroutine wavefunction(nm,l,v,r,D,lpoly)
  
          use numtype
          implicit none
          
          integer, intent(in) :: nm, l, D
          real(dp), intent(in) :: r, v
          real(dp), dimension(0:nm), intent(out) :: lpoly
          integer :: n
          real(dp) :: coeff

          do n = 0, nm 
            ! Multiply eaech Laguerre polynomial by the respective coefficient
            coeff = v**(0.25_dp)*(2*gamma(n+1._dp)/gamma(n+l+D/2._dp))**(0.5_dp)*exp(-v*r**2/2._dp)*(v*r**2)**(l/2._dp+(D-1)/4._dp)
            lpoly(n) = coeff*Laguerre_P(n,l+D/2._dp - 1,v*r**2)
          end do

        end subroutine wavefunction

        !This is a recursive function to calculate the Associated Laguerre polynomals
        recursive function Laguerre_P(n,a,x) result(s0)   
            use numtype
            implicit none

            integer, intent(in) :: n
            real(dp), intent(in) :: x, a
            real(dp) :: s0

            if(n < 0) then
                s0 = 0
            else if (n==0) then
                s0 = 1._dp
            else if (n==1) then
                s0 = 1 + a - x
            else 
                ! Calculation of the Laguerre polynomials for n>1
                s0 = ((2*n-1+a-x) * Laguerre_P(n-1,a,x) - (n-1+a) * Laguerre_P(n-2,a,x))/n
            end if    
    
        end function Laguerre_P

end program harmonic_oscillator_nd
