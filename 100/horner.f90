!A program to execute the horner scheme for a polynomial
program horner_test

    use numtype
    implicit none

    real(dp) :: coeff(0:100), x, y
    integer :: i, imax
    real(dp), external :: horner

    imax = 60 
    coeff(0) = 1
    do i = 1,imax
        coeff(i) = coeff(i-1)/i
    end do

    print *, 1/coeff(0:5)
    ! Exponential function

    x = 10.5_dp
    y = horner(x,imax,coeff)

    print *, x, exp(x), y

    ! log(1+x)
    imax = 60 
    coeff(0) = 0
    !coefficients of the log
    forall(i=1:imax) coeff(i) = (-1)**(i-1) * 1._dp/i

    print *, 1/coeff(1:5)

    x = 0.5_dp
    y = horner(x,imax,coeff)
    print *, x, log(1+x), y

    x = 1.5_dp
    y = horner(x,imax,coeff)
    print *, x, log(1+x), y
end program horner_test

function horner(x,nmax,an) result(y)

    use numtype !for double precision
    implicit none !dont assume declaration

    real(dp), intent(in) :: x !look up what intent means later
    integer, intent(in) :: nmax
    real(dp), dimension(0:nmax), intent(in) :: an !declare a real array from range 0 to nmax
    integer :: i 
    real(dp) :: y 

    y = an(nmax)
    do i = nmax - 1, 0, -1 !go backwards from nmax to 0
        y = an(i) + x * y



    end do

end function horner

