!A program to execute the horner scheme for a polynomial
program cftest

    use numtype
    use taylor_cf_approx
    implicit none

    real(dp) :: x, y
    integer :: i, imax, n

    imax = 40 
    t_coef(0) = 1
    do i = 1,imax
        t_coef(i) = t_coef(i-1)/i
    end do

    print *, 1/t_coef(0:5)
    ! Exponential function

    x = 10.5_dp

    print *, 'taylor : ',x, exp(x), horner(t_coef,imax,x), horner(t_coef,imax-1, x)

    n = imax
    !t_coef(0:n) = coeff(0:n)

    call taylor_cfrac(t_coef,n,cf_coef)  
    print *, 'c-frac : ', x, exp(x), evalcf(cf_coef,n,x), evalcf(cf_coef,n-1,x)

    !stop 
    print *, '-----------------log(1+x)------------------'
! log (1+x)
    imax = 60
    t_coef(0)=0
    forall(i=1:imax) t_coef(i) = (-1)**(i-1)*1._dp/i

    print *, 1/t_coef(1:5)

    x = 8.5_dp

    print *, 'taylor : ',x, log(1+x), horner(t_coef,imax,x), horner(t_coef,imax-1, x)

    !t_coef(0:n) = coeff(0:n)

    n = imax
    call taylor_cfrac(t_coef,n,cf_coef)  
    print *, 'c-frac : ', x, log(1+x), evalcf(cf_coef,n,x), evalcf(cf_coef,n-1,x)
   
end program cftest

!function horner(x,nmax,an) result(y)
!
!    use numtype !for double precision
!    implicit none !dont assume declaration
!
!    real(dp), intent(in) :: x !look up what intent means later
!    integer, intent(in) :: nmax
!    real(dp), dimension(0:nmax), intent(in) :: an !declare a real array from range 0 to nmax
!    integer :: i 
!    real(dp) :: y 
!
!    y = an(nmax)
!    do i = nmax - 1, 0, -1 !go backwards from nmax to 0
!        y = an(i) + x * y
!
!

 !   end do

!end function horner

