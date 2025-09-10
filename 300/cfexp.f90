! continued fraction code for exponential
program cfexp

  use numtype
  implicit none

  real(dp) :: x

  x = 2*pi/3

  ! exponentials
  print *, x, exp(x), exxp(x), eexp(x)

  x = 2.7_dp
  !logarithms
  print *, x, log(x), logg(x) 

  contains
    ! 1 way to calcualte exponential
    function logg(x) result(s)
    
      implicit none
      real(dp), intent(in) :: x
      real(dp) :: s, z
      integer :: i
      integer, parameter :: imax = 30
      
      z = x-1
      !continued fractions are built in here
      s = 2*imax + 1 
      do i = 2*imax + 1, 1, -2
        s = i + ((i+1)/2)**2*z / (i + 1 + (((i+1)/2)**2*z)/s) 
      end do
      s = z/s

    end function logg

    ! 1 way to calcualte logarithm(1+z)
    function eexp(z) result(s)
    
      implicit none
      real(dp), intent(in) :: z
      real(dp) :: x, s
      integer :: i
      integer, parameter :: imax = 30


      !continued fractions are built in here
      s = 0
      do i = imax, 1, -1
        s = (x/i)/(1+(x/i)-s)
      end do
      s= 1/(1-s)

    end function eexp



    ! 2nd way to calcualte exponential
    function exxp(x) result(s)
    
      implicit none
      real(dp), intent(in) :: x
      real(dp) :: s

      s = 1/(1-scf(1,x))

    end function exxp

    ! continued fraction loop
    recursive function scf(i, x) result(s) 
      implicit none
      real(dp) :: x, s
      integer :: i
      integer, parameter :: imax = 50

      if (i>= imax ) then 
        s = 0._dp
      else
        s = (x/i)/(1 + (x/i)-scf(i+1,x)) 
      end if  

    end function scf  



end program cfexp
