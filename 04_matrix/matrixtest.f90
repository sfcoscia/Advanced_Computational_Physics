! This is a program designed to use the linear algebra package for f90
program matrixtest
  use numtype
  implicit none
  integer , parameter :: ndim = 5, lwork = 5*ndim
  integer :: nn, i, info, j
  complex(dp), dimension(ndim,ndim) :: a, b, c, d, e ,f
  complex(dp) :: summe, work(lwork) !summe = sum
  real(dp), dimension(ndim) :: ww
  real(dp) :: rwork(3*ndim)

  nn = 3 !matrix dimension n x n
  ! Fortran stores matrices columnwise
  ! Take a = 7 0 0
  !          0 1 -i
  !          0 i -1
  !  NOTE : z0 represents 0
  !  NOTE:: one represents 1
  !  NOTE:: i1 represents complex i
  a(1:nn,1:nn) = reshape( (/7*one, z0, z0, z0, one, i1, z0, -i1, -one/) , (/ nn, nn/))!this is how we define arrays
  ! reshape sets it to an n x n matrix (nn x nn)

  ! Now for b = 1 0  3
  !             0 2i 0
  !             i 0 -5i
  b(1:nn,1:nn) = reshape( (/one, z0, i1, z0, 2*i1, z0, 3*one, z0, -5*i1/) , (/ nn, nn/))!this is how we define arrays

  print *, 'problem 2.9 zetilli'
  print *, 'matrix A & B'
  do i = 1, nn
    print '(10f8.2)', a(i,1:nn)               !10f means fixed format
  end do
  print *, '----------------------'
  do i = 1, nn
    print '(10f8.2)', b(i,1:nn)               !10f means fixed format
  end do


  print *, '-------Hermitian?------'
  
  c = conjg(transpose(a(1:nn,1:nn)))
  do i = 1, nn
  print '(10f8.2)', c(i,1:nn)-a(i,1:nn)               !iff A-C=0 is A hermitian
  end do
  print *, '----------------------'
  c = conjg(transpose(b(1:nn,1:nn)))
  do i = 1, nn
    print '(10f8.2)', c(i,1:nn)-b(i,1:nn)               
  end do

  print *, '------Commutator------'
  
  ! Calculate
  ! [A,B] (commutator)
  print *, '[A,B]'
  c(1:nn, 1:nn) = matmul(a(1:nn,1:nn),b(1:nn,1:nn)) - matmul(b(1:nn,1:nn),a(1:nn,1:nn))
  do i = 1, nn
    print '(10f8.2)', c(i,1:nn)
  end do


  print *, '------Trace------'
  summe = 0
  do i=1, nn
    summe = summe + c(i,i)
  end do
  print *, 'Trace([A,B])', summe, sum( (/(c(i,i),i=1, size(c,1))/)) !calculates 
  ! sum of the diagonals 
  !
  print *, '------Eigenvalues & Eigenvectors------'

  e(1:nn,1:nn) = a(1:nn, 1:nn)

  ! Call eigenvalue/eigenvector function for lapack
  call ZHEEV('v','u', nn, e, ndim, ww, WORK, LWORK, RWORK, INFO)
  do i = 1, nn
    print *, 'eigenvalues : ', ww(i)
    print '(10f8.2)', e(1:nn,i)
  end do

  d(1:nn,1:nn) = matmul(conjg(transpose(e(1:nn,1:nn))), e(1:nn, 1:nn))

  print *, '-----Orthogonality  : <phi_i | phi_j> = delta_ij -----'
  do i=1, nn
    do j=1, nn
      print *, i, j, dot_product(e(1:nn,i), e(1:nn,j)), d(i,j)
    end do
  end do

  print *, '-----Completeness : Sum(|phi_j> <phi_i) = Identity matrix -------'

  c(1:nn,1:nn) = matmul(conjg(transpose(e(1:nn,1:nn))), e(1:nn, 1:nn))
  do i = 1, nn
    print '(10f8.2)', c(i,1:nn)
  end do
end program matrixtest



