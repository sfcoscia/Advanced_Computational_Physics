!Set up the module 
module numtype

    !Save ? look up later
    save
    !Tell computer whether to use double or single precision
    !Natural double precision is 15,307
    integer, parameter :: dp = selected_real_kind(15,307)
    real(dp), parameter :: pi = 4 * atan(1._dp), tiny = 1.e-20_dp !_dp indicates 1 with 15 zeroes
    complex(dp), parameter :: i1 = (0._dp,1._dp), one = (1._dp,0._dp), z0 = (0._dp, 0._dp) ! Complex number i


end module numtype

