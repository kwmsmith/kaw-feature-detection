subroutine fft_r2c(rin, cout, M, N)
    use deriv_fft, only : fft_r2c_
    implicit none
    integer, intent(in) :: M, N
    real, intent(in) :: rin(M, N)
    complex, intent(out) :: cout(M/2+1, N)

    call fft_r2c_(rin, cout)
end subroutine fft_r2c

subroutine fft_c2r(cin, rout, M, N)
    use deriv_fft, only : fft_c2r_
    implicit none
    integer, intent(in) :: M, N
    complex, intent(inout) :: cin(M/2+1, N)
    real, intent(out) :: rout(M, N)

    call fft_c2r_(cin, rout)

end subroutine fft_c2r
