module deriv_fft
    implicit none
    include 'fftw3.f'

    contains

    subroutine fft_r2c_(rin, cout)
        implicit none
        real, intent(in) :: rin(:,:)
        complex, intent(out) :: cout(:,:)

        integer*8 :: r2c_plan
        integer :: mx, my

        mx = size(rin, 1)
        my = size(rin, 2)

        call sfftw_plan_dft_r2c_2d(r2c_plan, mx, my, rin, cout, FFTW_ESTIMATE)
        call sfftw_execute_dft_r2c(r2c_plan, rin, cout)

    end subroutine fft_r2c_

    subroutine fft_c2r_(cin, rout)
        implicit none
        complex, intent(inout) :: cin(:,:)
        real, intent(out) :: rout(:,:)

        integer*8 :: c2r_plan
        integer :: mx, my

        mx = size(rout, 1)
        my = size(rout, 2)

        call sfftw_plan_dft_c2r_2d(c2r_plan, mx, my, cin, rout, FFTW_ESTIMATE)
        call sfftw_execute_dft_c2r(c2r_plan, cin, rout)

    end subroutine fft_c2r_

end module deriv_fft
