module deriv_fft
    implicit none
    include 'fftw3.f'

    contains

    subroutine fft_r2c_(rin, cout)
        implicit none
        double precision, intent(in) :: rin(:,:)
        complex*16, intent(out) :: cout(:,:)

        integer*8 :: r2c_plan
        integer :: mx, my

        mx = size(rin, 1)
        my = size(rin, 2)

        call dfftw_plan_dft_r2c_2d(r2c_plan, mx, my, rin, cout, FFTW_ESTIMATE)
        call dfftw_execute_dft_r2c(r2c_plan, rin, cout)

    end subroutine fft_r2c_

    subroutine fft_c2r_(cin, rout)
        implicit none
        complex*16, intent(inout) :: cin(:,:)
        double precision, intent(out) :: rout(:,:)

        integer*8 :: c2r_plan
        integer :: mx, my

        mx = size(rout, 1)
        my = size(rout, 2)

        call dfftw_plan_dft_c2r_2d(c2r_plan, mx, my, cin, rout, FFTW_ESTIMATE)
        call dfftw_execute_dft_c2r(c2r_plan, cin, rout)

    end subroutine fft_c2r_

    subroutine deriv(rin, rout, axis, order)
        implicit none
        integer, intent(in) :: axis, order
        double precision, intent(in) :: rin(:,:)
        double precision, intent(out) :: rout(:,:)

        complex*16 :: carr(size(rin, 1)/2+1, size(rin, 2))
        integer :: i
        complex*16 :: ki

        call fft_r2c_(rin, carr)

        if(axis .eq. 1) then
            ! deriv along the 'short' dimension in k-space.
            do i = 1, size(carr, 1)
                ki = cmplx(0.0, i-1)
                carr(i,:) = carr(i,:)*(ki**order)
            enddo
        else if(axis .eq. 2) then
            ! deriv along the 'long' dimension in k-space.
            do i = 1, size(carr, 2)/2+1
                ki = cmplx(0.0, i-1)
                carr(:,i) = carr(:,i)*(ki**order)
            enddo
            do i = size(carr, 2)/2+2, size(carr, 2)
                ki = cmplx(0.0, i-size(carr, 2)-1)
                carr(:,i) = carr(:,i)*(ki**order)
            enddo
        endif
        call fft_c2r_(carr, rout)
    end subroutine deriv

end module deriv_fft
