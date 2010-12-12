subroutine fft_r2c(rin, cout, M, N)
    use deriv_fft, only : fft_r2c_
    implicit none
    integer, intent(in) :: M, N
    double precision, intent(in) :: rin(M, N)
    complex*16, intent(out) :: cout(M/2+1, N)

    call fft_r2c_(rin, cout)
end subroutine fft_r2c

subroutine fft_c2r(cin, rout, M, N)
    use deriv_fft, only : fft_c2r_
    implicit none
    integer, intent(in) :: M, N
    complex*16, intent(inout) :: cin(M/2+1, N)
    double precision, intent(out) :: rout(M, N)

    call fft_c2r_(cin, rout)

end subroutine fft_c2r

subroutine deriv(rin, rout, axis, order)
    use deriv_fft, only : deriv_ => deriv
    implicit none
    integer, intent(in) :: axis, order
    double precision, intent(in) :: rin(:,:)
    double precision, intent(out) :: rout(:,:)

    double precision :: norm

    call deriv_(rin, rout, axis, order)

    norm = 1.0D0 / (size(rout, 1) * size(rout, 2))

    rout = rout * norm

end subroutine deriv

subroutine vel_xy(rin, vx, vy)
    use deriv_fft, only : deriv
    implicit none
    double precision, intent(in) :: rin(:,:)
    double precision, intent(out) :: vx(:,:), vy(:,:)

    call deriv(rin, vx, axis=2, order=1)
    call deriv(rin, vy, axis=1, order=1)

    vy = - vy

end subroutine vel_xy

subroutine gamma_scal(u, v, gam_scl)
    use deriv_fft, only : deriv
    implicit none
    double precision, intent(in) :: u(:,:), v(:,:)
    double precision, intent(out) :: gam_scl(:,:)

    double precision :: u_x(size(u,1), size(u,2))
    double precision :: v_y(size(u,1), size(u,2))

    call deriv(u, u_x, axis=1, order=1)
    call deriv(v, v_y, axis=2, order=1)

    gam_scl = 0.5 * (u_x + v_y)

end subroutine gamma_scal

subroutine gamma_rotation(u, v, gam_rot)
    use deriv_fft, only : deriv
    implicit none
    double precision, intent(in) :: u(:,:), v(:,:)
    double precision, intent(out) :: gam_rot(:,:)

    double precision :: u_y(size(u,1), size(u,2))
    double precision :: v_x(size(u,1), size(u,2))

    call deriv(u, u_y, axis=2, order=1)
    call deriv(v, v_x, axis=1, order=1)

    gam_rot = 0.5 * (v_x - u_y)

end subroutine gamma_rotation

subroutine gamma_stretch(u, v, gam_str)
    use deriv_fft, only : deriv
    implicit none
    double precision, intent(in) :: u(:,:), v(:,:)
    double precision, intent(out) :: gam_str(:,:)

    double precision :: u_tmp(size(u,1), size(u,2))
    double precision :: v_tmp(size(v,1), size(v,2))

    call deriv(u, u_tmp, axis=1, order=1)
    call deriv(v, v_tmp, axis=2, order=1)

    gam_str = u_tmp - v_tmp
    gam_str = gam_str * gam_str

    call deriv(u, u_tmp, axis=2, order=1)
    call deriv(v, v_tmp, axis=1, order=1)

    u_tmp = u_tmp + v_tmp
    u_tmp = u_tmp * u_tmp

    gam_str = gam_str + u_tmp

    gam_str = 0.5 * sqrt(gam_str)

end subroutine gamma_stretch

subroutine make_gauss(rout, sigma)
    use deriv_fft, only : make_gauss_ => make_gauss, TWO_PI
    implicit none
    double precision, intent(out) :: rout(:,:)
    double precision, intent(in) :: sigma

    call make_gauss_(rout, sigma, TWO_PI)

end subroutine make_gauss

subroutine convolve_rr(rin1, rin2, rout)
    use deriv_fft, only : fft_r2c_, fft_c2r_
    implicit none
    double precision, intent(in) :: rin1(:,:), rin2(:,:)
    double precision, intent(out) :: rout(:,:)

    double complex :: c1(size(rin1, 1)/2+1, size(rin1, 2)), &
                      c2(size(rin2, 1)/2+1, size(rin2, 2))

    call fft_r2c_(rin1, c1)
    call fft_r2c_(rin2, c2)

    c1 = c1 * c2

    call fft_c2r_(c1, rout)

end subroutine convolve_rr
