module get_rhs
    use,intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
#ifdef USE_OPENMP
    use :: omp_lib
#endif 
    implicit none
    

    contains
    subroutine mom_rhs_compute(rhs_u,rhs_v,rhs_w,f,g)
        ! very important note:
        ! in rhs calculation there shouldnt be any pressure gradient.
        ! pressure effects must be added at intermediate velocity calculation
        type(grid_type)              :: g       
        type(field_type)             :: f
        real(C_DOUBLE)               :: rhs_u(0:g%nx+1,0:g%ny+1,0:g%nz+1)
        real(C_DOUBLE)               :: rhs_v(0:g%nx+1,1:g%ny+1,0:g%nz+1)
        real(C_DOUBLE)               :: rhs_w(0:g%nx+1,0:g%ny+1,0:g%nz+1)

        integer :: i,j,k,ip,im,kp,km,jp,jm

        real(C_DOUBLE) :: diff_ux,diff_uy,diff_uz
        real(C_DOUBLE) :: diff_vx,diff_vy,diff_vz
        real(C_DOUBLE) :: diff_wx,diff_wy,diff_wz

        real(C_DOUBLE) :: uu_p,uu_m,uv_p,uv_m,uw_p,uw_m
        real(C_DOUBLE) :: vu_p,vu_m,vv_p,vv_m,vw_p,vw_m
        real(C_DOUBLE) :: wu_p,wu_m,ww_p,ww_m,wv_p,wv_m

        real(C_DOUBLE) :: d_dx,d_dy,d_dz
        real(C_DOUBLE) :: d_ddx,d_ddy,d_ddz,div_re
        ! compute divisions and mult. before
        d_dx = 1.0d0/g%dx
        d_dy = 1.0d0/g%dy
        d_dz = 1.0d0/g%dz
        d_ddx = 1.0d0/g%dx**2
        d_ddy = 1.0d0/g%dy**2
        d_ddz = 1.0d0/g%dz**2
        div_re = 1.0d0/g%re
        !$omp do collapse(2) schedule(static) private(i,j,k)  
        do k = 1, g%nz
            do j = 1, g%ny
                do i = 1, g%nx                                                                                                                  
                    
                    ip = i+1
                    im = i-1
                    jp = j+1
                    jm = j-1
                    kp = k+1
                    km = k-1
                     
                    uu_p = 0.25d0*(f%un(ip,j,k)+f%un(i,j,k))**2
                    uu_m = 0.25d0*(f%un(i,j,k)+f%un(im,j,k))**2

                    uv_p = 0.25d0*(f%un(i,jp,k)+f%un(i,j,k))*(f%vn(i,jp,k)+f%vn(im,jp,k))
                    uv_m = 0.25d0*(f%un(i,j,k)+f%un(i,jm,k))*(f%vn(i,j,k)+f%vn(im,j,k))

                    uw_p = 0.25d0*(f%un(i,j,k)+f%un(i,j,kp))*(f%wn(i,j,kp)+f%wn(im,j,kp))
                    uw_m = 0.25d0*(f%un(i,j,k)+f%un(i,j,km))*(f%wn(i,j,k)+f%wn(im,j,k))

                    diff_ux = (f%un(im,j,k)-2.0d0*f%un(i,j,k)+f%un(ip,j,k))*d_ddx
                    diff_uy = (f%un(i,jm,k)-2.0d0*f%un(i,j,k)+f%un(i,jp,k))*d_ddy
                    diff_uz = (f%un(i,j,km)-2.0d0*f%un(i,j,k)+f%un(i,j,kp))*d_ddz

                    rhs_u(i,j,k) = &
                        -(uu_p-uu_m)*d_dx &
                        -(uv_p-uv_m)*d_dy &
                        -(uw_p-uw_m)*d_dz &
                        + div_re*(diff_ux + diff_uy + diff_uz) &
                        + f%b_x
                end do
            end do
        end do
        !$omp end do nowait
        ! y-momentum
        !$omp do collapse(2) schedule(static) private(i,j,k)
        do k = 1, g%nz
            do j = 2, g%ny
                do i = 1, g%nx

                    ip = i+1
                    im = i-1
                    jp = j+1
                    jm = j-1
                    kp = k+1
                    km = k-1

                    vu_p = 0.25d0*(f%vn(i,j,k)+f%vn(ip,j,k))*(f%un(ip,j,k)+f%un(ip,jm,k))
                    vu_m = 0.25d0*(f%vn(i,j,k)+f%vn(im,j,k))*(f%un(i,j,k)+f%un(i,jm,k))

                    vv_p = 0.25d0*(f%vn(i,j,k)+f%vn(i,jp,k))**2
                    vv_m = 0.25d0*(f%vn(i,j,k)+f%vn(i,jm,k))**2

                    vw_p = 0.25d0*(f%vn(i,j,kp)+f%vn(i,j,k))*(f%wn(i,j,kp)+f%wn(i,jm,kp))
                    vw_m = 0.25d0*(f%vn(i,j,km)+f%vn(i,j,k))*(f%wn(i,j,k)+f%wn(i,jm,k))

                    diff_vx = (f%vn(im,j,k)-2.0d0*f%vn(i,j,k)+f%vn(ip,j,k))*d_ddx
                    diff_vy = (f%vn(i,jm,k)-2.0d0*f%vn(i,j,k)+f%vn(i,jp,k))*d_ddy
                    diff_vz = (f%vn(i,j,km)-2.0d0*f%vn(i,j,k)+f%vn(i,j,kp))*d_ddz

                    rhs_v(i,j,k) = &
                        -(vu_p-vu_m)*d_dx &
                        -(vv_p-vv_m)*d_dy &
                        -(vw_p-vw_m)*d_dz &
                        + div_re*(diff_vx + diff_vy + diff_vz) & 
                        + f%b_y
                end do
            end do
        end do
        !$omp end do nowait
        ! z-momentum
        !$omp do collapse(2) schedule(static) private(i,j,k) 
        do k = 1, g%nz
            do j = 1, g%ny
                do i = 1, g%nx

                    ip = i+1
                    im = i-1
                    jp = j+1
                    jm = j-1
                    kp = k+1
                    km = k-1

                    wu_p = 0.25d0*(f%wn(i,j,k)+f%wn(ip,j,k))*(f%un(ip,j,k)+f%un(ip,j,km))
                    wu_m = 0.25d0*(f%wn(i,j,k)+f%wn(im,j,k))*(f%un(i,j,k)+f%un(i,j,km))

                    ww_p = 0.25d0*(f%wn(i,j,k)+f%wn(i,j,kp))**2
                    ww_m = 0.25d0*(f%wn(i,j,k)+f%wn(i,j,km))**2

                    wv_p = 0.25d0*(f%wn(i,j,k)+f%wn(i,jp,k))*(f%vn(i,jp,k)+f%vn(i,jp,km))
                    wv_m = 0.25d0*(f%wn(i,j,k)+f%wn(i,jm,k))*(f%vn(i,j,k)+f%vn(i,j,km))

                    diff_wx = (f%wn(im,j,k)-2.0d0*f%wn(i,j,k)+f%wn(ip,j,k))*d_ddx
                    diff_wy = (f%wn(i,jm,k)-2.0d0*f%wn(i,j,k)+f%wn(i,jp,k))*d_ddy
                    diff_wz = (f%wn(i,j,km)-2.0d0*f%wn(i,j,k)+f%wn(i,j,kp))*d_ddz

                    rhs_w(i,j,k) = &
                        -(wu_p-wu_m)*d_dx &
                        -(wv_p-wv_m)*d_dy &
                        -(ww_p-ww_m)*d_dz &
                        + div_re*(diff_wx + diff_wy + diff_wz) &
                        + f%b_z 

                end do
            end do
        end do
        !$omp end do 
    end subroutine mom_rhs_compute

end module get_rhs