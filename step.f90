!--------------------------!
!                          !
!       Time-stepper       !
!          module          !
!                          !
!--------------------------! 
! 
! authors: Dr.-Ing. Davide Gatti
!          B.Sc. Ahmet Cumhur
! 
! date:    28.04.26
! 

module step
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
#ifdef USE_IBM_G
    use :: ibmm, only: ibm_type
#endif
    implicit none

contains
#ifdef USE_IBM_G
    subroutine momentum(f, g,ibm)
#else
    subroutine momentum(f, g)
#endif
        implicit none
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g
#ifdef USE_IBM_G
        type(ibm_type), intent(in)      :: ibm
#endif

        integer :: i,j,k,ip,im,kp,km,jp,jm

        real(C_DOUBLE) :: diff_ux,diff_uy,diff_uz
        real(C_DOUBLE) :: diff_vx,diff_vy,diff_vz
        real(C_DOUBLE) :: diff_wx,diff_wy,diff_wz

        real(C_DOUBLE) :: uu_p,uu_m,uv_p,uv_m,uw_p,uw_m
        real(C_DOUBLE) :: vu_p,vu_m,vv_p,vv_m,vw_p,vw_m
        real(C_DOUBLE) :: wu_p,wu_m,ww_p,ww_m,wv_p,wv_m

        real(C_DOUBLE) :: dpx,dpy,dpz

        ! x-momentum
        do i = 1, g%nx
            do j = 1, g%ny
                do k = 1, g%nz

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

                    diff_ux = (f%un(im,j,k)-2.0d0*f%un(i,j,k)+f%un(ip,j,k))/g%dx**2
                    diff_uy = (f%un(i,jm,k)-2.0d0*f%un(i,j,k)+f%un(i,jp,k))/g%dy**2
                    diff_uz = (f%un(i,j,km)-2.0d0*f%un(i,j,k)+f%un(i,j,kp))/g%dz**2

                    dpx = (f%pn(i,j,k)-f%pn(im,j,k))/g%dx 
#ifdef USE_IBM_G
                    ! we apply the found coefficient as a coefficient to the eq.
                    f%us(i,j,k) = 1.0d0/(1.0d0-ibm%coef_u_lap(i,j,k)*g%dt*(1.0d0/g%re))*(f%un(i,j,k) + g%dt*( &
                        -(uu_p-uu_m)/g%dx &
                        -(uv_p-uv_m)/g%dy &
                        -(uw_p-uw_m)/g%dz &
                        - dpx  &
                        + f%b_x  &
                        + (1.0d0/g%re)*(diff_ux + diff_uy + diff_uz) ))
#else
                    f%us(i,j,k) = f%un(i,j,k) + g%dt*( &
                        -(uu_p-uu_m)/g%dx &
                        -(uv_p-uv_m)/g%dy &
                        -(uw_p-uw_m)/g%dz &
                        - dpx  &
                        + f%b_x  &
                        + (1.0d0/g%re)*(diff_ux + diff_uy + diff_uz) )
#endif
                end do
            end do
        end do

        ! y-momentum
        do i = 1, g%nx
            do j = 2, g%ny
                do k = 1, g%nz

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

                    diff_vx = (f%vn(im,j,k)-2.0d0*f%vn(i,j,k)+f%vn(ip,j,k))/g%dx**2
                    diff_vy = (f%vn(i,jm,k)-2.0d0*f%vn(i,j,k)+f%vn(i,jp,k))/g%dy**2
                    diff_vz = (f%vn(i,j,km)-2.0d0*f%vn(i,j,k)+f%vn(i,j,kp))/g%dz**2

                    dpy = (f%pn(i,j,k)-f%pn(i,jm,k))/g%dy
#ifdef USE_IBM_G
                    f%vs(i,j,k) = 1.0d0/(1.0d0-ibm%coef_v_lap(i,j,k)*g%dt*(1.0d0/g%re))*(f%vn(i,j,k) + g%dt*( &
                        -(vu_p-vu_m)/g%dx &
                        -(vv_p-vv_m)/g%dy &
                        -(vw_p-vw_m)/g%dz &
                        - dpy &
                        + f%b_y &
                        + (1.0d0/g%re)*(diff_vx + diff_vy + diff_vz) ))
#else 
                    f%vs(i,j,k) = f%vn(i,j,k) + g%dt*( &
                        -(vu_p-vu_m)/g%dx &
                        -(vv_p-vv_m)/g%dy &
                        -(vw_p-vw_m)/g%dz &
                        - dpy &
                        + f%b_y &
                        + (1.0d0/g%re)*(diff_vx + diff_vy + diff_vz) )
#endif
                end do
            end do
        end do

        ! z-momentum
        do i = 1, g%nx
            do j = 1, g%ny
                do k = 1, g%nz

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

                    diff_wx = (f%wn(im,j,k)-2.0d0*f%wn(i,j,k)+f%wn(ip,j,k))/g%dx**2
                    diff_wy = (f%wn(i,jm,k)-2.0d0*f%wn(i,j,k)+f%wn(i,jp,k))/g%dy**2
                    diff_wz = (f%wn(i,j,km)-2.0d0*f%wn(i,j,k)+f%wn(i,j,kp))/g%dz**2

                    dpz = (f%pn(i,j,k)-f%pn(i,j,km))/g%dz
#ifdef USE_IBM_G
                    f%ws(i,j,k) = 1.0d0/(1.0d0-ibm%coef_w_lap(i,j,k)*g%dt*(1.0d0/g%re))*(f%wn(i,j,k) + g%dt*( &
                        -(wu_p-wu_m)/g%dx &
                        -(wv_p-wv_m)/g%dy &
                        -(ww_p-ww_m)/g%dz &
                        - dpz &
                        + f%b_z &
                        + (1.0d0/g%re)*(diff_wx + diff_wy + diff_wz) ))
#else
                    f%ws(i,j,k) = f%wn(i,j,k) + g%dt*( &
                        -(wu_p-wu_m)/g%dx &
                        -(wv_p-wv_m)/g%dy &
                        -(ww_p-ww_m)/g%dz &
                        - dpz &
                        + f%b_z &
                        + (1.0d0/g%re)*(diff_wx + diff_wy + diff_wz) )
#endif
                end do
            end do
        end do

    end subroutine momentum


    subroutine corrector(f, g)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g

        integer :: i,j,k,im,km

        do i = 1, g%nx
            do j = 1, g%ny
                do k = 1, g%nz
                    f%un(i,j,k) = f%us(i,j,k) - g%dt*(f%pc(i,j,k)-f%pc(i-1,j,k))/g%dx
                end do
            end do
        end do

        do i = 1, g%nx
            do j = 2, g%ny
                do k = 1, g%nz
                    f%vn(i,j,k) = f%vs(i,j,k) - g%dt*(f%pc(i,j,k)-f%pc(i,j-1,k))/g%dy
                end do
            end do
        end do

        do i = 1, g%nx
            do j = 1, g%ny
                do k = 1, g%nz
                    f%wn(i,j,k) = f%ws(i,j,k) - g%dt*(f%pc(i,j,k)-f%pc(i,j,k-1))/g%dz
                end do
            end do
        end do
        
        f%pn = f%pn + f%pc

    end subroutine corrector


    subroutine divU(f, g)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g

        integer :: i,j,k

        do i = 1, g%nx
            do j = 1, g%ny
                do k = 1, g%nz

                    f%rhs(i,j,k) = ( &
                        (f%us(i+1,j,k)-f%us(i,j,k))/g%dx &
                      + (f%vs(i,j+1,k)-f%vs(i,j,k))/g%dy &
                      + (f%ws(i,j,k+1)-f%ws(i,j,k))/g%dz ) / g%dt

                end do
            end do
        end do

    end subroutine divU
    
    real function get_cfl(f,g)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g

        integer :: i,j,k
        get_cfl = max(maxval(abs(f%un/g%dx)),&
                      maxval(abs(f%vn/g%dy)),&
                      maxval(abs(f%wn/g%dz)))
    end function get_cfl

end module step
