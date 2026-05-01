# 1 "C:/Users/mehme/OneDrive/Desktop/thesis/nse/step.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "C:/Users/mehme/OneDrive/Desktop/thesis/nse/step.f90"
module step
    implicit none 
    
    contains
    subroutine moment(un,us,vn,vs,wn,ws,pn,re,dt,nx,ny,nz,dx,dy,dz)
        use, intrinsic :: iso_c_binding
        implicit none
        real(C_DOUBLE),intent(inout) :: un(nx,ny+2,nz),vn(nx,ny+1,nz),wn(nx,ny+2,nz),us(nx,ny+2,nz),vs(nx,ny+1,nz),ws(nx,ny+2,nz)
        real(C_DOUBLE),intent(inout) :: pn(nx,ny,nz)
        real(C_DOUBLE),intent(in) :: dx,dy,dz,re,dt
        integer, intent(in) :: nx,ny,nz
        integer :: i,j,k

        real(C_DOUBLE) :: diff_ux,diff_uy,diff_uz
        real(C_DOUBLE) :: diff_vx,diff_vy,diff_vz
        real(C_DOUBLE) :: diff_wx,diff_wy,diff_wz

        real(C_DOUBLE) :: uu_p,uu_m,uv_p,uv_m,uw_p,uw_m
        real(C_DOUBLE) :: vu_p,vu_m,vv_p,vv_m,vw_p,vw_m
        real(C_DOUBLE) :: wu_p,wu_m,ww_p,ww_m,wv_p,wv_m

        real(C_DOUBLE) :: dpx,dpy,dpz
        integer :: ip,im,kp,km
        ! x axis
        do i=1, nx
            do j= 2,ny-1
                do k=1,nz
                    ! we apply the periodic bc 
                    ip = i + 1
                    im = i - 1
                    kp = k + 1
                    km = k - 1

                    if (ip > nx) ip = 1
                    if (im < 1 ) im = nx

                    if (kp > nz) kp = 1
                    if (km < 1 ) km = nz

                    uu_p = 0.25d0*(un(ip,j,k)+un(i,j,k))*(un(ip,j,k)+un(i,j,k))
                    uu_m = 0.25d0*(un(i,j,k)+un(im,j,k))*(un(i,j,k)+un(im,j,k))
                    
                    uv_p = 0.25d0*(un(i,j+1,k)+un(i,j,k))*(vn(i,j+1,k)+vn(im,j+1,k))
                    uv_m = 0.25d0*(un(i,j,k)+un(i,j-1,k))*(vn(i,j,k)+vn(im,j,k))
                    
                    uw_p = 0.25d0*(un(i,j,k)+un(i,j,kp))*(wn(i,j,kp)+wn(im,j,kp))
                    uw_m = 0.25d0*(un(i,j,k)+un(i,j,km))*(wn(i,j,k)+wn(im,j,k))

                    diff_ux = (un(im,j,k)-2.0d0*un(i,j,k)+un(ip,j,k))/dx**2
                    diff_uy = (un(i,j-1,k)-2.0d0*un(i,j,k)+un(i,j+1,k))/dy**2
                    diff_uz = (un(i,j,km)-2.0d0*un(i,j,k)+un(i,j,kp))/dz**2

                    dpx = (pn(i,j,k)-pn(im,j,k))/dx

                    us(i,j,k) = un(i,j,k)+dt*(&
                                -(uu_p-uu_m)/dx+&
                                -(uv_p-uv_m)/dy+&
                                -(uw_p-uw_m)/dz+&
                                -dpx+&
                                (1.0d0/re)*(diff_ux + diff_uy + diff_uz))
                end do 
            end do 
        end do 
        
        ! y axis
        do i = 1,nx
            do j = 2, ny+1
                do k =1,nz

                    ip = i + 1
                    im = i - 1
                    kp = k + 1
                    km = k - 1

                    if (ip > nx) ip = 1
                    if (im < 1 ) im = nx

                    if (kp > nz) kp = 1
                    if (km < 1 ) km = nz


                    vu_p = 0.25d0*(vn(i,j,k)+vn(ip,j,k))*(un(ip,j,k)+un(ip,j-1,k))
                    vu_m = 0.25d0*(vn(i,j,k)+vn(im,j,k))*(un(i,j,k)+un(i,j-1,k))
                    
                    vv_p = 0.25d0*(vn(i,j,k)+vn(i,j+1,k))*(vn(i,j,k)+vn(i,j+1,k))
                    vv_m = 0.25d0*(vn(i,j,k)+vn(i,j-1,k))*(vn(i,j,k)+vn(i,j-1,k))
                    
                    vw_p = 0.25d0*(vn(i,j,kp)+vn(i,j,k))*(wn(i,j,kp)+wn(i,j-1,kp))
                    vw_m = 0.25d0*(vn(i,j,km)+vn(i,j,k))*(wn(i,j,k)+wn(i,j-1,k))

                    diff_vx = (vn(im,j,k)-2.0d0*vn(i,j,k)+vn(ip,j,k))/dx**2
                    diff_vy = (vn(i,j-1,k)-2.0d0*vn(i,j,k)+vn(i,j+1,k))/dy**2
                    diff_vz = (vn(i,j,km)-2.0d0*vn(i,j,k)+vn(i,j,kp))/dz**2

                    dpy = (pn(i,j,k)-pn(i,j-1,k))/dy

                    vs(i,j,k) = vn(i,j,k)+dt*(&
                                -(vu_p-vu_m)/dx+&
                                -(vv_p-vv_m)/dy+&
                                -(vw_p-vw_m)/dz+&
                                -dpy+&
                                (1.0d0/re)*(diff_vx+diff_vy+diff_vz))
                end do
            end do 
        end do 

        ! z axis
        do i = 1,nx
            do j = 2,ny+1
                do k = 1,nz
              
                    ip = i + 1
                    im = i - 1
                    kp = k + 1
                    km = k - 1

                    if (ip > nx) ip = 1
                    if (im < 1 ) im = nx

                    if (kp > nz) kp = 1
                    if (km < 1 ) km = nz
                    
                    wu_p = 0.25d0*(wn(i,j,k)+wn(ip,j,k))*(un(ip,j,k)+un(ip,j,km))
                    wu_m = 0.25d0*(wn(i,j,k)+wn(im,j,k))*(un(i,j,k)+un(i,j,km))
                    
                    ww_p = 0.25d0*(wn(i,j,k)+wn(i,j,kp))*(wn(i,j,k)+wn(i,j,kp))
                    ww_m = 0.25d0*(wn(i,j,k)+wn(i,j,km))*(wn(i,j,k)+wn(i,j,km))
                    
                    wv_p = 0.25d0*(wn(i,j,k)+wn(i,j+1,k))*(vn(i,j+1,k)+vn(i,j+1,km))
                    wv_m = 0.25d0*(wn(i,j,k)+wn(i,j-1,k))*(vn(i,j,k)+vn(i,j,km)) 
                    
                    diff_wx = (wn(im,j,k)-2.0d0*wn(i,j,k)+wn(ip,j,k))/dx**2
                    diff_wy = (wn(i,j-1,k)-2.0d0*wn(i,j,k)+wn(i,j+1,k))/dy**2
                    diff_wz = (wn(i,j,km)-2.0d0*wn(i,j,k)+wn(i,j,kp))/dz**2                     

                    dpz = (pn(i,j,k)-pn(i,j,km))/dz 
                    
                    ws(i,j,k) = wn(i,j,k)+dt*(&
                                -(wu_p-wu_m)/dx+&
                                -(wv_p-wv_m)/dy+&
                                -(ww_p-ww_m)/dz+&
                                -dpz+&
                                (1.0d0/re)*(diff_wx+diff_wy+diff_wz))
                end do 
            end do 
        end do


    end subroutine moment
    subroutine n_step (us,un,vs,vn,ws,wn,nx,ny,nz,pc,dt,dx,dy,dz)
        use, intrinsic :: iso_c_binding 
        implicit none
        real(C_DOUBLE),intent(in) :: dx,dy,dz
        integer,intent(in) :: nx,ny,nz
        integer :: i,j,k,km,kp,ip,im
        real(C_DOUBLE),intent(in) :: dt
        real(C_DOUBLE),intent(in) :: pc(nx,ny,nz)
        real(C_DOUBLE),intent(inout):: un(nx,ny+2,nz),vn(nx,ny+1,nz),wn(nx,ny+2,nz)
        real(C_DOUBLE),intent(in) :: us(nx,ny+2,nz),vs(nx,ny+1,nz),ws(nx,ny+2,nz)
        ! here we change the velocities to next time step
        do i = 1,nx
            do j = 2, ny-1
                do k = 1,nz
                    im = i - 1
                    if (im < 1 ) im = nx
                    un(i,j,k) = us(i,j,k)-dt*(pc(i,j,k)-pc(im,j,k))/dx
                end do 
            end do 
        end do
        do i = 1,nx
            do j = 2,ny
                do k = 1,nz
                    vn(i,j,k) = vs(i,j,k)-dt*(pc(i,j,k)-pc(i,j-1,k))/dy
                end do 
            end do 
        end do 

        do i = 1,nx
            do j = 2,ny-1
                do k=1 ,nz
                    km = k - 1
                    if (km < 1 ) km = nz
                    wn(i,j,k) = ws(i,j,k)-dt*(pc(i,j,k)-pc(i,j,km))/dz
                end do 
            end do
        end do 
        
    end subroutine n_step
    
    subroutine p_step(pn,pc,nx,ny,nz)
        use, intrinsic :: iso_c_binding
        implicit none
        integer :: i,j,k
        integer,intent(in) :: nx,ny,nz
        real(C_DOUBLE),intent(in) :: pc(nx,ny,nz)
        real(C_DOUBLE),intent(inout) :: pn(nx,ny,nz)
        do i = 1,nx
            do j = 1, ny
                do k = 1,nz
                    pn(i,j,k) = pn(i,j,k)+pc(i,j,k)
                end do 
            end do 
        end do 
    

    end subroutine p_step


    subroutine rhs_c (rhs,us,vs,ws,nx,ny,nz,dx,dy,dz,dt)
        use,intrinsic :: iso_c_binding
        implicit none
        integer,intent(in) :: nx,ny,nz
        integer :: i,j,k
        real(C_DOUBLE),intent(in) ::dx,dy,dz,dt 
        real(C_DOUBLE),intent(in):: us(nx,ny+2,nz),vs(nx,ny+1,nz),ws(nx,ny+2,nz)
        real(C_DOUBLE),intent(out)::rhs(nx,ny,nz)
        integer :: im,km
        ! here we take the divergence of the velocities
        do i = 1,nx
            do j = 1,ny
                do k = 1,nz

                im = i - 1
                km = k - 1
                if (im < 1 ) im = nx
                if (km < 1 ) km = nz

                rhs(i,j,k) = ((us(i,j,k)-us(im,j,k))/dx +&
                            (vs(i,j+1,k)-vs(i,j,k))/dy +&
                            (ws(i,j,k)-ws(i,j,km))/dz)&
                            /dt 
                end do 
            end do 
        end do 
    end subroutine rhs_c

end module step
