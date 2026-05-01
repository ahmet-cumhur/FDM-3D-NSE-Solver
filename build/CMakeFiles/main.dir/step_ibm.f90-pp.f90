# 1 "C:/Users/mehme/OneDrive/Desktop/thesis/nse/step_ibm.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "C:/Users/mehme/OneDrive/Desktop/thesis/nse/step_ibm.f90"
module step_ibm
    implicit none 
    
    contains
    ! this is the same as the old momentum part
    ! but i wanted a new one incase of problems
    subroutine moment_ibm(un,us,vn,vs,wn,ws,pn,re,dt,nx,ny,nz,dx,dy,dz,mask_u,mask_v,mask_w)
        use, intrinsic :: iso_c_binding
        use, intrinsic :: ieee_arithmetic
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
        ! i added the ghost boundaries too so i dont need to deal
        ! w/ their indices in the momentum loop
        integer,intent(in):: mask_u(nx,ny+2,nz)
        integer,intent(in):: mask_v(nx,ny+1,nz)
        integer,intent(in):: mask_w(nx,ny+2,nz)

        integer :: ip,im,kp,km,jp
        ! x axis
        do i=1, nx
            do j= 2,ny+1
                do k=1,nz
                    ! this is the body check
                    ! i deactived these cheks inside moment
                    ! since we already zero them in other function
                    ! @ apply_ibm_vel

                    !if(mask_u(i,j,k) == 1) then
                    !    us(i,j,k) = 0.0d0
                    !    un(i,j,k) = 0.0d0
                    !    cycle  
                    !end if

                    ! we apply the periodic bc 
                    ip = i + 1
                    im = i - 1
                    kp = k + 1
                    km = k - 1
                    jp = j-1
                    if (ip > nx) ip = 1
                    if (im < 1 ) im = nx

                    if (kp > nz) kp = 1
                    if (km < 1 ) km = nz


                    uu_p = 0.25d0*(un(ip,j,k)+un(i,j,k))*(un(ip,j,k)+un(i,j,k))
                    uu_m = 0.25d0*(un(i,j,k)+un(im,j,k))*(un(i,j,k)+un(im,j,k))
                    
                    uv_p = 0.25d0*(un(i,j+1,k)+un(i,j,k))*(vn(i,j,k)+vn(im,j,k))
                    uv_m = 0.25d0*(un(i,j,k)+un(i,j-1,k))*(vn(i,j-1,k)+vn(im,j-1,k))! might change
                    
                    uw_p = 0.25d0*(un(i,j,k)+un(i,j,kp))*(wn(i,j,kp)+wn(im,j,kp))
                    uw_m = 0.25d0*(un(i,j,k)+un(i,j,km))*(wn(i,j,k)+wn(im,j,k))

                    diff_ux = (un(im,j,k)-2.0d0*un(i,j,k)+un(ip,j,k))/dx**2
                    diff_uy = (un(i,j-1,k)-2.0d0*un(i,j,k)+un(i,j+1,k))/dy**2
                    diff_uz = (un(i,j,km)-2.0d0*un(i,j,k)+un(i,j,kp))/dz**2
                    ! we use jp here
                    dpx = (pn(i,jp,k)-pn(im,jp,k))/dx

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
            do j = 2, ny
                do k =1,nz
                    ! this is the body check
                    ! i deactived these cheks inside moment
                    ! since we already zero them in other function
                    ! @ apply_ibm_vel
                    
                    !if(mask_v(i,j,k) == 1) then
                    !    vs(i,j,k) = 0.0d0
                    !    vn(i,j,k) = 0.0d0
                    !    cycle  
                    !end if

                    ip = i + 1
                    im = i - 1
                    kp = k + 1
                    km = k - 1

                    if (ip > nx) ip = 1
                    if (im < 1 ) im = nx

                    if (kp > nz) kp = 1
                    if (km < 1 ) km = nz


                    vu_p = 0.25d0*(vn(i,j,k)+vn(ip,j,k))*(un(ip,j,k)+un(ip,j+1,k))
                    vu_m = 0.25d0*(vn(i,j,k)+vn(im,j,k))*(un(i,j,k)+un(i,j+1,k))
                    
                    vv_p = 0.25d0*(vn(i,j,k)+vn(i,j+1,k))*(vn(i,j,k)+vn(i,j+1,k))
                    vv_m = 0.25d0*(vn(i,j,k)+vn(i,j-1,k))*(vn(i,j,k)+vn(i,j-1,k))
                    
                    vw_p = 0.25d0*(vn(i,j,kp)+vn(i,j,k))*(wn(i,j,kp)+wn(i,j+1,kp))
                    vw_m = 0.25d0*(vn(i,j,km)+vn(i,j,k))*(wn(i,j,k)+wn(i,j+1,k))

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
                    ! this is the body check
                    ! i deactived these cheks inside moment
                    ! since we already zero them in other function
                    ! @ apply_ibm_vel

                    !if(mask_w(i,j,k) == 1) then
                    !    ws(i,j,k) = 0.0d0
                    !    wn(i,j,k) = 0.0d0
                    !    cycle  
                    !end if

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
                    
                    wv_p = 0.25d0*(wn(i,j,k)+wn(i,j+1,k))*(vn(i,j,k)+vn(i,j,km))
                    wv_m = 0.25d0*(wn(i,j,k)+wn(i,j-1,k))*(vn(i,j-1,k)+vn(i,j-1,km)) 
                    
                    diff_wx = (wn(im,j,k)-2.0d0*wn(i,j,k)+wn(ip,j,k))/dx**2
                    diff_wy = (wn(i,j-1,k)-2.0d0*wn(i,j,k)+wn(i,j+1,k))/dy**2
                    diff_wz = (wn(i,j,km)-2.0d0*wn(i,j,k)+wn(i,j,kp))/dz**2                     

                    jp = j-1
                    dpz = (pn(i,jp,k)-pn(i,jp,km))/dz 
                    
                    ws(i,j,k) = wn(i,j,k)+dt*(&
                                -(wu_p-wu_m)/dx+&
                                -(wv_p-wv_m)/dy+&
                                -(ww_p-ww_m)/dz+&
                                -dpz+&
                                (1.0d0/re)*(diff_wx+diff_wy+diff_wz))
                    
                end do 
            end do 
        end do


    end subroutine moment_ibm

    subroutine moment_ibm_parallel(un,us,vn,vs,wn,ws,pn,re,dt,nx,ny,nz,dx,dy,dz)
        use, intrinsic :: iso_c_binding
        use :: omp_lib
        implicit none
        integer,parameter :: NT = 6

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
        integer :: ip,im,kp,km,jp
        ! x axis
        ! $omp parallel default(private) shared(un,us,vn,vs,wn,ws,pn,re,dt,nx,ny,nz,dx,dy,dz) 
        ! $omp num_threads(NT) 
        ! $omp do
        do k=1, nz
            do j= 2,ny+1
                do i=1,nx
                    ! we apply the periodic bc 
                    ip = i + 1
                    im = i - 1
                    kp = k + 1
                    km = k - 1
                    jp = j-1
                    if (ip > nx) ip = 1
                    if (im < 1 ) im = nx

                    if (kp > nz) kp = 1
                    if (km < 1 ) km = nz


                    uu_p = 0.25d0*(un(ip,j,k)+un(i,j,k))*(un(ip,j,k)+un(i,j,k))
                    uu_m = 0.25d0*(un(i,j,k)+un(im,j,k))*(un(i,j,k)+un(im,j,k))
                    
                    uv_p = 0.25d0*(un(i,j+1,k)+un(i,j,k))*(vn(i,j,k)+vn(im,j,k))
                    uv_m = 0.25d0*(un(i,j,k)+un(i,j-1,k))*(vn(i,j-1,k)+vn(im,j-1,k))! might change
                    
                    uw_p = 0.25d0*(un(i,j,k)+un(i,j,kp))*(wn(i,j,kp)+wn(im,j,kp))
                    uw_m = 0.25d0*(un(i,j,k)+un(i,j,km))*(wn(i,j,k)+wn(im,j,k))

                    diff_ux = (un(im,j,k)-2.0d0*un(i,j,k)+un(ip,j,k))/dx**2
                    diff_uy = (un(i,j-1,k)-2.0d0*un(i,j,k)+un(i,j+1,k))/dy**2
                    diff_uz = (un(i,j,km)-2.0d0*un(i,j,k)+un(i,j,kp))/dz**2
                    ! we use jp here
                    dpx = (pn(i,jp,k)-pn(im,jp,k))/dx

                    us(i,j,k) = un(i,j,k)+dt*(&
                                -(uu_p-uu_m)/dx+&
                                -(uv_p-uv_m)/dy+&
                                -(uw_p-uw_m)/dz+&
                                -dpx+&
                                (1.0d0/re)*(diff_ux + diff_uy + diff_uz))

                end do 
            end do 
        end do
        ! $omp end do  
        
        ! y axis
        ! $omp do
        do k = 1,nz
            do j = 2, ny
                do i =1,nx

                    ip = i + 1
                    im = i - 1
                    kp = k + 1
                    km = k - 1

                    if (ip > nx) ip = 1
                    if (im < 1 ) im = nx

                    if (kp > nz) kp = 1
                    if (km < 1 ) km = nz


                    vu_p = 0.25d0*(vn(i,j,k)+vn(ip,j,k))*(un(ip,j,k)+un(ip,j+1,k))
                    vu_m = 0.25d0*(vn(i,j,k)+vn(im,j,k))*(un(i,j,k)+un(i,j+1,k))
                    
                    vv_p = 0.25d0*(vn(i,j,k)+vn(i,j+1,k))*(vn(i,j,k)+vn(i,j+1,k))
                    vv_m = 0.25d0*(vn(i,j,k)+vn(i,j-1,k))*(vn(i,j,k)+vn(i,j-1,k))
                    
                    vw_p = 0.25d0*(vn(i,j,kp)+vn(i,j,k))*(wn(i,j,kp)+wn(i,j+1,kp))
                    vw_m = 0.25d0*(vn(i,j,km)+vn(i,j,k))*(wn(i,j,k)+wn(i,j+1,k))

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
        ! $omp end do

        ! z axis
        ! $omp do
        do k = 1,nz
            do j = 2,ny+1
                do i = 1,nx

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
                    
                    wv_p = 0.25d0*(wn(i,j,k)+wn(i,j+1,k))*(vn(i,j,k)+vn(i,j,km))
                    wv_m = 0.25d0*(wn(i,j,k)+wn(i,j-1,k))*(vn(i,j-1,k)+vn(i,j-1,km)) 
                    
                    diff_wx = (wn(im,j,k)-2.0d0*wn(i,j,k)+wn(ip,j,k))/dx**2
                    diff_wy = (wn(i,j-1,k)-2.0d0*wn(i,j,k)+wn(i,j+1,k))/dy**2
                    diff_wz = (wn(i,j,km)-2.0d0*wn(i,j,k)+wn(i,j,kp))/dz**2                     

                    jp = j-1
                    dpz = (pn(i,jp,k)-pn(i,jp,km))/dz 
                    
                    ws(i,j,k) = wn(i,j,k)+dt*(&
                                -(wu_p-wu_m)/dx+&
                                -(wv_p-wv_m)/dy+&
                                -(ww_p-ww_m)/dz+&
                                -dpz+&
                                (1.0d0/re)*(diff_wx+diff_wy+diff_wz))
                    
                end do 
            end do 
        end do
        ! $omp end do


    end subroutine moment_ibm_parallel

    subroutine n_step (us,un,vs,vn,ws,wn,nx,ny,nz,pc,dt,dx,dy,dz)
        use, intrinsic :: iso_c_binding 
        implicit none
        real(C_DOUBLE),intent(in) :: dx,dy,dz
        integer,intent(in) :: nx,ny,nz
        integer :: i,j,k,km,kp,ip,im,jp
        real(C_DOUBLE),intent(in) :: dt
        real(C_DOUBLE),intent(in) :: pc(nx,ny,nz)
        real(C_DOUBLE),intent(inout):: un(nx,ny+2,nz),vn(nx,ny+1,nz),wn(nx,ny+2,nz)
        real(C_DOUBLE),intent(in) :: us(nx,ny+2,nz),vs(nx,ny+1,nz),ws(nx,ny+2,nz)
        ! here we change the velocities to next time step
        do i = 1,nx
            do j = 2, ny+1
                do k = 1,nz
                    im = i - 1
                    if (im < 1 ) im = nx
                    jp = j-1
                    un(i,j,k) = us(i,j,k)-dt*(pc(i,jp,k)-pc(im,jp,k))/dx
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
            do j = 2,ny+1
                do k=1 ,nz
                    km = k - 1
                    if (km < 1 ) km = nz
                    jp = j-1
                    wn(i,j,k) = ws(i,j,k)-dt*(pc(i,jp,k)-pc(i,jp,km))/dz
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
        integer :: im,km,ip,kp
        ! here we take the divergence of the velocities
        do i = 1,nx
            do j = 1,ny
                do k = 1,nz
                ip = i + 1 
                kp = k + 1
                im = i - 1
                km = k - 1
                if (ip > nx) ip = 1
                if (kp > nz) kp = 1
                if (im < 1 ) im = nx
                if (km < 1 ) km = nz

                rhs(i,j,k) = ((us(ip,j+1,k)-us(i,j+1,k))/dx +&
                            (vs(i,j+1,k)-vs(i,j,k))/dy +&
                            (ws(i,j+1,kp)-ws(i,j+1,k))/dz)&
                            /dt 
                end do 
            end do 
        end do 
    end subroutine rhs_c
    ! here we apply 0 velocities again
    subroutine apply_ibm_vel(mask_u,mask_v,mask_w,un,us,vn,vs,wn,ws,nx,ny,nz)
        use,intrinsic :: iso_c_binding
        implicit none
        integer,intent(in) :: nx,ny,nz
        real(C_DOUBLE),intent(inout) :: un(nx,ny+2,nz),us(nx,ny+2,nz)
        real(C_DOUBLE),intent(inout) :: vn(nx,ny+1,nz),vs(nx,ny+1,nz)
        real(C_DOUBLE),intent(inout) :: wn(nx,ny+2,nz),ws(nx,ny+2,nz)
        
        integer,intent(in) :: mask_u(nx,ny+2,nz)
        integer,intent(in) :: mask_v(nx,ny+1,nz)
        integer,intent(in) :: mask_w(nx,ny+2,nz)

        integer :: i,j,k
        do i = 1, nx
            do j = 2,ny+1
                do k = 1,nz
                    if(mask_u(i,j,k)== 1)then
                        un(i,j,k) = 0.0d0
                        us(i,j,k) = 0.0d0
                    end if
                end do 
            end do 
        end do 
        do i = 1, nx
            do j = 1,ny+1
                do k = 1,nz
                    if(mask_v(i,j,k)== 1)then
                        vn(i,j,k) = 0.0d0
                        vs(i,j,k) = 0.0d0
                    end if
                end do 
            end do 
        end do 
        do i = 1, nx
            do j = 2,ny+1
                do k = 1,nz
                    if(mask_w(i,j,k)== 1)then
                        wn(i,j,k) = 0.0d0
                        ws(i,j,k) = 0.0d0
                    end if
                end do 
            end do 
        end do 

        end subroutine apply_ibm_vel



end module step_ibm

