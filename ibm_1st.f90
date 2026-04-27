module apply_ibm_1st
    use, intrinsic :: iso_c_binding
    implicit none
    ! we apply a sinx shape here as staircase approximation 
    contains
        ! we get the shape of the wall 
        subroutine y_tot_sin(amp_x,n_wave_x,lx,phase_x,x,nx,amp_z,n_wave_z,lz,phase_z,z,nz,dy,y_wall)
            use,intrinsic :: iso_c_binding
            implicit none
            integer,intent(in) :: n_wave_x,nx,n_wave_z,nz
            real(C_DOUBLE),intent(in) :: amp_x,lx,phase_x,amp_z,lz,phase_z
            real(C_DOUBLE),intent(in) :: x(nx),z(nz),dy

            integer :: i,k
            real(C_DOUBLE),parameter :: pi = 3.141592653589793d0
            real(C_DOUBLE),intent(inout) :: y_wall(nx,nz) ! this might be problematic if nx and nz has diff. sizes
            real(C_DOUBLE) :: y0 ! some padding btwn lowest value of sine wave and y-boundary  
            
            
            ! and since we can have negative values we need to apply some 
            ! padding between y(1) and the lowest value here
            y0 = abs(amp_x) +abs(amp_z)+5*dy + 0.1d0! padding is here
            do i = 1,nx
                do k = 1,nz
                    y_wall(i,k) = y0 +amp_x*sin(2.0d0*pi*real(n_wave_x,kind=C_DOUBLE)*x(i)/lx+phase_x)+&
                    amp_z*sin(2.0d0*real(n_wave_z,kind=C_DOUBLE)*pi*z(k)/lz+phase_z)
                end do 
            end do
      
        end subroutine y_tot_sin

        subroutine get_vel_loc_u(x_u,y_u,z_u,dx,dy,dz,nx,ny,nz)
            use, intrinsic :: iso_c_binding
            implicit none
            integer,intent(in) :: nx,ny,nz
            real(C_DOUBLE),intent(in) :: dx,dy,dz
            real(C_DOUBLE),intent(inout):: x_u(nx),y_u(ny+2),z_u(nz)
            integer :: i,j,k
            do i=1,nx
                x_u(i) = real(i-1.0d0,kind=C_DOUBLE)*dx
            end do

            do j = 1, ny+2
                y_u(j) = real(j-1.5d0,kind=C_DOUBLE)*dy
            end do 

            do k = 1,nz
                z_u(k) = real(k-0.5d0,kind=C_DOUBLE)*dz
            end do 

        end subroutine get_vel_loc_u

        ! now we need to get the locations of the velocities and pressure
        subroutine get_mask_u(amp_x,amp_z,n_wave_x,n_wave_z,phase_x,phase_z,lx,lz,x_u,y_u,z_u,nx,ny,nz,dx,dy,dz,y_wall_u,mask_u)
            
            use,intrinsic :: iso_c_binding
            implicit none
            integer,intent(in) :: nx,ny,nz
            integer,intent(in) :: n_wave_x,n_wave_z
            real(C_DOUBLE),intent(in) :: amp_x,lx,phase_x,amp_z,lz,phase_z
            real(C_DOUBLE),intent(in) :: dy,dx,dz
            real(C_DOUBLE),intent(inout) :: x_u(nx),y_u(ny+2),z_u(nz)
            integer,intent(inout) :: mask_u(nx,ny+2,nz) 

            integer :: i,j,k
            real(C_DOUBLE),parameter :: pi = 3.141592653589793d0
            real(C_DOUBLE),intent(inout) :: y_wall_u(nx,nz) ! this might be problematic if nx and nz has diff. sizes
            real(C_DOUBLE) :: y0 ! some padding btwn lowest value of sine wave and y-boundary  
            
            
            ! and since we can have negative values we need to apply some 
            ! padding between y(1) and the lowest value here
            y0 = abs(amp_x) +abs(amp_z)+2.0d0*dy! padding is here


            do i = 1,nx
                do k = 1,nz
                    ! we changed the order of iterations cause we dont want to re-calculate the
                    ! height of the boundary every j iteration normal order is i-j-k
                    y_wall_u(i,k) = y0 +amp_x*sin(2.0d0*pi*real(n_wave_x,kind=C_DOUBLE)*x_u(i)/lx+phase_x)+&
                        amp_z*sin(2.0d0*real(n_wave_z,kind=C_DOUBLE)*pi*z_u(k)/lz+phase_z)
                    do j = 1,ny+2                         
                        ! this means the wall value is bigger than y_u
                        ! it means the velocity lies inside the solid
                        ! 1 means solid 
                        ! 0 means fluid
                        if (y_wall_u(i,k)>y_u(j)) then 
                            mask_u(i,j,k) = 1
                        else 
                            ! this means its in fluid
                            mask_u(i,j,k) = 0   
                        end if 


                    end do 
                end do 
            end do
      
        end subroutine get_mask_u

        subroutine get_vel_loc_v(x_v,y_v,z_v,dx,dy,dz,nx,ny,nz)
            use, intrinsic :: iso_c_binding
            implicit none
            integer,intent(in) :: nx,ny,nz
            real(C_DOUBLE),intent(in) :: dx,dy,dz
            real(C_DOUBLE),intent(inout):: x_v(nx),y_v(ny+1),z_v(nz)
            integer :: i,j,k
            do i=1,nx
                x_v(i) = real(i-0.5d0,kind=C_DOUBLE)*dx
            end do

            do j = 1, ny+1
                y_v(j) = real(j-1.0d0,kind=C_DOUBLE)*dy
            end do 

            do k = 1,nz
                z_v(k) = real(k-0.5d0,kind=C_DOUBLE)*dz
            end do 

        end subroutine get_vel_loc_v

        ! now we need to get the locations of the velocities and pressure
        subroutine get_mask_v(amp_x,amp_z,n_wave_x,n_wave_z,phase_x,phase_z,lx,lz,x_v,y_v,z_v,nx,ny,nz,dx,dy,dz,y_wall_v,mask_v)
            
            use,intrinsic :: iso_c_binding
            implicit none
            integer,intent(in) :: nx,ny,nz
            integer,intent(in) :: n_wave_x,n_wave_z
            real(C_DOUBLE),intent(in) :: amp_x,lx,phase_x,amp_z,lz,phase_z
            real(C_DOUBLE),intent(in) :: dy,dx,dz
            real(C_DOUBLE),intent(inout) :: x_v(nx),y_v(ny+1),z_v(nz)
            integer,intent(inout) :: mask_v(nx,ny+1,nz) 

            integer :: i,j,k
            real(C_DOUBLE),parameter :: pi = 3.141592653589793d0
            real(C_DOUBLE),intent(inout) :: y_wall_v(nx,nz) ! this might be problematic if nx and nz has diff. sizes
            real(C_DOUBLE) :: y0 ! some padding btwn lowest value of sine wave and y-boundary  
            
            
            ! and since we can have negative values we need to apply some 
            ! padding between y(1) and the lowest value here
            y0 = abs(amp_x) +abs(amp_z)+2.0d0*dy! padding is here


            do i = 1,nx
                do k = 1,nz
                    ! we changed the order of iterations cause we dont want to re-calculate the
                    ! height of the boundary every j iteration normal order is i-j-k
                    y_wall_v(i,k) = y0 +amp_x*sin(2.0d0*pi*real(n_wave_x,kind=C_DOUBLE)*x_v(i)/lx+phase_x)+&
                        amp_z*sin(2.0d0*real(n_wave_z,kind=C_DOUBLE)*pi*z_v(k)/lz+phase_z)
                    do j = 1,ny+1                         
                        ! this means the wall value is bigger than y_u
                        ! it means the velocity lies inside the solid
                        ! 1 means solid 
                        ! 0 means fluid
                        if (y_wall_v(i,k)>y_v(j)) then 
                            mask_v(i,j,k) = 1
                        else 
                            ! this means its in fluid
                            mask_v(i,j,k) = 0   
                        end if 


                    end do 
                end do 
            end do
      
        end subroutine get_mask_v

        subroutine get_vel_loc_w(x_w,y_w,z_w,dx,dy,dz,nx,ny,nz)
            use, intrinsic :: iso_c_binding
            implicit none
            integer,intent(in) :: nx,ny,nz
            real(C_DOUBLE),intent(in) :: dx,dy,dz
            real(C_DOUBLE),intent(inout):: x_w(nx),y_w(ny+2),z_w(nz)
            integer :: i,j,k
            do i=1,nx
                x_w(i) = real(i-0.5d0,kind=C_DOUBLE)*dx
            end do

            do j = 1, ny+2
                y_w(j) = real(j-1.5d0,kind=C_DOUBLE)*dy
            end do

            do k = 1,nz
                z_w(k) = real(k-1.0d0,kind=C_DOUBLE)*dz
            end do 

        end subroutine get_vel_loc_w

        ! now we need to get the locations of the velocities and pressure
        subroutine get_mask_w(amp_x,amp_z,n_wave_x,n_wave_z,phase_x,phase_z,lx,lz,x_w,y_w,z_w,nx,ny,nz,dx,dy,dz,y_wall_w,mask_w)
            
            use,intrinsic :: iso_c_binding
            implicit none
            integer,intent(in) :: nx,ny,nz
            integer,intent(in) :: n_wave_x,n_wave_z
            real(C_DOUBLE),intent(in) :: amp_x,lx,phase_x,amp_z,lz,phase_z
            real(C_DOUBLE),intent(in) :: dy,dx,dz
            real(C_DOUBLE),intent(inout) :: x_w(nx),y_w(ny+2),z_w(nz)
            integer,intent(inout) :: mask_w(nx,ny+2,nz) 

            integer :: i,j,k
            real(C_DOUBLE),parameter :: pi = 3.141592653589793d0
            real(C_DOUBLE),intent(inout) :: y_wall_w(nx,nz) ! this might be problematic if nx and nz has diff. sizes
            real(C_DOUBLE) :: y0 ! some padding btwn lowest value of sine wave and y-boundary  
            
            
            ! and since we can have negative values we need to apply some 
            ! padding between y(1) and the lowest value here
            y0 = abs(amp_x) +abs(amp_z)+2.0d0*dy! padding is here


            do i = 1,nx
                do k = 1,nz
                    ! we changed the order of iterations cause we dont want to re-calculate the
                    ! height of the boundary every j iteration normal order is i-j-k
                    y_wall_w(i,k) = y0 +amp_x*sin(2.0d0*pi*real(n_wave_x,kind=C_DOUBLE)*x_w(i)/lx+phase_x)+&
                        amp_z*sin(2.0d0*real(n_wave_z,kind=C_DOUBLE)*pi*z_w(k)/lz+phase_z)
                    do j = 1,ny+2                         
                        ! this means the wall value is bigger than y_u
                        ! it means the velocity lies inside the solid
                        ! 1 means solid 
                        ! 0 means fluid
                        if (y_wall_w(i,k)>y_w(j)) then 
                            mask_w(i,j,k) = 1
                        else 
                            ! this means its in fluid
                            mask_w(i,j,k) = 0   
                        end if 


                    end do 
                end do 
            end do
      
        end subroutine get_mask_w

        ! we combine all the mask subroutines in one 
        subroutine get_masks(amp_x,amp_z,n_wave_x,n_wave_z,phase_x,phase_z,lx,lz,x_u,y_u,z_u,x_v,y_v,z_v,x_w,y_w,z_w,nx,ny,nz,dx,dy,dz,&
            y_wall_u,y_wall_v,y_wall_w,mask_u,mask_v,mask_w)
            use,intrinsic:: iso_c_binding
            implicit none

            integer,intent(in) :: nx,ny,nz
            integer,intent(in) :: n_wave_x,n_wave_z
            real(C_DOUBLE),intent(in) :: dy,dx,dz
            real(C_DOUBLE),intent(in) :: amp_x,lx,phase_x,amp_z,lz,phase_z
            real(C_DOUBLE),intent(inout) :: x_u(nx),y_u(ny+2),z_u(nz)
            real(C_DOUBLE),intent(inout) :: x_v(nx),y_v(ny+1),z_v(nz)
            real(C_DOUBLE),intent(inout) :: x_w(nx),y_w(ny+2),z_w(nz)
            real(C_DOUBLE),intent(inout) :: y_wall_u(nx,nz)
            real(C_DOUBLE),intent(inout) :: y_wall_v(nx,nz)
            real(C_DOUBLE),intent(inout) :: y_wall_w(nx,nz)
            integer,intent(inout) :: mask_u(nx,ny+2,nz) 
            integer,intent(inout) :: mask_v(nx,ny+1,nz) 
            integer,intent(inout) :: mask_w(nx,ny+2,nz) 

            call get_vel_loc_u(x_u,y_u,z_u,dx,dy,dz,nx,ny,nz)
            call get_vel_loc_v(x_v,y_v,z_v,dx,dy,dz,nx,ny,nz)
            call get_vel_loc_w(x_w,y_w,z_w,dx,dy,dz,nx,ny,nz)
            call get_mask_u(amp_x,amp_z,n_wave_x,n_wave_z,phase_x,phase_z,lx,lz,x_u,y_u,z_u,nx,ny,nz,dx,dy,dz,y_wall_u,mask_u)
            call get_mask_v(amp_x,amp_z,n_wave_x,n_wave_z,phase_x,phase_z,lx,lz,x_v,y_v,z_v,nx,ny,nz,dx,dy,dz,y_wall_v,mask_v)
            call get_mask_w(amp_x,amp_z,n_wave_x,n_wave_z,phase_x,phase_z,lx,lz,x_w,y_w,z_w,nx,ny,nz,dx,dy,dz,y_wall_w,mask_w)

        end subroutine get_masks

    
end module apply_ibm_1st