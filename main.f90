program main
    use,intrinsic :: iso_c_binding
    use, intrinsic :: ieee_arithmetic

    use :: fftw_3d
    
    use :: initi
    use :: boundary
    use :: data_save
    !----------------------------------------!
    !if you want ibm then turn of these flags!                           
    !----------------------------------------!
    use :: apply_ibm_1st
    use :: step_ibm
    !use :: step

    implicit none

    integer :: i 
    real(C_DOUBLE) :: lx ,ly, lz 
    integer :: nx,ny,nz
    real(C_DOUBLE) :: re,dt,t_final,t_current
    real(C_DOUBLE) :: dx,dy,dz  
    real(C_DOUBLE), allocatable :: x(:),y(:),z(:)
    real(C_DOUBLE),allocatable :: un(:,:,:),vn(:,:,:),wn(:,:,:),us(:,:,:),vs(:,:,:),ws(:,:,:)
    real(C_DOUBLE),allocatable :: pn(:,:,:)
    real(C_DOUBLE),allocatable :: pc(:,:,:),rhs(:,:,:)
    real(C_DOUBLE),allocatable :: uc(:,:,:),vc(:,:,:),wc(:,:,:)
    character(len=256):: file_name


    !----------------------------------------!
    !ibm variables and arrays                !                           
    !----------------------------------------!

    ! they are initialized at initi file
    integer,allocatable :: mask_u(:,:,:),mask_v(:,:,:),mask_w(:,:,:)
    real(C_DOUBLE),allocatable :: y_wall_u(:,:),y_wall_w(:,:),y_wall_v(:,:)
    real(C_DOUBLE),allocatable :: x_u(:),y_u(:),z_u(:)
    real(C_DOUBLE),allocatable :: x_v(:),y_v(:),z_v(:)
    real(C_DOUBLE),allocatable :: x_w(:),y_w(:),z_w(:)
    real(C_DOUBLE) :: amp_x,phase_x,amp_z,phase_z
    integer :: n_wave_x,n_wave_z 


    !----------------------------------------!
    !ibm variables and arrays                !                           
    !----------------------------------------!

    
    print *, "initialzing fields and variables"
    ! initialize the variables
    call init_vars_ibm(lx,ly,lz,nx,ny,nz,dx,dy,dz,re,dt,t_final,t_current,n_wave_x,n_wave_z,amp_x,phase_x,amp_z,phase_z)
    print *, "variables initialized"
    ! initialize the fields
    call init_field_ibm(un,us,vn,vs,wn,ws,pn,pc,nx,ny,nz,x,y,z,dx,dy,dz,rhs,uc,vc,wc,mask_u,mask_v,mask_w,y_wall_u,y_wall_v,y_wall_w,&
        x_u,y_u,z_u,x_v,y_v,z_v,x_w,y_w,z_w)
    print *, "fields initialized" 
    
    !apply the initial condition
    call initi_c(un,us,vn,vs,wn,ws,nx,ny,nz)
    print *, "inital condition initialized"

    call get_masks(amp_x,amp_z,n_wave_x,n_wave_z,phase_x,phase_z,lx,lz,x_u,y_u,z_u,x_v,y_v,z_v,x_w,y_w,z_w,nx,ny,nz,dx,dy,dz,&
            y_wall_u,y_wall_v,y_wall_w,mask_u,mask_v,mask_w)
    call apply_ibm_vel(mask_u,mask_v,mask_w,un,us,vn,vs,wn,ws,nx,ny,nz)

    ! main time loop here
    print *, "main loop starting..."
    i = 0
    print*, i
    do while(t_current<t_final)
        t_current = t_current + dt
        i = i + 1
        ! check for problems 
        if (maxval(abs(un)) > 1.0d6 .or. &
            maxval(abs(vn)) > 1.0d6 .or. &
            maxval(abs(wn)) > 1.0d6 .or. &
            maxval(abs(pn)) > 1.0d6) then
            print *, "Velocity exploded at step", i
            print *, "un: ",maxval(abs(un)),"vn: ",maxval(abs(vn)),"wn: ",maxval(abs(wn)),"pn: ",maxval(abs(pn))
            stop
        end if
        !call moment(un,us,vn,vs,wn,ws,pn,re,dt,nx,ny,nz,dx,dy,dz)
        ! if we want to use ibm then use this one
        call apply_bc(un,us,vn,vs,wn,ws,nx,ny,nz)
        call moment_ibm_parallel(un,us,vn,vs,wn,ws,pn,re,dt,nx,ny,nz,dx,dy,dz)        
        call apply_bc(un,us,vn,vs,wn,ws,nx,ny,nz)
        call apply_ibm_vel(mask_u,mask_v,mask_w,un,us,vn,vs,wn,ws,nx,ny,nz)
        call rhs_c (rhs,us,vs,ws,nx,ny,nz,dx,dy,dz,dt)
        call poison_fft_3d (x,y,z,nx,ny,nz,lx,ly,lz,dx,dy,dz,rhs,pc)
        call n_step (us,un,vs,vn,ws,wn,nx,ny,nz,pc,dt,dx,dy,dz)        
        call p_step(pn,pc,nx,ny,nz)
        call apply_bc(un,us,vn,vs,wn,ws,nx,ny,nz)

        if (modulo(i,100) == 0)then
            print*, "current time step: ",i
            write(file_name,'("data_",I0,".vtk")') i
            call center_vel(un,vn,wn,uc,vc,wc,dx,dy,dz,nx,ny,nz)
            call data_output(uc,vc,wc,file_name,nx,ny,nz,dx,dy,dz,t_current)
        end if 
    end do 
    print *, "main loop ended..."
end program main
