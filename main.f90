program main
    use,intrinsic :: iso_c_binding
    
    use :: fftw_3d
    use :: step
    use :: initi
    use :: boundary
    use :: data_save

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

    ! initialize the variables
    call init_vars(lx,ly,lz,nx,ny,nz,dx,dy,dz,re,dt,t_final,t_current )
    ! initialize the fields
    call init_field(un,us,vn,vs,wn,ws,pn,pc,nx,ny,nz,x,y,z,dx,dy,dz,rhs,uc,vc,wc)
    
    !apply the initial condition
    call initi_c(un,us,vn,vs,wn,ws,nx,ny,nz)

    ! main time loop here
    i = 0
    do while(t_current<t_final)
        t_current = t_current + dt
        i = i + 1
        call apply_bc(un,us,vn,vs,wn,ws,nx,ny,nz)
        call moment(un,us,vn,vs,wn,ws,pn,re,dt,nx,ny,nz,dx,dy,dz)
        call apply_bc(un,us,vn,vs,wn,ws,nx,ny,nz)
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
end program main
