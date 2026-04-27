program test
    use :: apply_ibm_1st
    implicit none
    
    integer :: n_wave_x,nx,n_wave_z,nz,i,k
    real(C_DOUBLE) :: amp_x,lx,phase_x,amp_z,lz,phase_z
    real(C_DOUBLE),allocatable :: x(:),z(:),y_tot(:,:)
    real(C_DOUBLE) ::dx,dz,dy 
    
    

    n_wave_x = 1
    n_wave_z = 2
    nx = 10
    nz = 10
    allocate(x(nx),z(nz),y_tot(nx,nz))
    x(:) = 0.0d0
    z(:) = 0.0d0
    y_tot(:,:) = 0.0d0
    amp_x = 0.01d0
    amp_z = 0.01d0
    lx = 1.0d0
    lz = 1.0d0
    phase_x = 0.05d0
    phase_z = 0.05d0
    dx = lx/real(nx,kind=C_DOUBLE)
    dz = lz/real(nz,kind=C_DOUBLE)
    dy = 0.1
    do i = 1, nx
        x(i) = real(i,kind=C_DOUBLE)*dx 
    end do

    do k = 1, nz
        z(k) = real(k,kind=C_DOUBLE)*dz 
    end do

    call y_tot_sin(amp_x,n_wave_x,lx,phase_x,x,nx,amp_z,n_wave_z,lz,phase_z,z,nz,dy,y_tot)
    print*,x

end program test