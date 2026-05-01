# 1 "C:/Users/mehme/OneDrive/Desktop/thesis/nse/initi.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "C:/Users/mehme/OneDrive/Desktop/thesis/nse/initi.f90"
module initi
    use,intrinsic :: iso_c_binding
    implicit none

    contains
    ! this module contains variable/array initiation
    
    ! here we define variables 
    subroutine init_vars(lx,ly,lz,nx,ny,nz,dx,dy,dz,re,dt,t_final,t_current)
        use, intrinsic :: iso_c_binding
        implicit none
        real(C_DOUBLE),intent(out) ::lx,ly,lz,dx,dy,dz,re,dt,t_final,t_current 
        integer, intent(out) :: nx,ny,nz
        
        lx = 1.0d0
        ly = 1.0d0
        lz = 1.0d0
        
        nx = 32
        ny = 32
        nz = 32

        dx = lx/real(nx,C_DOUBLE)
        dy = ly/real(ny,C_DOUBLE)
        dz = lz/real(nz,C_DOUBLE)

        re = 100.0d0
        dt = 1e-3
        t_final = 1.0d0
        t_current = 0.0d0
    end subroutine init_vars
    ! here we define normal variables + ibm variables 
    subroutine init_vars_ibm(lx,ly,lz,nx,ny,nz,dx,dy,dz,re,dt,t_final,t_current,n_wave_x,n_wave_z,amp_x,phase_x,amp_z,phase_z)
        use, intrinsic :: iso_c_binding
        implicit none
        
        real(C_DOUBLE),intent(out) ::lx,ly,lz,dx,dy,dz,re,dt,t_final,t_current
        ! n_wave_,amp_,pahese_ is for getting the sin wave
        integer,intent(out) :: n_wave_x,n_wave_z 
        integer, intent(out) :: nx,ny,nz
        real(C_DOUBLE),intent(out) :: amp_x,phase_x,amp_z,phase_z
        
        print *, "entered init_vars_ibm"

        lx = 1.0d0
        ly = 1.0d0
        lz = 1.0d0
        
        nx = 100
        ny = 100
        nz = 100

        dx = lx/real(nx,C_DOUBLE)
        dy = ly/real(ny,C_DOUBLE)
        dz = lz/real(nz,C_DOUBLE)

        re = 100.0d0
        dt = 1e-3
        t_final = 2.0d0
        t_current = 0.0d0

        print *, "before ibm variables"
        ! ibm shape
        n_wave_x = 1
        n_wave_z = 1
        amp_x = 5.0d0*dy
        amp_z = 5.0d0*dy
        phase_x = 0.0d0
        phase_z = 0.0d0
        print *, "after ibm variables"
    end subroutine init_vars_ibm

    ! initiate the field 
    subroutine init_field(un,us,vn,vs,wn,ws,pn,pc,nx,ny,nz,x,y,z,dx,dy,dz,rhs,uc,vc,wc)
        use,intrinsic :: iso_c_binding
        implicit none
        integer :: i,j,k 
        integer, intent(in) :: nx,ny,nz
        real(C_DOUBLE), intent(in) :: dx,dy,dz
        real(C_DOUBLE),allocatable,intent(inout) :: un(:,:,:),vn(:,:,:),wn(:,:,:),us(:,:,:),vs(:,:,:),ws(:,:,:)
        real(C_DOUBLE),allocatable,intent(inout) :: pn(:,:,:),pc(:,:,:),x(:),y(:),z(:),rhs(:,:,:)
        real(C_DOUBLE),allocatable,intent(inout) :: uc(:,:,:),vc(:,:,:),wc(:,:,:)

        allocate(un(nx,ny+2,nz),us(nx,ny+2,nz),vn(nx,ny+1,nz),vs(nx,ny+1,nz))
        allocate(wn(nx,ny+2,nz),ws(nx,ny+2,nz),pn(nx,ny,nz),pc(nx,ny,nz),rhs(nx,ny,nz))
        allocate(x(nx),y(ny+1),z(nz))
        allocate(uc(nx,ny,nz),vc(nx,ny,nz),wc(nx,ny,nz))
        un(:,:,:) = 0.0d0
        us(:,:,:) = 0.0d0
        vn(:,:,:) = 0.0d0
        vs(:,:,:) = 0.0d0
        wn(:,:,:) = 0.0d0
        ws(:,:,:) = 0.0d0
        pc(:,:,:) = 0.0d0
        pn(:,:,:) = 0.0d0
        rhs(:,:,:) = 0.0d0
        uc(:,:,:) = 0.0d0
        vc(:,:,:) = 0.0d0
        wc(:,:,:) = 0.0d0

        do i = 1,nx
            x(i) = real(i-1, C_DOUBLE) * dx
        end do 
        do j = 1,ny+1
            y(j) = real(j-1,C_DOUBLE)* dy
        end do 
        do k = 1,nz
            z(k) = real(k-1, C_DOUBLE) * dz
        end do 

    end subroutine init_field

    ! initiate the field w/ ibm variables
    ! x_u/y/z, .... mask_u/y/z etc. 
    subroutine init_field_ibm(un,us,vn,vs,wn,ws,pn,pc,nx,ny,nz,x,y,z,dx,dy,dz,rhs,uc,vc,wc,mask_u,mask_v,mask_w,y_wall_u,y_wall_v,y_wall_w,&
        x_u,y_u,z_u,x_v,y_v,z_v,x_w,y_w,z_w)
        use,intrinsic :: iso_c_binding
        implicit none
        integer :: i,j,k 
        integer, intent(in) :: nx,ny,nz
        real(C_DOUBLE), intent(in) :: dx,dy,dz
        real(C_DOUBLE),allocatable,intent(inout) :: un(:,:,:),vn(:,:,:),wn(:,:,:),us(:,:,:),vs(:,:,:),ws(:,:,:)
        real(C_DOUBLE),allocatable,intent(inout) :: pn(:,:,:),pc(:,:,:),x(:),y(:),z(:),rhs(:,:,:)
        real(C_DOUBLE),allocatable,intent(inout) :: uc(:,:,:),vc(:,:,:),wc(:,:,:)
        ! ibm variables, arrays are here
        integer,allocatable,intent(inout) :: mask_u(:,:,:),mask_v(:,:,:),mask_w(:,:,:)
        real(C_DOUBLE),allocatable,intent(inout) :: y_wall_u(:,:),y_wall_w(:,:),y_wall_v(:,:)
        real(C_DOUBLE),allocatable,intent(inout) :: x_u(:),y_u(:),z_u(:)
        real(C_DOUBLE),allocatable,intent(inout) :: x_v(:),y_v(:),z_v(:)
        real(C_DOUBLE),allocatable,intent(inout) :: x_w(:),y_w(:),z_w(:)

        allocate(un(nx,ny+2,nz),us(nx,ny+2,nz),vn(nx,ny+1,nz),vs(nx,ny+1,nz))
        allocate(wn(nx,ny+2,nz),ws(nx,ny+2,nz),pn(nx,ny,nz),pc(nx,ny,nz),rhs(nx,ny,nz))
        allocate(x(nx),y(ny+1),z(nz))
        allocate(uc(nx,ny,nz),vc(nx,ny,nz),wc(nx,ny,nz))
        ! ibm variables, arrays are here
        allocate(mask_u(nx,ny+2,nz),mask_v(nx,ny+1,nz),mask_w(nx,ny+2,nz))
        allocate(y_wall_u(nx,nz),y_wall_v(nx,nz),y_wall_w(nx,nz))
        allocate(x_u(nx),y_u(ny+2),z_u(nz))
        allocate(x_v(nx),y_v(ny+1),z_v(nz))
        allocate(x_w(nx),y_w(ny+2),z_w(nz))
        
        un(:,:,:) = 0.0d0
        us(:,:,:) = 0.0d0
        vn(:,:,:) = 0.0d0
        vs(:,:,:) = 0.0d0
        wn(:,:,:) = 0.0d0
        ws(:,:,:) = 0.0d0
        pc(:,:,:) = 0.0d0
        pn(:,:,:) = 0.0d0
        rhs(:,:,:) = 0.0d0
        uc(:,:,:) = 0.0d0
        vc(:,:,:) = 0.0d0
        wc(:,:,:) = 0.0d0
        !ibm stuff
        mask_u(:,:,:) = 0
        mask_v(:,:,:) = 0
        mask_w(:,:,:) = 0

        y_wall_u(:,:) = 0.0d0
        y_wall_w(:,:) = 0.0d0
        y_wall_v(:,:) = 0.0d0

        x_u(:) = 0.0d0
        y_u(:) = 0.0d0
        z_u(:) = 0.0d0

        x_v(:) = 0.0d0
        y_v(:) = 0.0d0
        z_v(:) = 0.0d0
        
        x_w(:) = 0.0d0
        y_w(:) = 0.0d0
        z_w(:) = 0.0d0
        ! we dont use these 
        do i = 1,nx
            x(i) = real(i-1, C_DOUBLE) * dx
        end do 
        do j = 1,ny+1
            y(j) = real(j-1,C_DOUBLE)* dy
        end do 
        do k = 1,nz
            z(k) = real(k-1, C_DOUBLE) * dz
        end do 

    end subroutine init_field_ibm
end module initi
