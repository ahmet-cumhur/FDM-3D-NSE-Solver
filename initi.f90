module initi
    use,intrinsic :: iso_c_binding
    implicit none

    contains
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
        allocate(x(nx),y(ny),z(nz))
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
        do j = 1,ny
            y(j) = real(j-1,C_DOUBLE)* dy
        end do 
        do k = 1,nz
            z(k) = real(k-1, C_DOUBLE) * dz
        end do 

    end subroutine init_field
end module initi