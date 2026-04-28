module boundary
    implicit none 
    contains 
    
    subroutine apply_bc(un,us,vn,vs,wn,ws,nx,ny,nz)
        use, intrinsic :: iso_c_binding
        implicit none
        integer,intent(in) :: nx,ny,nz
        real(C_DOUBLE),intent(inout) :: un(nx,ny+2,nz),vn(nx,ny+1,nz),wn(nx,ny+2,nz),us(nx,ny+2,nz),vs(nx,ny+1,nz),ws(nx,ny+2,nz)
        
        vn(:,ny+1,:) = 0.0d0
        vn(:,1,:) = 0.0d0
        vs(:,ny+1,:) = 0.0d0
        vs(:,1,:) = 0.0d0 
        
        un(:,ny+2,:) = -un(:,ny+1,:)
        us(:,ny+2,:) = -us(:,ny+1,:) 
        un(:,1,:) = -un(:,2,:)
        us(:,1,:) = -us(:,2,:)

        wn(:,ny+2,:) = -wn(:,ny+1,:)
        ws(:,ny+2,:) = -ws(:,ny+1,:)
        wn(:,1,:) = -wn(:,2,:)
        ws(:,1,:) = -ws(:,2,:)

    end subroutine apply_bc
    subroutine initi_c(un,us,vn,vs,wn,ws,nx,ny,nz)
        use, intrinsic :: iso_c_binding
        implicit none
        integer,intent(in) :: nx,ny,nz
        real(C_DOUBLE),intent(inout) :: un(nx,ny+2,nz),vn(nx,ny+1,nz),wn(nx,ny+2,nz),us(nx,ny+2,nz),vs(nx,ny+1,nz),ws(nx,ny+2,nz)
        un(:,2:ny+1,:) = 1.0d0
        us(:,2:ny+1,:) = 1.0d0
    end subroutine initi_c
end module boundary
