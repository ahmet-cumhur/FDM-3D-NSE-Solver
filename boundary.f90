module boundary
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
    implicit none

contains

    subroutine apply_bc(f, g)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g

        ! Top wall: no slip
        ! -----------------
        f%un(:,g%ny+1,:) = -f%un(:,g%ny,:)	
        f%us(:,g%ny+1,:) = -f%us(:,g%ny,:)    
        
        f%vn(:,g%ny+1,:) = 0.d0
        f%vs(:,g%ny+1,:) = 0.d0
        
        f%wn(:,g%ny+1,:) = -f%wn(:,g%ny,:)	
        f%ws(:,g%ny+1,:) = -f%ws(:,g%ny,:)	
        
        ! Bottom wall: no slip
        ! -----------------
        f%un(:,0,:) = -f%un(:,1,:)	
        f%us(:,0,:) = -f%us(:,1,:)    
        
        f%vn(:,1,:) = 0.d0
        f%vs(:,1,:) = 0.d0
        
        f%wn(:,0,:) = -f%wn(:,1,:)	
        f%ws(:,0,:) = -f%ws(:,1,:)	    
        
        ! Inlet boundary: periodic ( ix=nx -> ix=0 )
        ! -----------------
        f%un(0,:,:) = f%un(g%nx,:,:)
        f%us(0,:,:) = f%us(g%nx,:,:)

        f%vn(0,:,:) = f%vn(g%nx,:,:)
        f%vs(0,:,:) = f%vs(g%nx,:,:)        
   
        f%wn(0,:,:) = f%wn(g%nx,:,:)
        f%ws(0,:,:) = f%ws(g%nx,:,:)      
        
        f%pn(0,:,:) = f%pn(g%nx,:,:)
        f%pc(0,:,:) = f%pc(g%nx,:,:)

        
        ! Outlet boundary: periodic ( ix=1 -> ix=nx+1 )
        ! -----------------
        f%un(g%nx+1,:,:) = f%un(1,:,:)
        f%us(g%nx+1,:,:) = f%us(1,:,:)

        f%vn(g%nx+1,:,:) = f%vn(1,:,:)
        f%vs(g%nx+1,:,:) = f%vs(1,:,:)        
   
        f%wn(g%nx+1,:,:) = f%wn(1,:,:)
        f%ws(g%nx+1,:,:) = f%ws(1,:,:)      
        
        f%pn(g%nx+1,:,:) = f%pn(1,:,:)
        f%pc(g%nx+1,:,:) = f%pc(1,:,:)
        
        
        ! Left boundary: periodic ( iz=nz -> iz=0 )
        ! -----------------
        f%un(:,:,0) = f%un(:,:,g%nz)
        f%us(:,:,0) = f%us(:,:,g%nz)

        f%vn(:,:,0) = f%vn(:,:,g%nz)
        f%vs(:,:,0) = f%vs(:,:,g%nz)        
   
        f%wn(:,:,0) = f%wn(:,:,g%nz)
        f%ws(:,:,0) = f%ws(:,:,g%nz)      
        
        f%pn(:,:,0) = f%pn(:,:,g%nz)
        f%pc(:,:,0) = f%pc(:,:,g%nz)

        
        ! Right boundary: periodic ( iz=1 -> iz=nz+1 )
        ! -----------------
        f%un(:,:,g%nz+1) = f%un(:,:,1)
        f%us(:,:,g%nz+1) = f%us(:,:,1)
        
        f%vn(:,:,g%nz+1) = f%vn(:,:,1)
        f%vs(:,:,g%nz+1) = f%vs(:,:,1)        
   
        f%wn(:,:,g%nz+1) = f%wn(:,:,1)
        f%ws(:,:,g%nz+1) = f%ws(:,:,1)      
        
        f%pn(:,:,g%nz+1) = f%pn(:,:,1)
        f%pc(:,:,g%nz+1) = f%pc(:,:,1)
        
        
    end subroutine apply_bc


end module boundary
