module data_save
    implicit none
    contains 
    subroutine center_vel(un,vn,wn,uc,vc,wc,dx,dy,dz,nx,ny,nz)
        use, intrinsic :: iso_c_binding 
        implicit none

        integer :: nx,ny,nz,i,j,k,kp,ip
        real(C_DOUBLE),intent(in) :: dx,dy,dz
        real(C_DOUBLE),intent(in) :: un(nx,ny+2,nz),vn(nx,ny+1,nz),wn(nx,ny+2,nz)
        real(C_DOUBLE),allocatable,intent(out) ::uc(:,:,:),vc(:,:,:),wc(:,:,:)
        allocate(uc(nx,ny,nz),vc(nx,ny,nz),wc(nx,ny,nz))        
        do i = 1, nx
            do j = 1,ny
                do k = 1,nz
                    ! periodic bc
                    ip = i + 1
                    kp = k + 1
                    if (ip > nx) ip = 1
                    if (kp > nz) kp = 1

                    uc(i,j,k) = 0.5d0 * (un(i ,j+1,k) + un(ip,j+1,k))
                    vc(i,j,k) = 0.5d0 * (vn(i,j,k) + vn(i,j+1,k))
                    wc(i,j,k) = 0.5d0 * (wn(i,j+1,k) + wn(i,j+1,kp))

                end do 
            end do 
        end do



    end subroutine center_vel

    subroutine data_output(uc,vc,wc,file_name,nx,ny,nz,dx,dy,dz,t_current)
        use, intrinsic :: iso_c_binding
        implicit none
        integer :: io,npts,i,j,k
        integer,intent(in) :: nx,ny,nz
        real(C_DOUBLE),intent(in) :: dx,dy,dz,t_current
        real(C_DOUBLE),intent(in) :: uc(nx,ny,nz),vc(nx,ny,nz),wc(nx,ny,nz)
        character(len=*), intent(in) :: file_name
        real(C_DOUBLE) :: x0,y0,z0
        x0 = 0.5d0*dx
        y0 = 0.5d0*dy
        z0 = 0.5d0*dz
        npts = nx*ny*nz

        open(newunit= io,file= trim(file_name),status="replace",action= "write" )
        write(io,'(A)') "# vtk DataFile Version 3.0"
        write(io,'(A)') "3D velocity field"
        write(io,'(A)') "ASCII"
        write(io,'(A)')  "DATASET STRUCTURED_POINTS"

        write(io,'(A,1X,I0,1X,I0,1X,I0)') "DIMENSIONS",nx+1,ny+1,nz+1 
        write(io,'(A,1X,3ES20.12)') "ORIGIN",0.0d0,0.0d0,0.0d0
        write(io,'(A,1X,3ES20.12)') "SPACING",dx,dy,dz

        write(io,'(A)') "FIELD FieldData 1"
        write(io,'(A)') "TIME 1 1 double"
        write(io,'(ES20.12)') t_current

        write(io,'(A,1X,I0)') "CELL_DATA",npts
        write(io,'(A)') "VECTORS U double"

        do k = 1,nz
            do j = 1,ny
                do i = 1,nx
                    write(io,'(3ES20.12)') uc(i,j,k),vc(i,j,k),wc(i,j,k)
                end do 
            end do 
        end do

        close(io)
    end subroutine data_output
end module data_save