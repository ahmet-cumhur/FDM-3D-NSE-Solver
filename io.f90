module io

    use :: init, only: grid_type,field_type
    implicit none
    contains 
    subroutine center_vel(f,g)
        use, intrinsic :: iso_c_binding 
        implicit none
        
        type(field_type), intent(inout) :: f
        type(grid_type), intent(in)     :: g
        integer :: i,j,k,ip,jp,kp
      
        do i = 1, g%nx
            do j = 1,g%ny
                do k = 1,g%nz
                    ip = i + 1
                    jp = j + 1
                    kp = k + 1
                    
                    f%uc(i,j,k) = 0.5d0 * (f%un(i,j,k) + f%un(ip,j,k))
                    f%vc(i,j,k) = 0.5d0 * (f%vn(i,j,k) + f%vn(i,jp,k))
                    f%wc(i,j,k) = 0.5d0 * (f%wn(i,j,k) + f%wn(i,j,kp))

                end do 
            end do 
        end do
    end subroutine center_vel

    subroutine data_output(f,g,file_name)
        use, intrinsic :: iso_c_binding
        implicit none

        type(field_type), intent(in) :: f
        type(grid_type), intent(in)     :: g
        
        integer :: io,npts,i,j,k
        character(len=*), intent(in) :: file_name
        real(C_DOUBLE) :: x0,y0,z0
        x0 = 0.5d0*g%dx
        y0 = 0.5d0*g%dy
        z0 = 0.5d0*g%dz
        npts = g%nx*g%ny*g%nz

        open(newunit= io,file= trim(file_name),status="replace",action= "write" )
        write(io,'(A)') "# vtk DataFile Version 3.0"
        write(io,'(A)') "3D velocity field"
        write(io,'(A)') "ASCII"
        write(io,'(A)')  "DATASET STRUCTURED_POINTS"

        write(io,'(A,1X,I0,1X,I0,1X,I0)') "DIMENSIONS",g%nx+1,g%ny+1,g%nz+1
        write(io,'(A,1X,3ES20.12)') "ORIGIN",0.0d0,0.0d0,0.0d0
        write(io,'(A,1X,3ES20.12)') "SPACING",g%dx,g%dy,g%dz

        write(io,'(A)') "FIELD FieldData 1"
        write(io,'(A)') "TIME 1 1 double"
        write(io,'(ES20.12)') g%t_current

        write(io,'(A,1X,I0)') "CELL_DATA",npts
        write(io,'(A)') "VECTORS U double"

        do k = 1,g%nz
            do j = 1,g%ny
                do i = 1,g%nx
 !                   print *, f%uc(i,j,k),f%vc(i,j,k),f%wc(i,j,k)
!                    write(io,*) f%uc(i,j,k),f%vc(i,j,k),f%wc(i,j,k)
                    write(io,'(3ES20.12)') f%uc(i,j,k),f%vc(i,j,k),f%wc(i,j,k)
                end do 
            end do 
        end do

        close(io)
    end subroutine data_output
end module io
