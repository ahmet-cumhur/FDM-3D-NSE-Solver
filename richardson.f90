module richardson
    use,intrinsic :: iso_c_binding
    use :: init, only: grid_type,field_type   
    implicit none
    contains
    ! is this a correct way? or do i need to add additioanlly check for any solid? 
    subroutine calc_mean_flow(f,g)
        implicit none
        type(grid_type)     :: g
        type(field_type)    :: f 
        integer             :: i,j,k
        real(C_DOUBLE)      :: A_x!,A_y_A_z
        A_x = g%dz*g%dy
        !A_y = g%dx*g%dz
        !A_z = g%dx*g%dy
        do k = 1, g%nz
            do j = 1, g%ny
                do i = 1,g%nx
                    f%q_x = f%q_x+f%un(i,j,k)*A_x
                end do 
            end do 
        end do
        f%q_x = f%q_x/g%nx
    end subroutine calc_mean_flow
    subroutine save_mean_flow(f,g,file_name)
        implicit none
        type(grid_type)     :: g
        type(field_type)    :: f
        integer :: io
        character(len=*), intent(in) :: file_name
        ! --------------------------------------------------------------------------
        call calc_mean_flow(f,g)
        open(newunit= io,file= trim(file_name),status="replace",action= "write" )
        
        write(io,'(A,1X,3ES20.12)') "current time[s]:       ", g%t_current
        write(io,'(A,1X,3ES20.12)') "time step size[s]:     ", g%dt
        write(io,'(A,1X,ES24.16)')  "mean flow:             ", f%q_x
        write(io,'(A,1X,I0)'      ) "nx:                    ", g%nx
        write(io,'(A,1X,I0)'      ) "ny:                    ", g%ny
        write(io,'(A,1X,I0)'      ) "nz:                    ", g%nz
        write(io,'(A,1X,3ES20.12)') "dx:                    ", g%dx
        write(io,'(A,1X,3ES20.12)') "dy:                    ", g%dy
        write(io,'(A,1X,3ES20.12)') "dz:                    ", g%dz
        
        close(io)

        
    end subroutine save_mean_flow
end module richardson