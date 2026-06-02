module richardson
    use,intrinsic :: iso_c_binding
    use :: init, only: grid_type,field_type 
#if defined(USE_IBM_G) || defined(USE_IBM) 
    use :: ibmm, only: ibm_type 
#endif
    implicit none
#define SOLID 1.0d30
    real(C_DOUBLE)::solid = 1.0d30
    contains
    ! is this a correct way? or do i need to add additioanlly check for any solid? 
#if defined(USE_IBM_G) || defined(USE_IBM)
    subroutine calc_mean_flow(f,g,ibm)
#else
    subroutine calc_mean_flow(f,g)
#endif
        implicit none
        type(grid_type)     :: g
        type(field_type)    :: f
#if defined(USE_IBM_G) || defined(USE_IBM)
        type(ibm_type)      :: ibm
#endif
        integer             :: i,j,k
        real(C_DOUBLE)      :: A_x
        A_x = g%dz*g%dy
        do k = 1, g%nz
            do j = 1, g%ny
                do i = 1,g%nx
#if defined(USE_IBM_G) || defined(USE_IBM)
                    if (ibm%coef_u(i,j,k)<0.5d0*SOLID)then
                        f%q_x = f%q_x+f%un(i,j,k)*A_x
                    end if
#else
                    f%q_x = f%q_x+f%un(i,j,k)*A_x
#endif
                end do 
            end do 
        end do
        f%q_x = f%q_x/g%nx
    end subroutine calc_mean_flow
#if defined(USE_IBM_G) || defined(USE_IBM)
    subroutine save_mean_flow(f,g,file_name,ibm)
#else
    subroutine save_mean_flow(f,g,file_name)
#endif
        implicit none
        type(grid_type)     :: g
        type(field_type)    :: f
        integer :: io
#if defined(USE_IBM_G) || defined(USE_IBM)
        type(ibm_type)      :: ibm
#endif
        character(len=*), intent(in) :: file_name
        ! --------------------------------------------------------------------------
#if defined(USE_IBM_G) || defined(USE_IBM)
        call calc_mean_flow(f,g,ibm)
#else
        call calc_mean_flow(f,g)
#endif 
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