!--------------------------!
!                          !
!       Time-stepper       !
!          module          !
!                          !
!--------------------------! 
! 
! authors: Dr.-Ing. Davide Gatti
!          B.Sc. Ahmet Cumhur
! 
! date:    28.04.26
! 

module step
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
    implicit none

contains


    subroutine corrector(f, g,dt_step)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g
        real(C_DOUBLE), intent(in)      :: dt_step

        integer :: i,j,k,im,km

        do i = 1, g%nx
            do j = 1, g%ny
                do k = 1, g%nz
                    f%un(i,j,k) = f%us(i,j,k) - dt_step*(f%pc(i,j,k)-f%pc(i-1,j,k))/g%dx
                end do
            end do
        end do

        do i = 1, g%nx
            do j = 2, g%ny
                do k = 1, g%nz
                    f%vn(i,j,k) = f%vs(i,j,k) - dt_step*(f%pc(i,j,k)-f%pc(i,j-1,k))/g%dy
                end do
            end do
        end do

        do i = 1, g%nx
            do j = 1, g%ny
                do k = 1, g%nz
                    f%wn(i,j,k) = f%ws(i,j,k) - dt_step*(f%pc(i,j,k)-f%pc(i,j,k-1))/g%dz
                end do
            end do
        end do
        
        f%pn = f%pn + f%pc

    end subroutine corrector


    subroutine divU(f, g,dt_step)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g
        real(C_DOUBLE), intent(in)      :: dt_step
        
        integer :: i,j,k

        do i = 1, g%nx
            do j = 1, g%ny
                do k = 1, g%nz

                    f%rhs(i,j,k) = ( &
                        (f%us(i+1,j,k)-f%us(i,j,k))/g%dx &
                      + (f%vs(i,j+1,k)-f%vs(i,j,k))/g%dy &
                      + (f%ws(i,j,k+1)-f%ws(i,j,k))/g%dz ) / dt_step

                end do
            end do
        end do

    end subroutine divU
    
    real function get_cfl(f,g)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g

        integer :: i,j,k
        get_cfl = max(maxval(abs(f%un/g%dx)),&
                      maxval(abs(f%vn/g%dy)),&
                      maxval(abs(f%wn/g%dz)))
    end function get_cfl

end module step
