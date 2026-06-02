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
#ifdef USE_IBM_G
    use :: ibmm, only: ibm_type
    use :: rk3_step_func, only: rk3_coeff
#elif USE_IBM
    use :: ibmm, only: ibm_type 
#endif
    implicit none
! LOCATION CHANGE !!! MOVED TO RK3 LOCATION
contains
#if defined(USE_IBM_G)
    subroutine corrector(f, g,rk3_c,n_loop)
#else
    subroutine corrector(f, g)
#endif
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g
#if defined(USE_IBM_G)
        type(rk3_coeff), intent(in)     :: rk3_c
        integer,intent(in)              :: n_loop
#endif
        integer :: i,j,k,im,km
        real(C_DOUBLE) :: dt_loop
#if defined(USE_IBM_G)
        dt_loop = rk3_c%c(n_loop)*g%dt
#else
        dt_loop = g%dt
#endif
        do i = 1, g%nx
            do j = 1, g%ny
                do k = 1, g%nz
                    
                    f%un(i,j,k) = f%us(i,j,k) - dt_loop*(f%pc(i,j,k)-f%pc(i-1,j,k))/g%dx
                end do
            end do
        end do

        do i = 1, g%nx
            do j = 2, g%ny
                do k = 1, g%nz
                    f%vn(i,j,k) = f%vs(i,j,k) - dt_loop*(f%pc(i,j,k)-f%pc(i,j-1,k))/g%dy
                end do
            end do
        end do

        do i = 1, g%nx
            do j = 1, g%ny
                do k = 1, g%nz
                    f%wn(i,j,k) = f%ws(i,j,k) - dt_loop*(f%pc(i,j,k)-f%pc(i,j,k-1))/g%dz
                end do
            end do
        end do
        
        f%pn = f%pn + f%pc

    end subroutine corrector


    
    real function get_cfl(f,g)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g

        integer :: i,j,k
        get_cfl = max(maxval(abs(f%un/g%dx)),&
                      maxval(abs(f%vn/g%dy)),&
                      maxval(abs(f%wn/g%dz)))
    end function get_cfl

end module step
