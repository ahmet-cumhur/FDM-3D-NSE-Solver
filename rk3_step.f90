module rk3_step_func
    use, intrinsic :: iso_c_binding
    use :: init, only : field_type,grid_type 
    use :: get_rhs
    use :: fftw_3d
    use :: boundary
    use :: step

#ifdef USE_IBM_G 
    use :: ibmm
#elif USE_IBM
    use :: ibmm
#endif 
    implicit none

    type rk3_coeff
    ! we have different A,B coefficients for each velocity!!
    real(C_DOUBLE)  ::   a(3) = [real(64.0d0/120.0d0,kind=C_DOUBLE),real(50.0d0/120.0d0,kind=C_DOUBLE),real(90.0d0/120.0d0,kind=C_DOUBLE)]
    real(C_DOUBLE)  ::   b(3) = [real(0.0d0/120.0d0,kind=C_DOUBLE),real(-34.0d0/120.0d0,kind=C_DOUBLE),real(-50.0d0/120.0d0,kind=C_DOUBLE)]
    real(C_DOUBLE)  ::   c(3) = [real(64.0d0/120.0d0,kind=C_DOUBLE),real(16.0d0/120.0d0,kind=C_DOUBLE),real(40.0d0/120.0d0,kind=C_DOUBLE)]
#ifdef USE_IBM_G   
    real(C_DOUBLE),allocatable :: A_rk3_u(:,:,:),B_rk3_u(:,:,:),A_rk3_v(:,:,:),B_rk3_v(:,:,:),A_rk3_w(:,:,:),B_rk3_w(:,:,:)
#endif     
    end type rk3_coeff

    contains
#ifdef USE_IBM_G
        subroutine init_rk3_arrays(g,rk3_c)
            implicit none
            type(rk3_coeff),intent(inout) :: rk3_c
            type(grid_type),intent(in)    :: g
            allocate(rk3_c%A_rk3_u(0:g%nx+1,0:g%ny+1,0:g%nz+1),rk3_c%B_rk3_u(0:g%nx+1,0:g%ny+1,0:g%nz+1))
            allocate(rk3_c%A_rk3_v(0:g%nx+1,1:g%ny+1,0:g%nz+1),rk3_c%B_rk3_v(0:g%nx+1,1:g%ny+1,0:g%nz+1))
            allocate(rk3_c%A_rk3_w(0:g%nx+1,0:g%ny+1,0:g%nz+1),rk3_c%B_rk3_w(0:g%nx+1,0:g%ny+1,0:g%nz+1))
            
            rk3_c%A_rk3_u = 1.0d0
            rk3_c%B_rk3_u = 1.0d0

            rk3_c%A_rk3_v = 1.0d0
            rk3_c%B_rk3_v = 1.0d0
            
            rk3_c%A_rk3_w = 1.0d0
            rk3_c%B_rk3_w = 1.0d0
        end subroutine init_rk3_arrays
#endif
#ifdef USE_IBM_G
            ! Here we calculte the A,B coefficients for rk3
            subroutine calc_a_b(g,ibm,rk3_c)
            implicit none
            ! we calculte in each real time step(non-rk3) 
            ! so we get a field of numbers
            type(ibm_type),intent(in)   :: ibm
            type(grid_type),intent(in)  :: g
            type(rk3_coeff),intent(inout) :: rk3_c
            integer :: i,j,k
            real(C_DOUBLE) :: ksi_u,ksi_v,ksi_w
            real(C_DOUBLE) :: e_u,e_v,e_w
            do k = 1,g%nz
                do j = 1,g%ny
                    do i = 1,g%nx
                        ! here we calculte the A and B
                        ! we need to add the re number!!
                        ksi_u = ibm%coef_u_lap(i,j,k)*g%dt / g%re
                        ksi_w = ibm%coef_w_lap(i,j,k)*g%dt / g%re
                        !check if exp part too small
                        e_u = exp(ksi_u)-1.0d0;e_w = exp(ksi_w)-1.0d0
                        ! get the B first ofc w/ taylor exp. if its too small
                        ! if we let ksi/ksi this causes NaN we need to simplfy
                        

                        !                   NEW QUESTIONS HERE                  !
                        ! this is where i calculate the A and B are they correct?
                        ! I used the way we talked about ---> approx. the e^x -1
                        ! if its too small ---> use taylor expansion
                        !-----------------------------------------------------!
                        ! is this a correct way to go? 
                        if (e_u >1e-5) then
                            rk3_c%B_rk3_u(i,j,k) =  (ksi_u)/(e_u)
                        else
                            rk3_c%B_rk3_u(i,j,k) =  1.0d0/(1.0d0+ksi_u/2.0d0+ksi_u**2/6.0d0+ksi_u**3/24.0d0)
                        end if

                        if (e_w >1e-5) then
                            rk3_c%B_rk3_w(i,j,k) =  (ksi_w)/(exp(ksi_w)-1)
                        else
                            rk3_c%B_rk3_w(i,j,k) =  1.0d0/(1.0d0+ksi_w/2.0d0+ksi_w**2/6.0d0+ksi_w**3/24.0d0)
                        end if
                        ! now get the A
                        rk3_c%A_rk3_u(i,j,k) = ksi_u+rk3_c%B_rk3_u(i,j,k)
                        rk3_c%A_rk3_w(i,j,k) = ksi_w+rk3_c%B_rk3_w(i,j,k) 
                    end do 
                end do
            end do            
            do k = 1,g%nz
                do j = 2, g%ny
                    do i = 1,g%nx
                        ! find the v dependent parameters here
                        ! since it has different shape than others
                        ksi_v = ibm%coef_v_lap(i,j,k)*g%dt / g%re
                        e_v = exp(ksi_v)-1.0d0
                        if (e_v >1e-5) then
                            rk3_c%B_rk3_v(i,j,k) =  (ksi_v)/(e_v)
                        else
                            rk3_c%B_rk3_v(i,j,k) =  1.0d0/(1.0d0+ksi_v/2.0d0+ksi_v**2/6.0d0+ksi_v**3/24.0d0)
                        end if
                        rk3_c%A_rk3_v(i,j,k) = ksi_v+rk3_c%B_rk3_v(i,j,k)
                    end do 
                end do 
            end do 

        end subroutine calc_a_b
#endif
        ! location changed
        ! should i add commented out part? 
        ! what about corrector step? 
        subroutine divU_rk3(f, g,rk3_c,n_loop)
            type(field_type),   intent(inout) :: f
            type(rk3_coeff),    intent(in)    :: rk3_c
            type(grid_type),    intent(in)    :: g
            integer,intent(in)                :: n_loop 
            integer :: i,j,k

            do i = 1, g%nx
                do j = 1, g%ny
                    do k = 1, g%nz

                        f%rhs(i,j,k) = ( &
                        (f%us(i+1,j,k)-f%us(i,j,k))/g%dx &
                        + (f%vs(i,j+1,k)-f%vs(i,j,k))/g%dy &         ! THIS PART
                        + (f%ws(i,j,k+1)-f%ws(i,j,k))/g%dz ) / (g%dt)!*rk3_c%c(n_loop))

                    end do
                end do
            end do

        end subroutine divU_rk3

        ! here we create a func. for each rk3 substep
        ! w/ its own projection steps
#ifdef USE_IBM
        subroutine main_loop(f,g,ibm,wss,rk3_c,n_loop)
#elif USE_IBM_G
        subroutine main_loop(f,g,ibm,wss,rk3_c,n_loop)
#else  
        subroutine main_loop(f,g,wss,rk3_c,n_loop)
#endif
            implicit none
            type(field_type) :: f
            type(grid_type) :: g
#ifdef USE_IBM_G
            type(ibm_type) :: ibm
#elif  USE_IBM
            type(ibm_type) :: ibm
#endif
            type(poisson_fft_workspace) :: wss
            type(rk3_coeff) :: rk3_c
            integer,intent(in) :: n_loop 
#ifdef USE_IBM
            call apply_ibm(f%us, ibm%coef_u, g)
            call apply_ibm(f%vs, ibm%coef_v, g)
            call apply_ibm(f%ws, ibm%coef_w, g)
#elif USE_IBM_G
            call apply_ibm(f%us, ibm%coef_u, g)
            call apply_ibm(f%vs, ibm%coef_v, g)
            call apply_ibm(f%ws, ibm%coef_w, g)
#endif
            call apply_bc(f,g)
            call divU_rk3(f, g,rk3_c,n_loop)
            call poisson(g, f, wss)
            call apply_bc(f,g)
            ! SHOULD I ADD N_LOOP?? 
            call corrector(f, g)
            call apply_bc(f,g)
        end subroutine main_loop

        ! first rk3 substep
        ! rk3 scheme is from luchini, gatti ibm paper
#ifdef USE_IBM_G
        subroutine rk3_first_st(f,g,rk3_c,ibm,wss)
#elif USE_IBM
        subroutine rk3_first_st(f,g,rk3_c,ibm,wss)
#else
        subroutine rk3_first_st(f,g,rk3_c,wss)
#endif
            implicit none
            type(field_type),intent(inout)            :: f
            type(grid_type),intent(in)                :: g
            type(rk3_coeff),intent(in)                :: rk3_c
#ifdef USE_IBM_G
            type(ibm_type),intent(in)                 :: ibm
#elif USE_IBM 
            type(ibm_type),intent(in)                 :: ibm
#endif
            type(poisson_fft_workspace),intent(inout) :: wss
            integer :: i,j,k
            call mom_rhs_compute(f%mom_rhs_u0,f%mom_rhs_v0,f%mom_rhs_w0,f,g)
            do k = 1, g%nz
                do j = 1, g%ny
                    do i = 1, g%nx
#ifdef USE_IBM_G
                    f%us(i,j,k) = (rk3_c%B_rk3_u(i,j,k) * f%un(i,j,k) + g%dt*rk3_c%a(1)*f%mom_rhs_u0(i,j,k))/rk3_c%A_rk3_u(i,j,k)
                    f%ws(i,j,k) = (rk3_c%B_rk3_w(i,j,k) * f%wn(i,j,k) + g%dt*rk3_c%a(1)*f%mom_rhs_w0(i,j,k))/rk3_c%A_rk3_w(i,j,k)
#else
                    f%us(i,j,k) = (f%un(i,j,k) + g%dt*rk3_c%a(1)*f%mom_rhs_u0(i,j,k))
                    f%ws(i,j,k) = (f%wn(i,j,k) + g%dt*rk3_c%a(1)*f%mom_rhs_w0(i,j,k))
#endif 
                    end do 
                end do 
            end do 
            ! we seperate the v and u,w because of their shape difference
            do k = 1,g%nz
                do  j = 2,g%ny
                    do i = 1, g%nx
#ifdef USE_IBM_G
                        f%vs(i,j,k) = (rk3_c%B_rk3_v(i,j,k) * f%vn(i,j,k) + g%dt*rk3_c%a(1)*f%mom_rhs_v0(i,j,k))/rk3_c%A_rk3_v(i,j,k)
#else
                        f%vs(i,j,k) = (f%vn(i,j,k) + g%dt*rk3_c%a(1)*f%mom_rhs_v0(i,j,k))
#endif
                    end do
                end do 
            end do 
            
            

#ifdef USE_IBM_G
            call main_loop(f,g,ibm,wss,rk3_c,1)
#elif USE_IBM
            call main_loop(f,g,ibm,wss,rk3_c,1)
#else 
            call main_loop(f,g,wss,rk3_c,1)
#endif
            ! now the us and un we found here are gonna be used in the next rk3 time step
        end subroutine rk3_first_st
        
        ! second rk3 substep
#ifdef USE_IBM_G
        subroutine rk3_second_st(f,g,rk3_c,ibm,wss)
#elif USE_IBM
        subroutine rk3_second_st(f,g,rk3_c,ibm,wss)
#else 
        subroutine rk3_second_st(f,g,rk3_c,wss)
#endif
            implicit none
            type(field_type),intent(inout)            :: f
            type(grid_type),intent(in)                :: g
            type(rk3_coeff),intent(in)                :: rk3_c
#ifdef USE_IBM_G
            type(ibm_type),intent(in)                 :: ibm
#elif  USE_IBM
            type(ibm_type),intent(in)                 :: ibm
#endif
            type(poisson_fft_workspace),intent(inout) :: wss
            
            integer :: i,j,k
            ! dont forget to update the mom_rhs_u0,v0,w0 in the main loop after this subroutine call
            call mom_rhs_compute(f%mom_rhs_u_int,f%mom_rhs_v_int,f%mom_rhs_w_int,f,g)
            do k = 1, g%nz
                do j = 1, g%ny
                    do i = 1, g%nx
#ifdef USE_IBM_G
                    f%us(i,j,k) = (rk3_c%B_rk3_u(i,j,k) * f%un(i,j,k) +&
                     g%dt*(rk3_c%a(2)*f%mom_rhs_u_int(i,j,k)+rk3_c%b(2)*f%mom_rhs_u0(i,j,k)))&
                    /rk3_c%A_rk3_u(i,j,k)
                    f%ws(i,j,k) = (rk3_c%B_rk3_w(i,j,k) * f%wn(i,j,k) +&
                    g%dt*(rk3_c%a(2)*f%mom_rhs_w_int(i,j,k)+rk3_c%b(2)*f%mom_rhs_w0(i,j,k)))&
                    /rk3_c%A_rk3_w(i,j,k)
#else
                    f%us(i,j,k) = (f%un(i,j,k) +&
                     g%dt*(rk3_c%a(2)*f%mom_rhs_u_int(i,j,k)+rk3_c%b(2)*f%mom_rhs_u0(i,j,k)))
                    f%ws(i,j,k) = (f%wn(i,j,k) +&
                    g%dt*(rk3_c%a(2)*f%mom_rhs_w_int(i,j,k)+rk3_c%b(2)*f%mom_rhs_w0(i,j,k)))
#endif    
                    end do 
                end do 
            end do 
            ! we  do the same seperation to the v here again
            do k = 1,g%nz
                do j = 2,g%ny
                    do i = 1, g%nx
#ifdef USE_IBM_G
                    f%vs(i,j,k) = (rk3_c%B_rk3_v(i,j,k) * f%vn(i,j,k) +&
                     g%dt*(rk3_c%a(2)*f%mom_rhs_v_int(i,j,k)+rk3_c%b(2)*f%mom_rhs_v0(i,j,k)))&
                    /rk3_c%A_rk3_v(i,j,k)
#else
                    f%vs(i,j,k) = (f%vn(i,j,k) +&
                     g%dt*(rk3_c%a(2)*f%mom_rhs_v_int(i,j,k)+rk3_c%b(2)*f%mom_rhs_v0(i,j,k)))
#endif
                    end do 
                end do 
            end do 



#ifdef USE_IBM_G
            call main_loop(f,g,ibm,wss,rk3_c,2)
#elif USE_IBM
            call main_loop(f,g,ibm,wss,rk3_c,2)
#else 
            call main_loop(f,g,wss,rk3_c,2)
#endif
        end subroutine rk3_second_st
        ! third rk3 substep
#ifdef USE_IBM_G
        subroutine rk3_third_st(f,g,rk3_c,ibm,wss)
#elif USE_IBM 
        subroutine rk3_third_st(f,g,rk3_c,ibm,wss)
#else
        subroutine rk3_third_st(f,g,rk3_c,wss)
#endif
            implicit none
            type(field_type),intent(inout)            :: f
            type(grid_type),intent(in)                :: g
            type(rk3_coeff),intent(in)                :: rk3_c
#ifdef USE_IBM_G
            type(ibm_type),intent(in)                 :: ibm
#elif  USE_IBM
            type(ibm_type),intent(in)                 :: ibm
#endif
            type(poisson_fft_workspace),intent(inout) :: wss
            
            integer :: i,j,k
            f%mom_rhs_u0 = f%mom_rhs_u_int
            f%mom_rhs_v0 = f%mom_rhs_v_int
            f%mom_rhs_w0 = f%mom_rhs_w_int
            call mom_rhs_compute(f%mom_rhs_u_int,f%mom_rhs_v_int,f%mom_rhs_w_int,f,g)
            do k = 1, g%nz
                do j = 1, g%ny
                    do i = 1, g%nx
#ifdef USE_IBM_G
                    f%us(i,j,k) = (rk3_c%B_rk3_u(i,j,k) * f%un(i,j,k) +&
                     g%dt*(rk3_c%a(3)*f%mom_rhs_u_int(i,j,k)+rk3_c%b(3)*f%mom_rhs_u0(i,j,k)))&
                    /rk3_c%A_rk3_u(i,j,k)
                    f%ws(i,j,k) = (rk3_c%B_rk3_w(i,j,k) * f%wn(i,j,k) +&
                    g%dt*(rk3_c%a(3)*f%mom_rhs_w_int(i,j,k)+rk3_c%b(3)*f%mom_rhs_w0(i,j,k)))&
                    /rk3_c%A_rk3_w(i,j,k)
#else
                    f%us(i,j,k) = (f%un(i,j,k) +&
                     g%dt*(rk3_c%a(3)*f%mom_rhs_u_int(i,j,k)+rk3_c%b(3)*f%mom_rhs_u0(i,j,k)))
                    f%ws(i,j,k) = (f%wn(i,j,k) +&
                    g%dt*(rk3_c%a(3)*f%mom_rhs_w_int(i,j,k)+rk3_c%b(3)*f%mom_rhs_w0(i,j,k))) 
#endif             
                    end do 
                end do 
            end do 
            do k = 1, g%nz
                do j = 2, g%ny
                    do i = 1,g%nx
#ifdef USE_IBM_G
                    f%vs(i,j,k) = (rk3_c%B_rk3_v(i,j,k) * f%vn(i,j,k) +&
                     g%dt*(rk3_c%a(3)*f%mom_rhs_v_int(i,j,k)+rk3_c%b(3)*f%mom_rhs_v0(i,j,k)))&
                    /rk3_c%A_rk3_v(i,j,k)
#else
                    f%vs(i,j,k) = (f%vn(i,j,k) +&
                     g%dt*(rk3_c%a(3)*f%mom_rhs_v_int(i,j,k)+rk3_c%b(3)*f%mom_rhs_v0(i,j,k)))
#endif
                    end do
                end do
            end do 

        
#ifdef USE_IBM_G
            call main_loop(f,g,ibm,wss,rk3_c,3)
#elif USE_IBM
            call main_loop(f,g,ibm,wss,rk3_c,3)
#else 
            call main_loop(f,g,wss,rk3_c,3)
#endif
        end subroutine rk3_third_st

end module rk3_step_func