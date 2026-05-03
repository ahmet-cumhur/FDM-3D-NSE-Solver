module rk3_binders
    use,intrinsic :: iso_c_binding
    use :: step
    use :: init, only: grid_type, field_type
    use :: fftw_3d
    use :: boundary
    use :: get_rhs, only: mom_rhs_compute         
#ifdef USE_IBM
    use :: ibmm            
#endif

    implicit none


    contains
        subroutine save_old_vel(f)
            ! at the start of the time loop we need to save
            ! un because we need that later
            type(field_type),intent(inout) :: f 
            f%u0(:,:,:) = f%un(:,:,:)
            f%v0(:,:,:) = f%vn(:,:,:)
            f%w0(:,:,:) = f%wn(:,:,:)
        end subroutine save_old_vel
        

        subroutine get_rhs0(f,g)
            type(field_type),intent(inout) :: f
            type(grid_type),intent(in)     :: g
            
            ! get the rhs0
            call mom_rhs_compute(f%rhs_0_u,f%rhs_0_v,f%rhs_0_w,f,g)
        end subroutine get_rhs0

#ifdef USE_IBM
        subroutine rk3_1st(f, g,ibm,ws)
#else   
        subroutine rk3_1st(f, g,ws)
#endif
            integer :: i,j,k
            ! here we apply the momentum equation normally only w/ dt/3
            type(grid_type)              :: g       
            type(field_type)             :: f
            real(C_DOUBLE)               :: dt_step
            type(poisson_fft_workspace)  :: ws       ! FFTW arrays, plans and parameters
#ifdef USE_IBM
            type(ibm_type)		         :: ibm      ! IBM arrays and parameters 
#endif    
            dt_step = real(g%dt/3,kind=C_DOUBLE) 
            call get_rhs0(f,g)
        
            do k = 1, g%nz
                do j = 1, g%ny
                    do i = 1, g%nx
                        f%us(i,j,k) = f%u0(i,j,k) +(dt_step)*f%rhs_0_u(i,j,k)
                        f%vs(i,j,k) = f%v0(i,j,k) +(dt_step)*f%rhs_0_v(i,j,k)
                        f%ws(i,j,k) = f%w0(i,j,k) +(dt_step)*f%rhs_0_w(i,j,k)
                    end do 
                end do 
            end do 
#ifdef USE_IBM
            call apply_ibm(f%us, ibm%coef_u, g)
            call apply_ibm(f%vs, ibm%coef_v, g)
            call apply_ibm(f%ws, ibm%coef_w, g)
#endif
            call apply_bc(f,g)
            call divU(f, g,dt_step)
            call poisson(g, f, ws)
            call apply_bc(f,g)
            call corrector(f, g,dt_step)
#ifdef USE_IBM
            call apply_ibm(f%un, ibm%coef_u, g)
            call apply_ibm(f%vn, ibm%coef_v, g)
            call apply_ibm(f%wn, ibm%coef_w, g)
#endif
            call apply_bc(f,g)    
        end subroutine rk3_1st

#ifdef USE_IBM
    subroutine rk3_2nd(f, g,ibm,ws)
#else
    subroutine rk3_2nd(f, g,ws)
#endif
        integer :: i,j,k
        ! here we apply the momentum equation normally only w/ dt/3
        type(grid_type)              :: g       
        type(field_type)             :: f
        real(C_DOUBLE)               :: dt_step
        type(poisson_fft_workspace)  :: ws       ! FFTW arrays, plans and parameters
#ifdef USE_IBM
        type(ibm_type)		         :: ibm      ! IBM arrays and parameters 
#endif    
        dt_step = real(2*g%dt/3,kind=C_DOUBLE) 
        call mom_rhs_compute(f%rhs_int_u,f%rhs_int_v,f%rhs_int_w,f,g)
        
        do k = 1, g%nz
            do j = 1, g%ny
                do i = 1, g%nx
                    f%us(i,j,k) = f%u0(i,j,k) +(dt_step)*f%rhs_int_u(i,j,k)
                    f%vs(i,j,k) = f%v0(i,j,k) +(dt_step)*f%rhs_int_v(i,j,k)
                    f%ws(i,j,k) = f%w0(i,j,k) +(dt_step)*f%rhs_int_w(i,j,k)
                end do 
            end do 
        end do 
#ifdef USE_IBM
        call apply_ibm(f%us, ibm%coef_u, g)
        call apply_ibm(f%vs, ibm%coef_v, g)
        call apply_ibm(f%ws, ibm%coef_w, g)
#endif
        call apply_bc(f,g)
        call divU(f, g,dt_step)
        call poisson(g, f, ws)
        call apply_bc(f,g)
        call corrector(f, g,dt_step)
#ifdef USE_IBM
        call apply_ibm(f%un, ibm%coef_u, g)
        call apply_ibm(f%vn, ibm%coef_v, g)
        call apply_ibm(f%wn, ibm%coef_w, g)
#endif
        call apply_bc(f,g)    
    end subroutine rk3_2nd

#ifdef USE_IBM
    subroutine rk3_3rd(f, g,ibm,ws)
#else
    subroutine rk3_3rd(f, g,ws)
#endif
        integer :: i,j,k
        type(grid_type)              :: g       
        type(field_type)             :: f
        real(C_DOUBLE)               :: dt_step
        type(poisson_fft_workspace)  :: ws       ! FFTW arrays, plans and parameters
#ifdef USE_IBM
        type(ibm_type)		         :: ibm      ! IBM arrays and parameters 
#endif

        dt_step = real(g%dt/4,kind=C_DOUBLE)
        call mom_rhs_compute(f%rhs_int_u,f%rhs_int_v,f%rhs_int_w,f,g)
         
        do k = 1, g%nz
            do j = 1, g%ny
                do i = 1, g%nx
                    f%us(i,j,k) = f%u0(i,j,k) +(dt_step)*(3*f%rhs_int_u(i,j,k)+f%rhs_0_u(i,j,k))
                    f%vs(i,j,k) = f%v0(i,j,k) +(dt_step)*(3*f%rhs_int_v(i,j,k)+f%rhs_0_v(i,j,k)) 
                    f%ws(i,j,k) = f%w0(i,j,k) +(dt_step)*(3*f%rhs_int_w(i,j,k)+f%rhs_0_w(i,j,k))
                end do 
            end do 
        end do 

#ifdef USE_IBM
            call apply_ibm(f%us, ibm%coef_u, g)
            call apply_ibm(f%vs, ibm%coef_v, g)
            call apply_ibm(f%ws, ibm%coef_w, g)
#endif
            call apply_bc(f,g)
            call divU(f, g,dt_step)
            call poisson(g, f, ws)
            call apply_bc(f,g)
            call corrector(f, g,dt_step)
#ifdef USE_IBM
            call apply_ibm(f%un, ibm%coef_u, g)
            call apply_ibm(f%vn, ibm%coef_v, g)
            call apply_ibm(f%wn, ibm%coef_w, g)
#endif
            call apply_bc(f,g)

    end subroutine rk3_3rd


end module rk3_binders