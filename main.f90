!--------------------------!
!                          !
!     Finite Difference    !
!      channel solver      !
!                          !
!--------------------------! 
! 
! authors: Dr.-Ing. Davide Gatti
!          B.Sc. Ahmet Cumhur
! 
! date:    28.04.26
! 

program main
    use,intrinsic :: iso_c_binding
    use, intrinsic :: ieee_arithmetic

    ! Use modules
    ! ------------------------
    use :: init   	   ! Initialisation
    use :: boundary	   ! Boundary conditions (BC)
    use :: io              ! I/O
    use :: step            ! Timestep
    use :: fftw_3d         ! FFT Poisson solver
#ifdef USE_IBM
    use :: ibmm            ! Immersed Boundary Method (IBM)
#elif USE_IBM_G 
    use :: ibmm
#endif
    ! rk3 is here
#ifdef USE_IBM_G
    use :: rk3_step_func, only: rk3_coeff, init_rk3_arrays, calc_a_b,rk3_first_st, rk3_second_st, rk3_third_st
#else
    use :: rk3_step_func, only: rk3_coeff,rk3_first_st, rk3_second_st, rk3_third_st
#endif
    use :: richardson,    only: calc_mean_flow,save_mean_flow
#if defined(USE_IBM_G) && defined(USE_IBM)
#error "Both Macros are on!"
#endif
#ifdef USE_OPENMP
    use :: omp_lib
#endif
    ! Define variables
    ! ------------------------
    integer(C_INT)               :: i        ! Index for the time loop 
    type(grid_type)              :: g        ! Grid and flow parameters (to be added)
    type(field_type)             :: f        ! Flow field
    type(poisson_fft_workspace)  :: wss       ! FFTW arrays, plans and parameters
    type(rk3_coeff)      :: rk3_c

#ifdef USE_IBM
    type(ibm_type)		 :: ibm      ! IBM arrays and parameters 
#elif USE_IBM_G 
    type(ibm_type)       :: ibm
#endif   
    character(len=256)           :: file_name
    character(len=256)           :: file_name_meanf
    real(C_DOUBLE)               :: t1,t2,t3,t_diff 

    ! Initialisations
    ! ------------------------
#ifdef USE_OPENMP
    print*,"working on number of threads: " ,omp_get_max_threads()
#endif 
    print *, "initialising grid..."
    call init_grid(g)
    print *, "initialising fields..."
    call init_field(f, g)
    print *, "initialising poisson solver..."
    call init_poisson_fft_workspace(wss, g)
#ifdef USE_IBM_G
    ! rk3 initialize
    print *,"initialize Runge-Kutta 3 variables..."
    call init_rk3_arrays(g,rk3_c)
#endif
#ifdef USE_IBM
    print *, "initialising IBM..."
    call init_ibm(ibm, g)
    call set_ibm_coeff(g, ibm, ibm%coef_u, 1, 0, 0,0,0,0)
    call set_ibm_coeff(g, ibm, ibm%coef_v, 0, 1, 0,0,1,0)
    call set_ibm_coeff(g, ibm, ibm%coef_w, 0, 0, 1,0,0,0)    
#endif    
#ifdef USE_IBM_G
    print *,"initializing IBM 2nd order..."
    call init_ibm(ibm, g)
    call set_ibm_coeff(g, ibm, ibm%coef_u, 1, 0, 0,0,0,0)
    call set_ibm_coeff(g, ibm, ibm%coef_v, 0, 1, 0,0,1,0)
    call set_ibm_coeff(g, ibm, ibm%coef_w, 0, 0, 1,0,0,0) 

    call set_ibm_coeff_2nd(g, ibm, ibm%coef_u, 1, 0, 0,ibm%coef_u_lap,0,0,0)
    call set_ibm_coeff_2nd(g, ibm, ibm%coef_v, 0, 1, 0,ibm%coef_v_lap,0,1,0)
    call set_ibm_coeff_2nd(g, ibm, ibm%coef_w, 0, 0, 1,ibm%coef_w_lap,0,0,0)
#endif
    ! Time loop
    ! ------------------------
#ifdef USE_OPENMP
t3 = omp_get_wtime()
#else
call cpu_time(t3)
#endif
    print *, "main loop starting..."
    i = 0
    do while(g%t_current<g%t_final)
        g%t_current = g%t_current + g%dt
        i = i + 1
        ! rk3 steps
#if defined(USE_IBM_G) || defined (USE_IBM)
        call rk3_first_st(f,g,rk3_c,ibm,wss)
        call rk3_second_st(f,g,rk3_c,ibm,wss)
        call rk3_third_st(f,g,rk3_c,ibm,wss)
#else 
        call rk3_first_st(f,g,rk3_c,wss)
        call rk3_second_st(f,g,rk3_c,wss)
        call rk3_third_st(f,g,rk3_c,wss)
#endif 
    
     
        if (modulo(i,100) == 0)then
            print*,"time step is at: ",i
            write(file_name,'("data_",I0,".vtk")') i
            print*, "current time step: ", i, "   filename: ", file_name, "   cfl:", g%cfl*g%dt
            call center_vel(f,g)
            call data_output(f,g,file_name)
        end if   
    end do 
    ! and we add a mean flow calculation at the end..
    ! calculation @ richardson.f90
#ifdef USE_OPENMP
    t2 = omp_get_wtime()
#else
    call cpu_time(t2)
#endif
    t_diff = t2-t3
    print*,"time spent for main loop: ",t_diff
    write(file_name_meanf,'("mf_data.txt")') 
#if defined(USE_IBM_G) || defined (USE_IBM)
    call save_mean_flow(f,g,file_name_meanf,ibm)
#else
    call save_mean_flow(f,g,file_name_meanf)
#endif
    !------------------------------------------------
    print *, "main loop ended..."
    call destroy_poisson_fft_workspace(wss)
end program main
