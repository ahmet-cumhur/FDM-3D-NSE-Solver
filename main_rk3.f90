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

program main_rk3
    use,intrinsic :: iso_c_binding
    use, intrinsic :: ieee_arithmetic

    ! Use modules
    ! ------------------------
    use :: init   	   ! Initialisation
    use :: boundary	   ! Boundary conditions (BC)
    use :: io              ! I/O
    use :: step            ! Timestep
    use :: fftw_3d         ! FFT Poisson solver
    use :: rk3_binders
    use :: get_rhs
#ifdef USE_IBM
    use :: ibmm            ! Immersed Boundary Method (IBM)
#endif

    ! Define variables
    ! ------------------------
    integer(C_INT)               :: i        ! Index for the time loop 
    type(grid_type)              :: g        ! Grid and flow parameters (to be added)
    type(field_type)             :: f        ! Flow field
    type(poisson_fft_workspace)  :: ws       ! FFTW arrays, plans and parameters
#ifdef USE_IBM
    type(ibm_type)		 :: ibm      ! IBM arrays and parameters 
#endif   
    character(len=256)           :: file_name

    ! Initialisations
    ! ------------------------
    print *, "initialising grid..."
    call init_grid(g)
    print *, "initialising fields..."
    call init_field(f, g)
    print *, "initialising poisson solver..."
    call init_poisson_fft_workspace(ws, g)
#ifdef USE_IBM
    print *, "initialising IBM..."
    call init_ibm(ibm, g)
    call set_ibm_coeff(g, ibm, ibm%coef_u, 1, 0, 0)
    call set_ibm_coeff(g, ibm, ibm%coef_v, 0, 1, 0)
    call set_ibm_coeff(g, ibm, ibm%coef_w, 0, 0, 1)    
#endif    

    ! Time loop
    ! ------------------------
    print *, "main loop starting..."
    i = 0
    do while(g%t_current<g%t_final)
        g%t_current = g%t_current + g%dt
        i = i + 1   
        call save_old_vel(f)
#ifdef USE_IBM
        call rk3_1st(f, g,ibm,ws)
        call rk3_2nd(f, g,ibm,ws)
        call rk3_3rd(f, g,ibm,ws)
#else
        call rk3_1st(f, g,ws)
        call rk3_2nd(f, g,ws)
        call rk3_3rd(f, g,ws)
#endif
        
        g%cfl=get_cfl(f,g)
        
        if (g%cflmax>0 .and. g%cfl>0) then
        	g%dt = min(g%cflmax/g%cfl,g%dtmax)
        end if
        

        if (modulo(i,100) == 0)then
            write(file_name,'("data_",I0,".vtk")') i
            print*, "current time step: ", i, "   filename: ", file_name, "   cfl:", g%cfl*g%dt
            call center_vel(f,g)
            call data_output(f,g,file_name)
        end if 
    end do 

    print *, "main loop ended..."
    call destroy_poisson_fft_workspace(ws)
end program main_rk3
