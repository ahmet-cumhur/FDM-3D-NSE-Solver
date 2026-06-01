!--------------------------!
!                          !
!    Immersed Boundary     !
!         Module           !
!                          !
!--------------------------! 
! 
! authors: Dr.-Ing. Davide Gatti
!          B.Sc. Ahmet Cumhur
! 
! date:    28.04.26
! 


module ibmm
    use, intrinsic :: iso_c_binding
    use, intrinsic :: ieee_arithmetic
    use :: init, only: grid_type,field_type
    implicit none


    !========================
    ! IBM TYPE
    !========================
    type :: ibm_type
        integer :: n_wave_x, n_wave_z
        real(C_DOUBLE) :: amp_x, phase_x
        real(C_DOUBLE) :: amp_z, phase_z

        real(C_DOUBLE), allocatable :: coef_u(:,:,:), coef_v(:,:,:), coef_w(:,:,:)
#ifdef USE_IBM_G
        real(C_DOUBLE) :: lambda
        real(C_DOUBLE), allocatable :: coef_u_lap(:,:,:), coef_v_lap(:,:,:), coef_w_lap(:,:,:)
#endif
    end type ibm_type

#define SOLID 1.0d30
        real(C_DOUBLE)::solid = 1.0d30
contains


!========================
! INITIALIZE IBM
!========================
subroutine init_ibm(ibm, g)
    type(ibm_type), intent(inout) :: ibm
    type(grid_type), intent(in)   :: g

    ibm%n_wave_x = 1
    ibm%n_wave_z = 1
    ibm%amp_x = 0.15d0
    ibm%amp_z = 0.10d0
    ibm%phase_x = 0.0d0
    ibm%phase_z = 0.0d0
#ifdef USE_IBM_G
    ibm%lambda = 0.0d0
#endif
    allocate(ibm%coef_u(0:g%nx+1,0:g%ny+1,0:g%nz+1))
    allocate(ibm%coef_v(0:g%nx+1,1:g%ny+1,0:g%nz+1))
    allocate(ibm%coef_w(0:g%nx+1,0:g%ny+1,0:g%nz+1))

#ifdef USE_IBM_G
    allocate(ibm%coef_u_lap(0:g%nx+1,0:g%ny+1,0:g%nz+1))
    allocate(ibm%coef_v_lap(0:g%nx+1,1:g%ny+1,0:g%nz+1))
    allocate(ibm%coef_w_lap(0:g%nx+1,0:g%ny+1,0:g%nz+1))
    ibm%coef_u_lap(:,:,:) = 0.0d0
    ibm%coef_v_lap(:,:,:) = 0.0d0
    ibm%coef_w_lap(:,:,:) = 0.0d0
#endif

end subroutine init_ibm


    logical function isInBody(x, y, z, ibm, g)
        implicit none

        real(C_DOUBLE), intent(in) :: x, y, z
        type(ibm_type), intent(in) :: ibm
        type(grid_type), intent(in) :: g

        real(C_DOUBLE), parameter :: pi = 3.141592653589793d0
        real(C_DOUBLE) :: y_body
        real(C_DOUBLE) :: y0
        y0 = 0.25d0
        y_body = y0+ibm%amp_x * 0.5d0 * &
                 (1.0d0 + sin(2.0d0*pi*real(ibm%n_wave_x,C_DOUBLE)*x/g%lx + ibm%phase_x))  + &
                 ibm%amp_z * 0.5d0 * &
                (1.0d0 + sin(2.0d0*pi*real(ibm%n_wave_z,C_DOUBLE)*z/g%lz + ibm%phase_z))

        isInBody = (y < y_body)

    end function isInBody


    subroutine set_ibm_coeff(g, ibm, coeff, dix, diy, diz,i0,j0,k0)
        implicit none

        type(grid_type), intent(in) :: g
        type(ibm_type), intent(in) :: ibm
        integer,intent(in)      :: i0,j0,k0
        real(C_DOUBLE), intent(inout) :: coeff(i0:,j0:,k0:)

        integer, intent(in) :: dix, diy, diz
        integer :: ix, iy, iz
        real(C_DOUBLE) :: x, y, z

        coeff = 0.0d0

        do iz = lbound(coeff,3), ubound(coeff,3)
            do iy = lbound(coeff,2), ubound(coeff,2)
                do ix = lbound(coeff,1), ubound(coeff,1)

                    x = (real(ix,C_DOUBLE) -0.5d0-real(dix,C_DOUBLE)*0.5d0        )*g%dx
                    y = (real(iy,C_DOUBLE) -0.5d0-real(diy,C_DOUBLE)*0.5d0        )*g%dy
                    z = (real(iz,C_DOUBLE) -0.5d0-real(diz,C_DOUBLE)*0.5d0        )*g%dz

                    if (isInBody(x, y, z, ibm, g)) then
                        coeff(ix,iy,iz) = SOLID
                    end if

                end do
            end do
        end do

    end subroutine set_ibm_coeff

#ifdef USE_IBM_G
    subroutine set_ibm_coeff_2nd(g, ibm, coeff, dix, diy, diz,coef_lap,i0,j0,k0)
        implicit none

        type(grid_type), intent(in) :: g
        type(ibm_type), intent(inout) :: ibm
        integer,intent(in)      :: i0,j0,k0
        real(C_DOUBLE), intent(inout) :: coeff(i0:,j0:,k0:)
        real(C_DOUBLE), intent(inout) :: coef_lap(i0:,j0:,k0:)

        integer, intent(in) :: dix, diy, diz
        integer :: ix, iy, iz
        real(C_DOUBLE) :: x, y, z
        ! these are for the checking the neigbours of the given point
        real(C_DOUBLE) :: x_ip,x_im, y_jp,y_jm, z_kp,z_km

        coeff = 0.0d0
        coef_lap = 0.0d0

        do iz = lbound(coeff,3), ubound(coeff,3)
            do iy = lbound(coeff,2), ubound(coeff,2)
                do ix = lbound(coeff,1), ubound(coeff,1)
                    ! we save the neigbours of the given point
                    
                    
                    x = (real(ix,C_DOUBLE) -0.5d0-real(dix,C_DOUBLE)*0.5d0        )*g%dx
                    y = (real(iy,C_DOUBLE) -0.5d0-real(diy,C_DOUBLE)*0.5d0        )*g%dy
                    z = (real(iz,C_DOUBLE) -0.5d0-real(diz,C_DOUBLE)*0.5d0        )*g%dz

                    x_ip= x + g%dx
                    x_im= x - g%dx 
                    y_jp= y + g%dy
                    y_jm= y - g%dy
                    z_kp= z + g%dz 
                    z_km= z - g%dz

                    ! then we check if the given point is in body.
                    if (isInBody(x, y, z, ibm, g)) then
                        coeff(ix,iy,iz) = SOLID
                        ! if it isnt in body and its neigbour is then; 
                        ! we add a coefficient to laplacian
                    end if
                        ! first we go for x dir
                    if (.not. isInBody(x,y,z,ibm,g) .and. isInBody(x_ip,y,z,ibm,g)) then
                        call find_btw_points_x(x,x_ip,y,z,ibm,g,coeff,ix,iy,iz,i0,j0,k0)
                        ! I added these checks if the point we are are at is solid then we skip directly
                        if (coeff(ix,iy,iz) == SOLID)then
                            coef_lap(ix,iy,iz) = 0.0d0
                            cycle
                        end if 
                        coef_lap(ix,iy,iz) =  coef_lap(ix,iy,iz) + ibm%lambda
                    end if
                    if (.not. isInBody(x,y,z,ibm,g) .and. isInBody(x_im,y,z,ibm,g)) then 
                        call find_btw_points_x(x,x_im,y,z,ibm,g,coeff,ix,iy,iz,i0,j0,k0)
                        if (coeff(ix,iy,iz) == SOLID)then
                            coef_lap(ix,iy,iz) = 0.0d0
                            cycle
                        end if
                        coef_lap(ix,iy,iz) =  coef_lap(ix,iy,iz) + ibm%lambda
                    end if
                        ! z dir    
                    if (.not. isInBody(x,y,z,ibm,g) .and. isInBody(x,y,z_kp,ibm,g)) then
                        call find_btw_points_z(z,z_kp,y,x,ibm,g,coeff,ix,iy,iz,i0,j0,k0)
                        if (coeff(ix,iy,iz) == SOLID)then
                            coef_lap(ix,iy,iz) = 0.0d0
                            cycle
                        end if
                        coef_lap(ix,iy,iz) =  coef_lap(ix,iy,iz) + ibm%lambda
                    end if
                    if (.not. isInBody(x,y,z,ibm,g) .and. isInBody(x,y,z_km,ibm,g)) then 
                        call find_btw_points_z(z,z_km,y,x,ibm,g,coeff,ix,iy,iz,i0,j0,k0)
                        if (coeff(ix,iy,iz) == SOLID)then
                            coef_lap(ix,iy,iz) = 0.0d0
                            cycle
                        end if
                        coef_lap(ix,iy,iz) =  coef_lap(ix,iy,iz) + ibm%lambda
                    end if
                        ! y dir
                    if (.not. isInBody(x,y,z,ibm,g) .and. isInBody(x,y_jp,z,ibm,g)) then
                        call find_btw_points_y(y,y_jp,x,z,ibm,g,coeff,ix,iy,iz,i0,j0,k0)
                        if (coeff(ix,iy,iz) == SOLID)then
                            coef_lap(ix,iy,iz) = 0.0d0
                            cycle
                        end if
                        coef_lap(ix,iy,iz) =  coef_lap(ix,iy,iz) + ibm%lambda
                    end if
                    if (.not. isInBody(x,y,z,ibm,g) .and. isInBody(x,y_jm,z,ibm,g)) then 
                        call find_btw_points_y(y,y_jm,x,z,ibm,g,coeff,ix,iy,iz,i0,j0,k0)
                        if (coeff(ix,iy,iz) == SOLID)then
                            coef_lap(ix,iy,iz) = 0.0d0
                            cycle
                        end if
                        coef_lap(ix,iy,iz) =  coef_lap(ix,iy,iz) + ibm%lambda
                    end if

                end do
            end do
        end do

    end subroutine set_ibm_coeff_2nd


    subroutine find_btw_points_x(x,x_n,y,z,ibm,g,coeff,ix,iy,iz,i0,j0,k0)
        implicit none
        real(C_DOUBLE), intent(in)       :: x, y, z,x_n
        type(grid_type),   intent(in)    :: g
        type(ibm_type), intent(inout)    :: ibm
        real(C_DOUBLE)                   :: x_fluid,x_solid,x_int,x_int_0
        integer                          :: i,n_iteration = 50
        real(C_DOUBLE)                   :: x_diff = 0.0d0
        integer,intent(in)               :: ix,iy,iz,i0,j0,k0
        real(C_DOUBLE), intent(inout)    :: coeff(i0:,j0:,k0:)
        real(C_DOUBLE)                   :: eps
        

        ! x is in fluid and x_n is in solid so we need to
        ! find their boundary step by step
        ! for this we will approximate the y_boundary first we will go by mid point
        ! the distance  between points are g%dL so it should be between 0 and dL. 
        ! we first checks if the mid point is in fluid or not
        x_fluid = x;x_solid = x_n
        x_int= 0.0d0; x_int_0 = 0.0d0
        ibm%lambda = 0.0d0
        do i = 1, n_iteration
            ! we save the middle value here
            x_int_0 = x_int
            x_int = real((x_fluid+x_solid)/2.0d0,kind=C_DOUBLE)

            if(isInBody(x_int,y,z,ibm,g))then
                x_solid = x_int
            else 
                x_fluid = x_int
            end if 

        end do
        ! we need to watch for the sign of the lambda
        ! it might amplify the velocities 
        
        ! -----------------!
        ! x diff goes very small numbers this causes to labmda to go inf
        x_diff = abs(x-x_int)   
        if (.not. ieee_is_finite(x_diff))then
            print *, "Warnig! x-IBM coefficient is not a finite number"
        end if
        ! I also added another check if the the distance between velocity point
        ! and boundary is too close it directly marks there as solid. This might be problem?
        ! otherwise either i have too big lambda or infinite lambda
        eps = 1.0d-10*g%dx
        if (x_diff<eps)then
            ibm%lambda= 0.0d0
            coeff(ix,iy,iz) = SOLID
        else
        ibm%lambda= real((1.0d0 / g%dx**2)*((g%dx/x_diff)-1.0d0),kind=C_DOUBLE)
        end if 
        ! we also need to add the 1/dx**2 to the lambda
    end subroutine find_btw_points_x

    ! same as x just names changed
    subroutine find_btw_points_z(z,z_n,y,x,ibm,g,coeff,ix,iy,iz,i0,j0,k0)
        implicit none
        real(C_DOUBLE), intent(in)       :: z, y, x,z_n
        type(grid_type),   intent(in)    :: g
        type(ibm_type), intent(inout)    :: ibm
        real(C_DOUBLE)                   :: z_fluid,z_solid,z_int,z_int_0
        integer                          :: i,n_iteration = 50
        real(C_DOUBLE)                   :: z_diff=0.0d0
        integer,intent(in)               :: ix,iy,iz,i0,j0,k0
        real(C_DOUBLE), intent(inout)    :: coeff(i0:,j0:,k0:)
        real(C_DOUBLE)                   :: eps

        z_fluid = z;z_solid = z_n
        z_int= 0.0d0; z_int_0 = 0.0d0
        ibm%lambda = 0.0d0
        do i = 1, n_iteration
            ! we save the middle value here
            z_int_0 = z_int
            z_int = real((z_fluid+z_solid)/2.0d0,kind=C_DOUBLE)

            if(isInBody(x,y,z_int,ibm,g))then
                z_solid = z_int
            else 
                z_fluid = z_int
            end if 

        end do
        z_diff = abs(z-z_int)
        if (.not.ieee_is_finite(z_diff))then
            print *, "Warnig! z-IBM coefficient is not a finite number"
        end if
        eps = 1.0d-10*g%dz
        if (z_diff<eps)then
            ibm%lambda= 0.0d0
            coeff(ix,iy,iz) = SOLID 
        else
        ibm%lambda= real((1.0d0 / g%dz**2)*((g%dz/z_diff)-1.0d0),kind=C_DOUBLE)
        end if
    end subroutine find_btw_points_z

    ! same as x just names changed
    subroutine find_btw_points_y(y,y_n,x,z,ibm,g,coeff,ix,iy,iz,i0,j0,k0)
        implicit none
        real(C_DOUBLE), intent(in)       :: y, z, x,y_n
        type(grid_type),   intent(in)    :: g
        type(ibm_type), intent(inout)    :: ibm
        real(C_DOUBLE)                   :: y_fluid,y_solid,y_int,y_int_0
        integer                          :: i,n_iteration = 50
        real(C_DOUBLE)                   :: y_diff= 0.0d0
        integer,intent(in)               :: ix,iy,iz,i0,j0,k0
        real(C_DOUBLE), intent(inout)    :: coeff(i0:,j0:,k0:)
        real(C_DOUBLE)                   :: eps

        y_fluid = y;y_solid = y_n
        y_int= 0.0d0; y_int_0 = 0.0d0
        ibm%lambda = 0.0d0
        do i = 1, n_iteration
            y_int_0 = y_int
            y_int = real((y_fluid+y_solid)/2.0d0,kind=C_DOUBLE)

            if(isInBody(x,y_int,z,ibm,g))then
                y_solid = y_int
            else 
                y_fluid = y_int
            end if 

        end do
        y_diff = abs(y-y_int)
        if (.not.ieee_is_finite(y_diff))then
            print *, "Warnig! y-IBM coefficient is not a finite number"
        end if
        eps = 1.0d-10*g%dy
        if (y_diff<eps)then
            ibm%lambda= 0.0d0
            coeff(ix,iy,iz) = SOLID
        else 
        ibm%lambda= real((1.0d0 / g%dy**2)*((g%dy/y_diff)-1.0d0),kind=C_DOUBLE)
        end if
    end subroutine find_btw_points_y

#endif 

    subroutine apply_ibm(field, coeff, g)
        implicit none

        real(C_DOUBLE), intent(inout)   :: field(:,:,:)
        real(C_DOUBLE), intent(in)      :: coeff(:,:,:)
        type(grid_type), intent(in)     :: g

        integer :: ix, iy, iz

        do iz = 1, size(field,3)
       	    do iy = 1, size(field,2)
	            do ix = 1, size(field,1)

                    field(ix,iy,iz) = field(ix,iy,iz) / &
                                    (1.0d0 + g%dt*coeff(ix,iy,iz))

	    end do
	        end do
                end do

end subroutine apply_ibm
end module ibmm
