!--------------------------!
!                          !
!     Initialisation       !
!         Module           !
!                          !
!--------------------------! 
! 
! authors: Dr.-Ing. Davide Gatti
!          B.Sc. Ahmet Cumhur
! 
! date:    28.04.26
! 

module init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Grid datatype
    ! ------------------------
    type :: grid_type
        integer :: nx, ny, nz
        real(C_DOUBLE) :: lx, ly, lz
        real(C_DOUBLE) :: dx, dy, dz
        real(C_DOUBLE) :: re, dt, t_final, t_current, cfl, cflmax, dtmax
    end type grid_type

    ! Flow field datatype
    ! ------------------------
    type :: field_type
        ! we add some body force
        real(C_DOUBLE) :: b_x,b_y,b_z
        real(C_DOUBLE), allocatable :: un(:,:,:), us(:,:,:)
        real(C_DOUBLE), allocatable :: vn(:,:,:), vs(:,:,:)
        real(C_DOUBLE), allocatable :: wn(:,:,:), ws(:,:,:)
        real(C_DOUBLE), allocatable :: pn(:,:,:), pc(:,:,:)
        real(C_DOUBLE), allocatable :: rhs(:,:,:)
        real(C_DOUBLE), allocatable :: uc(:,:,:), vc(:,:,:), wc(:,:,:)
    end type field_type

! -----------------------------
contains
! -----------------------------


! Grid and parameter initialisation
! -----------------------------
subroutine init_grid(g)
    type(grid_type), intent(inout) :: g
    integer :: i,j,k

#ifdef USE_IBM  
    g%nx = 50; g%ny = 50; g%nz = 20
#elif USE_IBM_G
    g%nx = 50; g%ny = 50; g%nz = 20
#else
    g%nx = 32; g%ny = 32; g%nz = 10
#endif

    g%lx = 1.0d0
    g%ly = 1.0d0
    g%lz = 1.0d0

    g%dx = g%lx / real(g%nx, C_DOUBLE)
    g%dy = g%ly / real(g%ny, C_DOUBLE)
    g%dz = g%lz / real(g%nz, C_DOUBLE)
    
    g%re = 100.0d0

#ifdef USE_IBM
    g%dt = 1.0d-4
#elif USE_IBM_G
    g%dt = 1.0d-3
#else
    g%dt = 1.0d-3
#endif
    g%cflmax = 0.1d0
    g%dtmax = 1.0d-3

    g%t_final = g%dt*10000
    g%t_current = 0.0d0

end subroutine init_grid

! Flow field initialisation
! -----------------------------
subroutine init_field(f, g)
    type(field_type), intent(inout) :: f
    type(grid_type), intent(in)     :: g

    allocate(f%un(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%us(0:g%nx+1,0:g%ny+1,0:g%nz+1))
    allocate(f%vn(0:g%nx+1,1:g%ny+1,0:g%nz+1), f%vs(0:g%nx+1,1:g%ny+1,0:g%nz+1))
    allocate(f%wn(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%ws(0:g%nx+1,0:g%ny+1,0:g%nz+1))

    allocate(f%pn(0:g%nx+1,1:g%ny,0:g%nz+1), f%pc(0:g%nx+1,1:g%ny,0:g%nz+1))
    allocate(f%rhs(1:g%nx,1:g%ny,1:g%nz))

    allocate(f%uc(1:g%nx,1:g%ny,1:g%nz), f%vc(1:g%nx,1:g%ny,1:g%nz), f%wc(1:g%nx,1:g%ny,1:g%nz))

    f%un = 0.0d0; f%us = 0.0d0
    f%vn = 0.0d0; f%vs = 0.0d0
    f%wn = 0.0d0; f%ws = 0.0d0
    f%pn = 0.0d0; f%pc = 0.0d0
    f%rhs = 0.0d0
    f%uc = 0.0d0; f%vc = 0.0d0; f%wc = 0.0d0
    ! we use body forces as driving force
    f%b_x = 1.0d0;f%b_y= 0.1d0;f%b_z= 0.0d0

end subroutine init_field


end module init
