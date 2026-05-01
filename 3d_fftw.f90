module fftw_3d
    use, intrinsic :: iso_c_binding
    use :: thomas
    use :: init, only: grid_type, field_type
    implicit none

    include "fftw3.f03"

    !============================================================
    ! Workspace for the FFT-based Poisson solver.
    !
    ! The workspace is allocated once at the beginning of the
    ! simulation and reused at every Poisson solve.
    !
    ! p_hat is used twice:
    !   1. first to store the Fourier-transformed RHS,
    !   2. then overwritten by the Fourier-space pressure solution.
    !============================================================
    type :: poisson_fft_workspace
        complex(C_DOUBLE_COMPLEX), allocatable :: p_hat(:,:,:)
        complex(C_DOUBLE_COMPLEX), allocatable :: plane_hat(:,:)

        real(C_DOUBLE), allocatable :: plane(:,:)

        complex(C_DOUBLE_COMPLEX), allocatable :: a(:), b(:), c(:)
        complex(C_DOUBLE_COMPLEX), allocatable :: d(:), sol(:)

        type(C_PTR) :: plan_fwd = C_NULL_PTR
        type(C_PTR) :: plan_bwd = C_NULL_PTR
    end type poisson_fft_workspace

contains

!============================================================
! Allocate work arrays and create FFTW plans.
! Call this once before the time loop.
!============================================================
subroutine init_poisson_fft_workspace(ws, g)
    type(poisson_fft_workspace), intent(inout) :: ws
    type(grid_type), intent(in) :: g

    allocate(ws%p_hat(g%nx/2+1, g%ny, g%nz))
    allocate(ws%plane_hat(g%nx/2+1, g%nz))
    allocate(ws%plane(g%nx, g%nz))

    allocate(ws%a(g%ny), ws%b(g%ny), ws%c(g%ny))
    allocate(ws%d(g%ny), ws%sol(g%ny))

    ws%plan_fwd = fftw_plan_dft_r2c_2d( &
        g%nz, g%nx, ws%plane, ws%plane_hat, FFTW_ESTIMATE)

    ws%plan_bwd = fftw_plan_dft_c2r_2d( &
        g%nz, g%nx, ws%plane_hat, ws%plane, FFTW_ESTIMATE)

end subroutine init_poisson_fft_workspace

!============================================================
! Solve the 3D Poisson equation with periodic BCs in x,z
! and Neumann BCs in y.
!
! Input:
!   f%rhs  = right-hand side
!
! Output:
!   f%pc   = pressure correction / Poisson solution
!
! The method:
!   1. FFT in x,z for each y-plane
!   2. solve one tridiagonal system in y for each Fourier mode
!   3. inverse FFT back to physical space
!============================================================
subroutine poisson(g, f, ws)
    type(grid_type), intent(in) :: g
    type(field_type), intent(inout) :: f
    type(poisson_fft_workspace), intent(inout) :: ws

    integer :: j, ikx, ikz
    integer :: kx, kz

    real(C_DOUBLE) :: kx_ph, kz_ph, k_tot
    real(C_DOUBLE), parameter :: pi = 3.141592653589793d0

    !--------------------------------------------------------
    ! Forward FFT of RHS.
    !
    ! Each x-z plane is transformed separately.
    ! The result is stored directly in p_hat.
    ! At this stage, p_hat contains rhs_hat.
    !--------------------------------------------------------
    do j = 1, g%ny
        ws%plane = f%rhs(:,j,:)

        call fftw_execute_dft_r2c( &
            ws%plan_fwd, ws%plane, ws%plane_hat)

        ws%p_hat(:,j,:) = ws%plane_hat
    end do
    


    !--------------------------------------------------------
    ! For each Fourier mode (kx,kz), solve the resulting
    ! one-dimensional Poisson problem in y.
    !--------------------------------------------------------
    do ikx = 1, g%nx/2 + 1

        kx = ikx - 1

        ! Modified wavenumber from second-order finite differences
        kx_ph = (4.0d0/g%dx**2) * &
                sin(pi*real(kx,C_DOUBLE)/real(g%nx,C_DOUBLE))**2

        do ikz = 1, g%nz

            ! Map FFT index to physical positive/negative mode number
            if (ikz <= g%nz/2 + 1) then
                kz = ikz - 1
            else
                kz = ikz - g%nz - 1
            end if

            kz_ph = (4.0d0/g%dz**2) * &
                    sin(pi*real(kz,C_DOUBLE)/real(g%nz,C_DOUBLE))**2

            k_tot = kx_ph + kz_ph

            ! Reset tridiagonal coefficients
            ws%a = cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)
            ws%b = cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)
            ws%c = cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)
            ws%d = cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)

            ! RHS for this Fourier mode
            ws%d(:) = ws%p_hat(ikx,:,ikz)

            if (.not. (kx == 0 .and. kz == 0)) then

                !------------------------------------------------
                ! Non-zero Fourier mode:
                !
                ! Solve:
                !   d2(p_hat)/dy2 - k_tot p_hat = rhs_hat
                !
                ! with Neumann BCs at bottom and top.
                !------------------------------------------------

                ! Bottom boundary
                ws%b(1) = cmplx(-1.0d0/g%dy**2 - k_tot, 0.0d0, kind=C_DOUBLE_COMPLEX)
                ws%c(1) = cmplx( 1.0d0/g%dy**2,         0.0d0, kind=C_DOUBLE_COMPLEX)

                ! Interior points
                ws%a(2:g%ny-1) = cmplx( 1.0d0/g%dy**2,        0.0d0, kind=C_DOUBLE_COMPLEX)
                ws%b(2:g%ny-1) = cmplx(-2.0d0/g%dy**2-k_tot, 0.0d0, kind=C_DOUBLE_COMPLEX)
                ws%c(2:g%ny-1) = cmplx( 1.0d0/g%dy**2,        0.0d0, kind=C_DOUBLE_COMPLEX)

                ! Top boundary
                ws%a(g%ny) = cmplx( 1.0d0/g%dy**2,        0.0d0, kind=C_DOUBLE_COMPLEX)
                ws%b(g%ny) = cmplx(-1.0d0/g%dy**2-k_tot, 0.0d0, kind=C_DOUBLE_COMPLEX)

            else

                !------------------------------------------------
                ! Zero Fourier mode:
                !
                ! The pure Neumann Poisson problem is singular.
                ! We fix one pressure value to remove the nullspace.
                !------------------------------------------------

                ws%b(1) = cmplx(1.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)
                ws%d(1) = cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)

                ! Interior points
                ws%a(2:g%ny-1) = cmplx( 1.0d0/g%dy**2, 0.0d0, kind=C_DOUBLE_COMPLEX)
                ws%b(2:g%ny-1) = cmplx(-2.0d0/g%dy**2, 0.0d0, kind=C_DOUBLE_COMPLEX)
                ws%c(2:g%ny-1) = cmplx( 1.0d0/g%dy**2, 0.0d0, kind=C_DOUBLE_COMPLEX)

                ! Top Neumann boundary
                ws%a(g%ny) = cmplx( 1.0d0/g%dy**2, 0.0d0, kind=C_DOUBLE_COMPLEX)
                ws%b(g%ny) = cmplx(-1.0d0/g%dy**2, 0.0d0, kind=C_DOUBLE_COMPLEX)

            end if

            ! Solve tridiagonal system for this Fourier mode
            call thomas_solver(ws%a, ws%b, ws%c, ws%d, ws%sol, g%ny)

            ! Overwrite RHS transform with pressure transform
            ws%p_hat(ikx,:,ikz) = ws%sol(:)

        end do
    end do


    !--------------------------------------------------------
    ! Inverse FFT.
    !
    ! Transform each y-plane of p_hat back to physical space.
    ! FFTW does not normalize inverse transforms, hence division
    ! by nx*nz.
    !--------------------------------------------------------
    do j = 1, g%ny
        ws%plane_hat = ws%p_hat(:,j,:)

        call fftw_execute_dft_c2r( &
            ws%plan_bwd, ws%plane_hat, ws%plane)

        f%pc(1:g%nx,j,1:g%nz) = ws%plane / real(g%nx*g%nz, C_DOUBLE)
    end do

end subroutine poisson

!============================================================
! Destroy FFTW plans and deallocate work arrays.
! Call this once after the time loop.
!============================================================
subroutine destroy_poisson_fft_workspace(ws)
    type(poisson_fft_workspace), intent(inout) :: ws

    if (c_associated(ws%plan_fwd)) call fftw_destroy_plan(ws%plan_fwd)
    if (c_associated(ws%plan_bwd)) call fftw_destroy_plan(ws%plan_bwd)

    if (allocated(ws%p_hat))    deallocate(ws%p_hat)
    if (allocated(ws%plane_hat)) deallocate(ws%plane_hat)
    if (allocated(ws%plane))    deallocate(ws%plane)

    if (allocated(ws%a))   deallocate(ws%a)
    if (allocated(ws%b))   deallocate(ws%b)
    if (allocated(ws%c))   deallocate(ws%c)
    if (allocated(ws%d))   deallocate(ws%d)
    if (allocated(ws%sol)) deallocate(ws%sol)

    ws%plan_fwd = C_NULL_PTR
    ws%plan_bwd = C_NULL_PTR

end subroutine destroy_poisson_fft_workspace

end module fftw_3d
