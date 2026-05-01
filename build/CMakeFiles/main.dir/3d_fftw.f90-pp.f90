# 1 "C:/Users/mehme/OneDrive/Desktop/thesis/nse/3d_fftw.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "C:/Users/mehme/OneDrive/Desktop/thesis/nse/3d_fftw.f90"
module fftw_3d
    use,intrinsic :: iso_c_binding
    use :: thomas
    implicit none
    include "fftw3.f03"
    ! here we solve the 3d poison equation w/ 2 periodic BC.
    contains
        subroutine poison_fft_3d (x,y,z,nx,ny,nz,lx,ly,lz,dx,dy,dz,rhs,pc)
            use, intrinsic :: iso_c_binding
            implicit none
            ! we get the grid data first
            real(C_DOUBLE), intent(in) :: x(nx),y(ny),z(nz) 
            integer, intent(in) :: nx,ny,nz
            real(C_DOUBLE), intent(in) :: lx,ly,lz 
            real(C_DOUBLE), intent(in) :: dx,dy,dz
            real(C_DOUBLE), intent(in) :: rhs(nx,ny,nz)
            real(C_DOUBLE),intent(out) :: pc(nx,ny,nz)
            
            integer :: i,j,k,ikx,ikz,kx,kz
            real(C_DOUBLE) :: kx_ph,kz_ph,k_tot
            real(C_DOUBLE), parameter :: pi = 3.141592653589793d0
            ! we put planes(x,z) into fft 
            complex(C_DOUBLE_COMPLEX) :: p_hat(nx/2+1,ny,nz),rhs_hat(nx/2+1,ny,nz),plane_hat(nx/2+1,nz)
            real(C_DOUBLE) :: plane_in(nx,nz),plane_out(nx,nz)
            ! we define the tridiagonal matrix coefficients here
            ! for each mode we solve a tridiagonal matrix sytem
            ! since its much easier than solving the matrix system direcly/iteratively
            ! this is where we gain time 
            complex(C_DOUBLE_COMPLEX) :: a(ny), b(ny), c(ny), d(ny), sol(ny)
            

            type(C_PTR) :: plan_fwd,plan_bwd

            plan_fwd = fftw_plan_dft_r2c_2d(nz,nx,plane_in,plane_hat,FFTW_ESTIMATE)
            plan_bwd = fftw_plan_dft_c2r_2d(nz,nx,plane_hat,plane_out,FFTW_ESTIMATE)
            
            ! lets do the forward transformation 
            do j = 1, ny
                plane_in = rhs(:,j,:)
                call fftw_execute_dft_r2c(plan_fwd,plane_in,plane_hat)
                rhs_hat(:,j,:) = plane_hat
            end do 
            ! we have the transformed rhs
            do ikx = 1, nx/2+1
                kx = ikx-1
                ! we find the actual kx and kz from the start
                ! the kx and kz definitions are got 
                kx_ph = (4.0d0/dx**2)*sin(pi*real(kx,C_DOUBLE)/real(nx,C_DOUBLE))**2
                do ikz = 1,nz
                    ! we need to do this in order to count for the negative modes
                    
                    ! careful here problem might occur
                    if (ikz <= nz/2+1) then
                        kz = ikz-1
                    else 
                        kz = ikz-nz-1
                    end if 
                    kz_ph = (4.0d0/dz**2)*sin(pi*real(kz,C_DOUBLE)/real(nz,C_DOUBLE))**2
                    ! we combine both kz and kx here
                    k_tot = kz_ph+kx_ph
                    
                    ! now lets define the a,b,c,d
                    a = (0.0d0,0.0d0)
                    b = (0.0d0,0.0d0)
                    c = (0.0d0,0.0d0)
                    d = (0.0d0,0.0d0)

                    do j = 1,ny
                        d(j) = rhs_hat(ikx,j,ikz)
                    end do

                    ! we need to check for 0 mode condition
                    if (.not.(kx == 0 .and. kz == 0)) then
                        ! lets add the bottom row first
                        ! we apply the  neuman BC here
                        b(1) = cmplx((-1.0d0/dy**2-k_tot),0.0d0,kind=C_DOUBLE_COMPLEX)
                        c(1) = cmplx((1.0d0/dy**2),0.0d0,kind=C_DOUBLE_COMPLEX)
                        
                        ! careful only interior points here
                        do j = 2, ny-1
                            a(j) = cmplx(1.0d0/dy**2,0.0d0,kind=C_DOUBLE_COMPLEX)
                            b(j) = cmplx(-2.0d0/dy**2-k_tot,0.0d0,kind=C_DOUBLE_COMPLEX)
                            c(j) = cmplx(1.0d0/dy**2,0.0d0,kind=C_DOUBLE_COMPLEX)
                        end do 
                        
                        !top row now
                        b(ny) = cmplx(-1.0d0/dy**2-k_tot,0.0d0)
                        a(ny) = cmplx(1.0d0/dy**2,0.0d0)

            
                    else 
                        ! fill here the 0 mode condition
                        ! we fix the pressure to zero in order to
                        ! get rid of singularity

                        b(1) = cmplx(1.0d0,0.0d0,kind=C_DOUBLE_COMPLEX)
                        d(1) = cmplx(0.0d0,0.0d0,kind=C_DOUBLE_COMPLEX)

                        do j = 2,ny-1
                            a(j) = cmplx( 1.0d0/dy**2, 0.0d0)
                            b(j) = cmplx(-2.0d0/dy**2, 0.0d0)
                            c(j) = cmplx( 1.0d0/dy**2, 0.0d0)
                        end do

                        ! we apply the  neuman BC here
                        a(ny) = cmplx( 1.0d0/dy**2, 0.0d0)
                        b(ny) = cmplx(-1.0d0/dy**2, 0.0d0)
                        
                    end if
                    ! now we solve the tridiagonal matrix
                    call thomas_solver(a,b,c,d,sol,ny)
                    do j = 1,ny
                        p_hat(ikx,j,ikz) = sol(j)
                    end do 
                end do
            end do
            do j = 1,ny
                plane_hat = p_hat(:,j,:)
                call fftw_execute_dft_c2r(plan_bwd,plane_hat,plane_out)
                pc(:,j,:) = plane_out/real(nx*nz,C_DOUBLE)
            end do    
            call fftw_destroy_plan(plan_fwd)
            call fftw_destroy_plan(plan_bwd)
        end subroutine poison_fft_3d
end module fftw_3d
