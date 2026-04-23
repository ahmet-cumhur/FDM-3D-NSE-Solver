module thomas
    use, intrinsic :: iso_c_binding
    implicit none 
    
    contains
    subroutine thomas_solver(a, b, c, d, x, n)
        implicit none
        integer, intent(in) :: n
        complex(C_DOUBLE_COMPLEX), intent(in)  :: a(n), b(n), c(n), d(n)
        complex(C_DOUBLE_COMPLEX), intent(out) :: x(n)

        complex(C_DOUBLE_COMPLEX) :: cp(n), dp(n), denom
        integer :: m

        cp(1) = c(1) / b(1)
        dp(1) = d(1) / b(1)

        do m = 2, n
            denom = b(m) - a(m)*cp(m-1)
            if (m < n) cp(m) = c(m) / denom
            dp(m) = (d(m) - a(m)*dp(m-1)) / denom
        end do

        x(n) = dp(n)
        do m = n-1, 1, -1
            x(m) = dp(m) - cp(m)*x(m+1)
        end do
    end subroutine thomas_solver
end module thomas