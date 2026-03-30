module Atmosphere
    use constants, only : dp
    implicit none

    real(dp), save :: rh0_val
    real(dp), save :: GroundLevel_val
    integer, save  :: n_layers
    real(dp), allocatable, save :: a_arr(:), b_arr(:), c_arr(:), h_top(:)

contains

    subroutine setup_atmosphere(filename, g_level)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: g_level
        integer :: unit_atm, i
        
        GroundLevel_val = g_level
        open(newunit=unit_atm, file=filename, status='old', action='read')
        read(unit_atm, *) rh0_val
        read(unit_atm, *) n_layers
        
        if (allocated(a_arr)) deallocate(a_arr, b_arr, c_arr, h_top)
        allocate(a_arr(n_layers), b_arr(n_layers), c_arr(n_layers), h_top(n_layers))
        
        do i = 1, n_layers
            read(unit_atm, *) a_arr(i), b_arr(i), c_arr(i), h_top(i)
        end do
        close(unit_atm)
    end subroutine setup_atmosphere

    subroutine AtmParams(h_in, a, b, c)
        real(dp), intent(in) :: h_in
        real(dp), intent(out) :: a, b, c
        integer :: j

        ! Default to the topmost layer immediately.
        ! If h_in is higher than all defined layers, it stays as the top layer.
        a = a_arr(n_layers)
        b = b_arr(n_layers)
        c = c_arr(n_layers)

        ! Loop through to find if it fits in a lower layer
        do j = 1, n_layers
            if (h_in <= h_top(j)) then
                a = a_arr(j)
                b = b_arr(j)
                c = c_arr(j)
                exit
            end if
        end do
    end subroutine AtmParams
    
end module Atmosphere
