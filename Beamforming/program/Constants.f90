     module constants
      integer, parameter:: dp=kind(0.d0)
      COMPLEX(dp), parameter :: CI=CMPLX(0._dp,1._dp)
      REAL(dp), parameter :: pi=2._dp *asin(1._dp)  ! This is the real constant pi
      REAL(dp), parameter :: c_l=0.299792458d0     ! speed of light in Vacuum [m/nsec]
      Real(dp), parameter :: R_Earth=6.3781d6 ! [m]
      real(dp), parameter :: rh0=0.000292d0  !  Refractivity ar ground level
     end module constants
!-----------------------------------------
!
subroutine AtmParams(height,a,b,c)
    use constants, only : dp,pi
    implicit none
    real(dp), intent(in) :: height
    real(dp), intent(out) :: a,b,c
    if (height.ge.10e3) then
         a = 0.61289; b = 1305.5948; c = 6361.4304
    elseif (height.ge.4e3) then
         a = -94.919d0; b = 1144.9069d0; c = 8781.5355d0
    else
         a = -186.555305d0; b = 1222.6562d0; c = 9941.8638d0
    endif
end subroutine AtmParams
