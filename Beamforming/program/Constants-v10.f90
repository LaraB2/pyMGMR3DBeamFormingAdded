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
