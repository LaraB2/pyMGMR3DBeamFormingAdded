!------------------------------
Module Atmosphere !------------------------
   use constants, only : dp, pi, rh0
   implicit none
   integer, parameter :: AtmHei_dim=2000
   real(dp), parameter ::  AtmHei_step=10.d0 ! [m]
   real(dp),save :: PenDepth(0:AtmHei_dim+1)
   real(dp),save :: xi(0:AtmHei_dim)
   real(dp),save :: dxi(0:AtmHei_dim)
   real(dp),save :: AirDensity(0:AtmHei_dim)
   real(dp),save :: TopAtmExpon
   Real(dp), parameter :: IndRefOrho=rh0/0.122980582373297d0  ! index of refractivity Over density; ground density=b/c=0.12298058237329704717942323852797
contains
!------------------------------
!------------------------------
Subroutine AtmosphereInit(ZenithAngle_rad, GroundLevel, HoriShwr)
! Adepted from MGMR3D_shower-v5.f90
! Generate refractivity Xi() and PenDepth()
   use constants, only : dp, R_Earth, pi
	implicit none
   Real(dp), intent(in) :: ZenithAngle_rad, GroundLevel
   Logical, intent(in) :: HoriShwr
   Real(dp) :: R_base, Z, Z_p, Height, Height_p, height_h, a,b,c
   integer :: i
   Real(dp) :: D1, D2, Lam, ddd, x_rh, Cos_Zenith, RPenDepth(0:AtmHei_dim+1)
   !Real(dp) :: D1, D2, Lam, ddd, x_rh, TopAtmExpon, ZenithAngle_shower, Cos_Zenith, RPenDepth(0:AtmHei_dim+1)
   R_base=R_Earth + GroundLevel
   RPenDepth(0)=0.  !  Relative penetration depth compared to that at ground-level
   Height=GroundLevel
   z=0.
   write(2,*) 'Generating atmosphere for shower zenith angle [degree]=', ZenithAngle_rad*180/pi
   Cos_Zenith = cos(ZenithAngle_rad)
   do i = 0, AtmHei_dim
     z_p=(i+1)*AtmHei_step  ! distance to ground along shower axis
     !height=z*Cos_Zenith + GroundLevel ! Vertical height above sea level [m]
     Height_p=sqrt(z_p*z_p + 2.d0*R_base*z_p*Cos_Zenith + R_base*R_base) - R_Earth ! Vertical height above sea level [m]
     ! This accounts for the curvature of Earth
     call AtmParams(height,a,b,c)
     AirDensity(i)=(b/c) * exp(-height/c)
     RPenDepth(i+1) = RPenDepth(i) + (b/c)*exp(-(Height+Height_p)/(2.*c)) * AtmHei_step ! assume density is about constant
     ! calculate  averaged refractivity per meter
     if(i.ge.1) xi(i) = IndRefOrho* RPenDepth(i)/z   ! should be independent of ZenithAngle_shower at ground
     if(i.eq.1) xi(0) = xi(1)
     if(i.ge.2) dxi(i-1)= (xi(i)-xi(i-2))/(2.*AtmHei_step)
     Height=Height_p
     z=z_p
   end do
   dxi(0)=dxi(1)
   dxi(AtmHei_dim)=dxi(AtmHei_dim-1)
   ! re-use value for z and Height from previous do-loop and calculate penetration depth for earlier layers
   ! Height=sqrt(z*z + 2.d0*R_base*z*Cos_Zenith + R_base*R_base) - R_Earth
   D1=z*sin(ZenithAngle_rad)  ! horizontal distance to core-on-ground
   D2=R_base+z*Cos_Zenith ! vertical discance to center Earth
   lam=atan2(D1,D2)
   TopAtmExpon=c/cos(ZenithAngle_rad-lam)
   ddd=z/10.   ! large step size for top atmosphere calculation
   X_rh=0 ! Penetration depth from infinity to top atmosphere that is calculated explicitly
   Do
     D1=z*sin(ZenithAngle_rad)  ! horizontal distance to core-on-ground
     D2=R_base+z*Cos_Zenith ! vertical discance to center Earth
     lam=atan2(D1,D2)
     !write(2,*) z/1000,height,lam,X_rh
     if(((ZenithAngle_rad-lam).lt.0.1) .or. (height.gt.1e5)) then ! include remaining tail of atmosphere
         call AtmParams(height,a,b,c)
         !height=sqrt(z*z + 2.d0*R_base*z*cos(ZenithAngle_shower-lam) + R_base*R_base) - R_Earth
         X_rh=X_rh+(a+b* exp(-height/c))/cos(ZenithAngle_rad-lam) ! Flat Earth assumption for remainder
         exit
     else
         z=z+ddd/2. ! to get at half distance
         !height=sqrt(z*z + 2.d0*R_base*z*cos(ZenithAngle_shower-lam) + R_base*R_base) - R_Earth  ! old form which I cannot reproduce
         height_h=sqrt(z*z + 2.d0*R_base*z*Cos_Zenith + R_base*R_base) - R_Earth  ! Height at half-way point
         z=z+ddd/2.   ! to be at top of this layer
         Height_p=sqrt(z*z + 2.d0*R_base*z*Cos_Zenith + R_base*R_base) - R_Earth
         call AtmParams(height_h,a,b,c)  ! assume exponential atmospher with constant coeff over interval
         X_rh=X_rh + b*(exp(-Height/c)-exp(-Height_p/c)) *ddd/(Height_p-Height) ! length times density at half point
         Height=Height_p
     endif
   enddo
   !write(2,*) '!Atmosphere Rpend:', RPenDepth(AtmHei_dim+1), X_rh, RPenDepth(0), RPenDepth(AtmHei_dim)
   If(HoriShwr) X_rh=0.
   do i = 0, AtmHei_dim
     PenDepth(i)=RPenDepth(AtmHei_dim+1)+X_rh - RPenDepth(i)
   enddo
   write(2,*) '!Atmosphere pendepth:', PenDepth(AtmHei_dim), PenDepth(0)
   !
   Return
End Subroutine AtmosphereInit
!----------------------------
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
!----------------------
End Module Atmosphere !------------------------
