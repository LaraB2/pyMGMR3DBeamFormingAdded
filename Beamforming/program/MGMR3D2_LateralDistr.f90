!------------------------------
Module LateralDistribution !------------------------
   use constants, only : dp, pi
   implicit none
   real(dp),save :: Moli_tc
   real(dp),save :: Norm_tc
contains
!------------------------------
Subroutine LateralInit(MoliereRadius, Xindx, AtmHei_step, alpha_tr)
   use constants, only : dp
   Use Pancakefunction, only : fieh
   implicit none
   Real(dp), intent(in) :: MoliereRadius
   Integer, intent(in) :: Xindx
   Real(dp), intent(in) :: AtmHei_step, alpha_tr(*)
   Real(dp) :: s, s_tc, alpha, lam, c, ddd, z, a
   Integer :: CD_i
   !
   Moli_tc=MoliereRadius
   If(Moli_tc.lt.0.) Moli_tc=-Moli_tc  !  ,MoliereRadius iXmx iImx*AtmHei_step/1000.
   If(MoliereRadius.eq.0.0) then
      s=Xindx*AtmHei_step/1000.  ! distance from ground to Xmax along shower axis
      Moli_tc=s*10.
      If(Moli_tc.gt.50) Moli_tc=50
      write(2,*) 'Moli_tc is set automatically to ', Moli_tc, ' m, distance from ground to X_max=',s,' km'
      Moli_tc=50.
      write(2,*) '***********Moli_tc is set back to the default ', Moli_tc, ' m, distance from ground to X_max=',s,' km'
   Endif
    Norm_tc=2.3/Moli_tc  ; s_tc=1.5
    s=0.        ! dummy
    alpha=alpha_tr(Xindx)  ! index pancake fie, at X_max
    lam=alpha_tr(10)      ! index pancake fie, near the ground
    c=Moli_tc/10.
    do CD_i = 1, 200
        ddd=(CD_i-0.5)*c    ! distance to core in [m]
        z=W_tc(ddd,a)
        s=s+z*c
        !b=ddd/1000.     ! distance to showerfront in [mm]
!        write(4,*) d,(1000.*2.5/(moli*2*pi))/((d/Moli+1.)**3.5), &
        !call fiehcore(b,1.d0,lam,nu,a,dfha)  ! @10*dz m
        !call fiehcore(b,1.d0,alpha,fh,dfhl,dfha)    ! @ X_max
        !dfha=dfha *(alpha_tr(iXmx+1)-alpha_tr(iXmx))/AtmHei_step
        ! if(test) write(4,*) ddd,z*100., nu,alpha,fh,dfhl+dfha,dfhl,dfha
    end do
    norm_tc=norm_tc/s   ! Norm for radial function w_tc(r)
    norm_tc=norm_tc/(2.*pi)
    write(2,*) '!Lateral; Moli_tc, norm_tc=', Moli_tc, norm_tc
    Return
End Subroutine LateralInit
!------------------------------
Real*8 Function W_tc(r,dwr)
!  W_tc == W = W_NKG *r
!  dwr == x d(W/x)/dx  !   Wprim= r {d\, w(r)/r \over dr}
!  -----   case  NKG=.true. : NKG for s=2
!    W_tc=Norm_tc*x**(s-1)*(x+1.)**(s-4.5)  ! r*NKG for s
!   for W/x=Norm* x**(s-2)*(x+1.)**(s-4.5)
!        gives dwr= x d(W/x)/dx= x Norm[(s-2)*x**(s-3)*(x+1.)**(s-4.5) + (s-4.5)*x**(s-2)*(x+1.)**(s-5.5)]=
!        thus dwr= x Norm[(s-2)*(x+1.) + (s-4.5)*x]x**(s-3)*(x+1.)**(s-5.5)= [(s-2)/x + (2s-6.5)] w_tc /(x+1.)
!  ----   case  NKG=.false. & tu1=.true. : tu form with t=6, u=1
!    W_tc=Norm_tc*W=x^u/(x+1)^t  ! r* non NKG form t, u
!        Max of W at x=u/(t-u)=a
!        thus set x=a* r/Moli   to have max independent of t and u, normalized for the value for (t=6, u=1) at r=moli
!   for W/x=Norm x^(u-1)/(x+1.)^t
!        gives dwr= x d(W/x)/dx= x Norm[(u-1)/x -t/(x+1)] x^(u-1)/(x+1.)^t= Norm [x(u-1-t)+u-1] w_tc /(x(x+1.))
   real*8 :: r,x,dwr,a
   !common /W_TrCurr/ Norm_tc,Moli_tc,s_tc,dalpha,Nrm_alpha(0:41)
   !real*8 :: Norm_tc,Moli_tc,s_tc,dalpha,Nrm_alpha ! ,alpha
   !Integer, parameter :: NKG=0     ! NKG, s=2 =(t=2.5,u=1)
   !Integer, parameter :: NKG=1      ! =(t=6,u=1)
   !Integer, parameter :: NKG=2      ! =(t=6,u=1.5)
   Integer, parameter :: NKG=3      ! =(t=6,u=2)
    !---------------------------------------------------
    select case (NKG)
    case (1)
      !    W_tc=Norm_tc*x/((x+1.)**4)  ! r* non NKG form t=4, u=1
      !    dwr= -4 *W_tc/(x+1.) !  for t=4, u=1
      a=1.5/(5.*Moli_tc) !  for t=6, u=1  ! factor 1.5 in a to keep the same value for Moli for the NKG (default) case
      x=r*a     !
      W_tc=Norm_tc*x/((x+1.)**6)  ! r* non NKG form t=6, u=1
      dwr= -6 *a*W_tc/(x+1.) !  for t=6, u=1
    case (2)
      a=1.5/(3.*Moli_tc)    ! for  t=6, u=1.5 ; *1.5 to normalize for NKG
      x=r*a
      W_tc=Norm_tc*x*sqrt(x)/((x+1.)**6)  ! r* non NKG form t=6, u=1.5  ! seems to need to increase J0Q
      dwr= a*(0.5/x - 5.5) *W_tc/(x+1.) !  for t=6, u=1.5
    case (3)
      a=1.5/(2.*Moli_tc) ! for  t=6, u=2 ; *1.5 to normalize for NKG
      x=r*a
      W_tc=Norm_tc*x*x/((x+1.)**6)  ! r* non NKG form t=6, u=2  ! seems to have a too small CE force near core
      dwr= a*(1./x-5) *W_tc/(x+1.) !  for t=6, u=2
    case default
      !    W_tc=Norm_tc*sqrt(x)/((x+1.)**3)  ! r*NKG for s=1.5
      !    dwr= (-0.5/x-3.5) *W_tc/(x+1.) !  for s=1.5
      a=1./Moli_tc !  for t=2.5, u=1 ; NKG, s=2
      x=r*a     !
      W_tc=Norm_tc*x/(x+1.)**2.5  ! r*NKG for s=2
      dwr= -2.5*a *W_tc/(x+1.) !  for s=2
      !    W_tc=Norm_tc*x*sqrt(x)/(x+1.)**2  ! r*NKG for s=2.5
      !    dwr= (0.5/x-1.5) *W_tc/(x+1.) !  for s=2.5
    End Select
   !
End function W_tc
!------------------------------
End Module LateralDistribution !------------------------
!-------------------------------------------------------------------------
