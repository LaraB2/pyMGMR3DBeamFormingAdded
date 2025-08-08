!---------------------------------------------------------
    Subroutine Initialize_shower
! v2: improved calculation of penetration depth for large zenith angles.
! v3: Changed rho dependence in CX and TC to improve norm factor
    use RFootPars, only : rh0,X_0,lamx,X_max,lam_tc,lam_100,XDepAlpha,J0t,J0Q, IntegrateCurrent
    use RFootPars, only : X_02, lamx2, X_max2, Energy_sh2
!            NPart=NPart+Energy_sh2*( (X_rh-X_02)/(X_max2-X_02))**((X_max2-X_02)/lamx2) * exp((X_max2-X_rh)/lamx2) !&
    use RFootPars, only : R_0,L_0, R_02,L_02, RL_param, MoliereRadius
    use RFootPars, only : alpha_frc0, alpha_vB, u0, F_over_beta, a_ChX, D_ESmooth, AlternativeSmooth
    use RFootPars, only : step,stpv, N_frc,h_frc,Force,alpha_frc,N_step_max, PenFacHeight
    use RFootPars, only : line,  h_frcL,ForceL,alpha_frcL, test
    use RFootPars, only : GroundLevel, PancakeIncField, vDrift2, SelectFh
    use BigArrays, only : IndRefOrho,TopAtmExpon
    use BigArrays, only : CoreDist_Dim, CoreDist_A, ObsDist_dim, ObsDist_Step, AtmHei_step, AtmHei_dim
    use BigArrays, only : PancakeThi,xi,dxi,PenDepth,Ix,Iy,IQ,IX_int,Iy_int,alpha_tr
    use constants, only : dp,pi
    use eventdata, only : Fitting,ZenithAngle_shower, F_max, Energy_sh, HoriShwr,FileShCurrent, ShPlane
    use eventdata, only : Voltages, Core_N, Core_E, Ant_N, Ant_E,Eoff,Noff, RelMx_N, RelMx_E, RelMx_U
    use eventdata, only : vBE,vBN,vBU,vvBE,vvBN,vvBU, vBxvvB, SN, SE, SU, sZS, Ang_Ux, alpha_Bz
   use RFootPars, only : FShift_x,FShift_y, D_IMax
    !use CrossProd,only : calc_alpha_vB
    !use CrossProd,only : calc_alpha_Bz
    implicit none
    real(dp) :: Force_x(0:AtmHei_dim),Force_y(0:AtmHei_dim),height, alpha_E, Force0, sin_alpha,cos_alpha
    Real(dp) :: Height_p, z_p, Height_h
    real(dp) :: RPenDepth(0:AtmHei_dim+1)
    real(dp) :: z,s,nu,a,b,c,ddd,X_rh,dfha,dfhl,fh,lam,r1,r2,dh ,D1,D2, NPart,NPart2,XmaxAve, ux,uy,BoF
    !real(dp) :: Eoff,Noff,vBE,vBN,vBU,vvBE,vvBN,vvBU,BoF
    Real(dp), save :: Xb_0=50, Xc_0=50     ! used in the parametrization of charge excess and drift velocity
    real(dp), parameter :: X_EM=36.7d0
    integer :: i,j,CD_i,iXmx, I1,I2
    common /W_TrCurr/ Norm_tc,Moli_tc,s_tc,dalpha,Nrm_alpha(0:41)
    real*8 :: Norm_tc, Moli_tc, s_tc, dalpha, Nrm_alpha, alpha, X_start, XXX, b_ground, H_s, NP_F
    real*8 :: W_BxvxB  ! weight factor for BxvxB Lorentz force
    real*8 :: W_tc, Cos_Zenith,AirDensity(0:AtmHei_dim),R_base,R_Earth=6.3781d6 ! [m]
    real(dp) :: Xmax_0,Isqr_max,Isqr, h_f(0:N_step_max),F(0:N_step_max),ang(0:N_step_max),dfx,dfy
    Real(dp) :: SigX0, SigX1, SigX2, FF_tc ! (FF=Fudge Factor); testing different longtidudinal functions
    integer :: iImx, indx(1:N_step_max), iQmx
    logical,save :: FiPa=.true.
    external W_tc
!
    write(2,201) rh0*10000,Moli_tc, J0Q, SelectFh, lam_tc, lam_100, XDepAlpha
201 format('Refractivity at sea level=',F8.3,'x 10^{-4}',&
      ', Moliere radius=',f6.1,'[m], J0Q=',f7.4 ,/, &
      'f(h) shape',i2,' with L=',f4.2,'[m]' , &
      ', L@100m=',f5.2,'[m]', ' index for X-dependence =',f4.2)
      XXX=16.
!
!       J0t[keV/m] = c[m/s]  * B[T} /1000 = 30 *B[G]
    sin_alpha = sin(alpha_vB*pi/180.d0)
    cos_alpha = cos(alpha_vB*pi/180.d0)
    Force0=J0t*sin_alpha ; alpha_frc0=0.d0               !   acts from max height down
    Cos_Zenith=cos(ZenithAngle_shower)
!    Force(1)=222.1 ; alpha_frc(1)=pi/2 ; h_frc(1)=4000 ; N_frc=0    !   acts from height(1) down
!
    dalpha=1/10.   ! Normalization of pancake function
    Nrm_alpha(:)=1.d0 !  just for set-up
    PancakeThi(:)=1.d0
    l_alpha: Do i=0,40
      alpha=1.d0+i * dalpha
      z=0. ; s=0. ; nu=0.
      l_idi: do j=1,10000
        a=j*alpha/40.d0    ! distance
        call fieh(a,1,alpha,ddd,b,c)  !    subroutine fieh(h,CoreDist_i,alpha,fh,dfhl,dfha)
        z=z+ddd*alpha/40.d0
      enddo l_idi
      xi(i)=1.d0/z
      if(test) write(2,*) 'alpha, nrm',alpha,1./z
    enddo l_alpha
    Do i=0,41
        Nrm_alpha(i)=xi(i)
    enddo
!
!   Initialize_shower; Radial dependence of pancake thickness
    b=(ObsDist_dim)*ObsDist_Step/100.
    c=lam_100*b/lam_tc
    a=c**(1./(CoreDist_Dim))
    b=(ObsDist_dim)*ObsDist_Step/c    ! =100. * lam_tc/lam_100  !
!    b=Line2Obsrv(ObsDist_dim)/c     as alternativenew def
    CoreDist_A(0)=b    ; PancakeThi(0)=lam_tc
!    write(2,*) 'b=',b
    Do CD_i=1,CoreDist_Dim  ! calculate pancake thickness-parameter as function of distance to shower-core
      PancakeThi(CD_i)=PancakeThi(CD_i-1)*a  ! Later this will be multiplied by alpha
      CoreDist_A(CD_i)=PancakeThi(CD_i)*b/lam_tc
    enddo
    if(CoreDist_A(CoreDist_Dim).lt.ObsDist_dim*ObsDist_Step) CoreDist_A(CoreDist_Dim)=ObsDist_dim*ObsDist_Step
    write(2,"('parametrization radial dependence of \lambda; @r=0:',F5.2,', @r=',f4.1,'m:',F5.2,', @r=',f6.1,'m:',F6.1)") &
        lam_tc,CoreDist_A(1),PancakeThi(1),ObsDist_dim*ObsDist_Step, c*lam_tc
    if(test) write(2,"('in ',I2,' steps of factor',f5.2,' in distance and',F5.2,' in lambda')") CoreDist_Dim,b,a
    !
    if(test) OPEN(UNIT=4,STATUS='unknown',FILE='plot/sh_lambda.dat')
    if(test) write(4,"('!',7x,'CoreDist [m],',18x,'lambda [m]')")
    Do CD_i=0,CoreDist_Dim  ! calculate pancake thickness-parameter as function of distance to shower-core
      If(.not. Fitting) write(2,*) 'CoreDistance=',CoreDist_A(CD_i),'[m], lambda=',PancakeThi(CD_i)
      if(test) write(4,*) CoreDist_A(CD_i),PancakeThi(CD_i)
    enddo
    if(test) close(unit=4)
    !
    !       Calculate things as function of height above ground @ GroundLevel
    ! xi(0)=rh0 If at sea level
    IndRefOrho=rh0/0.122980582373297d0  ! index of refraction Over density; ground density=b/c=0.12298058237329704717942323852797
    R_base=R_Earth + GroundLevel
    RPenDepth(0)=0.  !  Relative penetration depth compared to that at ground-level
    Height=GroundLevel
    z=0.
    do i = 0, AtmHei_dim
        !z=i*AtmHei_step  ! distance to ground along shower axis
        z_p=(i+1)*AtmHei_step  ! distance to ground along shower axis
        !height=z*Cos_Zenith + GroundLevel ! Vertical height above sea level [m]
        Height_p=sqrt(z_p*z_p + 2.d0*R_base*z_p*Cos_Zenith + R_base*R_base) - R_Earth ! Vertical height above sea level [m]
        ! This accounts for the curvature of Earth
        !write(2,*) height-(z*Cos_Zenith + GroundLevel)
        call AtmParams(height,a,b,c)
        AirDensity(i)=(b/c) * exp(-height/c)
        !If(HoriShwr) then
        !    ! AirDensity=(b/c) * exp(-height/c)
        !    X_rh = AirDensity(i) * (AtmHei_dim-i)*AtmHei_step
        !else
        !    X_rh=(a+b* exp(-height/c))/Cos_Zenith ! Penetration depth for flat Earth
        !endif
        !PenDepth(i)=X_rh
        !
        !D1=(z +R_base*Cos_Zenith)/(height+R_Earth)  ! = Cos_Zenith in the Earth curvature is ignored
        RPenDepth(i+1) = RPenDepth(i) + (b/c)*exp(-(Height+Height_p)/(2.*c)) * AtmHei_step ! assume density is about constant
        !write(2,*) height, a, PenDepth(0)-X_rh-RPenDepth(i)
        !
        ! calculate  averaged refractivity per meter
        if(i.ge.1) xi(i) = IndRefOrho* RPenDepth(i)/z   ! should be independent of ZenithAngle_shower at ground
        !if(i.ge.1) xi(i) = IndRefOrho* (PenDepth(0)-X_rh)/z   ! should be independent of ZenithAngle_shower at ground
        if(i.eq.1) xi(0) = xi(1)
        if(i.ge.2) dxi(i-1)= (xi(i)-xi(i-2))/(2.*AtmHei_step)
        Height=Height_p
        z=z_p
    end do
    dxi(0)=dxi(1)
    dxi(AtmHei_dim)=dxi(AtmHei_dim-1)
    ! re-use value for z and Height from previous do-loop and calculate penetration depth for earlier layers
    ! Height=sqrt(z*z + 2.d0*R_base*z*Cos_Zenith + R_base*R_base) - R_Earth
        D1=z*sin(ZenithAngle_shower)  ! horizontal distance to core-on-ground
        D2=R_base+z*cos(ZenithAngle_shower) ! vertical discance to center Earth
        lam=atan2(D1,D2)
    TopAtmExpon=c/cos(ZenithAngle_shower-lam)
    ddd=z/10.   ! step size for top atmosphere calculation
    X_rh=0 ! Penetration depth from infinity to top atmosphere that is calculated explicitly
    Do
        D1=z*sin(ZenithAngle_shower)  ! horizontal distance to core-on-ground
        D2=R_base+z*cos(ZenithAngle_shower) ! vertical discance to center Earth
        lam=atan2(D1,D2)
        !write(2,*) z/1000,height,lam,X_rh
        if(((ZenithAngle_shower-lam).lt.0.1) .or. (height.gt.1e5)) then ! include remaining tail of atmosphere
            call AtmParams(height,a,b,c)
            !height=sqrt(z*z + 2.d0*R_base*z*cos(ZenithAngle_shower-lam) + R_base*R_base) - R_Earth
            X_rh=X_rh+(a+b* exp(-height/c))/cos(ZenithAngle_shower-lam) ! Flat Earth assumption for remainder
            exit
        else
            z=z+ddd/2. ! to get at half distance
            !height=sqrt(z*z + 2.d0*R_base*z*cos(ZenithAngle_shower-lam) + R_base*R_base) - R_Earth  ! old form which I cannot reproduce
            height_h=sqrt(z*z + 2.d0*R_base*z*Cos_Zenith + R_base*R_base) - R_Earth
            z=z+ddd/2.   ! to be at top of this layer
            Height_p=sqrt(z*z + 2.d0*R_base*z*Cos_Zenith + R_base*R_base) - R_Earth
            call AtmParams(height_h,a,b,c)  ! assume exponential atmospher with constant coeff over interval
            X_rh=X_rh + b*(exp(-Height/c)-exp(-Height_p/c)) *ddd/(Height_p-Height) ! length times density at half point
            Height=Height_p
        endif
    enddo
   X_start=X_0 ! The X value where the shower has non-zero particle number
   If(RL_param) then
      L_0=X_0     !needed when fitting
      R_0=lamx
      X_start=X_max-L_0/R_0
!      L_0=R_0*(X_max+X_0)  !v3c:   X_0 is negative here!
      If(Energy_sh2.le.0.) Then
         write(2,*) 'Single bump LDF; Energy_sh, L_0, R_0, X_0=',Energy_sh, L_0, R_0, X_start
      Else
         write(2,*) 'Double bump LDF; Energy_sh, L_0, R_0, X_0=', Energy_sh, L_0, R_0, X_start, &
                ', second bump:', Energy_sh2, R_02,L_02, X_max2
      EndIf
   EndIf
    If(HoriShwr) then
        If(FiPa) write(2,*) 'AirDensity=',AirDensity(1),' , Estimated penetration depth offset at topheight=',X_rh,' is set to zero'
        X_rh=0.
    else
        If(FiPa) write(2,"('Estimated penetration depth offset at topheight=',g11.4,', air density at ground=',g11.4)") &
                    X_rh,AirDensity(1)
        If(X_start.lt.X_rh) then
            write(2,*) '!!! Warning, X_0=',X_start,' falls above topheight !!!'
        endif
    endif
    do i = 0, AtmHei_dim
        PenDepth(i)=RPenDepth(AtmHei_dim+1)+X_rh - RPenDepth(i)
    enddo
   !
    If(HoriShwr) then
        N_frc=0     ! no electric field for horizontal showers
    endif
    Force_x(:)=Force0*cos(alpha_frc0) ; Force_y(:)=Force0*sin(alpha_frc0)
    PenFacHeight=1. ! Penalty factor when two heights are too close
    b_ground=D_ESmooth*(PenDepth(0)-PenDepth(1))/(AtmHei_step*X_EM)  !
    !b_ground=0.
    If(N_frc.ge.1) then    ! skip when there is no electric field given
      if(step .or. stpv) then
        !call calc_alpha_Bz(alpha_Bz)
        If(FiPa) write(2,"('Lorentz force:',g11.4,'@ ',g11.4,'deg')") Force0,alpha_frc0
        F(0)=0.d0
        ang(0)=0.d0
        write(2,*) ' h_step[km], E_step[kV/m], a_step[deg];' &
                    ,' h_true[km], E_true[kV/m], a_true[deg]; ' &
                    ,'E_vz[kV/m], E_vvz[kV/m]; smoothing distance=',D_ESmooth,', u0=',u0
        if(step) then
            call step2stpv(h_frc,Force,alpha_frc,h_f,F,ang)
            Do i=1,N_frc
                write(2,"(2x,f8.3,4x,f9.3,6x,f7.1,3x,'; ',f8.3,4x,f9.3,6x,f7.1,3x,';',f9.3,4x,f9.3)") &
                  h_frc(i)/1000.,Force(i),alpha_frc(i)*180./pi , h_f(i)/1000.,F(i),ang(i)*180./pi &
                  , F(i)*cos(ang(i)-alpha_Bz), F(i)*sin(ang(i)-alpha_Bz)
                !write(2,"(f6.3,f8.3,f6.1,f6.3,f8.3,f6.1)") &
                !  h_frc(i)/1000.,Force(i),alpha_frc(i)*180./pi, h_f(i)/1000.,F(i),ang(i)*180./pi &
                !  , Force(i)*cos(alpha_frc(i)-alpha_Bz), Force(i)*sin(alpha_frc(i)-alpha_Bz)
            enddo
            Do i=1,N_frc
                h_f(i)=h_frc(i)
                F(i)= Force(i)
                ang(i)=alpha_frc(i)
            enddo
        else
            call stpv2step(h_frc,Force,alpha_frc,h_f,F,ang)
            Do i=1,N_frc
                !write(2,"(' h=',f6.3,'[km], E_step=',f8.3,'[kV/m] @',f7.1,'[deg] ; h=', &
                !  f6.3,'[km], E_true=',f8.3,'[kV/m] @',f7.1,'[deg]')") &
                write(2,"(2x,f8.3,4x,f9.3,6x,f7.1,3x,'; ',f8.3,4x,f9.3,6x,f7.1,3x,';',f9.3,4x,f9.3)") &
                  h_f(i)/1000.,F(i),ang(i)*180./pi , h_frc(i)/1000.,Force(i),alpha_frc(i)*180./pi &
                  , Force(i)*cos(alpha_frc(i)-alpha_Bz), Force(i)*sin(alpha_frc(i)-alpha_Bz)
            enddo
            Do i=1,N_frc-1
                if((h_frc(i)-h_frc(i+1)).lt.500 ) then
                    PenFacHeight=PenFacHeight*(1.d0 + 4.*((h_frc(i)-h_frc(i+1))/500.-1)**2)
                endif
            enddo
        endif
        ! h_f are heights above GroundLevel. Put penalty factor when this exceeds the top-level
        if((h_f(1)-AtmHei_dim*AtmHei_step).gt.1000 ) then
            PenFacHeight=PenFacHeight*(1.d0 + 4.*((h_f(1)-AtmHei_dim*AtmHei_step)/1000.-1)**2)
        endif
        if(PenFacHeight.gt.1.) write(2,*) 'penalty factor for E-field configuration=', PenFacHeight
        !
        ! AlternativeSmooth=.true.  !  ======
        H_s=700.  ! shift in height [m]
      !  H_s=000.  ! shift in height [m]
        write(2,*) 'H_s=',H_s,'[m] and D_ESmooth=', D_ESmooth
        Do j=1,N_frc
            z=abs(h_f(j))/Cos_Zenith ! calculated distance along shower axis, assuming flat Earth for simplicity
            i=z/AtmHei_step
            if(i.ge.AtmHei_dim-1) i=AtmHei_dim-1
            If(AlternativeSmooth) then ! Exponentially convoluted step in downward direction
                ! Alternative
                I1=i
                b=PenDepth(I1) -(z-I1*AtmHei_step)*AirDensity(I1)
                write(2,*) 'D_ESmooth=', D_ESmooth, &
                        ' with dx*e^{dx/DEs*X0} gives smoothing distance',D_ESmooth*X_EM/AirDensity(I1)
                Do i=I1,0,-1  ! calculate contribution of this step to all lower heights
                    height=i*AtmHei_step  ! distance along shower
                    !a=1.d0 -exp( -(PenDepth(i) - b )/(D_ESmooth*X_EM) ) ! folding with e^{dx/DEs*X0}
                    c=(PenDepth(i) - b )/(D_ESmooth*X_EM)
                    c=c*c  !  make it always positive
                    If(c.gt.5000.) Then
                        a=1.d0
                    Else
                       a=c/(1+c)
                    EndIf
                    !write(2,*) i,a,(PenDepth(i) - b )/(D_ESmooth*X_EM)
                    Force_x(i)=Force_x(i)+ a*F(j)*cos(ang(j))
                    Force_y(i)=Force_y(i)+ a*F(j)*sin(ang(j))
                enddo
            else  ! best fit, 2layr: D_ESmooth=2.5; 3layr: D_ESmooth=1.3
               b=(PenDepth(i)-PenDepth(i+1))/(AtmHei_step*X_EM)  ! convert AtmHei_step  to change in X  = AirDensity(i)
               b= D_ESmooth*b
               XXX=b*H_s
               !write(2,*) 'H_s=',H_s,'[m] and D_ESmooth=', D_ESmooth
               Do i=AtmHei_dim,0,-1    ! calculate transverse force
                 height=i*AtmHei_step  ! *Cos_Zenith for vertical height above ground
                 ! a=1.d0/(1.d0+exp(0.5*(height-abs(h_f(j)))/AtmHei_step))
                 If( (XXX+b*(height-z)) .gt. 1000.) then  ! distances measured along shower axis
                  a=0.
                 ElseIf( (XXX+b*(height-abs(h_f(j)))) .lt. -1000.) then
                  a=1.
                 Else
                  a=1.d0/(1.d0+exp(XXX+b*(height-abs(h_f(j)))) )  ! Smoothing of height dependence
                 EndIf
                 Force_x(i)=Force_x(i)+ a*F(j)*cos(ang(j))
                 Force_y(i)=Force_y(i)+ a*F(j)*sin(ang(j))
                enddo
            endif
        enddo
      elseif(line) then
        write(2,*) 'line definition for E-field:'
        Do j=0,N_frc
            write(2,*) h_frcL(j)/1000., ForceL(j), alpha_frcL(j)*180./pi
        enddo
        Do j=1,N_frc
            D1=abs(h_frcL(j-1)/Cos_Zenith)    ! convert to distance along shower axis
            D2=abs(h_frcL(j)/Cos_Zenith)
            I1=1 + D1/AtmHei_step
            I2=D2/AtmHei_step
            if(I2.ge.AtmHei_dim) I2=AtmHei_dim
            if(I2.lt.I1) then
                write(2,*) 'Heights for E-Field too close to distinguish', D1,D2,j,I1,I2
                stop
            endif
            a=(ForceL(j)*cos(alpha_frcL(j)) - ForceL(j-1)*cos(alpha_frcL(j-1)))/(D2-D1)  ! x-component
            b=(ForceL(j)*sin(alpha_frcL(j)) - ForceL(j-1)*sin(alpha_frcL(j-1)))/(D2-D1)  ! y=component
            Do i=I1,I2
                height=i*AtmHei_step  ! istance along shower axis
                Force_x(i)=Force_x(i)+ ForceL(j-1)*cos(alpha_frcL(j-1)) + a*(height-D1) !
                Force_y(i)=Force_y(i)+ ForceL(j-1)*sin(alpha_frcL(j-1)) + b*(height-D1) !
            enddo
        enddo
      endif
    endif    ! skip when there is no electric field given
!
!stop
   OPEN(UNIT=4,STATUS='unknown',FILE=trim(FileShCurrent)//'.dat')
   write(2,*) 'Initialize_shower:AtmHei_dim=',AtmHei_dim, X_start, PenDepth(AtmHei_dim), PenDepth(AtmHei_dim/2)
   write(4,"('!',7x,'z [km],',18x,'X [g/cm^2],',16x,'mean refractivity,',12x,'Ix,',23x,'Iy,',19x,'Ch Excess,',&
     19x,'dxi,',19x,'alpha_tr,',19x,'F_x,',23x,'F_y,', &
     ! 19x,'Ix/N,',23x,'IQ/N,', &
     19x,'|F| [keV/m],',19x,'phi')")
   Ix(AtmHei_dim)=0.  ; Iy(AtmHei_dim)=0.
   if(X_max .gt.2000.) X_max=2000.
   !if(x_0.gt. 0. .and. lamx .lt. 20.) lamx=20.
   !if(x_0.lt. 0. .and. lamx .gt. 1.) lamx=1./lamx
   If(HoriShwr) vDrift2=.false.
   alpha_tr(AtmHei_dim)=1.d0
   alpha_E=0.d0  ! increment in thickness due to the cumulative effect of E-field
   Isqr_max=0.0
   F_max=0.0
   Xmax_0 = 500 ! 540/Cos_Zenith changed June 24, 2018
   If(vDrift2) Xmax_0 = 540/Cos_Zenith   ! changed June 24, 2018 to =500
   !write(2,*) X_max, X_0, lamx, Energy_sh
   XmaxAve=X_max
   If(Energy_sh2.gt.0.) XmaxAve=(Energy_sh*X_max+Energy_sh2*X_max2)/(Energy_sh+Energy_sh2)
   Do i=AtmHei_dim-1,0,-1    ! calculate ixmx
      if(PenDepth(i).gt.XmaxAve) exit
   EndDo
   iXmx=i  ! corresponds to the index where Xmax is
   iImx=1
   iQmx=1
   W_BxvxB=0. ! Nov 2020: this gives good agreement for U/I, better than +/-1./10.
   !W_BxvxB=-1. ! Nov 2020: this W_BxvxB=(-1.) seems to explain the asymmetry in intensity along the vxvxB axis seen in CoREAS
                  ! but moves U/I off zero (negative) near the core.
   !W_BxvxB=+1. ! Nov 2020: Worse fit intensity than =-1, and moves U/I off zero (positive) near the core.
   !W_BxvxB=+3./10. ! Nov 2020: this gives good agreement for U/I
   !Xc_0=300
   SigX0=600.  ! used in transverse current
   SigX1=600.  ! used in charge excess
   Sigx2=0.
   !write(2,*) 'SigX0=',SigX0,', SigX1=',SigX1,', SigX2=',SigX2,', a_ChX=',a_ChX,', Xc_0=',Xc_0
   SigX0=SigX0*SigX0
   SigX1=SigX1*SigX1
   NP_F=0.15 ! 0.10
   FF_tc=1.0
   write(2,*) 'NP_F=',NP_F,' in "NPart*(1. +NP_F*(rho(1)/rho(i))*(Force_x(i)**2 + ..."'
   write(2,*) 'Shower-power:', ((X_max-X_0)/lamx), ((X_max2-X_02)/lamx2)
   flush(unit=2)
   Do i=AtmHei_dim-1,0,-1    ! calculate currents
      npart2=0.
        X_rh=PenDepth(i)
        If(X_rh.le.X_start) then
            Ix(i)=0.
            Iy(i)=0.
            IQ(i)=0.
            alpha_tr(i)=1.d0
            cycle
        endif
        z=i*AtmHei_step  !distance above ground along shower axis
        If(RL_param) then
         If(( 1 - R_0*(X_max-X_rh)/L_0).lt.0) then
            NPart=0
         else
            NPart=Energy_sh*( 1 - R_0*(X_max-X_rh)/L_0)**(1/(R_0*R_0)) * exp((X_max-X_rh)/(L_0*R_0)) !v3b and earlier
         endif
         If(( 1 - R_02*(X_max2-X_rh)/L_02).lt.0 .or. (Energy_sh2.le.0.)) then
            NPart2=0
         else
            NPart2=Energy_sh2*( 1 - R_02*(X_max2-X_rh)/L_02)**(1/(R_02*R_02)) * exp((X_max2-X_rh)/(L_02*R_02)) !v3b and earlier
            Npart=Npart+Npart2
         endif
        Else
         NPart=Energy_sh*( (X_rh-X_0)/(X_max-X_0))**((X_max-X_0)/lamx) * exp((X_max-X_rh)/lamx) !&
         If(Energy_sh2.gt.0. .and. ( (X_rh-X_02)/(X_max2-X_02) .gt. 0. ) ) Then
            NPart2=Energy_sh2*( (X_rh-X_02)/(X_max2-X_02))**((X_max2-X_02)/lamx2) * exp((X_max2-X_rh)/lamx2) !&
    !        write(2,*) 'Shower:', i,X_rh, Npart, Npart2  ! , &
    !           ( (X_rh-X_02)/(X_max2-X_02))**((X_max2-X_02)/lamx2), exp((X_max2-X_rh)/lamx2)
    !        write(2,*) 'Showerdetail:', i,X_rh,  &
    !           ( (X_rh-X_02)/(X_max2-X_02))**((X_max2-X_02)/lamx2), (X_rh-X_02)/(X_max2-X_02), ( (X_rh-X_0)/(X_max-X_0))
            If(NPart2.lt.0.) NPart2=0.
            Npart=Npart+Npart2
            ! Note that some shower-age type factors are calculated from the main X_Max
         EndIf
         !   * exp(AirDensity(iXmx)/AirDensity(i))     !v3c
        Endif
        NPart=NPart*(1. +NP_F*(AirDensity(1)/AirDensity(i))*(Force_x(i)**2 + Force_y(i)**2)/(100.*Force0*Force0)) ! to account for particle production with a strong force
        !----------------------------------
        ! FudgeFactor for testing different forms longitudinal function transverse current
        !FF_tc=exp(-(XmaxAve-X_rh)**2/SigX0)  !take out, Sept 2021
        !--------------------
        If(vDrift2) then
            BoF = 9.*X_rh*sqrt(XmaxAve*Xmax_0)/((XmaxAve + 2.*X_rh)*(XmaxAve + 2.*X_rh)) / F_over_beta
        else
            !v2!BoF = 4./((X_max-Xb_0)/(X_rh-Xb_0) + 3.)*sqrt(0.06/AirDensity(i)) / F_over_beta ! v2 Beta over force corrected for atm.
            !v3: BoF = 4./((X_max-Xb_0)/(X_rh-Xb_0) + 3.)/ F_over_beta ! v3 Beta over force corrected for atm.
            BoF = 4./((XmaxAve-Xb_0)/(X_rh-Xb_0) + 3.)*sqrt(AirDensity(iXmx)/AirDensity(i)) / F_over_beta ! v3b Beta over force corrected for atm.
            ! The problem with this BoF is that for double-humped showers the order in which sh1 and sh2 are given start to matter.
            ! see https://en.wikipedia.org/wiki/Drag_(physics) where the concept of 'terminal velocity'
            !  is discussed (final velocity under constant force). The terminal velocity is prop to \sqrt(force/density) for potato-shaped objects.
            !  From phenomenology we think we know that the velocity is prop to the magnetic field (no \sqrt)!
            !  Could still be that for the denisty dependence we need a sqrt(0.06/AirDensity(i))
        endif
        ux= Force_x(i) * BoF
        uy= Force_y(i) * BoF
        ! for the following Bxu_drift Lorentz force, only the magnetic field along the shower, B_z=J0t*cos(alpha), contributes to the transverse current
        uy = uy + ux*J0t*cos_alpha * BoF*W_BxvxB  ! Jan 2019, effect seems to be much too large for LOFAR, W_BxvxB=1/10.
        ux = ux - uy*J0t*cos_alpha * BoF*W_BxvxB  ! Nov 2020: this W_BxvxB=(-1.) seems to explain the asymmetry in intensity along the vxvxB axis seen in CoREAS
        s=u0/sqrt(u0*u0 + ux*ux + uy*uy)  ! u0=velocity, account for saturation effects
        !write(2,*) 's=',i,X_rh,vDrift2, NPart,npart2,Force_x(i), BoF, FF_tc
        Ix(i)=NPart * ux*s  * FF_tc
        Iy(i)=NPart * uy*s  * FF_tc
        Isqr=iy(i)*iy(i) + ix(i)*ix(i)
        if(Isqr.gt.Isqr_max) then
            Isqr_max=Isqr
            iImx=i
        endif
        D_IMax=iImx*AtmHei_step
        !
        If(vDrift2) then
            IQ(i)=J0Q* NPart * (1.+a_ChX)/(a_ChX + XmaxAve/X_rh)  ! Extra factor to model density dependent charge excess
        else
            b=(X_rh-XmaxAve)/(Xmax_0+X_rh-Xc_0)
            IQ(i)=J0Q* NPart* ((b/a_ChX) +1.)* (1.d0-exp(-5.d0*((X_rh-Xc_0)/(XmaxAve-Xc_0)))) &
               *sqrt(AirDensity(i)/.06) * AirDensity(iXmx)/.06 ! v3
        !       *(AirDensity(i)/.06) * AirDensity(iXmx)/.06 ! v3
            IQ(i)=IQ(i)*exp(-(XmaxAve-X_rh-SigX2)**2/SigX1)
        endif
        If(IQ(iQmx) .lt. IQ(i)) iQmx=i
        !
        !       Height dependence of pancake function may depend on X_rh instead of distance (in m) along shower axis
        alpha_E= (Force_x(i)**2 + Force_y(i)**2)*1.d-4  ! Effect of E=100 kV/m @ X_max (=600) is a 50% increase in lambda
        If(alpha_E.gt.F_max) F_max=alpha_E
        alpha_tr(i)=1.d0 + PancakeIncField*alpha_E  + XDepAlpha*sqrt(X_rh/X_max) ! the parameter that determines the width of the pancake function
        if(alpha_tr(i) .lt. 1.d0) alpha_tr(i)=1.d0
        if(alpha_tr(i) .gt. 5.) alpha_tr(i)=5.
        !if(x_rh.lt.x_max) iXmx=i
        write(4,"(20G14.5E3)") z/1000., X_rh,xi(i),Ix(i),Iy(i),IQ(i)/100., dxi(i),alpha_tr(i),Force_x(i),Force_y(i) &
            ! , Ix(i)/NPart, IQ(i)/NPart , X_rh-X_max, Iy(i)/NPart &
            , sqrt(alpha_E)*100., atan2(Force_y(i), Force_x(i))*180./pi
    end do
    Ix(AtmHei_dim)=0.d0 ;  Iy(AtmHei_dim)=0.d0   ; IQ(AtmHei_dim)=0.d0
    close(unit=4)
    F_max=sqrt(F_max)
    If(ShPlane) then
        write(2,"(1x,'X_max=',f7.2,', X_0=',f7.2,', D=',f7.3,'[km], H=',f7.3,'[km], E=',1pE9.2, &
            '[GeV], FootShift (x,y)=(',0pf8.2,',',f8.2,') [m], density@Xmax/.06=',f6.2,f9.5)") &
            X_max, X_start, iXmx*AtmHei_step/1000., iXmx*AtmHei_step*Cos_Zenith/1000.,Energy_sh, FShift_x,FShift_y &
               , (AirDensity(iXmx)/.06),AirDensity(AtmHei_dim)
    else
        !call calc_alpha_vB(vBE,vBN,vBU,vvBE,vvBN,vvBU)
        c = vBE*vvBN - vvBE*vBN
        Noff=(FShift_x*vvBE -FShift_y*vBE)/c
        Eoff=-(FShift_x*vvBN -FShift_y*vBN)/c
        write(2,"(1x,'X_max=',f7.2,', X_0=',f7.2,', D=',f7.3,'[km], H=',f7.3,'[km], E=',1pE9.2, &
            '[GeV], FootShift (x,y)=(',0pf8.2,',',f8.2,') or (E,N)=(',f8.2,',',f8.2,') [m]')") &
            X_max, X_start, iXmx*AtmHei_step/1000., iXmx*AtmHei_step*Cos_Zenith/1000.,Energy_sh, FShift_x,FShift_y, Eoff,Noff
        !write(2,"('Footprint shifted by (x,y)=(',f9.3,',',f9.3,') equivalent to (E,N)=(',f9.3,',',f9.3,') [m]')") &
        !    FShift_x,FShift_y
    endif
    If(Energy_sh2.gt.0.) Then
        write(2,"(1x,'X_max2=',f7.2,', X_02=',f7.2,', lamx2=',f7.3,', E2=',1pE9.2, '[GeV]')") &
            X_max2, X_02, lamx2, Energy_sh2
    EndIf
    write(2,"(A,F6.2,A,F6.1,A)") 'Charge excess max@ ',iQmx*AtmHei_step/1000.,'[km]=', PenDepth(iQmx),'[g/cm^2]'
    write(2,"(1x,'Curr_Xmax=',f7.2,', D=',f7.3,'[km], H=',f7.3,'[km], F_max=',F5.3,' x 100 keV/m, cos_Zen=',f4.2, &
        ' Peak current in x and y:',2g12.4)") &
        PenDepth(iImx), iImx*AtmHei_step/1000., iImx*AtmHei_step*Cos_Zenith/1000. ,F_max, Cos_Zenith, ix(iImx),iy(iImx)
!
!   Calculate magnitude of moving electric dipole through the integrated current
!
!    Ix_int(AtmHei_dim)=Ix(AtmHei_dim)  ; Iy_int(AtmHei_dim)=Iy(AtmHei_dim)
    Ix_int(AtmHei_dim)=0.0  ; Iy_int(AtmHei_dim)=0.0
!    IntegrateCurrent=0.2 ! Decay distance= s*X_0
!    IntegrateCurrent=0.02 ! Decay distance= s*X_0
    if(IntegrateCurrent.le.0.0003) then
        Ix_int(:)=0.d0
        Iy_int(:)=0.d0
        If(FiPa) write(2,"('No integrated charge dipole included')")
    else
        do i = AtmHei_dim-1, 0,-1
          z=i*AtmHei_step
          X_rh=(PenDepth(i+1)-PenDepth(i))/IntegrateCurrent !1.1  ! negative mass stepped through, scaled
          c=exp(X_rh/X_EM)
          !Ix_int(i)=Ix_int(i+1)*c+Ix(i)*AtmHei_step*(c-1.d0)*X_EM/X_rh
          !Iy_int(i)=Iy_int(i+1)*c+Iy(i)*AtmHei_step*(c-1.d0)*X_EM/X_rh
          Ix_int(i)=Ix_int(i+1)*c - Ix(i)*AtmHei_step*(c-1.d0)*X_EM/X_rh
          Iy_int(i)=Iy_int(i+1)*c - Iy(i)*AtmHei_step*(c-1.d0)*X_EM/X_rh
    !      Ix_int(i)=0.0  ; Iy_int(i)=0.0
        enddo
        If(FiPa) write(2,"('Integrated charge damped with',f5.3,' times X_EM=',f6.3)") IntegrateCurrent,X_EM
    endif
!
    if(test) OPEN(UNIT=4,STATUS='unknown',FILE='plot/sh_RadWei.dat')
    if(test) write(4,"('!',10x,'d[m/cm]',20x,'w(r) m',20x,'f(h) [cm]@100m',20x,'df/dh @100m' )")
   If(Moli_tc.lt.0.) Moli_tc=-Moli_tc  !  ,MoliereRadius iXmx iImx*AtmHei_step/1000.
   If(MoliereRadius.eq.0.0) then
      s=iXmx*AtmHei_step/1000.  ! distance from ground to Xmax along shower axis
      Moli_tc=s*10.
      If(Moli_tc.gt.50) Moli_tc=50
      write(2,*) 'Moli_tc is set automatically to ', Moli_tc, ' m, distance from ground to X_max=',s,' km'
   Endif
    Norm_tc=2.3/Moli_tc  ; s_tc=1.5
    s=0.        ! dummy
    alpha=alpha_tr(iXmx)  ! index pancake fie, at X_max
    lam=alpha_tr(10)      ! index pancake fie, near the ground
    c=Moli_tc/10.
    do CD_i = 1, 200
        ddd=(CD_i-0.5)*c    ! distance to core in [m]
        z=W_tc(ddd,a)
        s=s+z*c
        b=ddd/1000.     ! distance to showerfront in [mm]
!        write(4,*) d,(1000.*2.5/(moli*2*pi))/((d/Moli+1.)**3.5), &
        call fieh(b,1,lam,nu,a,dfha)  ! @10*dz m
        call fieh(b,1,alpha,fh,dfhl,dfha)    ! @ X_max
        dfha=dfha *(alpha_tr(iXmx+1)-alpha_tr(iXmx))/AtmHei_step
        if(test) write(4,*) ddd,z*100., nu,alpha,fh,dfhl+dfha,dfhl,dfha
    end do
    norm_tc=norm_tc/s   ! Norm for radial function w_tc(r)
!    write(2,214) norm_tc,1/s, moli_tc
!214 Format('transverse current: norm=',f7.3,', extra-norm=',f7.3,', radial parameter=',f7.2,'[m]')
    norm_tc=norm_tc/(2.*pi)
    if(test) close(unit=4)
!
!   Just for display purposes:
    if(test) then
        alpha=alpha_tr(iXmx)  ! index pancake fie, at X_max
        write(2,*) '@Xmax, alpha=',alpha, ', total force=',sqrt(Force_x(iXmx)**2 + Force_y(iXmx)**2)
        r1=0. ; r2=1.  ; dh=0.03  ! [m]
         call IntegFh(r1,r2,dh,alpha)
        r1=1. ; r2=4.  ; dh=0.03  ! [m]
         call IntegFh(r1,r2,dh,alpha)
        r1=4. ; r2=6.  ; dh=0.03  ! [m]
         call IntegFh(r1,r2,dh,alpha)
        r1=6. ; r2=10.  ; dh=0.05  ! [m]
         call IntegFh(r1,r2,dh,alpha)
        r1=10. ; r2=20.  ; dh=0.1  ! [m]
         call IntegFh(r1,r2,dh,alpha)
        r1=20. ; r2=50.  ; dh=0.1  ! [m]
         call IntegFh(r1,r2,dh,alpha)
        r1=50. ; r2=150.  ; dh=0.2  ! [m]
         call IntegFh(r1,r2,dh,alpha)
        r1=10. ; r2=100.  ; dh=0.1  ! [m]
         call IntegFh(r1,r2,dh,alpha)
!    r1=50. ; r2=70.  ; dh=0.2  ! [m]
!    call IntegFh(r1,r2,dh,alpha)
    endif
!
    If(FiPa) write(2,215) ZenithAngle_shower*180./pi !, moli_tc
215 Format('Zenith angle=',f7.3,' degree' ) ! ,', radial parameter=',f7.2,'[m]')
    FiPa=.false.
    Flush(unit=2)
!    stop
    return
    end
!
!===========================================================================================
subroutine IntegFh(r1,r2,dh,alpha)
!       calculate f(h) integrated over a range of core distances
    use constants, only : dp,pi
    implicit none
    real(dp), intent(in) :: r1,r2,dh,alpha
    integer :: ih,CD_i,ir
    integer, parameter :: NrInt=20
    real(dp) :: h,r,dwr,Nd, INd,IFh,fh,fh1,fh2,dfhl,dfha,CD_d,dr, accum,h25
    character*5 :: r1r2
    logical :: fiftypct,twfivpct
    external W_tc
    real*8 W_tc
!
    dr=(r2-r1)/NrInt
    accum=0.0
    fiftypct= .true.    ; twfivpct=.true.
    write(r1r2,"(I2.2,I3.3)") NINT(r1),NINT(r2)
    OPEN(UNIT=4,STATUS='unknown',FILE='plot/sh_IntFh-'//r1r2//'.csv') ! comma separated values
    write(4,"( '!   h   ',20x,',   f(h) ',20x,',  f(h)/N',20x,',  accumulated' )" )
!    write(2,*) r1r2,dh
    Do ih=1,500
        h=(ih-0.5)*dh
        IFh=0.  ; INd=0.
!        write(2,*) ih,h,dh
        Do ir=1,NrInt  ! integration over core distance
            r=r1+(ir-0.5)*dr
            Nd=W_tc(r,dwr)        ! ???
            INd=INd + Nd*dr
            call Find_CD_Interpol(r,CD_i,CD_d)
            call fieh(h,CD_i  ,alpha,fh1,dfhl,dfha)
            call fieh(h,CD_i+1,alpha,fh2,dfhl,dfha)
            Fh=(1.-CD_d)*Fh1 + CD_d*Fh2
            IFh=IFh + Fh*Nd*dr
        enddo
        accum=accum + IFh*dh/INd
        if(twfivpct .and. accum.gt.0.25) then
            twfivpct=.false.
            h25=h
        endif
        if(fiftypct .and. accum.gt.0.5) then
            fiftypct=.false.
            write(2,"( '25,50-percentile pancake at h=',F6.2,F6.2,' [m] between ',f5.1, f7.1,'[m]' )") h25,h,r1,r2
        endif
        if(IFh*dh/INd.gt. 10.d-12) write(4,*) h,' ,',IFh,' ,',IFh*dh/INd,' ,',accum,ih
!        write(2,*) ih,h,dh
    enddo
    close(unit=4)
!    write(2,*) 'fhint=',fhint
!
end subroutine IntegFh
