    Include 'Constants.f90'
!-----------------------------------------
    module RFootPars
      use constants, only : dp
!      real(dp), save :: J0y,lamy,X_may  ! Not used any more
      character (len=22), save :: release='August 5, 2018(Sept18)'
      integer, parameter :: N_FitPar_bas=10
      integer, parameter :: N_step_max=4    ! Used in the definition of the Atmospheric E-Fields
      integer, parameter :: N_line_max=6  ! N_line_max+1 should exceed N_step_max inorder not to have problems with the All_Params common block
      integer, parameter :: N_FitPar_max=N_FitPar_bas+3*(N_line_max+1)
      real(dp), parameter :: F_over_beta=250. ! [keV/m] from CONEX calculation by Krijn, corrected with cos(27)
      integer, save  :: SelectFh=4 ! option for pancake thickness parametrization
      real(dp) :: J0Q=0.25d0, lamx=100d0, XDepAlpha=0.0 ! Possible to fit
      real(dp), save :: lam_tc=0.05d0  ! july 202, moved out of fit common block
      real(dp), save :: lam_100=7.d0  ! july 202, moved out of fit common block
      real(dp), save :: rh0=0.000292d0, a_ChX=2.5
      real(dp) :: u0=60.d0/F_over_beta ! kV/m. determines saturation velocity, fit parameter per July 2020
      !real(dp), save :: F_over_beta=280. ! [keV/m] from CONEX calculation by Krijn
      real(dp), save :: PancakeIncField=0.41d0
      !integer, save :: NTo=100  ! Larger limit is needed for pulses at large distance, but value depends on TTRACE_STEP; works for  ObsDist_dim= 60
      integer, save :: NTo=150  ! Larger limit is needed for pulses at large distance, but value depends on TTRACE_STEP; FF_dim=16400=o(40020)
      integer, save :: Padding=5000  ! before & After padding-length in time-trace before down sampling
      real(dp) :: MoliereRadius=27.0d0  ![m] Possible to fit
      real(dp), save :: IntegrateCurrent=-0.0002 ! in negative no moving dipole is included
      integer, save :: StParRange  ! length of time-trace in number of down-sampled timestep used for calculating observables
      real(dp), save :: nu_min, nu_max, SamplingTime_dwn  ! Frequency window and sampling time for calculating observables
      real(dp), save :: GroundLevel=0.  ! [m]
      Logical, save :: RL_param=.false.
      logical, save :: AlternativeSmooth=.false.
      Real(dp), save :: R_0=0.28, L_0=215  ! Parametrization long. profile
      Real(dp), save :: R_02=0.28, L_02=215.  ! Parametrization long. profile for second bump
      real(dp)  :: X_0, X_max  ! Possible to fit
      real(dp)  :: X_02=100d0, lamx2=100d0, X_max2=600d0, Energy_sh2 ! second bump LDF
      real(dp), save :: D_IMax
      real(dp), save :: Zen_sh, Azi_sh  ! [deg], [North:Azi=90, East:Azi=0]
      real(dp), save :: alpha_vB=90., J0t, Zen_B, Azi_B  ! [deg], [North:Azi=90, East:Azi=0]
      real(dp)  :: h_frc(1:N_step_max),Force(1:N_step_max),alpha_frc(1:N_step_max),FShift_x=0.d0,FShift_y=0.d0
      real(dp)  :: h_frcL(0:N_line_max),ForceL(0:N_line_max),alpha_frcL(0:N_line_max)
!      real(dp)  :: h_frc_tr(1:N_line_max),Frc_tr(1:N_line_max),alpha_frc_tr(1:N_line_max)
      real(dp), save  ::  alpha_frc0, F_lim, PenFacHeight
      real(dp) ::  D_ESmooth=2.d0 ! smoothing E-field changes ; july 2020 moved to fit common block
      logical, save :: step=.false.,stpv=.false.,line=.false., test, Intensity_Weight=.false., vDrift2=.false.
      integer, save  :: N_frc
!      integer, parameter :: NAnt_max=1  ! 900
      Real(dp), save :: RnrmA=2.5, RnrmB=0.5
      real(dp) :: APar(N_FitPar_max)
      character (len=7), save :: AParMnm(N_FitPar_max)
      integer, save :: FitParam(N_FitPar_max), N_FitPar
      Logical, save :: Fit_StI=.false.   ! if =.true. only Stokes-I is fitted
      common / All_Params / J0Q,lamx,X_0,X_max,D_ESmooth,u0,XDepAlpha,MoliereRadius,FShift_x,FShift_y,ForceL,alpha_frcL,h_frcL
      character (len=7) :: AParMnmBas(N_FitPar_bas)=(/ &
          '    J0Q','   lamx','    X_0','  X_max','D_ESmth','     u0','XDepAlp','MollyRd','Shift_x','Shift_y'   /)
      equivalence (APar(1),J0Q)
      equivalence (APar(N_FitPar_bas + 1),h_frc(1))
      equivalence (Force(1),APar(N_FitPar_bas + 1+N_step_max))
      equivalence (alpha_frc(1),APar(N_FitPar_bas + 1+2*N_step_max))
      equivalence (AParMnmBas(1), AParMnm(1))
    end module RFootPars
!-----------------------------------------
     module eventdata
      use constants, only : dp
      integer, parameter :: N_ant_max=1000, N_scnt_max=30
      REAL(dp), save :: distance_antenna(N_ant_max), phi_antenna(N_ant_max),ZenithAngle_shower, F_max, Energy_sh
      REAL(dp), save :: St_I(N_ant_max), St_Q(N_ant_max), St_U(N_ant_max), St_V(N_ant_max)
      REAL(dp), save :: sigma_I(N_ant_max), sigma_Q(N_ant_max), sigma_U(N_ant_max), sigma_V(N_ant_max)
      REAL(dp), save :: StI_a(N_ant_max),StQ_a(N_ant_max),StU_a(N_ant_max),StV_a(N_ant_max)
      integer, save :: N_ant
      REAL(dp), save :: x_scnt(N_scnt_max), y_scnt(N_scnt_max), LORA(N_scnt_max), sigma_LORA(N_scnt_max), &
                         S_scnt(N_scnt_max)
      REAL(dp), save :: NoisePower=-1.  ! =\sigma_n^2
      REAL(dp), save :: Norm_I=1.  ! Normalisation factor for Stokes I to match MGMR with data
      character(len=30), save :: OFile  ! actual file used for writing grid data
      integer, save :: N_grid
      Real(dp), save ::  StI_max, d_grid
      integer, save :: N_scnt
    !use eventdata, only : Voltages, Core_N, Core_E, Ant_N, Ant_E,Eoff,Noff, RelMx_N, RelMx_E, RelMx_U
    !use eventdata, only : vBE,vBN,vBU,vvBE,vvBN,vvBU, vBxvvB, SN, SE, SU, sZS, Ang_Ux, alpha_Bz
      Real(dp), save ::  Core_N, Core_E, Eoff,Noff, RelMx_N, RelMx_E, RelMx_U
      Real(dp), save ::  vBE,vBN,vBU,vvBE,vvBN,vvBU, vBxvvB, SN, SE, SU, sZS, Ang_Ux, alpha_Bz
      Real(dp), save ::  Ant_N(N_ant_max), Ant_E(N_ant_max) !, ShPlane
      logical, save :: Voltages, Fitting, ReadInput, HoriShwr, ShPlane=.true. ! data read-in already projected on shower-plane
      character(len=40) :: FileFitResult='plot/FitResult'
      character(len=45) :: FileShCurrent='plot/sh_Current'
      character(len=30) :: FileGrid='plot/grid'
      character(len=10) :: OutFileLabel='v5'
     end module eventdata
!------------------------------
module CrossProd
    use RFootPars, only : Zen_Sh, Azi_Sh, Zen_B, Azi_B, alpha_vB
    use eventdata, only : ZenithAngle_shower, HoriShwr, alpha_Bz
    use eventdata, only : Voltages, vBE,vBN,vBU,vvBE,vvBN,vvBU, vBxvvB, SN, SE, SU, sZS, Ang_Ux  ! constants during fitting
    use constants, only : pi,dp
contains
    !subroutine calc_alpha_vB(vBE,vBN,vBU,vvBE,vvBN,vvBU)
    subroutine calc_alpha_vB()
    !Determine the kinematic shower parameters staying constant while fitting
    !    ! Azimuth convention: [ East:Azi=0, North:Azi=90, Up:Zen=0]
    !    sZS=sin(Zen_Sh*rad)  ;  sZB=sin(Zen_B*rad)
    !    SE=cos(Azi_Sh*rad)*sZS ; SN=sin(Azi_Sh*rad)*sZS ; SU=cos(Zen_Sh*rad) Components ov v_shower
    !    Antenna position w.r.t. showr core,  x=vxB
    !               ant_x=FShift_x + dist*cos(theta)
    !               ant_y=FShift_y + dist*sin(theta)
    !         Ant_N=(FShift_x*vvBE -FShift_y*vBE)/c  ! FShift_is the shift of the core in the antenna plane
    !         Ant_E=-(FShift_x*vvBN -FShift_y*vBN)/c
	 implicit none
    !real(dp),optional, intent(out) :: vBE,vBN,vBU,vvBE,vvBN,vvBU
    real(dp) :: sZB, cBz, sBz
    real(dp) :: BE,BN,BU
    real(dp) :: svB,cvB, rad
    real(dp) :: wBE,wBN,wBU, vzE, vzN, vzU
    logical, save :: first=.true.
    !
    rad=pi/180.d0
    HoriShwr = .false.
    If(Zen_sh.gt.70) then
        ! Zen_sh = 90.d0
        HoriShwr = .true.
        if(first) write(2,"(' Horizontal shower is assumed !!!!!!')")
    endif
    ZenithAngle_shower=Zen_sh*rad
    !
    If(Azi_B .eq. 0d0) then   ! alpha_vB has been entered on input
        if(first) write(2,"('alpha_vB=',F5.1,'deg, on input')") alpha_vB
        Zen_B=Zen_Sh-tan(alpha_vB)/rad
        Azi_B=Azi_Sh
    endif
    !
    sZS=sin(Zen_Sh*rad)  ;  sZB=sin(Zen_B*rad)
    SE=cos(Azi_Sh*rad)*sZS ; SN=sin(Azi_Sh*rad)*sZS ; SU=cos(Zen_Sh*rad)  ! points from the core to the shower
    BE=cos(Azi_B*rad)*sZB ; BN=sin(Azi_B*rad)*sZB ; BU=cos(Zen_B*rad)
    if(first) write(2,"('components vShower & B (East, North, UP)=(',3f7.3,') & (',3f7.3,')')") SE,SN,SU, BE,BN,BU
    ! direction of vXB
    cvB=(SE*BE+SN*BN+SU*BU)   ! =S.B  = cos(angle{S,B})
    wBE=SN*BU-SU*BN ;wBN=SU*BE-SE*BU ; wBU=SE*BN-SN*BE  ! =(SxB) . (E,N,U)  , not normalized
    svB=sqrt(wBE*wBE+wBN*wBN+wBU*wBU)   ! =S.B  = sin(angle{S,B})
    vBE=wBE/svB ;vBN=wBN/svB ; vBU=wBU/svB   ! =[SxB] . (E,N,U)  , normalized i.e. [] vs ()
    !write(2,"('components B (East, North, UP)=',3f7.3)") BE,BN,BU
    alpha_vB=atan2(svB,cvB)/rad  ! angle{S,B}
    !

    vzE=Sn ; vzN=-SE ; vzU=0.d0   ! =SxZ, no norm
    vvBE=SN*vBU-SU*vBN ;vvBN=SU*vBE-SE*vBU ; vvBU=SE*vBN-SN*vBE  ! = Sx[SxB]  no norm; should equal [S (S.B) - B]/svB
    cBz=(vzE*vBE+vzN*vBN) ! =(SxZ).[SxB]=cos(vxB & vxz)  apart from norm factor (SxZ)
    sBz=(vzE*vvBE+vzN*vvBN) ! =(SxZ).(Sx[SxB]) = sin(vxB & vxz)  apart from norm factor
    alpha_Bz=pi-atan2(sBz,cBz)  ! angle of vxz w.r.t. vxB; reason for the "pi -" is guessing, related to strange sign for vector v??
    !write(2,"('components vxB & vxvxB (East, North, UP)=(',3f7.3,') & (',3f7.3,')')") vBE,vBN,vBU, vvBE,vvBN,vvBU
    write(2,"('alpha_Bz=',F5.1,'deg, sin(alpha_Bs)=',4f5.2)") alpha_Bz/rad ! ,cvB,sin(alpha_vB*rad),cos(alpha_vB*rad)
    !  sBz=|(SxZ)x[SxB]|=Z.[SxB]=vBU ! use Ax(BxC)=B (A.C)-C (A.B)
    write(2,"('Sin_Bz-test=',2F10.5)") sBz, vBU ! =, indeed correct

      vBxvvB = vBE*vvBN - vvBE*vBN  ! =[SxB]x(Sx[SxB])_u = S_u  because of norn SxB, use Ax(BxC)=B (A.C)-C (A.B)
    write(2,"('vBxvvB_z-test=',2F10.5)") vBxvvB, SU ! =, indeed correct
      !If(first) write(2,*) FShift_x*vBN+FShift_y*vvBN,FShift_x*vBE+FShift_y*vvBE,FShift_x*vBU+FShift_y*vvBU
      !If(first) write(2,*) 'project up vector on ground:',-(FShift_x*vBU+FShift_y*vvBU)*SN/SU, -(FShift_x*vBU+FShift_y*vvBU)*SE/SU
      Ang_Ux = ATAN2(vvBU, vBU) ! Angle between Up (theta polarization) and x (=vxB)

    if(first) write(2,"('alpha_vB=',F5.1,'deg, sin(alpha_vB)=',4f5.2)") alpha_vB,svB ! ,cvB,sin(alpha_vB*rad),cos(alpha_vB*rad)
        vBE=wBE/svB ;vBN=wBN/svB ; vBU=wBU/svB
        ! direction of vx(vXB)
        vvBE=SN*vBU-SU*vBN ;vvBN=SU*vBE-SE*vBU ; vvBU=SE*vBN-SN*vBE
        if(first) write(2,"('components vxB & vxvxB (East, North, UP)=(',3f7.3,') & (',3f7.3,')')") vBE,vBN,vBU, vvBE,vvBN,vvBU
    !endif
    first=.false.
    !
    end subroutine calc_alpha_vB
!-----------------------------------------
    subroutine calc_alpha_Bz(alpha_Bz)
	implicit none
    ! Azimuth convention: [ East:Azi=0, North:Azi=90, Up:Zen=0]
    real(dp), intent(out) :: alpha_Bz
    real(dp) :: sZS,sZB
    real(dp) :: SE,SN,SU,BE,BN,BU,vBE,vBN,vBU,vvBE,vvBN,vvBU,vzE,vzN,vzU
    real(dp) :: sBz,cBz, rad
    !
    rad=pi/180.d0
    !
    If(Azi_B .eq. 0d0) then   ! alpha_vB has been entered on input
        return
    endif
    !
    sZS=sin(Zen_Sh*rad)  ;  sZB=sin(Zen_B*rad)
    SE=cos(Azi_Sh*rad)*sZS ; SN=sin(Azi_Sh*rad)*sZS ; SU=cos(Zen_Sh*rad)
    BE=cos(Azi_B*rad)*sZB ; BN=sin(Azi_B*rad)*sZB ; BU=cos(Zen_B*rad)
    ! direction of vXB
    vBE=SN*BU-SU*BN ;vBN=SU*BE-SE*BU ; vBU=SE*BN-SN*BE
    vzE=Sn ; vzN=-SE ; vzU=0.d0
    vvBE=SN*vBU-SU*vBN ;vvBN=SU*vBE-SE*vBU ; vvBU=SE*vBN-SN*vBE  ! vectors are not normalized
    cBz=(vzE*vBE+vzN*vBN) ! cos(vxB & vxz)  apart from norm factor
    sBz=(vzE*vvBE+vzN*vvBN) ! sin(vxB & vxz)  apart from norm factor
    !
    alpha_Bz=pi-atan2(sBz,cBz)  ! angle of vxz w.r.t. vxB; reason for the "pi -" is guessing, related to strange sign for vector v??
    !write(2,"('components vxB & vxvxB (East, North, UP)=(',3f7.3,') & (',3f7.3,')')") vBE,vBN,vBU, vvBE,vvBN,vvBU
    write(2,"('alpha_Bz=',F5.1,'deg, sin(alpha_Bs)=',4f5.2)") alpha_Bz/rad ! ,cvB,sin(alpha_vB*rad),cos(alpha_vB*rad)
    !
    end subroutine calc_alpha_Bz
end module CrossProd
!-----------------------------------------
