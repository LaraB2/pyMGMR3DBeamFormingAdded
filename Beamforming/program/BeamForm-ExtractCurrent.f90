    Include 'Constants.f90'
    Include 'MakeGLEplots.f90'
    Include 'MGMR3D_FFT.f90'
    Include 'MGMR3D2_Atmosphere.f90'
!    Include 'MGMR3D2_PancakeFunct.f90'
!    Include 'MGMR3D2_LateralDistr.f90'
    Include 'KernelFilter.f90'
    Include 'MGMR3D_spline.f90'
    Include 'Beam_DataRead-v20.f90'
!    Include 'GenerateKernel.f90'
    Include 'Beam_PSF-v20.f90'
!    Include 'Beam_Kernel.f90'
    Include 'Beam_DataBeaming.f90'
    !
    Include 'Beam_ShStruct-v10.f90'
    Include 'Beam_Folding.f90'
    !Include 'Beam_Edata-v10.f90'
    !
    !Include 'C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI\FORTRANsrc/nl2sol.f90' !Take this out when using the LINUX makefile
    Include '../LightningPros/nl2sol.f90' !Take this out when using the LINUX makefile
!------------------------------
!============================================
	program BeamFormingExtractCurrent
   use constants, only : pi, dp, c_l, ci
   use PSF, only : TestWeight, ReadKernel, Base
   Use PSF, only : hcnt_l, hcnt_u
   use DataBeam, only : DataBeaming
   !Use CurrentFolding, only : UnFoldJ
   Use Systemcommands, only : MakeGLEplots
   Use Atmosphere, only : AtmHei_dim, AtmHei_step
   Use Atmosphere, only : Xi, PenDepth
   Use ShowerData, only :  GroundLevel_sh, ZenithAngle_deg, ZenithAngle_rad, N_antennas, N_timesamples, SamplingTime, DistMax
   Use PSF, only : Grd_Kx, GrdPSFWghtx, GrdPwrWghtx, PSF_tracex
   Use PSF, only : Label2, Label3, PSF_max, N_PSF
   !Use PSF, only : GrdPSFWghtr, GrdPwrWghtr, PSF_tracer
   use Atmosphere, only : AtmosphereInit
	implicit none
   INTEGER DATE_T(8) !,WLength
   CHARACTER*12 REAL_C(3)
   character(Len=6)  :: date,time,dateId(100),timeId(100)
!
   Real(dp) :: Theta_zen, Cos_Zenith
   Integer :: nxx, k, i, i_d, N_d, i_t, N_t, nuTrace_dim, i_zf, N_zf, i_atm, i_tr, i_tr2, N_tr
   real(dp) :: z_low, Z_high, Del_z, Del_d, del_tr, d_tr  !  Del_t,
   real(dp), allocatable :: ttracex_zf(:,:), ttracey_zf(:,:)
   !Integer :: N_PSF
   Real(dp) :: dt_PSF, zeta_f, d_zeta
   Real(dp), allocatable :: powr_x(:), powr_y(:), PSF_pow_x(:), PSF_pow_y(:), Jx1PSF(:), Jx2PSF(:), Jx1Pwr(:), Jx2Pwr(:)
   Character(len=2) :: Label
   Character(len=30) :: PlotType, PlotName, PlotDataFile
   Integer :: i_r, i_r1, N_r, N_r1, Iresol
   Real(dp) :: del_r, del_r1, d_r1, d_r, zeta_K(1:PSF_max), rs_0, phi_0
   Real(dp) :: EnergyCR, X_max, R_0, L_0, qual
   integer :: N_fitpar, error, FitOptn, ih_min, ih_max, N_phi0
   logical RL_Switch
   Real(dp), allocatable :: RInt_PSF(:), RInt_Pwr(:), LongProfile_PSF(:), LongProfile_Pwr(:)
   Real(dp), allocatable :: RInt_Ker(:), LongProfile_Ker(:)
   Logical :: HoriShwr, AntennaLayout
   !Character(len=10) :: Label2
   Character(len=30) :: Base_pol, KerBase_pol
   Character(len=20), save :: KernelBase, DataBase
   character(Len=1) :: ext
   Real(dp) :: X,Y, z , R, S, W ! scratch
   Real(dp) :: Ix(0:AtmHei_dim), Iy(0:AtmHei_dim), IQ(0:AtmHei_dim), GoodQual, FracFold
   !============================
   !  Parameters for grid in Kernel calculation
   real(dp) :: d_h=0.125 ! =h_range/N_h; h_range=30.d0 ! h_range=60.d0 ! h_range=50.d0 ! h_range=100.d0 ! h_range=17.d0 ! h_range=200.d0 ! h_range=100.d0 ! h_range=20.d0 ! h_range=2.5d0  ! h_range=1.d0  !  [m]
   Integer :: N_h=120 ! N_h=240 ! N_h=200 ! N_h=400 ! N_h=50 ! N_h=20 ! N_h=100 ! N_h=26  ! N_h=80 ! number of steps in distance behind front
   Integer :: dN_tr=10 ! dN_tr=5 ! dN_tr=1 ! dN_tr=3 !  Step size in t_r in units of AtmHei_step
   ! fine N_phi grid important for power at lower altitudes
   Integer :: N_phir=5 ! N_phir=15 ! N_phir=3 ! N_phir=7 ! N_phir=17 !  Number of angles for integration over shower area
   ! fine N_ca grid important for quenching fluctuations
!   Integer, parameter ::  N_ca= 50 ! N_ca= 100 ! N_ca= 200 ! N_ca= 250 ! N_ca= 50 !messy) ! grid for antenna to shower-core distance
   ! finer than d_crMax=20.d0 grid not important for dmax=250
   real(dp) :: d_crMax=20.d0  ! d_crMax=10.d0  ! [m]  ! max increment in core-ray grid
   !  Notes on grid structure.
   !- h_range=50.d0; is certainly sufficient.
   !     Reducing it to 10[m] severely affects the Z_f<5km results, for higher altitudes the differences are much smaller.
   !     At present the distance limit is chosen symmetrical around zero, but the contribution from negative values is minor
   !     compared to that from the positive side. It may thus be better to take the limit asymmetrical.
   !- N_h=200; It is important that d_h=h_range/N_h is considerably smaller (factor 5) than the grid spacing imposed by
   !     the frequency filtering. For 30-80 MHz a value d_h=0.25[m] appears sufficient. If this is not small enough the
   !     gradual move of the peak in WeiFie(i_h,i_tr) as function of t_r is not properly reflected, resulting in jumps in
   !     the weight function after filtering. Going to an even smaller d_h appears to be useless.
   !- dN_tr=5; The resolution in height (=dN_tr*(AtmHei_step=10m) ) should be at least as good as the required precision
   !     for the currents to be extracted. 10 or 20 may thus be time-saving with the same results.
   !- N_phir=5; If taken smaller it will result in a very large scattering, in particular for those kinematic conditions
   !     where the Cherenkov angle just starts to be visible. If odd, phi_r=90, the angle where Ch becomes visible, is on
   !     the grid, and for too large grid spacing its contribution is strongly over-estimated.
   !====================================================
!
    NAMELIST /BeamPars/ Base, KernelBase, DataBase, AntennaLayout, del_z,  del_tr, &
      zeta_K, & ! rs_0, N_phi0, phi_0, & !      N_h, d_h, dN_tr, N_phir, d_crMax,
      DistMax, GoodQual, FracFold
!
   !=========================================================
   !  Run parameters
   Base='plot/'
   DataBase=''
   KernelBase=''
   AntennaLayout=.true.  ! use realistic antenna layout
   del_z =20. ![AtmHei_step] step size for focal-height grid
   del_tr = 20. ![AtmHei_step] step size for shower-axis distance grid
   !zeta_min = 500 ! Smallest focal height
   !zeta_max = 9000 ! Largest focal height
   !zeta_step = 1500 ! Step in focal height for Kernell calculation
   zeta_K(:)=0.
   rs_0 = 0.   !Distance focal point off shower axis
   N_phi0 = 1
   phi_0 = 0.  ! Angle of off-axis point w.r.t. antenna
   DistMax = 0.  ! will be set automatically is antenna layout is read-in (AntennaLayout=.true. )
   GoodQual=0.1
   FracFold=0.6
   read(*,NML = BeamPars)
   If(KernelBase .eq. '') KernelBase=Base
   If(DataBase .eq. '') DataBase=Base
   !
   !===================================================
   OPEN(UNIT=2,STATUS='unknown',FILE=TRIM(Base)//'BeamForm-ExCurr.out')
   CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T) !-----------------------------
   write(*,232) 'Initialize',(DATE_T(i),i=5,8)               !-----------------------------
   write(*,*) 'Work folder= "',TRIM(Base),'"'               !-----------------------------
   write(2,232) 'Initialize',(DATE_T(i),i=5,8)               !-----------------------------
232 FORMAT(3X,'start ',A,' @ ',I2,':',I2,':',I2,'.',I3,1X,25(1H-))
   Write(2,NML = BeamPars)
   Flush(unit=2)
   !=========================================================
   !  Determine grid in zeta
   del_z =del_z  * AtmHei_step
   del_tr=del_tr * AtmHei_step
   write(2,*) 'del_z, del_tr:',del_z, del_tr
   !
   !---------------------------------
   !  Reading E-field time traces
   !AntennaLayout=.false.
   !  - Read time traces for individual antennas where antenna function is already corrected for.
   !  - Fold with impulse response to find pulse-time. Note: identical to just find position of max in time trace.
   !  - Fit pulse times to reconstruct the core position.
   If(AntennaLayout) Then
      !   label='-test'
      !   OPEN(UNIT=4,STATUS='old',FILE=trim(Base)//trim(label)//'.dat',iostat=nxx) ! E13.4
      Call ReadData(DataBase)  ! true LOFAR antennas
   Else
      !   OPEN(UNIT=4,STATUS='unknown',FILE=trim(Base)//'th_000'//trim(label)//'.csv',iostat=nxx) ! E13.4
      Call ReadThetData(DataBase, DistMax) ! From theta grid
   EndIf
   write(2,*) 'DistMax=',DistMax,'[m]'
   !
   !-----------------------------------------------------
   ! Initialize penetration depth and index of refraction
   HoriShwr=.false.
   Call AtmosphereInit(ZenithAngle_rad, GroundLevel_sh, HoriShwr)
   !
   !------------------------------------
   !  Read Pulse Shape Function and kernel & weights
   !
   !Call PSF_setup(N_PSF, Label2, Label3)
   !write(2,*) 'N_PSF:', N_PSF
   !Flush(unit=2)
   !Call ReadPSF(N_PSF, dt_PSF,  z_low, del_z, del_tr,  N_zf, N_tr, &
   !    AtmHei_dim, AtmHei_step )
   ! Hack!!!!!!!!!!!!!!!!!!!!!!!!!
   !Label2='50_'   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !---------------------------
   Call PrepareKernel(Label2, Label3, PSF_max, N_PSF, zeta_K, rs_0, phi_0)
   KerBase_pol=TRIM(KernelBase)//'x-'    ! kernal is made for geomagnetic (== x pol) only
   Base_pol=TRIM(Base)//'x-'    ! kernal is made for geomagnetic (== x pol) only
   Call ReadKernel( N_PSF, dt_PSF, z_low, del_z, del_tr, N_zf, N_tr, AtmHei_dim, AtmHei_step &
      , Base_pol, KerBase_pol, Grd_Kx, GrdPSFWghtx, GrdPwrWghtx, PSF_tracex )
   ! First index in WSqG is retarded time (to be folded with current); second is zeta_f
   write(2,*) 'N_PSF, dt_PSF:', N_PSF, dt_PSF
   Flush(unit=2)
   !-----------------------------------------
   !  Calculate  beaming results using the antenna E-Fields
   !
   Call DataBeaming(Base, N_antennas, N_timesamples, ttracex_zf, ttracey_zf, ih_min, ih_max, z_low, del_z, N_zf, del_tr, N_tr)
   !
   !-----------------------------
   !  Testing
   !  Read currents Ix(t_r) , just for testing
   OPEN(UNIT=4,STATUS='old',FILE=TRIM(DataBase)//'sh_Current.dat',iostat=nxx) ! E13.4
   write(2,*) 'Reading file "',TRIM(DataBase)//'sh_Current.dat','"', AtmHei_step, AtmHei_dim
   If(nxx.ne.0) Then
      write(2,*) 'problems reading file:', TRIM(Base)//'sh_Current.dat'
      STOP 'shower profile file is problematic'
   EndIf
   Flush(Unit=2)
   Read(4,*)
   k=0
   Ix(:)=0.
   Iy(:)=0.
   IQ(:)=0.
   Do i=AtmHei_dim-1,1,-1
      Read(4,*,IOSTAT=nxx) z,x,y,R,S,W
      !write(2,*) '!Sh current:', k, z,x,y,R,S,W, nxx
      !Flush(Unit=2)
      k=NINT(z*1000./AtmHei_step)
      !IF(IS_IOSTAT_END(stat)) exit
      Ix(k)=R
      Iy(k)=S
      IQ(k)=W*100.  ! to compensate the /100. when writing   sh_Current.dat
      If(k.eq.0) exit  ! there may be fewer entries then AtmHei_dim when NPart=0 beyond certain heights
   EndDo
   Close(Unit=4)
   !
   Base_pol=TRIM(Base)//'x-'
   Call TestWeight(Base_pol, GrdPSFWghtx, GrdPwrWghtx, del_tr, N_tr, del_z, N_zf, z_low, Ix, AtmHei_step, AtmHei_dim, &
         Grd_Kx, hcnt_l, hcnt_u, ttracex_zf, ih_min, ih_max, PSF_tracex)
   !
   Base_pol=TRIM(Base)//'y-'
   Call TestWeight(Base_pol, GrdPSFWghtx, GrdPwrWghtx, del_tr, N_tr, del_z, N_zf, z_low, Iy, AtmHei_step, AtmHei_dim, &
         Grd_Kx, hcnt_l, hcnt_u, ttracey_zf, ih_min, ih_max, PSF_tracex)
   !
   !-----------------------------
   !
   !Del_t = 5*c_l  ! convert 5 ns to meters, read-in
   !Call GetEData(Base, N_d, dt_PSF, ttracex_zf, ttracer_zf, ih_min, ih_max, z_low, del_z, N_zf, del_tr, N_tr &
   !   , AtmHei_dim, AtmHei_step, Xi)
   ! ttrace_zf( ih_min:ih_max,0:N_zf) is now constructed
   write(2,*) 'h-dims(main):', ih_min, ih_max,';',hcnt_l,hcnt_u
   ih_min = MAX(ih_min, hcnt_l)  ! Should be negative; including  , -10  does not really change results, as expected
   ih_max = MIN(ih_max, hcnt_u)
   flush(unit=2)
   write(2,*) 'Timespan for beaming-power calculations:', ih_min*SamplingTime, ih_max*SamplingTime
   !======================================================
   !
   !-----------------------------------------------------
   ! Analyze traces
   !
   Allocate( powr_x(0:N_zf), PSF_pow_x(0:N_zf), powr_y(0:N_zf), PSF_pow_y(0:N_zf) )
   Do i_zf=0,N_zf ! Loop over focal height on the shower axis   hcnt_l:hcnt_u
      powr_x(i_zf) = sqrt(SUM( ttracex_zf(ih_min:ih_max, i_zf)**2 ))
      PSF_pow_x(i_zf)= SUM( ttracex_zf(ih_min:ih_max,i_zf)*PSF_tracex(ih_min:ih_max, i_zf))  ! R(z_f) = [E(h,z_f)*Psf(h)]
      powr_y(i_zf) = sqrt(SUM( ttracey_zf(ih_min:ih_max, i_zf)**2 ))
      PSF_pow_y(i_zf)= SUM( ttracey_zf(ih_min:ih_max,i_zf)*PSF_tracex(ih_min:ih_max, i_zf))  ! R(z_f) = [E(h,z_f)*Psf(h)]
     ! write(2,"(A,100G13.4)") '!powr,PSF-pow; x:', i_zf, i_zf*del_z+z_low &
     !    , SUM( PSF_tracex(ih_min:ih_max, i_zf)*PSF_tracex(ih_min:ih_max, i_zf)) &
     !    , SUM( PSF_tracex(hcnt_l:hcnt_u, i_zf)*PSF_tracex(hcnt_l:hcnt_u, i_zf))
      !   , SUM( ttracex_zf(-1:5,i_zf)*PSF_tracex(-1:5, i_zf)) &
      !   , SUM( ttracex_zf(-10:10,i_zf)*PSF_tracex(-10:10, i_zf)) &
      !   , SUM( ttracex_zf(-40:40,i_zf)*PSF_tracex(-40:40, i_zf)) &
      !   , SUM( ttracex_zf(-60:60,i_zf)*PSF_tracex(-60:60, i_zf)),PSF_pow_x(i_zf)
      !   powr_x(i_zf), PSF_pow_x(i_zf)
   EndDo ! i_zf=0,N_zf ! Loop over focal height on the shower axis
   i_zf=25
   !write(2,"(A,100G13.4)") '!powr,ttracex_zf:', i_zf, ttracex_zf(-1:15,i_zf)
   !write(2,"(A,100G13.4)") '!powr,PSF_tracex:', i_zf, PSF_tracex(-1:15,i_zf)
   !
   !-----------------------------------------------------
   ! Unfold InpPSFWght below height of 1 - 15 km.
   ! When gridding the current, the safest seems to be to take the grid points at the same heights as where PSF_Pow is calculated.
   ! This because the averaging of a sloping curve is the same as the value at the center.
   ! A grid distance of 1km may befine.
   ! J(t_r) = (WW)^-1 (PRW)
   ! WW = \sum_\zeta  N^2(\zeta)  *  W(t_r, \zeta)  *  W(t'_r, \zeta) * d\zeta
   ! PRW = \sum_\zeta  N(\zeta)  *  W(t_r, \zeta)  *  PSF_pow(\zeta) * d\zeta
   !
   ext='x'
   Call CurrentExtraction(Base, ext, GrdPwrWghtx, powr_x, GrdPSFWghtx, PSF_pow_x, N_tr,N_zf, del_tr, del_z, z_low, &
      GoodQual, FracFold)
   !
   If(  (SUM(powr_y(0:N_zf))/SUM(powr_x(0:N_zf)) ) .gt. 0.001) Then
      ext='y'
      Call CurrentExtraction(Base, ext, GrdPwrWghtx, powr_y, GrdPSFWghtx, PSF_pow_y, N_tr,N_zf, del_tr, del_z, z_low, &
          GoodQual, FracFold)
   EndIf
   !
   !ext='r'
   !Call CurrentExtraction(Base, ext, GrdPwrWghtr, powr_r, GrdPSFWghtr, PSF_pow_r, N_tr,N_zf, del_tr, del_z, z_low)
   !
   Call MakeGLEplots(Submit=.true.)
   !
End Program BeamFormingExtractCurrent
!------------------------
Subroutine CurrentExtraction(Base, ext, GrdPwrWght, powr, GrdPSFWght, PSF_pow, N_tr,N_zf, del_tr, del_z, z_low, GoodQual, Frac)
   use constants, only : dp
   !Use PSF, only : hcnt_l, hcnt_u
   !use EData, only : GetEData
   Use CurrentFolding, only : UnFoldJ
   Use Systemcommands, only : MakeGLEplots
   Use ShowerStruct, only : ParDim
   Use Atmosphere, only : AtmHei_dim, AtmHei_step
	implicit none
   !
   character(Len=*), intent(in)  :: Base
   character(Len=1), intent(in)  :: ext
   Integer, intent(in) :: N_tr, N_zf
   real(dp), intent(in) :: z_low, Del_z, del_tr
   Real(dp), intent(in) :: GrdPSFWght(0:N_tr,0:N_zf), GrdPwrWght(0:N_tr,0:N_zf), powr(0:N_zf), PSF_pow(0:N_zf)
   Real(dp), intent(in) :: GoodQual
   Real(dp), intent(in) :: Frac
   !real(dp), allocatable :: ttracex_zf(:,:), ttracer_zf(:,:)
   !Integer :: N_PSF
   !Integer :: nxx, idi, i, i_d, N_d, i_t, N_t, i_zf, i_atm, ih_max, ih_min, i_tr, i_tr2
   Real(dp) :: zeta_f, d_zeta
   !Real(dp), allocatable :: powr_x(:), powr_r(:), PSF_pow_x(:), PSF_pow_r(:), Jx1PSF(:), Jx2PSF(:), Jx1Pwr(:), Jx2Pwr(:)
   !Character(len=2) :: Label
   Character(len=30) :: PlotType, PlotName, PlotDataFile
   Integer :: i, i_zf, i_tr, i_r, i_r1, N_r, N_r1  ! , Iresol
   Real(dp) :: del_r, del_r1, d_r1, d_r, d_tr
   Real(dp) :: J_max(1:ParDim), X_max(1:ParDim), R_0(1:ParDim), L_0(1:ParDim), qual
   integer :: error, FitOptn
   logical RL_Switch
   Real(dp) :: RInt_Pwr(0:N_zf), RInt_PSF(0:N_zf), RInt_Ker(0:N_zf)
   Real(dp) :: LongProfile_Pwr(0:N_tr), LongProfile_PSF(0:N_tr), LongProfile_Ker(0:N_tr)
   Real(dp) :: Jx1(0:N_tr), Jx2(0:N_tr)
   Integer :: FitParams(12)
   Character(len=30) :: comment
   !
   !Iresol=4
   !N_r=N_tr/Iresol  ! Number of grid points defining the current, using a linear interpolation
   !Allocate( Jx2(0:N_tr) )
   Call UnFoldJ(trim(Base)//ext, GrdPSFWght, PSF_pow, N_tr,N_zf, del_tr, del_z, z_low, Frac, Jx2)
   !
   !Iresol=3  ! stepsize in unfolded current, del_r1=del_tr*Iresol
   !N_r=N_tr/Iresol  ! Number of grid points defining the current, using a linear interpolation
   !Allocate( Jx1(0:N_tr) )
   Jx1(:)=0.
   !Call UnFoldJ(GrdPwrWght, powr, N_tr,N_zf, del_tr, del_z, Iresol, Jx1, N_r1, del_r1)
   !
   !write(2,"(A,100G13.4)") 'powr(i_z):', powr(0:N_zf)
   OPEN(UNIT=4,STATUS='unknown',FILE=trim(Base)//ext//'-BeamPwr-a.dat')
   write(2,*) 'writing file "',trim(Base)//ext//'-BeamPwr-a.dat','" ',z_low, N_zf
   flush(unit=2)
   write(4,"(A)") '!    zeta_f,  powr_x,  PSF_pow_x,    J_x(unfolding)'
   Do i_zf=0,N_zf ! Loop over  Retarded time
      zeta_f=i_zf*del_z + z_low   ! focal height
      i_tr=i_zf
      If(i_tr.gt.N_tr) i_tr=N_tr
      !
      write(4,"(F10.4,21(1pG13.4))") zeta_f/1000., powr(i_zf), PSF_pow(i_zf), Jx1(i_tr), Jx2(i_tr) &
         , del_tr*SUM(GrdPSFWght(0:N_tr,i_zf)*Jx2(0:N_tr))
   EndDo ! i_zf=0,N_zf ! Loop over  Retarded time
   Close(Unit=4)
   !DeAllocate( Jx1, Jx2 )
   !--------------------------------------------------------
   !  Fit an analytic longitudinal profile
   write(2,*) '------------------- Power-weight fit ------------------'
   J_max(:)=0.
   J_max(1)= 0.114566E+09 ; X_max(1)=850. ; R_0(1)=0.2716 ; L_0(1)=200.4626 ; RL_Switch= .true. ;  FitOptn=0
   FitParams(1)=4 ;  FitParams(2)=7 ;  FitParams(3)=10 ;  FitParams(4)=-10    ; FitParams(5)=-1
   RInt_Pwr(:)=1.
   LongProfile_Pwr(:)=0.
   !Call Fit_ShSruct(FitOptn, powr, GrdPwrWght, N_tr, N_zf, del_tr, del_z, z_low, J_max, X_max, R_0, L_0, FitParams, &
   !      RL_Switch, qual, error, RInt_Pwr, LongProfile_Pwr)
   !
   write(2,*) '------------------- PSF fit ------------------'
   J_max(:)=0. ; X_max(:)=0
   !J_max(1)= 0.114566E+09 ; X_max(1)=850. ; R_0(1)=0.2716 ; L_0(1)=53.4626 ; RL_Switch= .true.  ;  FitOptn=0
   !FitParams(1)=1 ;  FitParams(2)=4 ;  FitParams(3)=7 ;  FitParams(4)=10    ; FitParams(5)=-1
   !FitParams(1)=41 ;  FitParams(2)=7 ;  FitParams(3)=10 ;  FitParams(4)=-10    ; FitParams(5)=-1
   J_max(1)= 0.114566e+09 ; X_max(1)=850. ; R_0(1)=0.2716 ; L_0(1)=200.4626 ; RL_Switch= .true.  ;  FitOptn=0
   ! 1 : not a fit parameter since average current is fixed
   ! 2,3 : ratio of current of 2nd or 3rd bump compared to that of 1sst
   ! 4,5,6 : X_max of each bump
   ! 7 :  R or X_0 value taken equal for all 3 bumps
   ! 8,9 : not fitparameters
   ! 10 :  lam or L_0 value taken equal for all 3 bumps
   ! 11,12 : not fitparameters
   !FitParams(1)=2 ;  FitParams(2)=3  ;  FitParams(3)=4 ;  FitParams(4)=5 ;  FitParams(5)=6 ;  FitParams(6)=7  ; FitParams(7)=10
   !FitParams(8)=-10
   !GoodQual=0.3
   !GoodQual=0.002
   comment='PSF fit - '//ext
   Call MultFit_ShSruct(FitOptn, PSF_pow, GrdPSFWght, N_tr, N_zf, del_tr, del_z, z_low, J_max, X_max, R_0, L_0,  &
         RL_Switch, qual, error, RInt_PSF, LongProfile_PSF,   GoodQual, comment)
   !
   write(2,*) '------------------- power-Kernel fit ------------------'
   J_max(:)=0. ; X_max(:)=0
   J_max(1)= 0.114566E+09 ; X_max(1)=850. ; R_0(1)=0.2716 ; L_0(1)=200 ; RL_Switch= .true.  ;  FitOptn=3
   If(ext.eq.'r') FitOptn=13
   !GoodQual=0.3
   comment='power-Kernel fit - '//ext
   Call MultFit_ShSruct(FitOptn, powr, GrdPwrWght, N_tr, N_zf, del_tr, del_z, z_low, J_max, X_max, R_0, L_0, &
         RL_Switch, qual, error, RInt_Ker, LongProfile_Ker,   GoodQual, comment)
   ! =====================================
   !  write results to file
   OPEN(UNIT=4,STATUS='unknown',FILE=trim(Base)//ext//'-BeamPwr-PwrFit.dat')
   write(2,*) 'writing file "',trim(Base)//ext//'-BeamPwr-PwrFit.dat','" ',z_low, N_zf
   flush(unit=2)
   write(4,"(A)") '!    zeta_f, powr,     RInt_Pwr,    J_x(Pwr),  PSF_pow,     RInt_PSF,     J_x(PSF)'
   Do i_zf=0,N_zf ! Loop over  Retarded time
      zeta_f=i_zf*del_z + z_low   ! focal height
      d_tr=i_zf*del_z/del_tr
      i_tr=INT(d_tr)
      If(i_tr.ge.N_tr) exit
      d_tr=d_tr-i_tr
      write(4,"(F10.4,21(1pG13.4))") zeta_f/1000. &
         , powr(i_zf), RInt_Pwr(i_zf), (LongProfile_Pwr(i_tr)*(1.-d_tr)+d_tr*LongProfile_Pwr(i_tr+1))  &
         , PSF_pow(i_zf), RInt_PSF(i_zf), (LongProfile_PSF(i_tr)*(1.-d_tr)+d_tr*LongProfile_PSF(i_tr+1)) &
         , RInt_Ker(i_zf), (LongProfile_Ker(i_tr)*(1.-d_tr)+d_tr*LongProfile_Ker(i_tr+1))
   EndDo ! i_tr=0,N_zf ! Loop over  Retarded time
   Close(Unit=4)
   !DeAllocate( RInt_Ker, RInt_PSF, RInt_Pwr, LongProfile_Ker, LongProfile_PSF, LongProfile_Pwr )
   !
   PlotType='BeamFormCurr'
   PlotName=trim(Base)//'currents-'//ext
   PlotDataFile=trim(Base)//ext//'-'
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Plot made with earlier GLE command
   !
   Return
End Subroutine CurrentExtraction
!----------------------------------------------
Subroutine MultFit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0, &
         RL_Switch, qual, error, RInt_I, LongProfile,   GoodQual, comment)
! Fit X_max, EnergyCR, [R_0, L_0  --or-- X_0, lamx]
   use constants, only : dp
   Use ShowerStruct, only : ParDim
   implicit none
!    integer, intent(in) :: N_ant
   !Integer, parameter :: ParDim=3
   Integer, intent(in) :: FitOption
   !  =1 fit PSF
   !  =2 fit with weights
   !  =3 fit with full kernel
   Real(dp), intent(in) :: R_Int_data(0:N_z)
   Real(dp), intent(in) :: R_GrdWght(0:N_tr,0:N_z)
   Integer, intent(in) :: N_tr, N_z
   Real(dp), intent(in) :: del_tr, del_z, z_low
   Real(dp), intent(inout) :: X_max(1:ParDim), J_max(1:ParDim), R_0(1:ParDim), L_0(1:ParDim)
   logical, intent(in) :: RL_Switch
   Real(dp), intent(out) :: qual
   integer, intent(out) ::  error
   Real(dp), intent(out) :: RInt_I(0:N_z), LongProfile(0:N_tr)
   Real(dp), intent(in) :: GoodQual
   Character(len=*), intent(in) :: comment
   !
   Integer :: i, N_fit, FitParams(12)
   Real(dp) :: KfitX(1:10,1:13)
   ! 1 : not a fit parameter since average current is fixed
   ! 2,3 : ratio of current of 2nd or 3rd bump compared to that of 1sst
   ! 4,5,6 : X_max of each bump
   ! 7 :  R or X_0 value taken equal for all 3 bumps
   ! 8,9 : not fitparameters
   ! 10 :  lam or L_0 value taken equal for all 3 bumps
   ! 11,12 : not fitparameters
   !
   FitParams(:)=0
   FitParams(1)=4 ;  FitParams(2)=10 ;  FitParams(3)=-10
   Call Fit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0, FitParams, &
         RL_Switch, qual, error, RInt_I, LongProfile)
   !FitParams(1)=4 ;  FitParams(2)=7 ;  FitParams(3)=10 ;  FitParams(4)=-10
   !Call Fit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0, FitParams, &
   !      RL_Switch, qual, error, RInt_I, LongProfile)
   N_fit=1  ; KfitX(N_fit,1)= qual ;
   KfitX(N_fit,2:4)=J_max(1:3) ; KfitX(N_fit,5:7)=X_max(1:3) ; KfitX(N_fit,8:10)=R_0(1:3) ; KfitX(N_fit,11:13)=L_0(1:3)
   If(qual.gt. GoodQual) Then
      write(2,*) '------- second round ------------ ',TRIM(comment),' ------------------'
      ! 1st Fit with 2
      X_max(2) = X_max(1)+100.  ; J_max(2)=J_max(1)
      FitParams(1)=2 ; FitParams(2)=4 ; FitParams(3)=5 ;  FitParams(4)=10 ;  FitParams(5)=-10 ;  FitParams(6)=-7
      Call Fit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0, FitParams, &
            RL_Switch, qual, error, RInt_I, LongProfile)
      N_fit=2  ;  KfitX(N_fit,1)= qual ;
      KfitX(N_fit,2:4)=J_max(1:3) ; KfitX(N_fit,5:7)=X_max(1:3) ; KfitX(N_fit,8:10)=R_0(1:3) ; KfitX(N_fit,11:13)=L_0(1:3)
      write(2,*) '------- third round ------------ ',TRIM(comment),' ------------------'
      ! 2st Fit with 2
      N_fit=1  ! restore first:
      J_max(1:3)=KfitX(N_fit,2:4) ; X_max(1:3)=KfitX(N_fit,5:7) ; R_0(1:3)=KfitX(N_fit,8:10) ; L_0(1:3)=KfitX(N_fit,11:13)
      X_max(2) = X_max(1)-100.  ; J_max(2)=J_max(1)   ! shift x-max to higher altitude
      !FitParams(1)=2 ; FitParams(2)=4 ; FitParams(3)=5 ;  FitParams(4)=7 ;  FitParams(5)=10 ;  FitParams(6)=-10
      FitParams(1)=2 ; FitParams(2)=4 ; FitParams(3)=5 ;  FitParams(4)=10 ;  FitParams(5)=-10 ;  FitParams(6)=-10
      Call Fit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0, FitParams, &
            RL_Switch, qual, error, RInt_I, LongProfile)
      N_fit=3  ; KfitX(N_fit,1)= qual ;
      KfitX(N_fit,2:4)=J_max(1:3) ; KfitX(N_fit,5:7)=X_max(1:3) ; KfitX(N_fit,8:10)=R_0(1:3) ; KfitX(N_fit,11:13)=L_0(1:3)
      !write(2,*) '------- best selected ------------ ',TRIM(comment),' ------------------'
      ! select best
      i=MINLOC(KfitX(1:3,1),1)  ; Qual=KfitX(i,1)
      If((i.ne. 1) .and. (Qual.gt. GoodQual) ) Then  ! try even more complicated form
         write(2,*) '------- fourth round ------------ ',TRIM(comment),' ------------------'
         J_max(1:3)=KfitX(i,2:4) ; X_max(1:3)=KfitX(i,5:7) ; R_0(1:3)=KfitX(i,8:10) ; L_0(1:3)=KfitX(i,11:13)
         If(i.eq.2) Then
            X_max(3) = X_max(2)+100.  ; J_max(3)=J_max(2)   ! shift x-max to higher altitude
         Else
            X_max(3) = X_max(2)-100.  ; J_max(3)=J_max(2)   ! shift x-max to higher altitude
         EndIf
         FitParams(1)=2 ; FitParams(2)=3 ; FitParams(3)=4 ; FitParams(4)=5 ; FitParams(5)=6 ; FitParams(6)=10 ; FitParams(7)=-10
         !FitParams(1)=2 ; FitParams(2)=3 ; FitParams(3)=4 ; FitParams(4)=5 ; FitParams(5)=6 ; FitParams(6)=7 ; FitParams(7)=10
         Call Fit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0,FitParams,&
               RL_Switch, qual, error, RInt_I, LongProfile)
         N_fit=4  ; KfitX(N_fit,1)= qual ;
         KfitX(N_fit,2:4)=J_max(1:3) ; KfitX(N_fit,5:7)=X_max(1:3) ; KfitX(N_fit,8:10)=R_0(1:3) ; KfitX(N_fit,11:13)=L_0(1:3)
         !
         write(2,*) '------- fifth round ------------ ',TRIM(comment),' ------------------'
         J_max(1:3)=KfitX(i,2:4) ; X_max(1:3)=KfitX(i,5:7) ; R_0(1:3)=KfitX(i,8:10) ; L_0(1:3)=KfitX(i,11:13)
         X_max(3) = (X_max(1)+X_max(2))/2.  ; J_max(3)=(J_max(1)+ J_max(2))/2.   ! shift x-max to inbetween
         !FitParams(1)=2 ; FitParams(2)=3 ; FitParams(3)=4 ; FitParams(4)=5 ; FitParams(5)=6 ; FitParams(6)=7 ; FitParams(7)=10
         FitParams(1)=2 ; FitParams(2)=3 ; FitParams(3)=4 ; FitParams(4)=5 ; FitParams(5)=6 ; FitParams(6)=10 ; FitParams(7)=-10
         Call Fit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0,FitParams,&
               RL_Switch, qual, error, RInt_I, LongProfile)
         N_fit=5  ; KfitX(N_fit,1)= qual ;
         KfitX(N_fit,2:4)=J_max(1:3) ; KfitX(N_fit,5:7)=X_max(1:3) ; KfitX(N_fit,8:10)=R_0(1:3) ; KfitX(N_fit,11:13)=L_0(1:3)
      EndIf
      i=MINLOC(KfitX(1:N_fit,1),1)   ! refinement for r parameter
      J_max(1:3)=KfitX(i,2:4) ; X_max(1:3)=KfitX(i,5:7) ; R_0(1:3)=KfitX(i,8:10) ; L_0(1:3)=KfitX(i,11:13)
      FitParams(:)=0 ; FitParams(1)=4  ; FitParams(2)=7 ; FitParams(3)=10
      If(i.gt.1) Then
         FitParams(4)=2 ; FitParams(5)=5
         If(i.gt.3) Then
            FitParams(6)=3 ; FitParams(7)=6
         EndIf
      EndIf
      i=MINLOC(KfitX(1:N_fit,1),1)  ; Qual=KfitX(i,1)
      If(Qual.gt.1.) Then
         write(2,*) '------- round #',N_fit+1,'------------ ',TRIM(comment),' ------------------ based on',i
         Call Fit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0,FitParams,&
               RL_Switch, qual, error, RInt_I, LongProfile)
         N_fit=N_fit+1  ; KfitX(N_fit,1)= qual ;
         KfitX(N_fit,2:4)=J_max(1:3) ; KfitX(N_fit,5:7)=X_max(1:3) ; KfitX(N_fit,8:10)=R_0(1:3) ; KfitX(N_fit,11:13)=L_0(1:3)
         !
         X_max(3)=0.; X_max(1)=850.; X_max(2)=800; J_max(3)=0.; J_max(1)=+1.e+9; J_max(2)=-1.e+9; R_0(1:3)=0.2716; L_0(1:3)=100.
         FitParams(1)=2 ; FitParams(2)=4 ; FitParams(3)=5 ;  FitParams(4)=10 ;  FitParams(5)=-10 ;  FitParams(6)=-10
         write(2,*) '------- round #',N_fit+1,'------------ ',TRIM(comment),' ------------------'
         Call Fit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0,FitParams,&
               RL_Switch, qual, error, RInt_I, LongProfile)
         N_fit=N_fit+1  ; KfitX(N_fit,1)= qual ;
         KfitX(N_fit,2:4)=J_max(1:3); KfitX(N_fit,5:7)=X_max(1:3); KfitX(N_fit,8:10)=R_0(1:3); KfitX(N_fit,11:13)=L_0(1:3)
         FitParams(1)=2 ; FitParams(2)=4 ; FitParams(3)=5 ;  FitParams(4)=7 ;  FitParams(5)=10 ;  FitParams(6)=-10
         write(2,*) '------- round #',N_fit+1,'------------ ',TRIM(comment),' ------------------'
         Call Fit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0,FitParams,&
               RL_Switch, qual, error, RInt_I, LongProfile)
         N_fit=N_fit+1  ; KfitX(N_fit,1)= qual ;
         KfitX(N_fit,2:4)=J_max(1:3) ; KfitX(N_fit,5:7)=X_max(1:3) ; KfitX(N_fit,8:10)=R_0(1:3) ; KfitX(N_fit,11:13)=L_0(1:3)
      EndIf
      !
      i=MINLOC(KfitX(1:N_fit,1),1)
      write(2,*) 'Best fit is #', i
      Qual=KfitX(i,1) ; J_max(1:3)=KfitX(i,2:4) ; X_max(1:3)=KfitX(i,5:7) ; R_0(1:3)=KfitX(i,8:10) ; L_0(1:3)=KfitX(i,11:13)
      FitParams(1)=-2
      Call Fit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0, FitParams, &
            RL_Switch, qual, error, RInt_I, LongProfile)
   EndIf
   write(2,*) '------- summary of all fit results for ',TRIM(comment),' ------------------'
      Do i=1,N_fit
         !write(2,"(20G13.4)") 'KfitX(N_fit,1:13)=', i, KfitX(i,1:13)
         write(2,"(A,I2,A,F8.3,A,3G13.4,A,3F7.1, A,F7.2, A,F7.2)") 'Fit#', i, ', qual=', KfitX(i,1), '%, I_m=',KfitX(i,2:4) &
            ,', X_m=',KfitX(i,5:7), '[g/cm^], R/lambda=',KfitX(i,8),', L/X_0=', KfitX(i,11)
      EndDo
   write(2,*) '===================== end fitting for ',TRIM(comment),' ==================='
   Return
End Subroutine MultFit_ShSruct
