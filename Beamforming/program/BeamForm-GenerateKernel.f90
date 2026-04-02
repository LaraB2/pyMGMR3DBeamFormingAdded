    Include 'Constants.f90'
    Include 'MakeGLEplots.f90'
    Include 'MGMR3D_FFT.f90'
    Include 'MGMR3D2_Atmosphere.f90'
    Include 'MGMR3D2_PancakeFunct.f90'
    Include 'MGMR3D2_LateralDistr.f90'
    Include 'KernelFilter.f90'
    Include 'MGMR3D_spline.f90'
    Include 'Beam_DataRead-v20.f90'
    Include 'GenerateKernel.f90'
!    Include 'Beam_PSF-v20.f90'
    Include 'Beam_Kernel.f90'
!    Include 'Beam_DataBeaming.f90'
    !
!    Include 'Beam_ShStruct-v10.f90'
!    Include 'Beam_Folding.f90'
    !Include 'Beam_Edata-v10.f90'
    !
    !Include 'C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI\FORTRANsrc/nl2sol.f90' !Take this out when using the LINUX makefile
    Include '../LightningPros/nl2sol.f90' !Take this out when using the LINUX makefile
!------------------------------
!============================================
	program BeamFormingGenerateKernel
   use constants, only : pi, dp, c_l, ci
   !Use CurrentFolding, only : UnFoldJ
   Use Systemcommands, only : MakeGLEplots
   Use Atmosphere, only : AtmHei_dim, AtmHei_step
   Use Atmosphere, only : Xi, PenDepth
   Use ShowerData, only :  GroundLevel_sh, ZenithAngle_deg, ZenithAngle_rad, N_antennas, N_timesamples, SamplingTime, DistMax
   !Use PSF, only : GrdPSFWghtr, GrdPwrWghtr, PSF_tracer
   use Atmosphere, only : AtmosphereInit
	implicit none
   Integer, parameter :: PSF_max=25
   INTEGER DATE_T(8) !,WLength
   CHARACTER*12 REAL_C(3)
   character(Len=6)  :: date,time,dateId(100),timeId(100)
!
   Real(dp) :: Theta_zen, Cos_Zenith
   Integer :: nxx, k, i, i_d, N_d, i_t, N_t, nuTrace_dim, i_zf, N_zf, i_atm, i_tr, i_tr2, N_tr
   real(dp) :: z_low, Z_high, Del_z, Del_d, del_tr, d_tr  !  Del_t,
   real(dp), allocatable :: ttracex_zf(:,:), ttracey_zf(:,:)
   Real(dp) :: dt_PSF, zeta_f, d_zeta
   Real(dp), allocatable :: powr_x(:), powr_y(:), PSF_pow_x(:), PSF_pow_y(:), Jx1PSF(:), Jx2PSF(:), Jx1Pwr(:), Jx2Pwr(:)
   Character(len=2) :: Label
   Character(len=30) :: PlotType, PlotName, PlotDataFile
   Integer :: i_r, i_r1, N_r, N_r1, Iresol
   Real(dp) :: del_r, del_r1, d_r1, d_r, zeta_K(1:PSF_max), rs_0, phi_0, theta_0
   Real(dp) :: EnergyCR, X_max, R_0, L_0, qual
   integer :: N_fitpar, error, FitOptn, ih_min, ih_max, N_phi0
   logical RL_Switch
   Real(dp), allocatable :: RInt_PSF(:), RInt_Pwr(:), LongProfile_PSF(:), LongProfile_Pwr(:)
   Real(dp), allocatable :: RInt_Ker(:), LongProfile_Ker(:)
   Logical :: HoriShwr, AntennaLayout, ReCreate
   !Character(len=10) :: Label2
   Character(len=30) :: Base_pol
   character(Len=1) :: ext
   Real(dp) :: X,Y, z , R, S, W ! scratch
   Integer :: N_PSF
   Character(len=10), save :: Label2
   Character(len=13), save :: Label3(1:PSF_max)
   Character(len=20), save :: Base, DataBase
   character(len=50), save :: NamdeHDF5File
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
    NAMELIST /BeamPars/ Base, DataBase, AntennaLayout, del_z,  del_tr, &
      zeta_K,ZenithAngle_deg,NamdeHDF5File, &  ! rs_0, N_phi0, phi_0, &
      N_h, d_h, dN_tr, N_phir, d_crMax, DistMax
!
   !=========================================================
   !  Run parameters
   Base='plot/'
   DataBase=''
   !AntennaLayout=.true.  ! use realistic antenna layout
   del_z =20. ![AtmHei_step] step size for focal-height grid
   del_tr = 20. ![AtmHei_step] step size for shower-axis distance grid
!   ReCreate = .false. ! recreate kernell
   zeta_K(:)=0.
   rs_0 = 0.   !Distance focal point off shower axis
   theta_0 = 0. !Angle where focal point off shower axis, measured as angle from the positive V x B axis.
   write(*,*) "Taking RS_0 =", rs_0
   N_phi0 = 1
   phi_0 = 0.  ! Angle of off-axis point w.r.t. antenna
   DistMax = 0.  ! will be set automatically is antenna layout is read-in (AntennaLayout=.true. )
   read(*,NML = BeamPars)   !Reads all the parameters in for the beam from the beamforming .in file.
   If(DataBase .eq. '') DataBase=Base
   !
   !===================================================
   OPEN(UNIT=2,STATUS='unknown',FILE=TRIM(Base)//'BeamForm-GenK.out')
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
      !   OPEN(UNIT=4,STATUS='old',FILE=trim(Base)//trim(label)//'.dat',iostat=nxx) ! E13.4    Zen_sh = 43.
        ZenithAngle_rad=ZenithAngle_deg*pi/180.

      Call CorePosDistances()  ! true LOFAR antennas
   Else
      !   OPEN(UNIT=4,STATUS='unknown',FILE=trim(Base)//'th_000'//trim(label)//'.csv',iostat=nxx) ! E13.4
      Call ReadThetData(DataBase, DistMax) ! From theta grid
   EndIf
   write(2,*) 'DistMax=',DistMax,'[m]'
   !
   !-----------------------------------------------------
   ! Initialize penetration depth and index of refraction
   HoriShwr=.false.
   write(*,*) "Zenith angle =", ZenithAngle_rad
   Call AtmosphereInit(ZenithAngle_rad, GroundLevel_sh, HoriShwr)
   !
   !
   !------------------------------------------
   ! When (ReCreate=.true.) generate kernel for selected focal heights and write to file
   !  Use Shower angle and antenna layout as used in the data
   !
   ReCreate=.true.
   Call GenerateKernel(Base, Label2, Label3, PSF_max, N_PSF, ReCreate, zeta_K, rs_0, N_phi0, phi_0, &
      N_h, d_h, dN_tr, N_phir, d_crMax, theta_0,AntennaLayout)
   !------------------------------------
   !
End Program BeamFormingGenerateKernel
!------------------------
