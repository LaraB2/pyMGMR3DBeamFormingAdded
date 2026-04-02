!    Include 'MGMR3D_BA-v4.f90'
! ----------------------------------------------------------------
Subroutine GenerateKernel(PlotBase, Label2, Label3, PSF_max, N_PSF, ReCreate, zeta_K, rs_0, N_phi0, phi_0, &
       N_h, d_h, dN_tr, N_phir, d_crMax,theta_0,AntennaLayout)
!
   use constants, only : rh0,dp,pi
   Use Systemcommands, only : MakeGLEplots
!   use eventdata, only : PlotBase
   use Atmosphere, only : AtmHei_step, AtmHei_dim, xi,dxi,PenDepth  !,Ix,Iy,IQ,IX_int,Iy_int,alpha_tr
   Use ShowerData, only :  ZenithAngle_deg
   Use Pancakefunction, only : CoreDist_Dim, PancakeThi, CoreDist_A, alpha_tr
   use Pancakefunction, only : PancakeInit, SelectFh
   use LateralDistribution, only : LateralInit
   implicit none
   character(Len=*), intent(in)  :: PlotBase
   Character(len=10), intent(out) :: Label2
   Character(len=13), intent(out) :: Label3(1:PSF_max)
   Integer, intent(in) :: PSF_max
   Integer, intent(out) :: N_PSF
   real(dp), intent(inout) :: zeta_K(1:PSF_max)
    real(dp), intent(inout) :: phi_0, rs_0, theta_0
    integer, intent(inout) :: N_phi0
   Integer, intent(in) :: N_h
   Real(dp), intent(in) :: d_h
   Integer, intent(in) :: dN_tr
   Integer, intent(in) :: N_phir
   Real(dp), intent(in) :: d_crMax
!    logical first
!    integer, parameter :: NAnt_max=1  ! 900
      real(dp), parameter :: XDepAlpha=0.0d0
      real(dp), parameter :: MoliereRadius = 0.  !  Use default
      real(dp), parameter :: lam_tc=0.05d0  ! july 202, moved out of fit common block
      real(dp), parameter :: lam_100=7.d0  ! july 202, moved out of fit common block, important
      Integer, parameter :: ObsDist_dim=75
      real(dp), parameter :: ObsDist_Step=10. ! [m]
!
   real(dp) :: zeta_f, phi_0_init !, N_phi0, phi_0
   real(dp) :: zeta_max, zeta_min, zeta_m
   Logical :: ReCreate, AntennaLayout
    INTEGER DATE_T(8),i,k ,nxx,idi, CD_i ! ,inu,ith,Nth, Ntr
    CHARACTER*12 REAL_C (3)
    character*80 line
    real(dp) :: r,z,T_o,t_r,nu
    real(dp) :: zeta
	 Character*80 :: lname
	 Character*20 :: label
!
!===========================================================================================
    CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T) !-----------------------------
    write(*,232) 'Initialize',(DATE_T(i),i=5,8)               !-----------------------------
    write(2,232) 'Initialize',(DATE_T(i),i=5,8)               !-----------------------------
232 FORMAT(3X,'start ',A,' @ ',I2,':',I2,':',I2,'.',I3,1X,25(1H-))
    Flush(unit=2)
!===========================================================================================
    ! ------------------------------------------------------------------------------------------------------
!
    write(2,201) rh0*10000,MoliereRadius, SelectFh, lam_tc, lam_100
201 format('Refractivity at sea level=',F8.3,'x 10^{-4}',&
      ', Moliere radius=',f6.1,'[m], f(h) shape',i2,' with L=',f4.2,'[m]' , &
      ', L@100m=',f5.2,'[m]')
!
   !-----------------------------------------------
   !
   !   Initialize Radial dependence of pancake thickness & normalisation
   Call PancakeInit(XDepAlpha, PenDepth, AtmHei_dim, ObsDist_dim, ObsDist_Step, lam_100, lam_tc)
   Write(2,"('parametrization radial dependence of \lambda; @r=0:',F5.2,', @r=',f4.1,'m:',F5.2,', @r=',f6.1,'m:',F6.1)") &
        lam_tc,CoreDist_A(1),PancakeThi(1),ObsDist_dim*ObsDist_Step, lam_tc
   write(2,*) 'Alpha_tr set equal to ',alpha_tr(AtmHei_dim/2)
   !-----------------------------------------------
   !
   !   Initiate the lateral distribution function
   k = AtmHei_dim*2/3
   Call  LateralInit(MoliereRadius, k, AtmHei_step, alpha_tr)
   !-----------------------------------------------
   !
   ! ------------------------------------------------------------------------------------------------------
   ! nu_min=30.  ; nu_max= 80. ; SamplingTime=5 ! [ns]
   ! StParRange=11  ! Stokes parameter sampling
   !
   !  Calculate Beamforming kernel for selected fucus distances
   Call PrepareKernel(Label2, Label3, PSF_max, N_PSF, zeta_K, rs_0, phi_0)
   !
   write(2,*) '!ReCreate?', ReCreate, rs_0, N_phi0, phi_0
   write(2,"(A,30F7.3)") '!zeta_K=',zeta_K(:)
      flush(unit=2)
   !If( zeta_min .lt. 100.) zeta_min= 100.
   !N_psf = 2+ (zeta_max-zeta_min)/zeta_step
   !If(N_psf .gt. PSF_max) N_psf = PSF_max
   !zeta_max = zeta_min + (N_psf-0.5)*zeta_step
   zeta_min=zeta_K(1)*1000. -10.
   zeta_max=maxval(zeta_K(:))*1000. +50.
   phi_0_init = phi_0
   If(N_phi0.lt.1) N_phi0=1
   !
   !$OMP PARALLEL PRIVATE(zeta_m, zeta_f,phi_0,label,REAL_C,DATE_T)
   !$OMP DO
   Do k=1,N_PSF
      zeta_f = zeta_K(k)*1000.
      If(k+10.gt.PSF_max) Then
         zeta_m=zeta_max
      Else
         zeta_m=zeta_K(k+10)*1000. +50.
         If(zeta_m .lt. zeta_f) zeta_m=zeta_max
      EndIf
      !flush(unit=2)
      If(ReCreate) Then
         CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T)  !-----------------------------
         write(*,233) 'Generate Kernel @', zeta_f, (DATE_T(i),i=5,8)  !-----------------------------
   233 FORMAT(3X,'start ',A,F6.0,' @ ',I2,':',I2,':',I2,'.',I3,1X,25(1H-))
         write(2,232) 'Generate Kernel',(DATE_T(i),i=5,8)          !-----------------------------
         phi_0=phi_0_init*pi**k
         label=TRIM(label2)//TRIM(Label3(k))
         write(2,*) 'SourceKernel,zeta_f,rs:',zeta_f,rs_0, N_phi0, phi_0, label
         if (AntennaLayout) Then
            Call KernelCreateIndivid(PlotBase, zeta_f, rs_0, zeta_min, zeta_m, N_phi0, phi_0, label, &
                N_h, d_h, dN_tr, N_phir, d_crMax, theta_0)
        else
            Call KernelCreate(PlotBase, zeta_f, rs_0, zeta_min, zeta_m, N_phi0, phi_0, label, &
                N_h, d_h, dN_tr, N_phir, d_crMax)
            EndIf
      EndIf
   EndDo ! while (nxx.eq.0)
     !$OMP END PARALLEL
   Call MakeGLEplots(Submit=.true.)
   If(zeta_K(2).lt.0.) Then
   ! Just for testing purposes
   !
      stop 'zeta_K(2)<0'
   EndIf
!===========================================================================================
End Subroutine GenerateKernel
