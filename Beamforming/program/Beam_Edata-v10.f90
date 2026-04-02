Module EData
   use constants, only : pi,dp, ci
contains
!-------------------------------------
Subroutine GetEData(Base, N_d, dt_PSF, ttracex_zf, ttracer_zf, ih_min, ih_max, z_low, del_z, N_zf, del_tr, N_tr &
      , AtmHei_dim, AtmHei_step, Xi)
   use FFT, only : FFTransform_su, DAssignFFT, RFTransform_CF, RFTransform_CF2RT
   Use Systemcommands, only : MakeGLEplots
   Implicit none
   Integer, intent(inout) :: N_d ! number of distances or E-files
   character(Len=*), intent(in)  :: Base
   Real(dp), intent(in) :: dt_PSF  !  Assumed sampling time [m]
   real(dp), allocatable, intent(out) :: ttracex_zf(:,:), ttracer_zf(:,:)
   Integer, intent(out) :: ih_min, ih_max
   Real(dp), intent(in) :: z_low
   Real(dp), intent(in) :: del_z
   Real(dp), intent(in) :: del_tr
   Integer, intent(in) :: N_zf
   Integer, intent(in) :: N_tr
   Integer, intent(in) :: AtmHei_dim
   real(dp), intent(in) :: AtmHei_step ! [m]
   Real(dp), intent(in) :: Xi(0:AtmHei_dim)
   !
   real(dp), allocatable :: Ex_ant(:,:), Ey_ant(:,:), Er_ant(:,:), d_ant(:), RNorm(:)
   Complex(dp), allocatable :: Cnux_f(:), Cnur_f(:)
   Complex(dp), allocatable :: Cnux_ant(:,:), Cnur_ant(:,:)
   real(dp) :: Sample_offset, dsample_offset  ! time offset in [sample]= t_offset/Del_t
   Integer :: isample_offset
   Integer :: i_d, N_t, nutrace_dim, i, i_t, i_atm, i_zf, nxx
   Character(len=3) :: Label
   Real(dp) :: Ant_d, del_t, t_offset, NN_f
   Integer :: i_nu_ini, i_nu_max
   Real(dp) :: EX_max(100)
   Complex(dp) :: dph_shift, ph_shift
   Complex(dp), parameter ::    ipi=ci*pi  ! This is the real constant i*pi
   real(dp) :: dd, dt, t_shift
   Integer :: i_phi, N_Phi0
   Real(dp) :: zeta_f, Phi_f, d_phi0, rs_f, Factr, weight,  R_f, rc_a, d_zeta, rc_a_l, rc_a_u
   Character(len=30) :: PlotType, PlotName, PlotDataFile
   !-----------------------------
   !  Reading E-field time traces
   !
   !Del_t = 5*c_l  ! convert 5 ns to meters, read-in
   !N_d=50  !  Number of distances
   nxx=0
   Do i_d=1,N_d
      Write(label,"(I2.2)") i_d
      OPEN(UNIT=4,STATUS='unknown',FILE=trim(Base)//'th_000'//trim(label)//'.csv',iostat=nxx) ! E13.4
      !write(2,*) 'reading file: "', trim(Base)//'th_000'//trim(label)//'.csv'
      If(i_d.eq.1) Then
         If(nxx.ne.0) Then
            write(2,*) 'problems reading file:', trim(Base)//'th_00001.csv'
            STOP TRIM(Base)//'th_00001.csv is problematic'
         EndIf
         Read(4,*) label, Ant_d, N_t, del_t, t_offset, i_nu_ini, i_nu_max
         write(2,*) '*** check dt_PSF,del_t', N_t, dt_PSF, del_t, dt_PSF-del_t
         nutrace_dim=N_t/2
         allocate(Ex_ant(1:N_t,1:N_d), Ey_ant(1:N_t,1:N_d), Er_ant(1:N_t,1:N_d) )
         allocate( Cnux_ant(0:nutrace_dim,1:N_d), Cnur_ant(0:nutrace_dim,1:N_d) )
         allocate( d_ant(0:N_d+1) )
         d_ant(0)=0
      Else
         If(nxx.ne.0) exit
         Read(4,*) label, Ant_d
      EndIf
      d_ant(i_d)=Ant_d
      Do i_t=1,N_t
         Read(4,*) i, Ex_ant(i_t,i_d), Ey_ant(i_t,i_d), Er_ant(i_t,i_d)
      End Do ! i_t=1,N_t
      Close(Unit=4)
      EX_Max(i_d)=MaxLoc(ABS(Ex_ant(:,i_d)),1)
   End Do ! i_d=1,N_d
   If(i_d .lt. N_d ) N_d=i_d
   d_ant(N_d+1)=d_ant(N_d)
   !
   !------------------------------------
   !  FFT, Go to frequency space
   Call FFTransform_su(N_t)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !
   Do i_d=1,N_d
      Call RFTransform_CF(Ex_ant(1,i_d), Cnux_ant(0,i_d))  ! transform to frequency
      Call RFTransform_CF(Er_ant(1,i_d), Cnur_ant(0,i_d))  ! transform to frequency
      ! Call RFTransform_CF2RT(K_nu(0, i_tr),K_h(1-nudim_new, i_tr)) !transform to real-time trace
      !write(2,*) Sample_offset,i_d, Ex_ant(iSample_offset-1:iSample_offset+5,i_d)*1.d-5
   End Do ! i_d=1,N_d
   Deallocate(Ex_ant, Ey_ant, Er_ant )
   !
   !------------------------------------------
   !
   N_phi0 =1
   nutrace_dim=N_t/2
   !Sample_offset = t_offset/del_t - 0.5 ! t=0 point at sample t_offset/del_t+1 while in PSF t=0 point at sample +0.5
   !
   !Sample_offset = t_offset/del_t -1.0 ! t=0 point at sample t_offset/del_t+1 while in PSF t=0 point at sample +0.5
   Sample_offset = t_offset/del_t -1.25 ! t=0 point at sample t_offset/del_t+1 while in PSF t=0 point at sample +0.5
   isample_offset=NINT(Sample_offset) ! t=0 point at sample t_offset/del_t+1 while in PSF t=0 point at sample +0.5
   ! a more negative constant moves the data to the right w.r.t. the PSF model
   ! Thus let ttrace_zf run from [(-isample_offset+0.5)*del_t : N_t-isample_offset-0.5)*del_t] in steps of [del_t]
   dsample_offset=Sample_offset - isample_offset
   Allocate( Cnux_f(0:nutrace_dim), ttracex_zf( -isample_offset:N_t-isample_offset-1,0:N_zf), RNorm(0:N_zf) )
   Allocate( Cnur_f(0:nutrace_dim), ttracer_zf( -isample_offset:N_t-isample_offset-1,0:N_zf) )
   ! calculate norm consistent with that of the Kernel
   rc_a_l=145 ;  rc_a_u=155  ! S15/50
   !rc_a_l=147.5 ;  rc_a_u=152.5  ! S30/100
   !rc_a_l=220 ;  rc_a_u=230  ! S45/100

!   write(2,*) 'rc_a limits:', rc_a_l, rc_a_u, Sample_offset
   RNorm(:)=0.
   Do i_zf=0,N_zf
      zeta_f=i_zf*del_z + z_low   ! focal height
      Do i_d=1, N_d
         rc_a=d_ant(i_d) ! Distance of antenna to shower axis
!         If((rc_a.lt.rc_a_l) .or. (rc_a.gt.rc_a_u)) cycle
!         write(2,*) 'rc_a:',i_d, rc_a, (d_ant(i_d+1)-d_ant(i_d-1))/2.
         R_f=sqrt(rc_a*rc_a +zeta_f*zeta_f)
         weight=2.*rc_a*(d_ant(i_d+1)-d_ant(i_d-1))/2.   ! Integration area ring
         RNorm(i_zf)=RNorm(i_zf) + weight/R_f  !   /(R_0*R_0)
      EndDo
   EndDo ! idi=1, idi_N
   ! Calculate R from 'measured' fields
   Do i_zf=0,N_zf ! Loop over focal height on the shower axis
      zeta_f=i_zf*del_z + z_low   ! focal height
      ! ttrace is calculated in relatively small steps of del_z=200 [m]
      !write(2,*) 'making ttrace for zeta_f:',i_zf,zeta_f,N_d, del_z, z_low
      !Flush(unit=2)
      Cnux_f(:)=0.
      Cnur_f(:)=0.
      i_atm= zeta_f/AtmHei_step
      NN_f=1.+ xi(i_atm)
      !Do i_d=20, 20  ! sum over all antenna distances
      Do i_d=1, N_d  ! sum over all antenna distances
         rc_a=d_ant(i_d) ! Distance of antenna to shower axis
!         If((rc_a.lt.rc_a_l) .or. (rc_a.gt.rc_a_u)) cycle
         phi_f=0. ; i_phi=0  ; d_phi0=pi ; rs_f=0. ! limit to shower axis
         weight=2.*rc_a*d_phi0*(d_ant(i_d+1)-d_ant(i_d-1))/2.   ! Integration area ring
         !Do i_phi=1,N_phi0  ! angle between ObsDist and rs_f
         !   phi_f=(i_phi-0.5)*d_phi0
         !   If(N_phi0.eq.1) phi_f=0.
         !R_f=sqrt(zeta_f*zeta_f+ rc_a*rc_a + rs_f*rs_f - 2*rs_f*rc_a*cos(phi_f)) ! Distance between antenna and source in shower
         R_f=sqrt(zeta_f*zeta_f+ rc_a*rc_a ) ! Distance between antenna and source in shower
         t_shift=(NN_f*R_f -zeta_f)/del_t + dsample_offset !+ 5. !  + h_f ! time offset in antenna trace [m]
         !write(2,*) '!t_shift diff:', (NN_f*R_f -zeta_f)/del_t
         Factr=weight/(NN_f*R_f*RNorm(i_zf))    ! Normalisation?????? maybe the same as for weights
         !
         ! Shift antenna trace and add (in frequency space)
         dph_shift=exp(-ipi*t_shift/nuTrace_dim)
         ph_shift= exp(-ipi*t_shift*i_nu_ini/nuTrace_dim)  ! factor -CI to match sign template
         DO i=i_nu_ini,i_nu_max
            Cnux_f(i)=Cnux_f(i) + Cnux_ant(i,i_d) * Factr * ph_shift
            Cnur_f(i)=Cnur_f(i) + Cnur_ant(i,i_d) * Factr * ph_shift
            ph_shift=ph_shift*dph_shift
         ENDDO
         !EndDo ! i_phi=1,N_phi0  ! angle between ObsDist and rs_f
      EndDo ! i_d=1, N_d  ! sum over all antenna distances
      ! At this point the beamformed pulse is determined for height zeta_f
      ! transform to time trace
      Call RFTransform_CF2RT(Cnux_f(0),ttracex_zf( -isample_offset,i_zf)) !transform to real-time trace, Should have t=0 point at isample_offset
      Call RFTransform_CF2RT(Cnur_f(0),ttracer_zf( -isample_offset,i_zf)) !transform to real-time trace, Should have t=0 point at isample_offset
      !If(zeta_f.ge.5400. .and. zeta_f.le.5600) Then
      !   write(2,"(2F9.2,A,10g13.3)") zeta_f, (NN_f*R_f -zeta_f)/del_t, &
      !   ', ttrace_zf(-1,+2):', ttrace_zf( -1:2,i_zf)
      !EndIf
   EndDo ! i_tr=0,N_tr ! Loop over focal height on the shower axis
   Call DAssignFFT()      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   DeAllocate( Cnux_f, Cnur_f )
   !
   ih_min = -isample_offset
   ih_max = N_t-isample_offset-1
   !write(2,*) 'h-dims(EData):', ih_min, ih_max, N_t
   !flush(unit=2)
   !------------------------------------------------
   !  Make data file for plotting
   OPEN(UNIT=4,STATUS='unknown',FILE=trim(Base)//'x-BeamPwr-tt.dat')
   write(2,*) 'writing file "',trim(Base)//'x-BeamPwr-tt.dat','"'
   flush(unit=2)
   write(4,"(A,17I4)") '!    i_h,  ttracex_zf(i_h,0:N_zf) i_zf=      ', (NINT((i_d*1000.-z_low)/del_z), i_d=1,17)
   write(4,"(A,17I4)") '!    i_h,  ttracex_zf(i_h,0:N_zf) zeta_f[km]=',(i_d, i_d=1,17)
   Do  i=ih_min,ih_max
      write(4,"(22(1pG13.4))") (i-0.5)*del_t, (ttracex_zf(i,NINT((i_d*1000.-z_low)/del_z)), i_d=1,17)!, PSF_trace(i)
   EndDo ! i_zf=0,N_zf
   Close(Unit=4)
   PlotType='Test_kernel_1-3'
   PlotName=trim(Base)//'testKernel_x'
   PlotDataFile=trim(Base)//'x-'
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Submit, CleanPlotFile
   !
   PlotType='W_PSF'
   PlotName=trim(Base)//'WPSF_x'
   PlotDataFile=trim(Base)//'x-'
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Submit, CleanPlotFile
   !
   OPEN(UNIT=4,STATUS='unknown',FILE=trim(Base)//'r-BeamPwr-tt.dat')
   write(2,*) 'writing file "',trim(Base)//'r-BeamPwr-tt.dat','"'
   flush(unit=2)
   write(4,"(A,17I4)") '!    i_h,  ttracer_zf(i_h,0:N_zf) i_zf=      ', (NINT((i_d*1000.-z_low)/del_z), i_d=1,17)
   write(4,"(A,17I4)") '!    i_h,  ttracer_zf(i_h,0:N_zf) zeta_f[km]=',(i_d, i_d=1,17)
   Do  i=ih_min,ih_max
      write(4,"(22(1pG13.4))") (i-0.5)*del_t, (ttracer_zf(i,NINT((i_d*1000.-z_low)/del_z)), i_d=1,17)!, PSF_trace(i)
   EndDo ! i_zf=0,N_zf
   Close(Unit=4)
   PlotType='Test_kernel_1-3'
   PlotName=trim(Base)//'testKernel_r'
   PlotDataFile=trim(Base)//'r-'
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Submit, CleanPlotFile
   !
   PlotType='W_PSF'
   PlotName=trim(Base)//'WPSF_r'
   PlotDataFile=trim(Base)//'r-'
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Submit, CleanPlotFile
   !
   Return
End Subroutine GetEData
!-------------------------------
!------------------------
End Module EData

!       w3IfZYkpeQLscX4WkfAS
