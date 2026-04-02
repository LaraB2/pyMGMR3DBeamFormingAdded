Module DataBeam
contains
Subroutine DataBeaming(Base, N_antennas, N_timesamples, ttracex_zf, ttracey_zf, ih_min, ih_max, z_low, del_z, N_zf, del_tr, N_tr)
!  Calculate time traces for different heights along the shower axis for beam-forming using all antennas
!     for a realistic antenna layout.
   use constants, only : dp, pi, ci, c_l
   Use ShowerData, only : core_pos, d_ant
   Use ShowerData, only : DistMax, Time_maxx, PulseTime_x
   Use ShowerData, only : ant_pos, ant_toffset, ant_xtrace, ant_ytrace, Sample_Offset
   Use ShowerData, only : SamplingTime, nu_min, nu_max
   use FFT, only : FFTransform_su, DAssignFFT, RFTransform_CF, RFTransform_CF2RT
   use Atmosphere, only : AtmHei_step, AtmHei_dim, xi,dxi,PenDepth
   Use Systemcommands, only : MakeGLEplots
   Implicit none
   Integer, intent(in) :: N_antennas ! number of antennas
   Integer, intent(in) :: N_timesamples ! Length of trace
   character(Len=*), intent(in)  :: Base
   !Real(dp)  :: dt_PSF  !  Assumed sampling time [m]
   real(dp), allocatable, intent(out) :: ttracex_zf(:,:), ttracey_zf(:,:)
   Integer, intent(out) :: ih_min, ih_max
   Real(dp), intent(in) :: z_low
   Real(dp), intent(in) :: del_z
   Real(dp), intent(in) :: del_tr
   Integer, intent(in) :: N_zf
   Integer, intent(in) :: N_tr
   !
   real(dp), allocatable :: Ex_ant(:,:), Ey_ant(:,:), Er_ant(:,:)
   !
   Complex(dp) :: Cnux_f(0:N_timesamples/2), Cnuy_f(0:N_timesamples/2)
   Complex(dp) :: Cnux_ant(0:N_timesamples/2,1:N_antennas), Cnuy_ant(0:N_timesamples/2,1:N_antennas)
   real(dp) :: RNorm ! time offset in [sample]= t_offset/Del_t
   !Integer :: isample_offset
   Integer :: i_ant, nutrace_dim, i, k, i_t, i_atm, i_zf, nxx
   Character(len=3) :: Label
   Real(dp) :: del_t, t_offset, NN_f, NuSample
   Integer :: i_nu_ini, i_nu_max
   Real(dp) :: EX_max(100)
   Complex(dp) :: dph_shift, ph_shift
   Complex(dp), parameter ::    ipi=ci*pi  ! This is the real constant i*pi
   real(dp) :: dd, dt, t_shift, Zeta_step
   Integer :: i_phi, N_Phi0
   Real(dp) :: zeta_f, Phi_f, d_phi0, rs_f, Factr, weight,  R_f, rc_a, d_zeta, rc_a_l, rc_a_u
   Character(len=30) :: PlotType, PlotName, PlotDataFile
   !
   !------------------------------------
   !  FFT, Go to frequency space
   Call FFTransform_su(N_timesamples)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !
   Do i_ant=1,N_antennas
      Call RFTransform_CF(ant_xtrace(1,i_ant), Cnux_ant(0,i_ant))  ! transform to frequency
      Call RFTransform_CF(ant_ytrace(1,i_ant), Cnuy_ant(0,i_ant))  ! transform to frequency
   End Do !
   !
   !nutrace_dim=N_timesamples/2
   !NuSample=1000.*c_l/(nutrace_dim*SamplingTime*2.)    ! in [MHz]
   !write(2,*) '!Frequency spectra'
   !Do i=0,nutrace_dim
   !   write(2,"(I4,3g13.4,20g13.4)") i,500.*i/nutrace_dim,i*NuSample &
   !      , (ABS(Cnux_ant(i,i_ant)),i_ant=1,N_antennas,30)
   !EndDo
   !------------------------------------------
   !
   N_phi0 =1
   d_phi0= pi/N_phi0
   del_t = SamplingTime
   nutrace_dim=N_timesamples/2
   NuSample=1000.*c_l/(nutrace_dim*SamplingTime*2.)    ! in [MHz]
   i_nu_ini= nu_min/NuSample  ; i_nu_max=nu_max/NuSample+0.5
   if(i_nu_ini .lt. 0) i_nu_ini=0
   if(i_nu_max .ge. nutrace_dim) i_nu_max=nutrace_dim
   !Sample_offset = t_offset/del_t - 0.5 ! t=0 point at sample t_offset/del_t+1 while in PSF t=0 point at sample +0.5
   !
   !Sample_offset = t_offset/del_t -1.0 ! t=0 point at sample t_offset/del_t+1 while in PSF t=0 point at sample +0.5
   !Sample_offset = t_offset/del_t -1.25 ! t=0 point at sample t_offset/del_t+1 while in PSF t=0 point at sample +0.5
   !isample_offset=NINT(Sample_offset) ! t=0 point at sample t_offset/del_t+1 while in PSF t=0 point at sample +0.5
   ! a more negative constant moves the data to the right w.r.t. the PSF model
   ! Thus let ttrace_zf run from [(-isample_offset+0.5)*del_t : N_t-isample_offset-0.5)*del_t] in steps of [del_t]
   !dsample_offset=Sample_offset - isample_offset
   ih_min = -sample_offset
   ih_max = N_timesamples-sample_offset-1
   Allocate(  ttracex_zf( ih_min:ih_max,0:N_zf) )
   Allocate(  ttracey_zf( ih_min:ih_max,0:N_zf) )
   ! calculate norm consistent with that of the Kernel

!   write(2,*) 'rc_a limits:', rc_a_l, rc_a_u, Sample_offset ! ant_toffset
   ! Calculate R from 'measured' fields
   ttracex_zf(:,:) =0.
   ttracey_zf(:,:) =0.
   Do i_zf=0,N_zf ! Loop over focal height on the shower axis
      zeta_f=i_zf*del_z + z_low   ! focal height
      ! ttrace is calculated in relatively small steps of del_z=200 [m]
      !write(2,*) 'making ttrace for zeta_f:',i_zf,zeta_f,N_d, del_z, z_low
      !Flush(unit=2)
      Cnux_f(:)=0.
      Cnuy_f(:)=0.
      i_atm= zeta_f/AtmHei_step
      NN_f=1.+ xi(i_atm)
      RNorm=0.
      Do i_ant=1,N_antennas
         rc_a = d_ant(i_ant)    ! Distance of antenna to shower axis
         R_f=NN_f*sqrt(rc_a*rc_a +zeta_f*zeta_f)
         weight=2.*rc_a   ! Integration area ring of some arbitrary width
         RNorm=RNorm + weight/R_f  !   /(R_0*R_0)
      EndDo ! i_ant=1,N_antennas
      !RNorm=RNorm * 10.*1000.
      !write(2,*) '!RNorm:', i_zf, RNorm
      !Flush(Unit=2)
      Do i_ant=1,N_antennas  ! sum over all antenna distances
         rc_a=d_ant(i_ant) ! Distance of antenna to shower axis
         !phi_f=0. ; i_phi=0  ; d_phi0=pi ; rs_f=0. ! limit to shower axis
         weight=2.*rc_a *d_phi0      ! Integration area ring of some arbitrary width
         !Do i_phi=1,N_phi0  ! angle between ObsDist and rs_f
         !   phi_f=(i_phi-0.5)*d_phi0
         !   If(N_phi0.eq.1) phi_f=0.
         !R_f=sqrt(zeta_f*zeta_f+ rc_a*rc_a + rs_f*rs_f - 2*rs_f*rc_a*cos(phi_f)) ! Distance between antenna and source in shower
         R_f=NN_f*sqrt(zeta_f*zeta_f+ rc_a*rc_a ) ! Distance between antenna and source in shower
         !  Move each spectrum to the left (lowering time) by this much:
         t_shift=  ((R_f -zeta_f - ant_toffset(i_ant))/del_t - Sample_Offset ) !  + h_f ! time offset in antenna trace [samples]
         !write(2,*) '!t_shift diff:', (NN_f*R_f -zeta_f)/del_t
         Factr=weight/(R_f*RNorm)    ! Normalisation
         !
         ! Shift antenna trace and add (in frequency space)
         dph_shift=exp(-ipi*t_shift/nuTrace_dim)
         ph_shift= exp(-ipi*t_shift*i_nu_ini/nuTrace_dim)  ! factor -CI to match sign template
         !write(2,*) '!Factr', i_ant, Factr, ph_shift, dph_shift, i_nu_ini,i_nu_max
         !If(zeta_f.ge.6000. .and. zeta_f.le.6500) Then
         !   !k = Time_maxx(i_ant) -t_shift -(ant_toffset(i_ant)-ant_toffset(1))/del_t+0.5
         !   k = NINT((R_f -zeta_f -ant_toffset(i_ant))/del_t)
         !   i=k-3 ; if(i.lt.1) i=1
         !   write(2,"(2F9.2,2i4, 3F7.3, 20g13.3)") zeta_f, d_ant(i_ant), i_ant, k &
         !      , t_shift &
         !      , ((R_f -zeta_f) -ant_toffset(i_ant))/del_t - k &
         !      , PulseTime_x(i_ant)-(R_f -zeta_f) & !+ (NN_f-1)*zeta_f & ! should be close to zero and positive
         !      !Time_maxx(i_ant) -t_shift -(ant_toffset(i_ant)/del_t - Sample_Offset+1)-k, &
         !      !Time_maxx(i_ant) -t_shift -(ant_toffset(i_ant)-ant_toffset(1))/del_t-k, &
         !      !, Time_maxx(i_ant) -t_shift -k &
         !      , ant_xtrace(i:k+4,i_ant)
         !EndIf
         !Flush(Unit=2)
         DO i=i_nu_ini,i_nu_max
            Cnux_f(i)=Cnux_f(i) + Cnux_ant(i,i_ant) * Factr * ph_shift
            Cnuy_f(i)=Cnuy_f(i) + Cnuy_ant(i,i_ant) * Factr * ph_shift
            ph_shift=ph_shift*dph_shift
         ENDDO
         !EndDo ! i_phi=1,N_phi0  ! angle between ObsDist and rs_f
      EndDo ! i_ant  ! sum over all antenna distances
      ! At this point the beamformed pulse is determined for height zeta_f
      ! transform to time trace
      Call RFTransform_CF2RT(Cnux_f(0),ttracex_zf( -Sample_Offset,i_zf)) !transform to real-time trace, Should have t=0 point at isample_offset
      Call RFTransform_CF2RT(Cnuy_f(0),ttracey_zf( -Sample_Offset,i_zf)) !transform to real-time trace, Should have t=0 point at isample_offset
      !If(zeta_f.ge.6000. .and. zeta_f.le.6500) Then
      !   write(2,"(2g10.5,A,10g13.3)") zeta_f, Sample_Offset &
      !   , ', ttrace_zf(-3,+4):', ttracex_zf( -3:4,i_zf)
      !EndIf
   EndDo ! i_tr=0,N_tr ! Loop over focal height on the shower axis
   Call DAssignFFT()      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !DeAllocate( Cnux_f, Cnuy_f )
   !
   !write(2,*) 'h-dims(EData):', ih_min, ih_max, N_t
   !flush(unit=2)
   !stop 'DataBeaming'
   !------------------------------------------------
   !  Make data file for plotting
   !i_zf=25
   !write(2,*) '!Beamed E:',ih_min,ih_max, sqrt( SUM(ttracex_zf( ih_min:ih_max,i_zf)*ttracex_zf( ih_min:ih_max,i_zf)) )&
   !   , sqrt( SUM(ttracex_zf( -20:+20,i_zf)*ttracex_zf( -20:+20,i_zf)) )
   OPEN(UNIT=4,STATUS='unknown',FILE=trim(Base)//'x-BeamPwr-tt.dat')
   write(2,*) 'writing file "',trim(Base)//'x-BeamPwr-tt.dat','"'
   flush(unit=2)
   Zeta_step = N_zf * del_z / 10.
   !write(4,"(A,17I4)") '!    i_h,  ttracex_zf(i_h,0:N_zf) i_zf=      ', (NINT((k*1000.-z_low)/del_z), k=1,17)
   write(4,"(A,17G13.4)") '!    i_h,  ttracex_zf(i_h,zeta_f[km]=',(k*Zeta_step, k=3,10)
   Do  i=ih_min,ih_max
      write(4,"(22(1pG13.4))") (i-0.5)*SamplingTime, (ttracex_zf(i,NINT(k*Zeta_step/del_z)), k=1,10)  !, PSF_trace(i)
   EndDo ! i_zf=0,N_zf
   Close(Unit=4)
   PlotType='Test_kernel_1-3'
   PlotName=trim(Base)//'testKernel_x'
   PlotDataFile=trim(Base)//'x-'
   !Call MakeGLEplots(PlotType, PlotName, PlotDataFile, mode='ExC') !, Submit, CleanPlotFile
   !
   PlotType='W_PSF'
   PlotName=trim(Base)//'WPSF_x'
   PlotDataFile=trim(Base)//'x-'
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Submit, CleanPlotFile
   !
   OPEN(UNIT=4,STATUS='unknown',FILE=trim(Base)//'y-BeamPwr-tt.dat')
   write(2,*) 'writing file "',trim(Base)//'y-BeamPwr-tt.dat','"'
   flush(unit=2)
   write(4,"(A,17I4)") '!    i_h,  ttracer_zf(i_h,0:N_zf) i_zf=      ', (NINT((k*1000.-z_low)/del_z), k=1,17)
   write(4,"(A,17I4)") '!    i_h,  ttracer_zf(i_h,0:N_zf) zeta_f[km]=',(k, k=1,17)
   Do  i=ih_min,ih_max
      write(4,"(22(1pG13.4))") (i-0.5)*del_t, (ttracey_zf(i,NINT(k*Zeta_step/del_z)), k=1,10)  !, PSF_trace(i)
   EndDo ! i_zf=0,N_zf
   Close(Unit=4)
   PlotType='Test_kernel_1-3'
   PlotName=trim(Base)//'testKernel_y'
   PlotDataFile=trim(Base)//'y-'
   !Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Submit, CleanPlotFile
   !
   PlotType='W_PSF'
   PlotName=trim(Base)//'WPSF_y'
   PlotDataFile=trim(Base)//'y-'
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Submit, CleanPlotFile
   !
   !
   Return
End Subroutine DataBeaming
End Module DataBeam
