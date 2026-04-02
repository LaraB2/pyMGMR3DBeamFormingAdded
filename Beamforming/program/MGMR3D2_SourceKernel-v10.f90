!------------------------------
!------------------------------
!---------------------------------------------------------
      module KernelFilter
      contains
Subroutine FreqFilterKernel(K_in, NuDim_in, N_tr, SampleTime_old, padding, SampleTime_new, nu_min, nu_max, &
      K_nu, K_h, NuDim_new)
   use constants, only : dp, pi, c_l
   use FFT, only : FFTransform_su,DAssignFFT, RFTransform_CF, RFTransform_CF2RT
   implicit none
   Integer, intent(in) :: N_tr
   Integer, intent(in) :: NuDim_in ! half the length of input trace [samples]
   REAL(dp), intent(in) :: K_in(1-NuDim_in:NuDim_in,0:N_tr)
   real(dp), intent(in) :: SampleTime_old ! step size before FreqFilt [m]
   integer, intent(in) :: padding ! number of zero samples added at begin and end input trace
   real(dp), intent(in) :: SampleTime_new  ! sampling time trace after FreqFilt [ns]; typically 5 for LOFAR
   real(dp), intent(in) :: nu_min, nu_max ! block filter frequencies [MHz]
   complex(dp), allocatable, intent(out) :: K_nu(:,:)
   Real(dp), allocatable, intent(out) :: K_h(:,:)
   Integer, intent(out) :: NuDim_new ! half the length of output trace [samples]
   integer :: i_tr
   integer :: tDim_in
   real(dp), allocatable :: tTrace(:)  ! scratch array
   integer :: FF_dim, i_nu_ini, i_nu_max  !, tTrace_Offset_dwn
   real(dp) :: NuSample_new
   integer :: tdim_new ! length of input trace [samples]
   complex(dp), allocatable :: filt(:)
   complex(dp), allocatable :: Cnu(:)
    !
    tDim_in=NuDim_in*2
    FF_dim=(NuDim_in+padding)*2 ! dimension time trace before frequency filtering, padding 'padding' number of zeros at begin and end
    !
    if(SampleTime_new.lt.SampleTime_old) then ! new should be down sampled from input, not up
        write(2,*) 'SampleTime_new should not be smaller than SampleTime_old=',SampleTime_new ,SampleTime_old
        Stop 'FreqFilterKernel:SampleTime_new'
    endif
    tdim_new=FF_dim*SampleTime_old/SampleTime_new  ! cover the same length in time [m] for the two traces
    nudim_new=tdim_new/2  ; tdim_new=2*NuDim_new
    NuSample_new=c_l/(tdim_new*SampleTime_new)    ! in [GHz]
    write(2,"(A,I4,A,F6.3,A,F7.2,A,F7.2)") 'padding=',padding, ', sampling time=', SampleTime_old, '[m], input trace length=', &
      tDim_in*SampleTime_old, '[m] after padding:', FF_dim*SampleTime_old
    i_nu_ini= nu_min/(1000.*NuSample_new)  ; i_nu_max=nu_max/(1000.*NuSample_new)+0.5
    if(i_nu_ini .lt. 0) i_nu_ini=0
    if(i_nu_max .ge. NuDim_new) i_nu_max=NuDim_new
    write(2,"(A,F5.3,A,F7.2,A,F7.2,A,F5.1,A,F7.1,A,F5.1,A)") 'Down-sampling to samples of ',SampleTime_new/c_l,&
      '[ns] and trace length of ',tdim_new*SampleTime_new,'[m]=',tdim_new*SampleTime_new/c_l, &
      '[ns], frequency filter from ',i_nu_ini*NuSample_new*1000.,' till ',i_nu_max*NuSample_new*1000., &
      '[MHz], steps of ',NuSample_new*1000.,'[MHz]'
   !
   allocate(filt(0:nudim_new))
   allocate(tTrace(1:FF_dim))
   allocate(Cnu(0:NuDim_in+padding))
   allocate(K_nu(0:nudim_new,0:N_tr))
   allocate(K_h(1-nudim_new:nudim_new,0:N_tr))
   write(2,*) '!nudim_new:', nudim_new,nudim_in, ',  i_nu_ini:', i_nu_ini,i_nu_max
   !
   call FFTransform_su(FF_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   filt(:)=0.      ! Frequency filter
   filt(i_nu_ini: i_nu_max)=1.
   !
!   write(2,*) 'FreqFilterKernel:', FF_dim, N_tr, SampleTime_old, ' nu_max_old=',0.5*c_l/SampleTime_old
!   write(2,*) 'FreqFilterKernel:', tdim_new*SampleTime_new, FF_dim*SampleTime_old
!   write(2,*) (Real(NuSample_new*i_tr),i_tr =0,nudim_new)
   Do i_tr =0,N_tr
      tTrace(1:padding)=0.0d0
      tTrace(padding+1:padding+tDim_in)=K_in(1-NuDim_in:NuDim_in,i_tr)
      tTrace(tDim_in+padding+1:FF_dim)=0.0d0
      Call RFTransform_CF(tTrace,Cnu(0))  ! transform to frequency
      K_nu(0:nudim_new, i_tr)=filt(0:nudim_new)*Cnu(0:nudim_new) ! apply filter and reduce nu-trace length
!      If(i_tr.lt.5) Then
!         write(2,*) real(abs(Cnu(0:nudim_new)))
!      EndIf
   EndDo !  i_tr =1,ObsDist_dim
   Call DAssignFFT()      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
!   write(2,*) 'FreqFilterKernel:', tdim_new, NuSample_new, SampleTime_new, ' nu_max_new=',0.5*c_l/SampleTime_new
   Call FFTransform_su(tdim_new)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   Do i_tr =0,N_tr
      Call RFTransform_CF2RT(K_nu(0, i_tr),K_h(1-nudim_new, i_tr)) !transform to real-time trace
   EndDo !  i_tr =1,ObsDist_dim
   Call DAssignFFT()      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
   Deallocate(tTrace,filt,Cnu)
!   Flush(Unit=2)
!   stop
   Return
   End subroutine FreqFilterKernel
end module KernelFilter
!-------------------------------------------------------------------------
!------------------------------
Subroutine WideSource_Kernel(zeta_f, rs_0, Z_min, Z_max, N_phi0, phi_0, DistMax, label)
!   Calculate the source @ height=\zeta as function of tmh=(t-h) for a fixed observerdistance to this shower
!   A backtrace integration is performed over observer position
!   CD_i (CoreDist_A(CD_i)) determines the position of the pencil shower w.r.t. the core of the real shower which determines the pancake thickness.
!   zeta= height in the atmosphere of the back-trace voxel
!   Shower time varies with h, distance behind the shower front.
!
    use RFootPars, only : nu_min, nu_max
    use Atmosphere, only : xi,dxi, AtmHei_dim, AtmHei_step
   use Pancakefunction, only : PancakeThi, CoreDist_Dim, alpha_tr
    ! Ex_to(i=1,tTrace_dim_o;idi=1,ObsDist_dim) Field as function of time for different observer positions
    !    where t_to(i)=i*tTrace_step   !  starts after shouwer touching ground.
   !
   use BigArrays, only : Line2Core
   use BigArrays, only : ObsDist_dim, ObsDist_Step ! , tTrace_dim_b, tTrace_dim_o, tTrace_step
    use BigArrays, only : tTrace_step !, tTrace_dim_b, tTrace_dim_o, tTrace_Offset, nuTrace_step, nuTrace_dim,
    use BigArrays, only : Ix,Iy,IQ
   use eventdata, only : PlotBase
   use constants, only : dp, pi, c_l
   use KernelFilter, only : FreqFilterKernel
   Use Systemcommands, only : MakeGLEplots
   Use Pancakefunction, only : GetLambda
   Use Pancakefunction, only : fiehLam
   use LateralDistribution, only : W_tc
   implicit none
!    real*8 :: Ex(0:CoreDist_Dim),AxD(0:CoreDist_Dim),Ey(0:CoreDist_Dim),AyD(0:CoreDist_Dim),Ar(0:CoreDist_Dim),nt_r,t_o
   real(dp), intent(in) :: zeta_f, rs_0, Z_min, Z_max  ! zeta_f is the distance from ground-core, rs_0 is distance from shower axis of focus point
   real(dp), intent(inout) :: phi_0
   Integer, intent(in) :: N_phi0
   !real(dp), parameter :: DistMax=500.d0 ! DistMax=500.d0 ! DistMax=700.d0 !! maximum for antenna to shower-core distance
   real(dp), intent(in) :: DistMax !! maximum distance of antenna to shower-core
	Character*10, intent(in) :: label
!PlotBase="plot00-up-L040" ! h_range=50.d0; N_h=200; N_phir=5; N_ca=50
!PlotBase="plot00-up-L040R" ! h_range=50.d0; N_h=50; N_phir=15
!PlotBase="plot00-up-L040ER" ! h_range=50.d0; N_h=200; N_phir=15
!PlotBase="plot00-up-L040W" ! h_range=100.d0; N_h=400
   real(dp), parameter :: h_range=20.d0 ! h_range=30.d0 ! h_range=50.d0 ! h_range=100.d0 ! h_range=17.d0 ! h_range=200.d0 ! h_range=100.d0 ! h_range=20.d0 ! h_range=2.5d0  ! h_range=1.d0  !  [m]
   Integer, parameter :: N_h=100 ! N_h=100 !N_h=120 !  N_h=200 ! N_h=400 ! N_h=50 ! N_h=20 ! N_h=100 ! N_h=26  ! N_h=80 ! number of steps in distance behind front
   Integer, parameter :: dN_tr=5 ! dN_tr=10 ! dN_tr=1 ! dN_tr=3 !  Step size in t_r in units of AtmHei_step
   ! fine N_phi grid important for power at lower altitudes
   Integer, parameter :: N_phir=5 ! N_phir=15 ! N_phir=3 ! N_phir=7 ! N_phir=17 !  Number of angles for integration over shower area
   ! fine N_ca grid important for quenching fluctuations
   Integer, parameter ::  N_ca= 25 ! N_ca= 50 ! N_ca= 100 ! N_ca= 200 ! N_ca= 250 ! N_ca= 50 !messy) ! grid for antenna to shower-core distance
   ! finer than d_crMax=20.d0 grid not important for dmax=250
   real(dp), parameter :: d_crMax=20.d0  ! d_crMax=10.d0  ! [m]  ! max increment in core-ray grid
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
   ! Default: d_h, N_ca, d_crMax:   0.25       100   20.00
   real(dp) :: h_0,R_0, NN_0,dr_0, d_phi0, rsa_0 ! all _0 mark the quantities at the focus point of retraced rays
   real(dp) :: rc_a, t_a, d_ca ! all _0 mark the quantities at the focus point of retraced rays
   real(dp) ::  zeta_s, h_s, NN_s, drefrac_s, tr_s
   Real(dp) :: rc_r, dc_r, xc_r, rr_a, dr_a, xr_a, phi_r, d_pr, c_pr, lim_ca, R_r, h_r, CD_d
   Integer :: i_pr, ic_r, NC_r, ir_a
   Integer, parameter :: CoreRaySteps_dim=500
   real(dp) :: CoreRaySteps(0:CoreRaySteps_dim)
   real(dp) :: di, weight, norm, nrm, wei, Factr, d_h
   real(dp) :: Kx_ca, Kx, dKx_h, Sx_h, Rx_h
   real(dp) :: Kr_ca, Kr, dKr_h, Sr_h, Rr_h
   Real(dp), allocatable :: WeiFie_x(:,:)
   Real(dp), allocatable :: WeiFie_r(:,:)
   Complex(dp), allocatable :: K_nu(:,:)
   Real(dp), allocatable :: K_h(:,:)
   Real(dp) :: MaxWF, R_z
   integer :: i, i_ca, i_t, i_h, i_tr, i_phi, CD_i, it_min, it_max, N_tr, i_atm
   real(dp) :: lambda, DerLambda, Alpha, DerAlpha, dwr, rcr_max, fh, dfh, dfhdr, D_s, X_tst
   !integer, save :: il_start = 1, CD_i_prev=-1
   Integer :: trhMax(1:2),trMax(1-N_h-52:N_h+52) !, trhMin(1:2),trMin(1-N_h-50:N_h+50)
   Integer :: padding, NuDim_new, hcnt_l, hcnt_u  ! , base_h, hCount
   Real(dp) :: SampleTime_new ! , nu_min, nu_max
   character(len=10) :: txt
   character(len=45) :: FMTSK
   Integer :: i_zf, i_z0, N_ca_l, N_ca_u
   Real(dp) :: K_IncPwr(0:AtmHei_dim), K_trace(1-N_h-52:N_h+52), K_PSF(0:AtmHei_dim)
   !Real(dp) :: K_IncPwrA(0:AtmHei_dim),  K_IncPwrB(0:AtmHei_dim), K_IncPwrC(0:AtmHei_dim), K_IncPwrD(0:AtmHei_dim)
   Real(dp) :: K_CohPwr(-AtmHei_dim/2:AtmHei_dim/2), dK_CohPwr(-AtmHei_dim/2:AtmHei_dim/2)
   real(dp)  delKX_h
   Character(len=50) :: PlotType, PlotName, PlotDataFile
   !
   !write(2,*) 'SourceKernel, zeta_f, rs_f:', zeta_f, rs_0, ', N_phif, phi_f:', N_phi0, phi_0*180/pi
   write(2,"('input parameters: zeta_f=',F7.0,' rs_0=',F7.1,' Z_min=',F7.0,' Z_max=',F7.0,' phi_0=',F6.0, &
      ' N_phi0=',I3,' label=',A)") &
      zeta_f, rs_0, Z_min, Z_max, phi_0*180/pi, N_phi0, label
   write(2,"('hard set parameters: N_h=',I4.0,' h_range=',F5.1,' dN_tr=',I3,' N_phir=',I3,' N_ca=',I4,' DistMax=',F7.0)") &
       N_h, h_range, dN_tr, N_phir, N_ca, DistMax
   Flush(unit=2)
   !
   d_h=h_range/N_h
   d_phi0=pi/N_phi0
   If((ObsDist_dim-1)*ObsDist_Step .lt. DistMax) Then
      Write(2,*) '*** error in WideSource_Kernel; DistMax toolarge:', DistMax, (ObsDist_dim-1)*ObsDist_Step
      stop 'error in WideSource_Kernel'
   EndIf
   d_ca=DistMax/N_ca
   !
   If(z_max.le. zeta_f) stop 'WideSource_Kernel: too low z_max'
   it_min=z_min/(dN_tr*AtmHei_step)
   If(it_min.lt.1) it_min=1
   it_max=z_max/(dN_tr*AtmHei_step)
   If(it_max.ge. AtmHei_dim/dN_tr) it_max=AtmHei_dim/dN_tr - 1
   N_tr=(it_max-it_min)
   !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   If(N_tr.lt. 51) N_tr=51  ! otherwise plotting the distributions give problems
   !allocate( WeiFie(1-N_h:N_h,it_min:it_max) )
   allocate( WeiFie_x(1-N_h:N_h,0:N_tr) )
   allocate( WeiFie_r(1-N_h:N_h,0:N_tr) )
   !
   NC_r=CoreDist_Dim*2+1  ! two less than in the definition of Line2Core
   rcr_max= Line2Core(NC_r)
   If(CoreRaySteps_dim.le.NC_r) Then
      write(2,*) 'Error ***** CoreRaySteps_dim too small:',CoreRaySteps_dim,NC_r
      Stop 'Error-CoreRaySteps_dim'
   EndIf
   CoreRaySteps(0:2)=Line2Core(0:2)
   Do ic_r=3, NC_r+1
      CoreRaySteps(ic_r)=Line2Core(ic_r)
      dc_r=(Line2Core(ic_r+1)-Line2Core(ic_r-1))/2.
      If(dc_r .gt. d_crMax) exit    ! step to (ic_r+1) too large
   EndDo
   !write(2,*) 'ic_r:', ic_r, Line2Core(ic_r), dc_r, d_crMax
   If(ic_r.ge.NC_r+1) then  ! looped till the end
      CoreRaySteps(NC_r+1)=Line2Core(NC_r+1)
   Else     ! refine final step sizes
      NC_r=ic_r + (Line2Core(NC_r)-Line2Core(ic_r))/d_crMax
      If(NC_r .gt. CoreRaySteps_dim) NC_r=CoreRaySteps_dim
      Do i=ic_r+1,NC_r+1
         CoreRaySteps(i)=CoreRaySteps(i-1) + d_crMax
      EndDo
   EndIf
   !write(2,*) 'NC_r:', NC_r, CoreRaySteps(NC_r)
   !
   nrm=1.e6
   norm=0.
   !N_ca_l=15 ; N_ca_u=15  ! nominally from 1 to N_ca
   N_ca_l=1 ; N_ca_u=N_ca  ! nominally from 1 to N_ca
   write(2,*) 'Consider antennas in band width', (N_ca_l-0.5)*d_ca, ' to', (N_ca_u+0.5)*d_ca
   Do i_ca=N_ca_l, N_ca_u
      rc_a=i_ca*d_ca !
      R_0=sqrt(rc_a*rc_a +zeta_f*zeta_f)
      weight=2.*rc_a*d_ca
      norm=norm+weight/R_0  !   /(R_0*R_0)
   EndDo ! idi=1, idi_N
   !norm=1.
   !norm=norm*nrm
   !normc=normc*nrm
   WeiFie_x(:,:) = 0.
   WeiFie_r(:,:) = 0.
   R_z = 0.
   !
   !N_tr=2
   write(2,"(A,I4,I3, F6.3,3I4,F6.2)") 'N_tr, dN_tr, d_h, N_ca_l,_u, N_phir, d_crMax:', &
         N_tr, dN_tr, d_h, N_ca_l, N_ca_u, N_phir, d_crMax
   !write(2,*) '!(N_tr+it_min)*dN_tr=', (N_tr+it_min)*dN_tr, N_tr,it_min,dN_tr, it_max
   !flush(unit=2)
   !
   write(2,*) '!Alpha_tr before zeroing:'
   write(2,"('!',50g13.6)") alpha_tr(:)
   alpha_tr(:)=1.
   write(2,*) '!!!!!!!!!!!!! Hack: alpha_tr(:)=1**********'
   !
   Do i_h=1-N_h,N_h
   !Do i_h=10, 10
      h_0=(i_h-0.5)*d_h
      !If(h_0.lt.40.5 .or. h_0.gt.41.) cycle
      !tr_0=h_0-zeta_f ! should be negative, all quantities at the source are given
      !write(2,*) 'Source_zeta; z,h=',zeta_f,h
      !Flush(unit=2)
      i=(zeta_f/AtmHei_step)
      if(i.lt.1) return
      di=zeta_f/AtmHei_step-i
      NN_0=1.+ di*xi(i+1) + (1.-di)*xi(i)
      dr_0=   di*dxi(i+1) + (1.-di)*dxi(i)
      !
      !write(2,*) 'dto:', 'dto, ObsDist, Rx/norm, Ex/(NN*R), Ex, cald '
      Sx_h=0.
      Kx=0.
      Sr_h=0.
      Kr=0.
      !it_min=500 ; it_max=610
      !Do i_tr=1,1 ! Loop over source time for obtaining weighting factors for shower current
      Do i_tr=0,N_tr ! Loop over source time for obtaining weighting factors for shower current
         i_atm=(i_tr+it_min)*dN_tr
         tr_s=-i_atm*AtmHei_step
         !
         !JQ= IQ(i_atm)
         !Jx= Ix(i_atm)
         !Jy= Iy(i_atm)
         Alpha=alpha_tr(i_atm) ! the parameter that determines the width of the pancake function
         DerAlpha=(alpha_tr(i_atm+1)-alpha_tr(i_atm-1))/(2*AtmHei_step) ! its derivative
         !write(2,*) 'Alpha,DerAlpha:', Alpha,DerAlpha
         !Write(2,*) 'h_0,tr_s:', h_0, tr_s, '============================='
         !Flush(unit=2)
         !
         Do i_ca=N_ca_l, N_ca_u ! integrate over antenna position and time, keeping the focus point given by zeta_f, rs_0 & phi_0
         !Do i_ca=40,40 ! 50
            Kx_ca =0
            Kr_ca =0
            !write(2,*) 'i_h, i_atm, idi',i_h, i_atm, idi
            !     External geometry, not related to shower structure
            ! zeta_f = focus height (input)
            ! rs_0   = distance focus point to shower axis, where we will consider the focal-ray (input)
            ! h_0    = distance behind shower front of focus point = subtracted retarded time at focus point (loop)
            ! rc_a   = distance antenna from shower axis (integrated over)
            ! phi_0  = angle between rc_a and rs_0 (may be integrated)
            !  Derived quantities:
            ! rsa_0  = distance antenna to ray at rs_0
            ! R_0    = distance to focus point on focal-ray
            ! t_a    = arrival time of signal from focus point on focal-ray
            !
            rc_a=i_ca*d_ca ! Distance of antenna to shower axis
            !write(2,*) 'rc_a:',i_ca, rc_a
            weight=2.*rc_a*d_ca*d_phi0   !*(DistMax-ObsDist)
            Do i_phi=1,N_phi0  ! angle between ObsDist and rs_0
               If(N_phi0.gt.1) phi_0=(i_phi-0.5)*d_phi0
               rsa_0=sqrt(rc_a*rc_a + rs_0*rs_0 - 2*rs_0*rc_a*cos(phi_0)) ! distance to shower
               R_0=sqrt(zeta_f*zeta_f+ rc_a*rc_a + rs_0*rs_0 - 2*rs_0*rc_a*cos(phi_0)) ! Distance between antenna and source in shower
               t_a=NN_0*R_0 + h_0-zeta_f ! time in antenna is derived quantity
     ! If(i_h.eq.0 .and. i_tr.eq.10) write(2,*) 't_shift:', NN_0*R_0-zeta_f, zeta_f, rc_a
               Factr=1.d0/(NN_0*R_0)
               !Factr=Factr*2.*(DistMax-rc_a)/DistMax   ! linear decrease towards rim
               If(t_a.lt.0.) cycle
               !
               !     Internal geometry, related to the extent of the source that may give signal in the antenna
               ! tr_s   = retarded time of the shower source, outermost loop since this determines the current
               ! rc_r   = distance ray to core ; ic_r, dc_r, xc_r (pojection along rc_a direction)
               ! rr_a   = distance antenna to shower ray ; ir_a, dr_a, xr_a
               ! phi_r  = angle between rc_a and rc_r or rc_a and rr_a depending on ray integration ; i_pr, d_pr
               !  Derived quantities:
               ! h_s    = Distance behind shower front of source
               ! R_r    = distance from antenna to position on ray at height zeta_s
               ! lim_ca = rc_a/2. half way point (determines what to take as integration variable for ray position)
               lim_ca=rc_a/2.
               d_pr=pi/N_phir
               If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.10 .and. NINT(rc_a/10).eq.8) &
                  write(2,"(30g13.5)") '!Kx_c', zeta_f,i_h,i_tr, i_pr, rc_a &
                  , NN_s, drefrac_s, t_a, tr_s, CoreRaySteps(0:5)
               Do i_pr=1,N_phir       ! angle of point off shower axis with antenna-axis direction
                  phi_r=(i_pr-0.5)*d_pr
                  c_pr=cos(phi_r)
                  dKx_h=0.
                  dKr_h=0.
                  !write(2,"(T26,'dKx_h',T38,' dfh',T50,' wei',T62,'  h_r',T74,' rc_r',T86,' rr_a',2I4)") i_pr,i_ca           !
                  Do ic_r=1,Nc_r     ! distance ray to shower axis, use predefined grid
                     rc_r=CoreRaySteps(ic_r)
                     dc_r=(CoreRaySteps(ic_r+1)-CoreRaySteps(ic_r-1))/2.
                     xc_r=rc_r*c_pr
                     !write(2,*) '!,ic_r:', xc_r, lim_ca, CoreRaySteps(ic_r+1)*c_pr, rr_a
                     If(xc_r.gt. lim_ca) exit ! getting closer to antenna than core
                     If(CoreRaySteps(ic_r+1)*c_pr .gt. lim_ca) then ! next step would move out of range, take care of proper area (approx)
                        dc_r=lim_ca/c_pr-(rc_r-dc_r/2.)  ! remaining distance along the center line betwee observer and core
                     endif
                     rr_a=sqrt(rc_a*rc_a + rc_r*rc_r - 2*xc_r*rc_a)
                     Call Findhs(t_a, rr_a, tr_s, h_r, NN_s, drefrac_s)
!  "Findhs" solves   NN_s^2 [zeta_s^2 + r_a^2] = (t_a-tr_s)^2  with  zeta_s = h_r-tr_s
                     If(h_r.le.0.) cycle
                     If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.10 .and. NINT(rc_a/10).eq.8) &
                        write(2,"(30g13.5)") '!Findhs:', t_a, rr_a, tr_s, h_r, NN_s, drefrac_s
                     D_s=NN_s* NN_s*(h_r-tr_s) + drefrac_s*(t_a-tr_s)*(t_a-tr_s)/NN_s
                     Call GetLambda(rc_r, Lambda, DerLambda, CD_i, CD_d)
                     If(h_r/(Lambda*Alpha) .gt. 20. ) cycle
                     If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.10 .and. NINT(rc_a/10).eq.8) &
                        write(2,"(30g13.5)") '!GetLambda:', rc_r, Lambda, DerLambda, CD_i, CD_d
                     Call fiehLam(h_r,Lambda, DerLambda, Alpha,DerAlpha,fh,dfh,dfhdr)
                     If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.10 .and. NINT(rc_a/10).eq.8) &
                        write(2,"(30g13.5)") '!fiehLam:',  h_r,Lambda, DerLambda, Alpha,DerAlpha,fh,dfh
                     wei=W_tc(rc_r,dwr)
                     If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.10 .and. NINT(rc_a/10).eq.8) &
                        write(2,"(30g13.5)") '!W_tc:',  wei, rc_r,dwr
! x-polarisation:
!        Ex =Ex  + weigh * (dfh*Jx -fh *dJx) ! still needs to be times w(r), as done in following line
!         Ex_to(:,idi)=Ex_to(:,idi) + gth* ( wei *DEx(:) - DAxD(:)*weid ) !  (from LateralInt)
                     !dKx_h = dKx_h+ weight*Factr*dfh*wei*dc_r*d_pr/D_s
                     dKx_h = dKx_h + weight*Factr*dfh*wei * d_pr*dc_r/D_s
! radial polarisation:
!         Ar =Ar  -(fh *JQ) * weigh    ! still needs to be times dw/dr
!         Erh =Erh  -(dfhdr *JQ) * weigh    ! still needs to be times w(r)
!         Er_to(:,idi)=Er_to(:,idi) + gth*(DAr(:)*weid + DErh(:)*wei)*c_pr !  (from LateralInt)
                     dKr_h = dKr_h - weight*Factr*(dfhdr*wei + fh*dwr)*c_pr *  d_pr*dc_r/D_s
               !If(rc_r.lt.1. .and. rc_a.gt.195. .and. rc_a.lt.205.) &
               !   Write(2,"(20g13.3)") h_0, rc_a, rr_a, rc_r, tr_s,'a-Ex:', dfh,wei, ', Er:', fh,dwr, dfhdr, c_pr &
               !      ,' ;',h_r, lambda
                     If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.10 .and. NINT(rc_a/10).eq.8) &
                        write(2,"(30g13.5)") '!Kx_', zeta_f,i_h,i_tr, i_pr, rc_a &
                        , dKx_h/norm,  ic_r  ,dfh, wei, h_r, rc_r,rr_a, h_r/Lambda !h_r-tr_s, rc_a!,'fh', fh,dfh
                     !   Write(2,"(A, 1pg12.4, A,9(1pg12.4))") '! delKX_h=weight*Factr*dfh*wei*dc_r/D_s', delKX_h, ';' &
                     !      ,weight, Factr, dfh, wei, dc_r, d_pr, 1./D_s,  lim_ca, c_pr, CoreRaySteps(ic_r+1)
                        !Write(2,*) '! h_r,Lambda, DerLambda ;  fh,dfh,dfhdr', &
                        !   h_r,Lambda, DerLambda, ';',fh,dfh,dfhdr
                     !EndIf
                  EndDo ! ic_r=0,CoreDist_Dim-1CD_i
                  !write(2,"(A, 3I5, 1pg12.4, A,9(1pg12.4))") '!dKx_h1,ic_r:',i_ca,ic_r,i_pr, dKx_h, ';' &
                  !    ,dfh, wei, h_r, rc_r,rr_a, h_r/Lambda, Kx_ca, dKx_h !h_r-tr_s, rc_a!,'fh', fh,dfh
                  !    Flush(unit=2)
                  Kx_ca = Kx_ca + 2*dKx_h
                  Kr_ca = Kr_ca + 2*dKr_h
                  If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.10 .and. NINT(rc_a/10).eq.8) &
                     write(2,"(30g13.5)") '!Kx_ca', zeta_f,i_h,i_tr, i_pr, rc_a &
                     , Kx_ca/norm, dKx_h/norm, wei, rc_r, Lambda,rc_r, D_s
                  !write(2,"(20g13.3)") '!aKx_h,ic_r:',i_atm, i_ca, i_pr, i_h,  dKx_h, dKr_h, tr_s
                  !
                  dr_a=ObsDist_Step ! to be able to take care of rs=0
                  dKx_h=0.
                  dKr_h=0.
                  Do ir_a=1, ObsDist_dim      ! distance ray to antenna, use regular grid
                     rr_a=(ir_a-0.5)*ObsDist_Step
                     xr_a=rr_a*c_pr
                     If(xr_a.gt. lim_ca) exit ! getting closer to core than antenna
                     If((ir_a+0.5)*ObsDist_Step*c_pr .gt. lim_ca) then ! next step would move out of range, take care of proper area (approx)
                        dr_a=lim_ca/c_pr-(ir_a-1)*ObsDist_Step  ! remaining distance along the center line betwee observer and core
                     endif
                     rc_r=sqrt(rc_a*rc_a + rr_a*rr_a - 2*xr_a*rc_a)  !
                     If(rc_r .gt. rcr_max) cycle
                     Call Findhs(t_a, rr_a, tr_s, h_r, NN_s, drefrac_s) !  "Findhs" solves   NN_s^2 [zeta_s^2 + r_a^2] = (t_a-tr_s)^2  with  zeta_s = h_r-tr_s
                     If(h_r.le.0.) cycle
                     D_s=NN_s* NN_s*(h_r-tr_s) + drefrac_s*(t_a-tr_s)*(t_a-tr_s)/NN_s
                     Call GetLambda(rc_r, Lambda, DerLambda, CD_i, CD_d)
                     If(h_r/(Lambda*Alpha) .gt. 20. ) cycle
                     Call fiehLam(h_r,Lambda, DerLambda, Alpha,DerAlpha,fh,dfh,dfhdr)
                     wei=W_tc(rc_r,dwr)*rr_a/rc_r
                     dwr=dwr*rr_a/rc_r
!        wei=dr_a*d_pr*W_tc(rc_r,dwr)*rr_a/rc_r  ! !  (from LateralInt)
!        weid=dr_a*d_pr*dwr*rr_a/rc_r          ! !  (from LateralInt)
! x-polarisation:
!        Ex =Ex  + weigh * (dfh*Jx -fh *dJx) ! still needs to be times w(r)
!        Ex_to(:,idi)=Ex_to(:,idi) + gth*wei *DEx(:) - gth*DAxD(:)*weid  !  (from LateralInt)
                     dKx_h = dKx_h+ weight*Factr*dfh*wei * dr_a*d_pr/D_s
!
! radial polarisation:
!         Ar =Ar  -(fh *JQ) * weigh    ! still needs to be times dw/dr
!         Erh =Erh  -(dfhdr *JQ) * weigh    ! still needs to be times w(r)
!        Er_to(:,idi)=Er_to(:,idi) + gth*(DAr(:)*weid + DErh(:)*wei)*((rc_a-xr_a)/rc_r)   !  (from LateralInt)
                     dKr_h = dKr_h - weight*Factr*(dfhdr*wei + fh*dwr)*((rc_a-xr_a)/rc_r) * dr_a*d_pr/D_s
               !If(i_atm.eq.60 .and. i_ca.eq.20 ) &
                  If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.10 .and. NINT(rc_a/10).eq.8) &
                     write(2,"(30g13.5)") '!Kx_', zeta_f,i_h,i_tr, i_pr, rc_a &
                     , dKx_h/norm,  ic_r  ,dfh, wei, h_r, rc_r,rr_a, h_r/Lambda !h_r-tr_s, rc_a!,'fh', fh,dfh
                  EndDo ! ir_a=1, ObsDist_dim      ! distance ray to antenna, use regular grid
                  Kx_ca = Kx_ca + 2*dKx_h
                  Kr_ca = Kr_ca + 2*dKr_h
                  !If(abs(10.*dKx_h*Ix(i_atm)/normc).gt. 1.0)
                  !if(i_ca.eq.30) write(2,*) '!dKx_h,ir_a:',i_atm, i_pr, dKx_h, tr_s, h_r
                  !If(abs(dKx_h).gt.1.e-6 .or. abs(dKr_h).gt.1.e-6 ) &
                  !write(2,"(20g13.3)") '!bKx_h,ic_r:',i_atm, i_ca, i_pr, i_h,  dKx_h, dKr_h, tr_s
                  If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.10 .and. NINT(rc_a/10).eq.8) &
                     write(2,"(30g13.5)") '!Kx_cb', zeta_f,i_h,i_tr, i_pr, rc_a &
                     , Kx_ca/norm, dKx_h/norm, wei, rc_r, Lambda, rc_r, D_s
               EndDo ! i_pr=1,      ! angle with antenna direction
            EndDo ! i_phi=1,N_phi0
            WeiFie_x(i_h,i_tr) = WeiFie_x(i_h,i_tr) + Kx_ca
            WeiFie_r(i_h,i_tr) = WeiFie_r(i_h,i_tr) + Kr_ca
            !write(2,*) 'cumu-WeiFie(i_h,i_tr)',i_tr,i_h,WeiFie(i_h,i_tr), Kx_ca, i_ca
         EndDo ! i_ca=1,N_ca
         Sx_h=Sx_h +WeiFie_x(i_h,i_tr)*AtmHei_step*Ix(i_atm)*dN_tr/(norm*nrm)
         Sr_h=Sr_h +WeiFie_r(i_h,i_tr)*AtmHei_step*IQ(i_atm)*dN_tr/(norm*nrm)
         !write(2,*) 'nrm-WeiFie(i_h,i_tr)',i_tr,i_h,WeiFie(i_h,i_tr)/normc
         !If( (abs(WeiFie(i_h,i_tr))/norm .gt. 4.e-6) ) Then
         !   write(2,*) 'WideSource_Kernel,i_h, i_tr:',h_0, i_tr, WeiFie(i_h,i_tr)/norm
         !   flush(unit=2)
         !EndIf
      EndDo ! i_tr=0,N_tr  !
      !write(2,"(A,F8.3,3(1pG13.4))") 'Unfiltered Kernel*I; @h=',h_0, Sx_h ! , dfh*Jx*nrm, fh *dJx*nrm
      !write(2,*) '!Rx_h, Sx_h:', rc_a, i_phi,t_a,';', Rx_h/normc, Sx_h/normc,', normc=',normc
      !flush(unit=2)
      R_z = R_z + Sx_h*Sx_h*d_h
      !
      !If(i_h .eq. -N_h/2) stop 'Wource_Kernel-i_h'
   EndDo ! i_h=1-N_h,N_h
   write(2,"(A,2(1pG13.4),A,2(1pG13.4))") 'sq-integrated gives R^2, R:', R_z, sqrt(R_z), ', dt=',d_h
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !stop 'Wource_Kernel-end'
   !Return !================================================================
   !11111111111111111111111111111!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   If(norm.gt.0.) Then
      WeiFie_x(:,:) = WeiFie_x(:,:)/norm
      WeiFie_r(:,:) = WeiFie_r(:,:)/norm
   EndIf
   !write(2,*) 'UnFil@i_h=8, i_tr=15-25', WeiFie(8,15:25)
   !
   !-------------------------------------------------------------
   write(txt,"(I3.3,i3.3,I2.2)") NINT(zeta_f/100.),NINT(rs_0),NINT(phi_0*10/pi)
   !
   trhMax(:)=MaxLoc(ABS(WeiFie_x(:,0:N_tr)))
   trhMax(2)=trhMax(2)-1
   trhMax(1)=trhMax(1)-N_h
   !trhMin(:)=MinLoc(WeiFie(:,it_min:it_max))
   !trhMin(2)=trhMin(2)+it_min-1
   !trhMin(1)=trhMin(1)-N_h
   write(2,*) 'trhMax(1:2)',trhMax(1:2) ,', N_h=', N_h ! , ', trhMin(1:2)',trhMin(1:2)
   Flush(unit=2)
   If(trhMax(1).lt.-N_h+5) trhMax(1)=-N_h+5
   If(trhMax(1).gt.N_h-4) trhMax(1)=N_h-4
   !
   i_h=2000./(dN_tr*AtmHei_step)
   i_t=zeta_f/(dN_tr*AtmHei_step)-it_min  !  count corresponding to t_r=zeta_f
   i_h=min(i_h,N_h,i_t-1,N_tr-i_t-1,30)
   write(2,*) 'i_h,N_h,i_t-1,N_tr-i_t-1:',i_h,N_h,i_t-1,N_tr-i_t-1
   flush(unit=2)
   trMax(-i_h:i_h)=MAXLOC(ABS(WeiFie_x(1-N_h:N_h,i_t-i_h:i_t+i_h)),1)-N_h  ! determined the h where a max occurs for fixed t_r with (i_tr=i_t-i_h:i_t+i_h)
   write(2,"(A,61F10.1)") '  t_r [m]:',((i+it_min)*dN_tr*AtmHei_step,i=i_t-i_h,i_t+i_h)
   write(2,"(A,61F10.3)") 'h_max [m]:',((trMax(i)-0.5)*d_h,i=-i_h,+i_h)
   write(2,"(A,61g10.2)") 'Value@max:',(WeiFie_x(i,i_t+i),i=-i_h,+i_h)
   !
   trMax(1-N_h:N_h)=MAXLOC(ABS(WeiFie_x(1-N_h:N_h,0:N_tr)),2)-1  ! determined the t_r where a max occurs for fixed h
   !trMin(1-N_h:N_h)=MINLOC(WeiFie_x(1-N_h:N_h,it_min:it_max),2)+it_min-1
   write(2,*) 'trMax(:)',trMax(1-N_h:N_h)
   MaxWF=MAXVAL(ABS(WeiFie_x))
   !MinWF=MINVAL(WeiFie_x)
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)//'_'//trim(txt)//'.grd')
   write(2,*) 'writing file "',TRIM(PlotBase)//'x-SrcKrnl-'//trim(txt)//'.grd','"'
   write(4,"(2(1pG13.4),6(1pG13.4), (1pG13.4),0pF6.0,2I5,A)") zeta_f, rs_0, 1.-N_h, 1.*N_h, d_h, &
         it_min, it_min+N_tr, &
         dN_tr*AtmHei_step, MaxWF, phi_0*180/pi, trhMax(1), trhMax(2), ', ! WideSource_Kernel'
   !freadln inchan z0 rs0 ih1 ih2 delh it1 it2 dt maxK phi_f h_max itmax
   write(4,"('! Generating parameters, zeta_f=',F7.0,' rs_0=',F7.1,' Z_min=',F7.0,' Z_max=',F7.0,' phi_0=',F6.0, &
      ' N_phi0=',I3,' label=',A,' N_h=',I4.0,' h_range=',F5.1,' N_phir=',I3,' N_ca=',I4,' dN_tr=',I3,' DistMax=',F7.0)") &
      zeta_f, rs_0, Z_min, Z_max, phi_0, N_phi0, label, N_h, h_range, N_phir, N_ca, dN_tr, DistMax
   If((trhMax(2)+25) .gt. N_tr) trhMax(2)=N_tr-25
   If((trhMax(2)-25) .lt. 0) trhMax(2)=25
   write(4,"(A, 21(F6.0,7x))") '!    h,  t_r(max)  W@t_r=',(((trhMax(2)+i*5)+it_min)*dN_tr*AtmHei_step,i=-5,+5)
   Do  i_h=1-N_h,N_h
      write(4,"(F10.4,F8.0,21(1pG13.4))") (i_h-0.5)*d_h, &
         (trMax(i_h)+it_min)*dN_tr*AtmHei_step, (WeiFie_x(i_h,trhMax(2)+i*5),i=-5,+5)
   EndDo ! i_h=1-N_h,N_h
   Close(Unit=4)
   write(2,"(A,2(1pG13.4))") ' MaxWF,', MaxWF
   !
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'r-SrcKrnl-'//trim(label)//'_'//trim(txt)//'.grd')
   write(2,*) 'writing file "',TRIM(PlotBase)//'r-SrcKrnl-'//trim(txt)//'.grd','"'
   MaxWF=MAXVAL(ABS(WeiFie_r))
   write(4,"(2(1pG13.4),6(1pG13.4), (1pG13.4),0pF6.0,2I5,A)") zeta_f, rs_0, 1.-N_h, 1.*N_h, d_h, &
         it_min, it_min+N_tr, &
         dN_tr*AtmHei_step, MaxWF, phi_0*180/pi, trhMax(1), trhMax(2), ', ! WideSource_Kernel'
   !freadln inchan z0 rs0 ih1 ih2 delh it1 it2 dt maxK phi_f h_max itmax
   write(4,"('! Generating parameters, zeta_f=',F7.0,' rs_0=',F7.1,' Z_min=',F7.0,' Z_max=',F7.0,' phi_0=',F6.0, &
      ' N_phi0=',I3,' label=',A,' N_h=',I4.0,' h_range=',F5.1,' N_phir=',I3,' N_ca=',I4,' dN_tr=',I3,' DistMax=',F7.0)") &
      zeta_f, rs_0, Z_min, Z_max, phi_0, N_phi0, label, N_h, h_range, N_phir, N_ca, dN_tr, DistMax
   If((trhMax(2)+25) .gt. N_tr) trhMax(2)=N_tr-25
   If((trhMax(2)-25) .lt. 0) trhMax(2)=25
   write(4,"(A, 21(F6.0,7x))") '!    h,  t_r(max)  W@t_r=',(((trhMax(2)+i*5)+it_min)*dN_tr*AtmHei_step,i=-5,+5)
   Do  i_h=1-N_h,N_h
      write(4,"(F10.4,F8.0,21(1pG13.4))") (i_h-0.5)*d_h, &
         (trMax(i_h)+it_min)*dN_tr*AtmHei_step, (WeiFie_r(i_h,trhMax(2)+i*5),i=-5,+5)
   EndDo ! i_h=1-N_h,N_h
   Close(Unit=4)
   write(2,"(A,2(1pG13.4))") ' MaxWF,', MaxWF
   !
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)//'_'//trim(txt)//'.z') ! E13.4
   write(2,*) 'writing file "',TRIM(PlotBase)//'x-SrcKrnl-'//trim(txt)//'.z','"'
   write(4,"('! nx ',I4,' ny ',I5,' xmin ',F6.1,' xmax ',F6.1,' ymin ',F7.0,' ymax ',F7.0)") &
      2*N_h, N_tr+1, (0.5-N_h)*d_h, (N_h-0.5)*d_h, it_min*dN_tr*AtmHei_step, (N_tr+it_min)*dN_tr*AtmHei_step
   write(FMTSK,"(A,I5,A)") '("!",3I4,I5,I4,2G13.6,A,', 2*N_h, 'F13.4)'  !  "(A,200F13.4)"
   write(4,FMTSK) N_tr, dN_tr, it_min, 1-N_h, N_h, zeta_f, d_h, & ! i_atm=(i_tr+it_min)*dN_tr
         ' h-columns:',((i-0.5)*d_h,i=5-N_h,N_h)
   write(FMTSK,"(A,I4,A)") '(',2*N_h,'(1pG13.4))'  !  '(',2*N_h,'E13.4)'
   Do i_tr=0,N_tr  !   i_atm=(i_tr+it_min)*dN_tr
      write(4,FMTSK) WeiFie_x(1-N_h:N_h,i_tr)
   EndDo
   Close(Unit=4)
   i_tr=20
   write(2,*) i_tr,it_min,dN_tr, (i_tr)*dN_tr
   Write(2,"(30g13.4)") '!WeiFie_x(-4:5, 20)=', WeiFie_x(-4:5, i_tr-it_min)
   !
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'r-SrcKrnl-'//trim(label)//'_'//trim(txt)//'.z') ! E13.4
   write(2,*) 'writing file "',TRIM(PlotBase)//'r-SrcKrnl-'//trim(txt)//'.z','"'
   write(4,"('! nx ',I4,' ny ',I5,' xmin ',F6.1,' xmax ',F6.1,' ymin ',F7.0,' ymax ',F7.0)") &
      2*N_h, N_tr+1, (0.5-N_h)*d_h, (N_h-0.5)*d_h, it_min*dN_tr*AtmHei_step, (N_tr+it_min)*dN_tr*AtmHei_step
   write(FMTSK,"(A,I5,A)") '("!",3I4,I5,I4,2G13.6,A,', 2*N_h, 'F13.4)'  !  "(A,200F13.4)"
   write(4,FMTSK) N_tr, dN_tr, it_min, 1-N_h, N_h, zeta_f, d_h, & ! i_atm=(i_tr+it_min)*dN_tr
         ' h-columns:',((i-0.5)*SampleTime_new,i=5-N_h,N_h)
   write(FMTSK,"(A,I4,A)") '(',2*N_h,'(1pG13.4))'  !  '(',2*N_h,'E13.4)'
   Do i_tr=0,N_tr  !   i_atm=(i_tr+it_min)*dN_tr
      write(4,FMTSK) WeiFie_r(1-N_h:N_h,i_tr)
   EndDo
   Close(Unit=4)
   !
   PlotType='SrcKrnlFiltMap'
   PlotName=trim(PlotBase)//'SrcKrnlMap_x-'//trim(label)//'_'//trim(txt)
   PlotDataFile=TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)//'_'//trim(txt)
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Submit, CleanPlotFile
   !PlotType=TRIM(PlotBase)//'SrcKrnl-'//trim(label)//'_'//trim(txt)//'.z'
   !Call MakeGLEplots(CleanPlotFile=PlotType) !, Submit, CleanPlotFile
   !
   PlotType='SrcKrnlFiltMap'
   PlotName=trim(PlotBase)//'SrcKrnlMap_r-'//trim(label)//'_'//trim(txt)
   PlotDataFile=TRIM(PlotBase)//'r-SrcKrnl-'//trim(label)//'_'//trim(txt)
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Submit, CleanPlotFile
   !
   !Call MakeGLEplots(Submit=.true.)
   !stop 'Wource_Kernel-end'
   !
   ! Remainder of this routine is obsolete because it is moved to "BeamForm"
   !
   !-------------------------------------------------------------
   ! Perform frequency filtering
   SampleTime_new=5*c_l  ! 5~ns converted to units [m]
   !SampleTime_new=d_h  ! do not change the sampling [m]
   !
   Padding=1+50*SampleTime_new/d_h  ! allow for at least 50 new samples left and right of the region of interest.
   !Padding=0  ! allow for at least 50 new samples left and right of the region of interest.
   write(2,*) 'SampleTime_old, SampleTime_new:', d_h, SampleTime_new, &
      '[m], nu_min, nu_max:', nu_min, nu_max,'[MHz]'
   Call FreqFilterKernel(WeiFie_x, N_h, N_tr, d_h, padding, SampleTime_new, nu_min, nu_max, &
      K_nu, K_h, NuDim_new)
   !
   write(2,*) 'DownSampled; txt=',txt, NuDim_new, label
   ! ---------------------------------
   i_h=2000./(dN_tr*AtmHei_step)
   i_t=(zeta_f/(AtmHei_step*dN_tr)-it_min)   !  count corresponding to t_r=zeta_f
   i_h=min(i_h,N_h,i_t-1,N_tr-i_t-1,30)
   write(2,*) 'i_h,N_h,i_t-1,N_tr-i_t-1:',i_h,N_h,i_t-1,N_tr-i_t-1
   flush(unit=2)
   trMax(-i_h:i_h)=MAXLOC(ABS(WeiFie_x(1-N_h:N_h,i_t-i_h:i_t+i_h)),1)-N_h  ! determined the h where a max occurs for fixed t_r with (i_tr=i_t-i_h:i_t+i_h)
   write(2,"(A,61F10.1)") '  t_r [m]:',((i+it_min)*dN_tr*AtmHei_step,i=i_t-i_h,i_t+i_h)
   write(2,"(A,61F10.3)") 'h_max [m]:',((trMax(i)-0.5)*SampleTime_new,i=-i_h,+i_h)
   write(2,"(A,61g10.2)") 'Value@max:',(WeiFie_x(i,i_t+i),i=-i_h,+i_h)
   !
   !
   trhMax(:)=MaxLoc(ABS(K_h(1-nudim_new:nudim_new,0:N_tr)))
   trhMax(2)=trhMax(2)-1
   trhMax(1)=trhMax(1)-NuDim_new
   write(2,*) 'Filtered, trhMax(1:2)',trhMax(1:2), ', max=', K_h(trhMax(1),trhMax(2)) ! , ', trhMin(1:2)',trhMin(1:2)
   trMax(1-NuDim_new:NuDim_new)=MAXLOC(K_h(1-NuDim_new:NuDim_new,0:N_tr),2)-1
   write(2,*) 'Filtered, trMax(1-NuDim_new:NuDim_new)',trMax(1-NuDim_new:NuDim_new)
   write(2,*) 'Filtered, Max(1-NuDim_new:0)',K_h(1-NuDim_new:0,trMax(1-NuDim_new:0))
   flush(unit=2)
   MaxWF=MAXVAL(K_h)
   !base_h=45
   !hCount=NuDim_new-base_h
   !hcnt_l=-(NuDim_new-45 -1)
   !hcnt_u=NuDim_new-45
   hcnt_l = -40./SampleTime_new+0.5  ! corresponds to h=-40 [m]
   hcnt_u = +50./SampleTime_new+0.5  ! corresponds to h=+50 [m]
   write(2,*) 'writing file "',TRIM(PlotBase)//'SrcKFilt-'//trim(label)//'_'//trim(txt)//'.grd','"'
   flush(unit=2)
   !
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'SrcKFilt-'//trim(label)//'_'//trim(txt)//'.grd')
   write(4,"(2(1pG13.4),2(2I5,1pG13.4), (1pG13.4),0pF6.0,2I5,A)") &
         zeta_f, rs_0, hcnt_l, hcnt_u, SampleTime_new, it_min, it_min+N_tr, &
         dN_tr*AtmHei_step,  MaxWF, phi_0*180/pi, trhMax(1), (trhMax(2)+it_min), ', ! WideSource_Kernel Filtered'
   write(4,"('! Generating parameters, zeta_f=',F7.0,' rs_0=',F7.1,' Z_min=',F7.0,' Z_max=',F7.0,' phi_0=',F6.0, &
      ' N_phi0=',I3,' label=',A,' N_h=',I4.0,' d_h=',F5.3,' N_phir=',I3,' N_ca=',I4,' dN_tr=',I3,' DistMax=',F7.0)") &
      zeta_f, rs_0, Z_min, Z_max, phi_0, N_phi0, label, NuDim_new, SampleTime_new, N_phir, N_ca, dN_tr, DistMax
   write(4,"(A, 21(F6.0,7x))") '!    h,  t_r(max)  W@t_r=',(((trhMax(2)+i*5)+it_min)*dN_tr*AtmHei_step,i=-5,+5)
   Flush(unit=4)
   write(2,*) 'NuDim_new:',NuDim_new, trhMax(2)
   If((trhMax(2)+25) .gt. N_tr) trhMax(2)=N_tr-25
   If((trhMax(2)-25) .lt. 0) trhMax(2)=25
   write(2,*) 'NuDim_new:',NuDim_new, trhMax(2), N_tr
   Flush(Unit=2)
   R_z=0.
   X_tst=0.
   Do  i_h=hcnt_l, hcnt_u
      write(4,"(F10.4,F8.0,21(1pG13.4))") (i_h-0.5)*SampleTime_new, &
         (trMax(i_h)+it_min)*dN_tr*AtmHei_step, (K_h(i_h,trhMax(2)+i*5),i=-5,+5)
      !write(4,"(F10.4,F8.0,21(1pG13.4))") (i_h-0.5)*SampleTime_new, &
      !   trMax(i_h)*AtmHei_step, (K_h(i_h,trhMax(2)+i*5-it_min+1),i=-10,+10)
      Sx_h=0.
      Rx_h=0.
      Do i_tr=0,N_tr ! Loop over source time for obtaining weighting factors for shower current
         i_atm=(i_tr+it_min)*dN_tr
         Sx_h=Sx_h +K_h(i_h,i_tr)*AtmHei_step*dN_tr*Ix(i_atm)/nrm  ! after summing = E(h,zeta_f)
      EndDo
      R_z = R_z + Sx_h*Sx_h  ! should equal Sum_h[ E^2(h,zeta_f) ]
      !Write(2,*) '!(Filtered kernel* current) @h=', (i_h-0.5)*SampleTime_new, Sx_h ! , Rx_h
      If(ABS((i_h-0.5)*SampleTime_new) .lt. 30.) X_tst=X_tst + Sx_h*Sx_h
   EndDo ! i_h=1-N_h,N_h
   Close(Unit=4)
   write(2,"(A,(1pG13.4),A,(1pG13.4),A,(1pG13.4),A,(1pG13.4))") 'MaxWF:', MaxWF, &
      ', sqrt(E^2, h-integrated) = R(E) =', sqrt(R_z), ', dt=',SampleTime_new ,', h from -30 to +30 gives R(E)=',sqrt(X_tst)
!MaxWF:   1.7293E-05, sq-integrated gives R^2, R:    15.91        3.989    , dt=    1.499
! writing file "plot00-up-L040/SrcKFilt-13_05500000.z"
!integrated after folding with PSF gives R:    3.216        3.938    , dt=    1.499       1.9354E-05
   !
   !Write(2,*) 'NuDim_new:', NuDim_new, SampleTime_new
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'SrcKFilt-'//trim(label)//'_'//trim(txt)//'.z') ! E13.4
   write(2,*) 'writing file "',TRIM(PlotBase)//'SrcKFilt-'//trim(label)//'_'//trim(txt)//'.z','"'
   flush(unit=2)
   !write(FMTSK,"(A,I3,A)") '(',2*hCount,'(1pG13.4))'  !  '(',2*N_h,'E13.4)'
   write(FMTSK,"(A,I3,A)") '(', hcnt_u-hcnt_l+1, '(1pG13.4))'  !  '(',2*N_h,'E13.4)'
   write(4,"('! nx ',I3,' ny ',I4,' xmin ',F6.1,' xmax ',F6.1,' ymin ',F7.0,' ymax ',F7.0,g13.6)") &
     ! 2*hCount, N_tr+1, (0.5-hCount)*SampleTime_new, (hCount-0.5)*SampleTime_new, &
      (hcnt_u-hcnt_l+1), N_tr+1, (hcnt_l-0.5)*SampleTime_new, (hcnt_u-0.5)*SampleTime_new, &
      it_min*dN_tr*AtmHei_step, (N_tr+it_min)*dN_tr*AtmHei_step
   write(4,"('!',3I4,I4,I3,2G13.6,A,200F13.4)") N_tr, dN_tr, it_min, hcnt_l, hcnt_u, zeta_f, SampleTime_new, & ! i_atm=(i_tr+it_min)*dN_tr
         ' h-columns:',((i-0.5)*SampleTime_new,i=5+hcnt_l,hcnt_u)
   Do i_tr=0,N_tr  !   i_atm=(i_tr+it_min)*dN_tr
      write(4,FMTSK) K_h(hcnt_l:hcnt_u, i_tr)
   EndDo
   Close(Unit=4)
   !
   PlotType='SrcKrnlFiltMap'
   PlotName=trim(PlotBase)//'SrcKrnlFiltMap-'//trim(label)//'_'//trim(txt)
   PlotDataFile=TRIM(PlotBase)//'SrcKFilt-'//trim(label)//'_'//trim(txt)
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Submit, CleanPlotFile
   PlotType=TRIM(PlotBase)//'*cdata.dat'
   Call MakeGLEplots(CleanPlotFile=PlotType) !, Submit, CleanPlotFile
   PlotType=TRIM(PlotBase)//'*clabels.dat'
   Call MakeGLEplots(CleanPlotFile=PlotType) !, Submit, CleanPlotFile
   PlotType=TRIM(PlotBase)//'*cvalues.dat'
   Call MakeGLEplots(CleanPlotFile=PlotType) !, Submit, CleanPlotFile
   !PlotType=TRIM(PlotBase)//'SrcKFilt-'//trim(label)//'_'//trim(txt)//'.z'
   !Call MakeGLEplots(CleanPlotFile=PlotType) !, Submit, CleanPlotFile
   !-------------------------------------------------------
   !  Power analysis filtered kernel
   !
   !
   !     coherent power
   i_zf=zeta_f/AtmHei_step
   i_z0=(i_zf/dN_tr-it_min)
   !write(2,*) 'i_z0:', i_z0, i_zf, zeta_f, it_min, dN_tr
   K_trace(1-nudim_new:nudim_new)=K_h(1-nudim_new:nudim_new, i_z0)
   K_CohPwr(:) = 0.
   K_CohPwr(0)= SUM(K_trace(hcnt_l:hcnt_u)**2)
   dK_CohPwr(0)=K_CohPwr(0)
   Do i=1,N_tr-i_z0
      If(i_z0-i .le. 0) exit
      K_trace(1-nudim_new:nudim_new)=K_trace(1-nudim_new:nudim_new)+ K_h(1-nudim_new:nudim_new, i_z0-i)  ! SUM(K_h(1-hCount:hCount, i_z0-i:i_z0+i-1))
      !write(2,*) 'K_trace-i:',i,K_trace(1-hCount:hCount)
      !write(2,*) 'K_trace-0i:',Sum(K_h(0, i_z0-i:i_z0+i-1)),K_h(0, i_z0-i),K_h(0, i_z0+i-1),i_z0-i,i_z0+i-1
      K_CohPwr(-i)= SUM(K_trace(1-nudim_new:nudim_new)**2)
      dK_CohPwr(-i)= 0.5*(K_CohPwr(-i) - K_CohPwr(i-1))/(2.*i)
      K_trace(1-nudim_new:nudim_new)=K_trace(1-nudim_new:nudim_new)+ K_h(1-nudim_new:nudim_new, i_z0+i)  !  =SUM(K_h(1-hCount:hCount, i_z0-i:i_z0+i))
      K_CohPwr(i)= SUM(K_trace(1-nudim_new:nudim_new)**2)
      dK_CohPwr(i)= 0.5*(K_CohPwr(i) - K_CohPwr(-i))/(2.*i+1.)
      If((i_z0-i).le.1) exit
   End   Do ! i=0,it_max-i_zf
   !     incoherent power
   Do i_tr=0,N_tr  !   i_atm=(i_tr+it_min)*dN_tr
      K_IncPwr(i_tr)= SUM(K_h(1-nudim_new:nudim_new, i_tr)**2)
   EndDo
   dwr=K_IncPwr(i_z0)
   K_IncPwr(:)=sqrt(K_IncPwr(:)/dwr)
   !Write(2,*) 'dwr:', dwr, i_z0, N_tr
   !     Pulse shape folding
   !If( K_h(1, i_z0) .lt.0.) Then  ! make sure the PSF(t_r=z) is positive
   !   K_h(:, :)=-K_h(:, :)
   !EndIf
   K_trace(1-nudim_new:nudim_new)=K_h(1-nudim_new:nudim_new, i_z0)  !  PSF
   R_z=0.
   X_tst=0.
   Do i_tr=0,N_tr  !   i_atm=(i_tr+it_min)*dN_tr
      K_PSF(i_tr)= SUM(K_h(1-nudim_new:nudim_new, i_tr)*K_trace(1-nudim_new:nudim_new))/dwr
      i_atm=(i_tr+it_min)*dN_tr
      R_z = R_z + K_PSF(i_tr)*AtmHei_step*dN_tr*(Ix(i_atm)/nrm) ! *SampleTime_new
      X_tst=X_tst + K_IncPwr(i_tr)*AtmHei_step*dN_tr*(Ix(i_atm)/nrm)
   EndDo
   write(2,"(A,(1pG13.4),A,(1pG13.4), A,(1pG13.4))") 'Folding with PSF gives R_PSF=', &
      R_z*sqrt(dwr), ', norm=', sqrt(dwr), ', folding with Weights gives R_Pwr=', X_tst*sqrt(dwr)
   !
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'KPSF-'//trim(label)//'_'//trim(txt)//'.dat') ! E13.4
   write(2,*) 'writing file "',TRIM(PlotBase)//'KPSF-'//trim(txt)//'.dat','"'
   write(4,*) '!', sqrt(dwr), nudim_new, SampleTime_new, hcnt_l, hcnt_u, zeta_f
   flush(unit=2)
   Do i_h=1-nudim_new,nudim_new  !   h= (i_h-0.5)*SampleTime_new
      write(4,"(I5,3(1pG13.4))") i_h, K_trace(i_h)/sqrt(dwr)
   EndDo
   Close(Unit=4)
   !
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'KPwr-'//trim(label)//'_'//trim(txt)//'.dat') ! E13.4
   write(2,*) 'writing file "',TRIM(PlotBase)//'KPwr-'//trim(txt)//'.dat','"'
   write(4,*) '!', sqrt(dwr), N_tr, zeta_f
   flush(unit=2)
   Do i_tr=0,N_tr  !   i_atm=(i_tr+it_min)*dN_tr  ! dK_CohPwr not really used anywhere
      write(4,"(F8.1,3(1pG13.4))") (i_tr+it_min)*dN_tr*AtmHei_step, K_IncPwr(i_tr), dK_CohPwr(i_tr-i_z0), K_PSF(i_tr)
   EndDo
   Close(Unit=4)
   !
   Deallocate( K_nu, K_h, WeiFie_x, WeiFie_r )
   !
   !Call MakeGLEplots(Submit=.true.)
   !stop 'Wource_Kernel-end'
!    EndDo
    return
End Subroutine WideSource_Kernel
