!------------------------------
Subroutine KernelCreate(PlotBase, zeta_f, rs_0, Z_min, Z_max, N_phi0, phi_0, label, &
      N_h, d_h, dN_tr, N_phir, d_crMax)
!   Calculate the source @ height=\zeta as function of tmh=(t-h) for a fixed observerdistance to this shower
!   A backtrace integration is performed over observer position
!   CD_i (CoreDist_A(CD_i)) determines the position of the pencil shower w.r.t. the core of the real shower which determines the pancake thickness.
!   zeta= height in the atmosphere of the back-trace voxel
!   Shower time varies with h, distance behind the shower front.
!
   use Atmosphere, only : xi,dxi, AtmHei_dim, AtmHei_step
   Use ShowerData, only : d_ant, DistMax, N_antennas
   !use eventdata, only : PlotBase
   use constants, only : dp, pi, c_l
   Use Pancakefunction, only : CoreDist_Dim, CoreDist_A
   Use Systemcommands, only : MakeGLEplots
   Use Pancakefunction, only : GetLambda
   Use Pancakefunction, only : fiehLam, alpha_tr
   use LateralDistribution, only : W_tc
   implicit none
   character(Len=*), intent(in)  :: PlotBase
   real(dp), intent(in) :: zeta_f, rs_0, Z_min, Z_max  ! zeta_f is the distance from ground-core, rs_0 is distance from shower axis of focus point
   real(dp), intent(inout) :: phi_0
   Integer, intent(in) :: N_phi0
   Character(len=*), intent(in) :: label
   Integer, intent(in) :: N_h
   Real(dp), intent(in) :: d_h
   Integer, intent(in) :: dN_tr
   Integer, intent(in) :: N_phir
   Real(dp), intent(in) :: d_crMax
   !
   !real(dp), parameter :: h_range=30.d0 ! h_range=60.d0 ! h_range=50.d0 ! h_range=100.d0 ! h_range=17.d0 ! h_range=200.d0 ! h_range=100.d0 ! h_range=20.d0 ! h_range=2.5d0  ! h_range=1.d0  !  [m]
   !Integer, parameter :: N_h=120 ! N_h=240 ! N_h=200 ! N_h=400 ! N_h=50 ! N_h=20 ! N_h=100 ! N_h=26  ! N_h=80 ! number of steps in distance behind front
   !Integer, parameter :: dN_tr=10 ! dN_tr=5 ! dN_tr=1 ! dN_tr=3 !  Step size in t_r in units of AtmHei_step
   ! fine N_phi grid important for power at lower altitudes
   !Integer, parameter :: N_phir=5 ! N_phir=15 ! N_phir=3 ! N_phir=7 ! N_phir=17 !  Number of angles for integration over shower area
   ! fine N_ca grid important for quenching fluctuations
   !   !Integer, parameter ::  N_ca= 50 ! N_ca= 100 ! N_ca= 200 ! N_ca= 250 ! N_ca= 50 !messy) ! grid for antenna to shower-core distance
   ! finer than d_crMax=20.d0 grid not important for dmax=250
   !real(dp), parameter :: d_crMax=20.d0  ! d_crMax=10.d0  ! [m]  ! max increment in core-ray grid
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
   real(dp) :: rc_a, t_a ! all _0 mark the quantities at the focus point of retraced rays
   real(dp) ::  zeta_s, h_s, NN_s, drefrac_s, tr_s
   Real(dp) :: rc_r, dc_r, xc_r, rr_a, dr_a, xr_a, phi_r, d_pr, c_pr, lim_ca, R_r, h_r, CD_d
   Integer :: i_pr, ic_r, NC_r, ir_a
   Integer, parameter :: CoreRaySteps_dim=500
   real(dp) :: CoreRaySteps(0:CoreRaySteps_dim)
   real(dp) :: di, weight, norm, nrm, wei, Factr
   real(dp) :: Kx_ca, Kx, dKx_h, Sx_h, Rx_h
   Real(dp), allocatable :: WeiFie_x(:,:)
   Real(dp) :: MaxWF !, R_z
   integer :: i, i_t, i_h, i_tr, i_phi, CD_i, it_min, it_max, N_tr, i_atm, i_ant
   real(dp) :: lambda, DerLambda, Alpha, DerAlpha, dwr, rcr_max, fh, dfh, dfhdr, D_s
   Integer :: trhMax(1:2),trMax(1-N_h-52:N_h+52) !,
   character(len=10) :: txt
   character(len=45) :: FMTSK
   Integer :: i_zf, i_z0, N_ra
   Real(dp) :: K_CohPwr(-AtmHei_dim/2:AtmHei_dim/2), dK_CohPwr(-AtmHei_dim/2:AtmHei_dim/2)
   Character(len=50) :: PlotType, PlotName, PlotDataFile
   real(dp) :: h_range
   !
   h_range=d_h*N_h
   d_phi0=pi/N_phi0
   N_ra = 2.*DistMax/ d_crMax
   !
   write(2,*) '!KernelCreate, zeta_f, rs_0, Z_min, Z_max, N_phi0, phi_0:', zeta_f, rs_0, Z_min, Z_max, N_phi0, phi_0
   write(2,*) '!KernelCreate, N_h, d_h, dN_tr, N_phir, d_crMax:', N_h, d_h, dN_tr, N_phir, d_crMax
   write(2,"('input parameters: zeta_f=',F7.0,', rs_0=',F7.1,', Z_min=',F7.0,', Z_max=',F7.0,', phi_0=',F6.0, &
      ', N_phi0=',I3,', DistMax=',F7.0,', label=',A)") &
      zeta_f, rs_0, Z_min, Z_max, phi_0*180/pi, N_phi0, DistMax, label
   write(2,"('hard set parameters: N_h=',I4.0,', h_range=',F5.1,', dN_tr=',I3,', N_phir=',I3)") &
       N_h, h_range, dN_tr, N_phir
   Flush(unit=2)
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
   !allocate( WeiFie_r(1-N_h:N_h,0:N_tr) )
   !
   !-----------------------------------------------
   !
   NC_r=2*CoreDist_Dim-1  ! one less than in the definition of CoreDist_A
   rcr_max= CoreDist_A(CoreDist_Dim)
   ! write(2,*) '!CoreDist_A', CoreDist_A(0:4)
   CoreRaySteps(0)=0.
   CoreRaySteps(1)=CoreDist_A(0)/2.
   CoreRaySteps(2)=CoreDist_A(0)
   !CoreRaySteps(0:2)=CoreDist_A(0:2)
   Do ic_r=1, NC_r/2-1
      !CoreRaySteps(ic_r)=CoreDist_A(ic_r)
      !dc_r=(CoreDist_A(ic_r+1)-CoreDist_A(ic_r-1))/2.
      CoreRaySteps(2*ic_r+1)=(CoreDist_A(ic_r-1)+CoreDist_A(ic_r))/2.
      CoreRaySteps(2*ic_r+2)=CoreDist_A(ic_r)
      dc_r=(CoreDist_A(ic_r+1)-CoreDist_A(ic_r))/2.
      If(dc_r .gt. d_crMax) exit    ! step to (ic_r+1) too large
   EndDo
   If(ic_r.ge.NC_r/2-1) then  ! looped till the end =  CoreRaySteps_dim
      CoreRaySteps(NC_r+1)=CoreDist_A(NC_r/2)
   Else     ! refine final step sizes
      NC_r=(2*ic_r+2) + (CoreDist_A(CoreDist_Dim)-CoreDist_A(ic_r))/d_crMax
      If(NC_r .ge. CoreRaySteps_dim) NC_r=CoreRaySteps_dim-1
      Do i=2*ic_r+3,NC_r+1
         CoreRaySteps(i)=CoreRaySteps(i-1) + d_crMax
      EndDo
   EndIf
   !write(2,*) 'NC_r:', NC_r, CoreRaySteps(2*ic_r+1:2*ic_r+4),';', CoreRaySteps(NC_r), CoreRaySteps(NC_r+1)
   !
   nrm=1.e6
   norm=0.
   Do i_ant=1,N_antennas
      rc_a=d_ant(i_ant)
      R_0=sqrt(rc_a*rc_a +zeta_f*zeta_f)
      weight=2.*rc_a
      norm=norm+weight/R_0  !   /(R_0*R_0)
   EndDo ! idi=1, idi_N
   WeiFie_x(:,:) = 0.
   !R_z = 0.
   !
   !N_tr=2
   !write(2,"(A,I4,2F6.3,I4,F6.2)") 'N_tr, dN_tr, d_h, N_phir, d_crMax:', &
   !      N_tr, dN_tr, d_h, N_phir, d_crMax
   !write(2,*) '!(N_tr+it_min)*dN_tr=', (N_tr+it_min)*dN_tr, N_tr,it_min,dN_tr, it_max
   !flush(unit=2)
   Do i_h=1-N_h,N_h
      h_0=(i_h-0.5)*d_h
      !tr_0=h_0-zeta_f ! should be negative, all quantities at the source are given
      !write(2,*) 'Source_zeta; z,h=',zeta_f,h
      !Flush(unit=2)
      i=(zeta_f/AtmHei_step)
      if(i.lt.1) return
      If(i.lt. AtmHei_dim) Then
         di=zeta_f/AtmHei_step-i
         NN_0=1.+ di*xi(i+1) + (1.-di)*xi(i)
         dr_0=   di*dxi(i+1) + (1.-di)*dxi(i)
      Else
         NN_0=1.+ xi(AtmHei_dim)
         dr_0=   dxi(AtmHei_dim)
      EndIf
      !
      !write(2,*) 'dto:', 'dto, ObsDist, Rx/norm, Ex/(NN*R), Ex, cald '
      Sx_h=0.
      Kx=0.
      Do i_tr=0,N_tr ! Loop over source time for obtaining weighting factors for shower current
         i_atm=(i_tr+it_min)*dN_tr
         tr_s=-i_atm*AtmHei_step
         !
         If(i.lt. AtmHei_dim) Then
            Alpha=alpha_tr(i_atm) ! the parameter that determines the width of the pancake function
            DerAlpha=(alpha_tr(i_atm+1)-alpha_tr(i_atm-1))/(2*AtmHei_step) ! its derivative
         Else
            Alpha=alpha_tr(AtmHei_dim) ! the parameter that determines the width of the pancake function
            DerAlpha=(alpha_tr(AtmHei_dim)-alpha_tr(AtmHei_dim-1))/(AtmHei_step) ! its derivative
         EndIf
         !
         Do i_ant=1,N_antennas
            rc_a = d_ant(i_ant)
            Kx_ca =0
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
            weight=2.*rc_a*d_phi0   !*(DistMax-ObsDist)
            Do i_phi=1,N_phi0  ! angle between ObsDist and rs_0
               If(N_phi0.gt.1) phi_0=(i_phi-0.5)*d_phi0
               rsa_0=sqrt(rc_a*rc_a + rs_0*rs_0 - 2*rs_0*rc_a*cos(phi_0)) ! distance to shower
               R_0=sqrt(zeta_f*zeta_f+ rc_a*rc_a + rs_0*rs_0 - 2*rs_0*rc_a*cos(phi_0)) ! Distance between antenna and source in shower
               t_a=NN_0*R_0 + h_0-zeta_f ! time in antenna is derived quantity
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
               !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
               !   write(2,"(30g13.5)") '!Kx_c', zeta_f,i_h,i_tr, i_pr, rc_a &
               !   , NN_s, drefrac_s, t_a, tr_s, CoreRaySteps(0:5)
               Do i_pr=1,N_phir       ! angle of point off shower axis with antenna-axis direction
                  phi_r=(i_pr-0.5)*d_pr
                  c_pr=cos(phi_r)
                  dKx_h=0.
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
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !write(2,"(30g13.5)") '!Findhs:', t_a, rr_a, tr_s, h_r, NN_s, drefrac_s
                     D_s=NN_s* NN_s*(h_r-tr_s) + drefrac_s*(t_a-tr_s)*(t_a-tr_s)/NN_s
                     Call GetLambda(rc_r, Lambda, DerLambda, CD_i, CD_d)
                     If(h_r/(Lambda*Alpha) .gt. 20. ) cycle
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !   write(2,"(30g13.5)") '!GetLambda:', rc_r, Lambda, DerLambda, CD_i, CD_d
                     Call fiehLam(h_r,Lambda, DerLambda, Alpha,DerAlpha,fh,dfh,dfhdr)
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !   write(2,"(30g13.5)") '!fiehLam:',  h_r,Lambda, DerLambda, Alpha,DerAlpha,fh,dfh
                     wei=W_tc(rc_r,dwr)
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !   write(2,"(30g13.5)") '!W_tc:',  wei, rc_r,dwr
! x-polarisation:
!        Ex =Ex  + weigh * (dfh*Jx -fh *dJx) ! still needs to be times w(r), as done in following line
!         Ex_to(:,idi)=Ex_to(:,idi) + gth* ( wei *DEx(:) - DAxD(:)*weid ) !  (from LateralInt)
                     !dKx_h = dKx_h+ weight*Factr*dfh*wei*dc_r*d_pr/D_s
                     dKx_h = dKx_h + weight*Factr*dfh*wei * d_pr*dc_r/D_s
                     !   write(2,"(20g13.5)") '!dKx_h:', dKx_h, weight, Factr, dfh, wei,  d_pr, dc_r, D_s, ic_r, c_pr, rc_r &
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !   write(2,"(30g13.5)") '!Kx_', zeta_f,i_h,i_tr, i_pr, rc_a &
                     !   , dKx_h/norm,  ic_r  ,dfh, wei, h_r, rc_r,rr_a, h_r/Lambda !h_r-tr_s, rc_a!,'fh', fh,dfh
                  EndDo ! ic_r=0,CoreDist_Dim-1CD_i
                  !If(dKx_h.ne.0. .and. i_h .eq. -57) write(2,*) '!Kx_ca + 2*dKx_h:', i_h, i_tr, i_ant, Kx_ca, dKx_h
                  !Flush(unit=2)
                  Kx_ca = Kx_ca + 2*dKx_h
                  !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                  !   write(2,"(30g13.5)") '!Kx_ca', zeta_f,i_h,i_tr, i_pr, rc_a &
                  !   , Kx_ca/norm, dKx_h/norm, wei, rc_r, Lambda, rc_r, D_s
                  !
                  dr_a = d_crMax ! to be able to take care of rs=0
                  dKx_h=0.
                  Do ir_a=1, N_ra      ! distance ray to antenna, use regular grid
                     rr_a=(ir_a-0.5)*d_crMax
                     xr_a=rr_a*c_pr
                     If(xr_a.gt. lim_ca) exit ! getting closer to core than antenna
                     If((ir_a+0.5)* d_crMax *c_pr .gt. lim_ca) then ! next step would move out of range, take care of proper area (approx)
                        dr_a=lim_ca/c_pr-(ir_a-1)* d_crMax  ! remaining distance along the center line betwee observer and core
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
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !   write(2,"(30g13.5)") '!Kx_', zeta_f,i_h,i_tr, i_pr, rc_a &
                     !   , dKx_h/norm,  ic_r  ,dfh, wei, h_r, rc_r,rr_a, h_r/Lambda !h_r-tr_s, rc_a!,'fh', fh,dfh
                  EndDo ! ir_a=1, N_ra     ! distance ray to antenna, use regular grid
                  Kx_ca = Kx_ca + 2*dKx_h
                  !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                  !   write(2,"(30g13.5)") '!Kx_cb', zeta_f,i_h,i_tr, i_pr, rc_a &
                  !   , Kx_ca/norm, dKx_h/norm, wei, rc_r, Lambda, rc_r, D_s
               EndDo ! i_pr=1,      ! angle with antenna direction
            EndDo ! i_phi=1,N_phi0
            WeiFie_x(i_h,i_tr) = WeiFie_x(i_h,i_tr) + Kx_ca
         EndDo ! i_ca=1,N_ca
      EndDo ! i_tr=0,N_tr  !
      !
      !If(i_h .eq. -N_h/2) stop 'Wource_Kernel-i_h'
   EndDo ! i_h=1-N_h,N_h
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !stop 'Wource_Kernel-end'
   !Return !================================================================
   !11111111111111111111111111111!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   If(norm.gt.0.) Then
      WeiFie_x(:,:) = WeiFie_x(:,:)/norm
   EndIf
   !write(2,*) 'UnFil@i_h=8, i_tr=15-25', WeiFie(8,15:25)
   !
   !-------------------------------------------------------------
   write(txt,"(I3.3,i3.3,I2.2)") NINT(zeta_f/100.),NINT(rs_0),NINT(phi_0*10/pi)
   !
   trhMax(:)=MaxLoc(ABS(WeiFie_x(:,0:N_tr)))
   trhMax(2)=trhMax(2)-1
   trhMax(1)=trhMax(1)-N_h
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
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)//'.grd')
   write(2,*) 'writing file "',TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)//'.grd','"'
   write(4,"(2(1pG13.4),6(1pG13.4), (1pG13.4),0pF6.0,2I5,A)") zeta_f, rs_0, 1.-N_h, 1.*N_h, d_h, &
         it_min, it_min+N_tr, &
         dN_tr*AtmHei_step, MaxWF, phi_0*180/pi, trhMax(1), trhMax(2), ', ! WideSource_Kernel'
   !freadln inchan z0 rs0 ih1 ih2 delh it1 it2 dt maxK phi_f h_max itmax
   write(4,"('! Generating parameters, zeta_f=',F7.0,' rs_0=',F7.1,' Z_min=',F7.0,' Z_max=',F7.0,' phi_0=',F6.0, &
      ' N_phi0=',I3,' label=',A,' N_h=',I4.0,' h_range=',F5.1,' N_phir=',I3,'  dN_tr=',I3,' DistMax=',F7.0)") &
      zeta_f, rs_0, Z_min, Z_max, phi_0, N_phi0, label, N_h, h_range, N_phir, dN_tr, DistMax
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
   !
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)//'.z') ! E13.4
   write(2,*) 'writing file "',TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)//'.z','"'
   write(4,"('! nx ',I4,' ny ',I5,' xmin ',F6.1,' xmax ',F6.1,' ymin ',F6.2,' ymax ',F6.2)") &
      2*N_h, N_tr+1, (0.5-N_h)*d_h/c_l, (N_h-0.5)*d_h/c_l, &
      it_min*dN_tr*AtmHei_step/1000., (N_tr+it_min)*dN_tr*AtmHei_step/1000.
   write(4,"('! Generating parameters, zeta_f=',F7.0,' rs_0=',F7.1,' Z_min=',F7.0,' Z_max=',F7.0,' phi_0=',F6.0, &
      ' N_phi0=',I3,' label=',A,' N_h=',I4.0,' h_range=',F5.1,' N_phir=',I3,'  dN_tr=',I3,' DistMax=',F7.0)") &
      zeta_f, rs_0, Z_min, Z_max, phi_0, N_phi0, label, N_h, h_range, N_phir, dN_tr, DistMax
   write(FMTSK,"(A,I5,A)") '("!",3I4,I5,I4,2G13.6,A,', 2*N_h, 'F13.4)'  !  "(A,200F13.4)"
   write(4,FMTSK) N_tr, dN_tr, it_min, 1-N_h, N_h, zeta_f, d_h, & ! i_atm=(i_tr+it_min)*dN_tr
         ' h-columns:',((i-0.5)*d_h ,i=5-N_h,N_h)
   write(FMTSK,"(A,I4,A)") '(',2*N_h,'(1pG13.4))'  !  '(',2*N_h,'E13.4)'
   Do i_tr=0,N_tr  !   i_atm=(i_tr+it_min)*dN_tr
      write(4,FMTSK) WeiFie_x(1-N_h:N_h,i_tr)
   EndDo
   Close(Unit=4)
   i_tr=20
   !write(2,*) i_tr,it_min,dN_tr, (i_tr)*dN_tr
   !Write(2,"(30g13.4)") '!WeiFie_x(-4:5, 20)=', WeiFie_x(-4:5, i_tr-it_min)
   !
   !
   PlotType='SrcKrnlFiltMap'
   PlotName=trim(PlotBase)//'SrcKrnlMap_x-'//trim(label)
   PlotDataFile=TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile, mode='GenK') !, Submit, CleanPlotFile
   !PlotType=TRIM(PlotBase)//'SrcKrnl-'//trim(label)//'_'//trim(txt)//'.z'
   !Call MakeGLEplots(CleanPlotFile=TRIM(PlotDataFile)//'-c*.dat') !, Submit, CleanPlotFile
   !
   !
   !Call MakeGLEplots(Submit=.true.)
   !stop 'Wource_Kernel-end'
   !
   !Deallocate( K_nu, K_h, WeiFie_x )
   Deallocate(  WeiFie_x )
   !
   !Call MakeGLEplots(Submit=.true.)
   !stop 'Wource_Kernel-end'
!    EndDo
    return
End Subroutine KernelCreate


!------------------------------
Subroutine KernelCreateIndivid(PlotBase, zeta_f, rs_0, Z_min, Z_max, N_phi0, phi_0, label, &
      N_h, d_h, dN_tr, N_phir, d_crMax, theta_0)
!   Calculate the source @ height=\zeta as function of tmh=(t-h) for a fixed observerdistance to this shower
!   A backtrace integration is performed over observer position
!   CD_i (CoreDist_A(CD_i)) determines the position of the pencil shower w.r.t. the core of the real shower which determines the pancake thickness.
!   zeta= height in the atmosphere of the back-trace voxel
!   Shower time varies with h, distance behind the shower front.
!
   use Atmosphere, only : xi,dxi, AtmHei_dim, AtmHei_step
   Use ShowerData, only : d_ant, DistMax, N_antennas, theta_ant, Weights_ant
   !use eventdata, only : PlotBase
   use constants, only : dp, pi, c_l
   Use Pancakefunction, only : CoreDist_Dim, CoreDist_A
   Use Systemcommands, only : MakeGLEplots
   Use Pancakefunction, only : GetLambda
   Use Pancakefunction, only : fiehLam, alpha_tr
   use LateralDistribution, only : W_tc
   implicit none
   character(Len=*), intent(in)  :: PlotBase
   real(dp), intent(in) :: zeta_f, rs_0, Z_min, Z_max, theta_0  ! zeta_f is the distance from ground-core, rs_0 is distance from shower axis of focus point
   real(dp), intent(inout) :: phi_0
   Integer, intent(in) :: N_phi0
   Character(len=*), intent(in) :: label
   Integer, intent(in) :: N_h
   Real(dp), intent(in) :: d_h
   Integer, intent(in) :: dN_tr
   Integer, intent(in) :: N_phir
   Real(dp), intent(in) :: d_crMax
   !
   !real(dp), parameter :: h_range=30.d0 ! h_range=60.d0 ! h_range=50.d0 ! h_range=100.d0 ! h_range=17.d0 ! h_range=200.d0 ! h_range=100.d0 ! h_range=20.d0 ! h_range=2.5d0  ! h_range=1.d0  !  [m]
   !Integer, parameter :: N_h=120 ! N_h=240 ! N_h=200 ! N_h=400 ! N_h=50 ! N_h=20 ! N_h=100 ! N_h=26  ! N_h=80 ! number of steps in distance behind front
   !Integer, parameter :: dN_tr=10 ! dN_tr=5 ! dN_tr=1 ! dN_tr=3 !  Step size in t_r in units of AtmHei_step
   ! fine N_phi grid important for power at lower altitudes
   !Integer, parameter :: N_phir=5 ! N_phir=15 ! N_phir=3 ! N_phir=7 ! N_phir=17 !  Number of angles for integration over shower area
   ! fine N_ca grid important for quenching fluctuations
   !   !Integer, parameter ::  N_ca= 50 ! N_ca= 100 ! N_ca= 200 ! N_ca= 250 ! N_ca= 50 !messy) ! grid for antenna to shower-core distance
   ! finer than d_crMax=20.d0 grid not important for dmax=250
   !real(dp), parameter :: d_crMax=20.d0  ! d_crMax=10.d0  ! [m]  ! max increment in core-ray grid
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
   real(dp) :: rc_a, t_a, theta_a ! all _0 mark the quantities at the focus point of retraced rays
   real(dp) ::  zeta_s, h_s, NN_s, drefrac_s, tr_s
   Real(dp) :: rc_r, dc_r, xc_r, rr_a, dr_a, xr_a, phi_r, d_pr, c_pr, lim_ca, R_r, h_r, CD_d
   Integer :: i_pr, ic_r, NC_r, ir_a
   Integer, parameter :: CoreRaySteps_dim=500
   real(dp) :: CoreRaySteps(0:CoreRaySteps_dim)
   real(dp) :: di, weight, norm, nrm, wei, Factr
   real(dp) :: Kx_ca, Kx, dKx_h, Sx_h, Rx_h
   Real(dp), allocatable :: WeiFie_x(:,:)
   Real(dp) :: MaxWF !, R_z
   integer :: i, i_t, i_h, i_tr, i_phi, CD_i, it_min, it_max, N_tr, i_atm, i_ant
   real(dp) :: lambda, DerLambda, Alpha, DerAlpha, dwr, rcr_max, fh, dfh, dfhdr, D_s
   Integer :: trhMax(1:2),trMax(1-N_h-52:N_h+52) !,
   character(len=10) :: txt
   character(len=45) :: FMTSK
   Integer :: i_zf, i_z0, N_ra
   Real(dp) :: K_CohPwr(-AtmHei_dim/2:AtmHei_dim/2), dK_CohPwr(-AtmHei_dim/2:AtmHei_dim/2)
   Character(len=50) :: PlotType, PlotName, PlotDataFile
   real(dp) :: h_range
   !
   h_range=d_h*N_h
   d_phi0=pi/N_phi0
   N_ra = 2.*DistMax/ d_crMax
   !
   write(2,*) '!KernelCreate, zeta_f, rs_0, Z_min, Z_max, N_phi0, phi_0:', zeta_f, rs_0, Z_min, Z_max, N_phi0, phi_0
   write(2,*) '!KernelCreate, N_h, d_h, dN_tr, N_phir, d_crMax:', N_h, d_h, dN_tr, N_phir, d_crMax
   write(2,"('input parameters: zeta_f=',F7.0,', rs_0=',F7.1,', Z_min=',F7.0,', Z_max=',F7.0,', phi_0=',F6.0, &
      ', N_phi0=',I3,', DistMax=',F7.0,', label=',A)") &
      zeta_f, rs_0, Z_min, Z_max, phi_0*180/pi, N_phi0, DistMax, label
   write(2,"('hard set parameters: N_h=',I4.0,', h_range=',F5.1,', dN_tr=',I3,', N_phir=',I3)") &
       N_h, h_range, dN_tr, N_phir
   Flush(unit=2)
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
   !allocate( WeiFie_r(1-N_h:N_h,0:N_tr) )
   !
   !-----------------------------------------------
   !
   NC_r=2*CoreDist_Dim-1  ! one less than in the definition of CoreDist_A
   rcr_max= CoreDist_A(CoreDist_Dim)
   ! write(2,*) '!CoreDist_A', CoreDist_A(0:4)
   CoreRaySteps(0)=0.
   CoreRaySteps(1)=CoreDist_A(0)/2.
   CoreRaySteps(2)=CoreDist_A(0)
   !CoreRaySteps(0:2)=CoreDist_A(0:2)
   Do ic_r=1, NC_r/2-1
      !CoreRaySteps(ic_r)=CoreDist_A(ic_r)
      !dc_r=(CoreDist_A(ic_r+1)-CoreDist_A(ic_r-1))/2.
      CoreRaySteps(2*ic_r+1)=(CoreDist_A(ic_r-1)+CoreDist_A(ic_r))/2.
      CoreRaySteps(2*ic_r+2)=CoreDist_A(ic_r)
      dc_r=(CoreDist_A(ic_r+1)-CoreDist_A(ic_r))/2.
      If(dc_r .gt. d_crMax) exit    ! step to (ic_r+1) too large
   EndDo
   If(ic_r.ge.NC_r/2-1) then  ! looped till the end =  CoreRaySteps_dim
      CoreRaySteps(NC_r+1)=CoreDist_A(NC_r/2)
   Else     ! refine final step sizes
      NC_r=(2*ic_r+2) + (CoreDist_A(CoreDist_Dim)-CoreDist_A(ic_r))/d_crMax
      If(NC_r .ge. CoreRaySteps_dim) NC_r=CoreRaySteps_dim-1
      Do i=2*ic_r+3,NC_r+1
         CoreRaySteps(i)=CoreRaySteps(i-1) + d_crMax
      EndDo
   EndIf
   !write(2,*) 'NC_r:', NC_r, CoreRaySteps(2*ic_r+1:2*ic_r+4),';', CoreRaySteps(NC_r), CoreRaySteps(NC_r+1)
   !
   nrm=1.e6
   norm=0
   Do i_ant=1,N_antennas
      rc_a=d_ant(i_ant)
      R_0=sqrt(rc_a*rc_a +zeta_f*zeta_f)
      weight=Weights_ant(i_ant)
      norm=norm+weight !   /(R_0*R_0)
   EndDo ! idi=1, idi_N
   WeiFie_x(:,:) = 0.
   !R_z = 0.
   !
   !N_tr=2
   !write(2,"(A,I4,2F6.3,I4,F6.2)") 'N_tr, dN_tr, d_h, N_phir, d_crMax:', &
   !      N_tr, dN_tr, d_h, N_phir, d_crMax
   !write(2,*) '!(N_tr+it_min)*dN_tr=', (N_tr+it_min)*dN_tr, N_tr,it_min,dN_tr, it_max
   !flush(unit=2)
   Do i_h=1-N_h,N_h
      h_0=(i_h-0.5)*d_h
      !tr_0=h_0-zeta_f ! should be negative, all quantities at the source are given
      !write(2,*) 'Source_zeta; z,h=',zeta_f,h
      !Flush(unit=2)
      i=(zeta_f/AtmHei_step)
      if(i.lt.1) return
      If(i.lt. AtmHei_dim) Then
         di=zeta_f/AtmHei_step-i
         NN_0=1.+ di*xi(i+1) + (1.-di)*xi(i)
         dr_0=   di*dxi(i+1) + (1.-di)*dxi(i)
      Else
         NN_0=1.+ xi(AtmHei_dim)
         dr_0=   dxi(AtmHei_dim)
      EndIf
      !
      !write(2,*) 'dto:', 'dto, ObsDist, Rx/norm, Ex/(NN*R), Ex, cald '
      Sx_h=0.
      Kx=0.
      Do i_tr=0,N_tr ! Loop over source time for obtaining weighting factors for shower current
         i_atm=(i_tr+it_min)*dN_tr
         tr_s=-i_atm*AtmHei_step
         !
         If(i.lt. AtmHei_dim) Then
            Alpha=alpha_tr(i_atm) ! the parameter that determines the width of the pancake function
            DerAlpha=(alpha_tr(i_atm+1)-alpha_tr(i_atm-1))/(2*AtmHei_step) ! its derivative
         Else
            Alpha=alpha_tr(AtmHei_dim) ! the parameter that determines the width of the pancake function
            DerAlpha=(alpha_tr(AtmHei_dim)-alpha_tr(AtmHei_dim-1))/(AtmHei_step) ! its derivative
         EndIf
         !
         Do i_ant=1,N_antennas
            rc_a = d_ant(i_ant)
            theta_a = theta_ant(i_ant)
            Kx_ca =0
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
            weight=Weights_ant(i_ant)   !*(DistMax-ObsDist)
               rsa_0=sqrt(rc_a*rc_a + rs_0*rs_0 - 2*rs_0*rc_a*cos(theta_a-theta_0)) ! distance to shower
               R_0=sqrt(zeta_f*zeta_f+ rc_a*rc_a + rs_0*rs_0 - 2*rs_0*rc_a*cos(theta_a-theta_0)) ! Distance between antenna and source in shower
               t_a=NN_0*R_0 + h_0-zeta_f ! time in antenna is derived quantity
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
               ! h_s    = Distance behind showe            Do i_phi=1,N_phi0  ! angle between ObsDist and rs_0r front of source
               ! R_r    = distance from antenna to position on ray at height zeta_s
               ! lim_ca = rc_a/2. half way point (determines what to take as integration variable for ray position)
               lim_ca=rc_a/2.
               d_pr=pi/N_phir
               !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
               !   write(2,"(30g13.5)") '!Kx_c', zeta_f,i_h,i_tr, i_pr, rc_a &
               !   , NN_s, drefrac_s, t_a, tr_s, CoreRaySteps(0:5)
               Do i_pr=1,N_phir       ! angle of point off shower axis with antenna-axis direction
                  phi_r=(i_pr-0.5)*d_pr
                  c_pr=cos(phi_r)
                  dKx_h=0.
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
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !write(2,"(30g13.5)") '!Findhs:', t_a, rr_a, tr_s, h_r, NN_s, drefrac_s
                     D_s=NN_s* NN_s*(h_r-tr_s) + drefrac_s*(t_a-tr_s)*(t_a-tr_s)/NN_s
                     Call GetLambda(rc_r, Lambda, DerLambda, CD_i, CD_d)
                     If(h_r/(Lambda*Alpha) .gt. 20. ) cycle
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !   write(2,"(30g13.5)") '!GetLambda:', rc_r, Lambda, DerLambda, CD_i, CD_d
                     Call fiehLam(h_r,Lambda, DerLambda, Alpha,DerAlpha,fh,dfh,dfhdr)
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !   write(2,"(30g13.5)") '!fiehLam:',  h_r,Lambda, DerLambda, Alpha,DerAlpha,fh,dfh
                     wei=W_tc(rc_r,dwr)
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !   write(2,"(30g13.5)") '!W_tc:',  wei, rc_r,dwr
! x-polarisation:
!        Ex =Ex  + weigh * (dfh*Jx -fh *dJx) ! still needs to be times w(r), as done in following line
!         Ex_to(:,idi)=Ex_to(:,idi) + gth* ( wei *DEx(:) - DAxD(:)*weid ) !  (from LateralInt)
                     !dKx_h = dKx_h+ weight*Factr*dfh*wei*dc_r*d_pr/D_s
                     dKx_h = dKx_h + weight*Factr*dfh*wei * d_pr*dc_r/D_s
                     !   write(2,"(20g13.5)") '!dKx_h:', dKx_h, weight, Factr, dfh, wei,  d_pr, dc_r, D_s, ic_r, c_pr, rc_r &
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !   write(2,"(30g13.5)") '!Kx_', zeta_f,i_h,i_tr, i_pr, rc_a &
                     !   , dKx_h/norm,  ic_r  ,dfh, wei, h_r, rc_r,rr_a, h_r/Lambda !h_r-tr_s, rc_a!,'fh', fh,dfh
                  EndDo ! ic_r=0,CoreDist_Dim-1CD_i
                  !If(dKx_h.ne.0. .and. i_h .eq. -57) write(2,*) '!Kx_ca + 2*dKx_h:', i_h, i_tr, i_ant, Kx_ca, dKx_h
                  !Flush(unit=2)
                  Kx_ca = Kx_ca + 2*dKx_h
                  !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                  !   write(2,"(30g13.5)") '!Kx_ca', zeta_f,i_h,i_tr, i_pr, rc_a &
                  !   , Kx_ca/norm, dKx_h/norm, wei, rc_r, Lambda, rc_r, D_s
                  !
                  dr_a = d_crMax ! to be able to take care of rs=0
                  dKx_h=0.
                  Do ir_a=1, N_ra      ! distance ray to antenna, use regular grid
                     rr_a=(ir_a-0.5)*d_crMax
                     xr_a=rr_a*c_pr
                     If(xr_a.gt. lim_ca) exit ! getting closer to core than antenna
                     If((ir_a+0.5)* d_crMax *c_pr .gt. lim_ca) then ! next step would move out of range, take care of proper area (approx)
                        dr_a=lim_ca/c_pr-(ir_a-1)* d_crMax  ! remaining distance along the center line betwee observer and core
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
                     !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                     !   write(2,"(30g13.5)") '!Kx_', zeta_f,i_h,i_tr, i_pr, rc_a &
                     !   , dKx_h/norm,  ic_r  ,dfh, wei, h_r, rc_r,rr_a, h_r/Lambda !h_r-tr_s, rc_a!,'fh', fh,dfh
                  EndDo ! ir_a=1, N_ra     ! distance ray to antenna, use regular grid
                  Kx_ca = Kx_ca + 2*dKx_h
                  !If(NINT(zeta_f/500.).eq.5 .and. i_h.eq.1 .and. i_tr.eq.11 .and. NINT(rc_a/10).eq.8) &
                  !   write(2,"(30g13.5)") '!Kx_cb', zeta_f,i_h,i_tr, i_pr, rc_a &
                  !   , Kx_ca/norm, dKx_h/norm, wei, rc_r, Lambda, rc_r, D_s
               EndDo ! i_pr=1,      ! angle with antenna direction
            WeiFie_x(i_h,i_tr) = WeiFie_x(i_h,i_tr) + Kx_ca
         EndDo ! i_ca=1,N_ca
      EndDo ! i_tr=0,N_tr  !
      !
      !If(i_h .eq. -N_h/2) stop 'Wource_Kernel-i_h'
   EndDo ! i_h=1-N_h,N_h
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !stop 'Wource_Kernel-end'
   !Return !================================================================
   !11111111111111111111111111111!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   If(norm.gt.0.) Then
      WeiFie_x(:,:) = WeiFie_x(:,:)/norm
   EndIf
   !write(2,*) 'UnFil@i_h=8, i_tr=15-25', WeiFie(8,15:25)
   !
   !-------------------------------------------------------------
   write(txt,"(I3.3,i3.3,I2.2)") NINT(zeta_f/100.),NINT(rs_0),NINT(phi_0*10/pi)
   !
   trhMax(:)=MaxLoc(ABS(WeiFie_x(:,0:N_tr)))
   trhMax(2)=trhMax(2)-1
   trhMax(1)=trhMax(1)-N_h
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
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)//'.grd')
   write(2,*) 'writing file "',TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)//'.grd','"'
   write(4,"(2(1pG13.4),6(1pG13.4), (1pG13.4),0pF6.0,2I5,A)") zeta_f, rs_0, 1.-N_h, 1.*N_h, d_h, &
         it_min, it_min+N_tr, &
         dN_tr*AtmHei_step, MaxWF, phi_0*180/pi, trhMax(1), trhMax(2), ', ! WideSource_Kernel'
   !freadln inchan z0 rs0 ih1 ih2 delh it1 it2 dt maxK phi_f h_max itmax
   write(4,"('! Generating parameters, zeta_f=',F7.0,' rs_0=',F7.1,' Z_min=',F7.0,' Z_max=',F7.0,' phi_0=',F6.0, &
      ' N_phi0=',I3,' label=',A,' N_h=',I4.0,' h_range=',F5.1,' N_phir=',I3,'  dN_tr=',I3,' DistMax=',F7.0)") &
      zeta_f, rs_0, Z_min, Z_max, phi_0, N_phi0, label, N_h, h_range, N_phir, dN_tr, DistMax
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
   !
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)//'.z') ! E13.4
   write(2,*) 'writing file "',TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)//'.z','"'
   write(4,"('! nx ',I4,' ny ',I5,' xmin ',F6.1,' xmax ',F6.1,' ymin ',F6.2,' ymax ',F6.2)") &
      2*N_h, N_tr+1, (0.5-N_h)*d_h/c_l, (N_h-0.5)*d_h/c_l, &
      it_min*dN_tr*AtmHei_step/1000., (N_tr+it_min)*dN_tr*AtmHei_step/1000.
   write(4,"('! Generating parameters, zeta_f=',F7.0,' rs_0=',F7.1,' Z_min=',F7.0,' Z_max=',F7.0,' phi_0=',F6.0, &
      ' N_phi0=',I3,' label=',A,' N_h=',I4.0,' h_range=',F5.1,' N_phir=',I3,'  dN_tr=',I3,' DistMax=',F7.0)") &
      zeta_f, rs_0, Z_min, Z_max, phi_0, N_phi0, label, N_h, h_range, N_phir, dN_tr, DistMax
   write(FMTSK,"(A,I5,A)") '("!",3I4,I5,I4,2G13.6,A,', 2*N_h, 'F13.4)'  !  "(A,200F13.4)"
   write(4,FMTSK) N_tr, dN_tr, it_min, 1-N_h, N_h, zeta_f, d_h, & ! i_atm=(i_tr+it_min)*dN_tr
         ' h-columns:',((i-0.5)*d_h ,i=5-N_h,N_h)
   write(FMTSK,"(A,I4,A)") '(',2*N_h,'(1pG13.4))'  !  '(',2*N_h,'E13.4)'
   Do i_tr=0,N_tr  !   i_atm=(i_tr+it_min)*dN_tr
      write(4,FMTSK) WeiFie_x(1-N_h:N_h,i_tr)
   EndDo
   Close(Unit=4)
   i_tr=20
   !write(2,*) i_tr,it_min,dN_tr, (i_tr)*dN_tr
   !Write(2,"(30g13.4)") '!WeiFie_x(-4:5, 20)=', WeiFie_x(-4:5, i_tr-it_min)
   !
   !
   PlotType='SrcKrnlFiltMap'
   PlotName=trim(PlotBase)//'SrcKrnlMap_x-'//trim(label)
   PlotDataFile=TRIM(PlotBase)//'x-SrcKrnl-'//trim(label)
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile, mode='GenK') !, Submit, CleanPlotFile
   !PlotType=TRIM(PlotBase)//'SrcKrnl-'//trim(label)//'_'//trim(txt)//'.z'
   !Call MakeGLEplots(CleanPlotFile=TRIM(PlotDataFile)//'-c*.dat') !, Submit, CleanPlotFile
   !
   !
   !Call MakeGLEplots(Submit=.true.)
   !stop 'Wource_Kernel-end'
   !
   !Deallocate( K_nu, K_h, WeiFie_x )
   Deallocate(  WeiFie_x )
   !
   !Call MakeGLEplots(Submit=.true.)
   !stop 'Wource_Kernel-end'
!    EndDo
    return
End Subroutine KernelCreateIndivid
