Module ShowerData
   use constants, only : dp
   Integer :: N_antennas, N_timesamples, N_radial
   Real(dp) :: GroundLevel_sh, ZenithAngle_deg,ZenithAngle_rad, Azimuth_deg, EstCore_N, EstCore_E
   Real(dp) :: core_pos(1:2)  ! core position, determines in FitCorePos
   Real(dp), allocatable :: ant_pos(:,:)     ! [m](x/y,ant#); antenna position in (x,y) (whatever they are, in the shower plane)
   Real(dp), allocatable :: ant_toffset(:)   ! [ns](ant#); time offset first sample in time trace
   Real(dp), allocatable :: ant_xtrace(:,:), ant_ytrace(:,:)   ! [arb unit](sample,ant#); recorded time traces [arb unit](sample,ant#)
   Real(dp), allocatable :: Time_maxx(:), Val_maxx(:), Time_maxy(:), Val_maxy(:) !
   Real(dp), allocatable :: Time_errx(:), Time_erry(:)     ! [/sample^2]
   Real(dp), allocatable :: PulseTime_x(:), PulseTime_y(:)  ! [ns](ant#); pulse time, corrected for offset
   Real(dp), allocatable :: tx_fit(:), ty_fit(:), tr_fit(:)  ! [ns](radial); fit to pulse times
   Real(dp), allocatable :: d_ant(:),theta_ant(:),Weights_ant(:)    ! Distance of antenna to shower axis  (in shower plane), angle from x axis to antenna location in shower axis
   Real(dp),save :: DistMax    ! max Distance of antenna to shower core (in shower plane)
   Real(dp),save :: SamplingTime    ! sampling time [m]
   Integer, save :: Sample_Offset   ! sample in time trace that corresponds to t=0 when beamforming
   Real(dp),save :: nu_min, nu_max    ! [MHz] used block frequency filter
   Character(len=50), allocatable :: ant_file(:)
!Contains
!---------------------------
End Module ShowerData
!============================================
Subroutine ReadData(Base)
!  - Read time traces for individual antennas where antenna function is already corrected for.
!  - Fold with impulse response to find pulse-time. Note: identical to just find position of max in time trace.
!  - Fit pulse times to reconstruct the core position.
   use constants, only : dp, pi, c_l
   Use ShowerData, only : N_antennas, N_timesamples, N_radial, core_pos, d_ant, DistMax
   Use ShowerData, only : ant_pos, ant_toffset, ant_xtrace, ant_ytrace, ant_file, Sample_Offset
   Use ShowerData, only : SamplingTime, nu_min, nu_max
   Use ShowerData, only : Time_maxx, Time_errx, Val_maxx, Time_maxy, Time_erry, Val_maxy
   Use ShowerData, only : PulseTime_x, PulseTime_y, tx_fit, ty_fit, tr_fit
   Use ShowerData, only :  GroundLevel_sh, ZenithAngle_deg, ZenithAngle_rad, Azimuth_deg, EstCore_N, EstCore_E
   !use Atmosphere, only : AtmosphereInit
	implicit none
   !Integer, intent(out) :: N_antennas
   Character(len=*) :: Base
   Integer :: nxx, i, i_min, i_max, i_ant, N_samples
   Character(len=10) :: Label
   Character(len=20) :: Labela, labelb
   Real(dp) :: t_offsetx, t_offsety, NN_f, tTrace_Offset_dwn, D_max, D_time
   Real(dp) :: a,b,c,d , ant_x, ant_y ! scratch variables
   Real(dp), allocatable :: Shape_x(:), Shape_x_pp(:), Shape_y(:), Shape_y_pp(:)
   Character(len=200) :: LineIn
   !
   label='-test'
   OPEN(UNIT=4,STATUS='old',FILE=trim(Base)//trim(label)//'.dat',iostat=nxx) ! E13.4
   write(2,*) 'reading file: "', trim(Base)//trim(label)//'.dat'
   If(nxx.ne.0) Then
      write(2,*) 'problems reading file:', trim(Base)//trim(label)//'.dat'
      STOP  trim(Base)//trim(label)//'.dat is problematic'
   EndIf
!    write(4,"(3I5 )") N_ant, StParRange, GroundLevel, Zen_sh, Azi_sh, Core_N, Core_E, tTrace_dim_dwn
   Read(4,*) N_antennas, N_timesamples, GroundLevel_sh, ZenithAngle_deg, Azimuth_deg, EstCore_N, EstCore_E, nu_min, nu_max
   write(2,*) 'Use actual antenna positions, N_antennas, N_timesamples', N_antennas, N_timesamples, GroundLevel_sh &
         , ZenithAngle_deg, Azimuth_deg, EstCore_N, EstCore_E, nu_min, nu_max
   flush(unit=2)
   ZenithAngle_rad=ZenithAngle_deg*pi/180.
   !
   !-----------------------------------------------
   Allocate( ant_file(N_antennas), ant_pos(1:2,N_antennas), ant_toffset(N_antennas))
   Allocate( Shape_x(N_antennas), Shape_x_pp(N_antennas))
   Allocate( Shape_y(N_antennas), Shape_y_pp(N_antennas))
   Allocate( ant_xtrace(1:N_timesamples,N_antennas), ant_ytrace(1:N_timesamples,N_antennas) )
   Allocate( PulseTime_x(N_antennas), PulseTime_y(N_antennas) )
   Allocate( Time_maxx(N_antennas), Time_errx(N_antennas), Val_maxx(N_antennas) )
   Allocate( Time_maxy(N_antennas), Time_erry(N_antennas), Val_maxy(N_antennas) )
   Allocate( d_ant(N_antennas) )
   !
   Do i_ant=1,N_antennas
      !Read(4,"(A200)") LineIn
      !write(2,*) LineIn
      Read(4,*) ant_file(i_ant), ant_pos(1:2,i_ant), i_min, i_max, SamplingTime, tTrace_Offset_dwn
      !  tTrace_Offset_dwn = time corresponding to the first sample that is read-in
      write(2,*) i_ant,TRIM(ant_file(i_ant)),';', &
            ant_pos(1:2,i_ant), i_min, i_max, SamplingTime, tTrace_Offset_dwn
      !flush(unit=2)
      ant_toffset(i_ant)=i_min*SamplingTime-tTrace_Offset_dwn  !  [m]
!                 TRIM(PrntFile), Ant_North, Ant_East, &
!               i_min, MIN(i_min+StParRange-1,tTrace_dim_dwn), SamplingTime, tTrace_Offset_dwn ! t=(i*SamplingTime-tTrace_Offset_dwn)/c      Read(4,*) ant_xtrace(1:N_timesamples, i)
      OPEN(UNIT=14,STATUS='old',FILE=trim(ant_file(i_ant))//'.csv',iostat=nxx) ! E13.4
      !write(2,*) 'reading file: "', trim(ant_file(i_ant))//'.csv'
      If(nxx.ne.0) Then
         write(2,*) 'problems reading file:', trim(ant_file(i_ant))//'.csv'
         STOP  trim(ant_file(i_ant))//'.csv is problematic'
      EndIf
      N_samples=N_timesamples
      read(14,*) labela
      Do i=1,N_timesamples
         Read(14,*,iostat=nxx) a, ant_xtrace(i, i_ant), b, ant_ytrace(i, i_ant), c
         If(nxx.ne.0) Then
            write(2,*) '********* stop reading at i=',i
            ant_xtrace(i:N_timesamples, i_ant)=0.
            ant_ytrace(i:N_timesamples, i_ant)=0.
            Shape_x(i:N_timesamples)=0.
            Shape_y(i:N_timesamples)=0.
            N_samples=i-1
            exit
         EndIf
         Shape_x(i)=ant_xtrace(i, i_ant)*ant_xtrace(i, i_ant) + b*b
         Shape_y(i)=ant_ytrace(i, i_ant)*ant_ytrace(i, i_ant) + c*c
      EndDo
      Close(unit=14)
      !
      Call FindPulseMax(N_samples, Shape_x, Time_maxx(i_ant), Time_errx(i_ant), Val_maxx(i_ant) )
      PulseTime_x(i_ant)=Time_maxx(i_ant)*SamplingTime + ant_toffset(i_ant)
      i=Time_maxx(i_ant)
      If( Val_maxx(i_ant) .lt. (Shape_x(i)+Shape_x(i+1))/2.) &
         write(2,*) '*** Non-max in Hilbert_x@',trim(Base)//trim(label)//'.dat',';',Val_maxx(i_ant),i,Shape_x(i-1:i+2)
      !write(2,*) 'FindPulseMax-x:', Time_maxx(i_ant), SamplingTime, tTrace_Offset_dwn, Time_errx(i_ant)
      !write(2,*) 'FindPulseMax-x:', Time_maxx(i_ant), SamplingTime, tTrace_Offset_dwn, Time_errx(i_ant)
      Call FindPulseMax(N_samples, Shape_y, Time_maxy(i_ant), Time_erry(i_ant), Val_maxy(i_ant) )
      PulseTime_y(i_ant)=Time_maxy(i_ant)*SamplingTime + ant_toffset(i_ant)
      If( Val_maxy(i_ant) .lt. (Shape_y(i)+Shape_y(i+1))/2.) &
         write(2,*) '*** Non-max in Hilbert_y@',trim(Base)//trim(label)//'.dat',';',Val_maxy(i_ant),i,Shape_y(i-1:i+2)
      !write(2,*) 'FindPulseMax-y:', N_samples,Time_maxy(i_ant), Time_erry(i_ant)
      !
   EndDo
   write(2,*) 'SampleTime:', SamplingTime,'[m]=',SamplingTime/c_l,'[ns], nu_min, nu_max:', nu_min, nu_max,'[MHz]'
   !---------------------
   !  Find shower arrival time and correct toffset & derived quantities
   Sample_Offset=Time_maxx(1)+1  ! [samples] rather arbitrary, should correspond to t=0 while beamforming
   a =MINVAL( PulseTime_x(:),1 ) ! [m]
   b =MINVAL( PulseTime_y(:),1 )
   ! the more positive, the more back-shifted:
   c = +.125 -0.4* SamplingTime  ! 3.5 determined by adjusting the PSF position; reason unknown
   !ant_toffset(:) = ant_toffset(:) -0.25 +1.5* SamplingTime    ! [m]  Determined by fixing the peak in the PSF 1+1.5* SamplingTime=3
   write(2,*) 'Timing off-set correction=', -1.5 +0.0* SamplingTime  ,'[m]', +0.5 +1.0* SamplingTime,-0.25 +1.5* SamplingTime
   write(2,*) 'toffset:', Sample_Offset,'[samples],', Sample_Offset*SamplingTime,a,b, c,'[m]',SamplingTime
   ant_toffset(:) = ant_toffset(:) - c    ! [m]
   PulseTime_x(:) = PulseTime_x(:) - c
   PulseTime_y(:) = PulseTime_y(:) - c
   write(2,*) 'PulseTime(1)_:', PulseTime_x(1), PulseTime_y(1), ant_toffset(1),'[m]'
   !  ----------------------------
   !  Fit pulse times
   !
   D_time=20.     ! [m]
   D_max=400.
   N_radial=1+D_max/D_time
   core_pos(1)= 0.
   core_pos(2)= 0.
   Allocate( tx_fit(0:N_radial), ty_fit(0:N_radial), tr_fit(0:N_radial) )
   !Call FitPulseTimes(N_antennas, PulseTime_x, Time_errx, PulseTime_y, Time_erry, &
   !   ant_pos, core_pos, D_time, N_radial, tx_fit, ty_fit, tr_fit )
   Call FitCorePos(core_pos )
   !
   !core_pos(1)= 20.
   !core_pos(2)= -24.
   !write(2,*) '************ core position adjusted to:', core_pos
   Do i_ant=1,N_antennas
      ant_x= ant_pos(1,i_ant) - core_pos(1)
      ant_y= ant_pos(2,i_ant) - core_pos(2)
      d_ant(i_ant)= sqrt(ant_x**2 + ant_y**2)    ! Distance of antenna to shower axis
   EndDo ! i_ant=1,N_antennas
   DistMax = MaxVal( d_ant(:), 1 )
   !
   !stop 'ReadData'
   Return
End Subroutine ReadData
!----------------------------------------------------------
Subroutine FitCorePos(core_pos )
   use constants, only : dp, pi
   Use ShowerData, only : N_antennas, N_radial
   Use ShowerData, only : ant_pos
   Use ShowerData, only : PulseTime_x, PulseTime_y, Time_errx, Time_erry, tx_fit, ty_fit, tr_fit
   Implicit none
   Real(dp), intent(inout) :: core_pos(1:2)
   !
   Integer :: NF,j,k,m, error
   Real(dp) :: Qual
   integer ( kind = 4 ) :: meqn  ! Number of data points
   integer ( kind = 4 ) :: v_dim  ! dim=93 + n*p + 3*n + p*(3*p+33)/2 ,n-meqn, p=nvar
   integer ( kind = 4 ) :: nvar    ! number of parameters
   integer ( kind = 4 ) iv(65)   ! 60+nvar
   external FitPulseTimes  ! subroutine that compares FitFunction to Data
   external ufparm  ! Dummy external routine
   integer ( kind = 4 ) :: uiparm(1)    ! Not really used
   real ( kind = 8 ) :: urparm(1)    ! Not really used
   real ( kind = 8 ), allocatable :: v(:) ! dim=93 + n*p + 3*n + p*(3*p+33)/2 ,n-meqn, p=nvar
   real ( kind = 8 ) :: x(2) ! parameters that are optimized
   !
   meqn = 2*N_antennas
   nvar = 2
   X(1:2) = core_pos(1:2)
   !
   v_dim=93 + Meqn*nvar + 3*Meqn + nvar*(3*nvar+33)/2
   write(2,*) 'Fitting; v_dim are:', v_dim
   Allocate( v(v_dim) )

    call dfault( iv, v)
    iv(1) = 12 ! 12= do not call dfault again
    iv(14) = 1 ! 1: means print a covariance matrix at the solution.
    iv(15) = 2 ! if = 1 or 2, then a finite difference hessian approximation h is obtained.
               ! if positive: with step sizes determined using v(delta0=44), a multiplicative factor)
               ! If negative: then only function values are used with step sizes determined using v(dltfdc=40)
    iv(19) = 0 ! controls the number and length of iteration summary lines printed
    iv(21) = 2 ! is the output unit number on which all printing is done.
    iv(22) = 1 ! print out the value of x returned (as well as the corresponding gradient and scale vector d).
    iv(23) = 1 ! print summary statistics upon returning.
    iv(24) = 0 ! print the initial x and scale vector d
    v(32) =v(32)*1.d+4 ! is the relative function convergence tolerance
    v(36) =v(36)*1.d+3  ! step size for derivatives
    v(40) =v(40)*1.d+4 ! the step size used when computing the covariance matrix when iv(covreq=15) = -1 or -2, step size = v(dltfdc=40) * max(abs(x(i)), 1/d(i))
    v(44) =v(44)*1.d+4 ! the factor used in choosing the finite difference step size used in computing the covariance matrix when
                !    iv(covreq=15) = 1 or 2, step size = v(delta0=44) * max(abs(x(i)), 1/d(i)) * sign(x(i))
!
! iv(nfcall)... iv(6) is the number of calls so far made on calcr (i.e.,
!             function evaluations, including those used in computing
!             the covariance).
! iv(mxfcal)... iv(17) gives the maximum number of function evaluations
!             (calls on calcr, excluding those used to compute the co-
!             variance matrix) allowed.  if this number does not suf-
!             fice, then nl2sol returns with iv(1) = 9.  default = 200.
! iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
!             it also indirectly limits the number of gradient evalua-
!             tions (calls on calcj, excluding those used to compute
!             the covariance matrix) to iv(mxiter) + 1.  if iv(mxiter)
!             iterations do not suffice, then nl2sol returns with
!             iv(1) = 10.  default = 150.
! Of particular interest are the entries in v that limit the length of
!                  the first step attempted (lmax0), specify conver-
!                  gence tolerances (afctol, rfctol, xctol, xftol),
!                  and help choose the step size used in computing the
!                  covariance matrix (delta0).  see the section on
!                  (selected) v input values below.!    iv(18)=1
! v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the
!             very first step that nl2sol attempts.  it is also used
!             in testing for singular convergence -- if the function
!             reduction predicted for a step of length bounded by
!             v(lmax0) is at most v(rfctol) * abs(f0), where  f0  is
!             the function value at the start of the current iteration,
!             and if nl2sol does not return with iv(1) = 3, 4, 5, or 6,
!             then it returns with iv(1) = 7.    default = 100.
! iv(nfcov).... iv(40) is the number of calls made on calcr when
!             trying to compute covariance matrices.
! iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
!             calcj) so far done (including those used for computing
!             the covariance).
! iv(ngcov).... iv(41) is the number of calls made on calcj when
!             trying to compute covariance matrices.
! iv(niter).... iv(31) is the number of iterations performed.
!
   If(meqn .gt. nvar) then   ! otherwise the system is underdetermined
      write(2,*) 'nl2sol diagnostics, iv 1,15,16',IV(1),IV(15),IV(16)
       !V(44)=1.d+0
       !V(40)=1.d+0
      write(2,*) 'nl2sol diagnostics, v 32,36,40,44', V(32), V(36), V(40), V(44)
     call nl2sno ( meqn, nvar, x, FitPulseTimes, iv, v, uiparm, urparm, ufparm, error )
     !
     write(2,"('Result, chi^2/ndf=',g15.5)") 2*v(10)/(meqn-nvar)
     write(2,"('# fie calls=',i3,' and # of iterations=',i3)") iv(6),iv(31)
   endif
   !
   write(2,*) 'Best fit recalculation'
   write(*,*) 'Best fit recalculation:',x(1:nvar)
   !
   NF=-1
   call FitPulseTimes( meqn, nvar, x, NF, v, uiparm, urparm, ufparm )
   qual=sum(V(1:Meqn)*V(1:Meqn))
   core_pos(1:2)=x(1:2)
   !
   write(2,*) 'chi^2/df=', qual/(meqn-nvar), ', core @', core_pos(1:2)
   DeAllocate( v )
   return
End Subroutine FitCorePos
!----------------------------------------------------
Subroutine FitPulseTimes( meqn, nvar, x, NFit, R, uiparm, urparm, ufparm )
!     meqn = 2*N_antennas
!     nvar = 2
!     X(1:2) = core_pos(1:2)
!  Fit pulse arrival times with a linearly interpolated form, wher values at distances D_time are determined
!  Angular dependence fitted as _x + _r cos(phi) ; _y + _r sin(phi)
   Use ShowerData, only : N_antennas, N_radial
   Use ShowerData, only : ant_pos
   Use ShowerData, only : PulseTime_x, PulseTime_y, Time_errx, Time_erry, tx_fit, ty_fit, tr_fit
   use constants, only : dp, pi
   Implicit none
   Real, parameter :: Drad_min=20.  ! [m]
   Integer, intent(in) :: meqn, nvar, NFit
   Real(dp), intent(in) :: X(*)
   real ( kind = 8 ), intent(out) :: r(meqn)
   integer ( kind = 4 ), intent(in) :: uiparm(*)
   real ( kind = 8 ), intent(in) :: urparm(*)
   external ufparm
   !
   Real(dp) :: D_radial, core_pos(1:2)
   Integer :: Numb(0:N_radial), Njd(1:N_radial+1)
   Real(dp) :: A_x(0:N_radial), A_r(0:N_radial), A_y(0:N_radial)
   Real(dp) :: B_xx(0:N_radial,0:N_radial), B_xr(0:N_radial,0:N_radial), B_rr(0:N_radial,0:N_radial)
   Real(dp) :: B_yr(0:N_radial,0:N_radial), B_yy(0:N_radial,0:N_radial)
   Integer :: i, i_d, j_ant, j_d, N_points, NDoF, j_d_prev
   Real(dp) :: Ant_x, ant_y, Dist, SinTh, CosTh, del_j, Rx, Ry, t_x, t_y, t_r, theta, D_max, chisq
   !
   Character(len=1) :: UPLO
   Integer :: N, NRHS, LWork, INFO, Numb_lim
   Integer, allocatable :: IPIV(:)
   Real(dp), allocatable :: Res(:,:)
   Real(dp), allocatable :: ReWork(:), B(:,:)
   Real(dp) :: D, Dm(1)
   !
   D_max=0.
   If(meqn.ne.2*N_antennas) Then
      write(2,*) '***** error in parameter counting', meqn, 2*N_antennas
      Stop 'parameter counting'
   EndIf
   core_pos(1:2)=x(1:2)
   Do j_ant=1,N_antennas
      ant_x= ant_pos(1,j_ant) - core_pos(1)
      ant_y= ant_pos(2,j_ant) - core_pos(2)
      dist = sqrt(ant_x**2 + ant_y**2)
      If( dist .gt. D_max) D_max=dist
   EndDo
   D_radial=D_max/(N_radial-1)
   If(D_radial.lt.Drad_min) D_radial=Drad_min
   !
   Numb(:)=0
   A_x(:)=0.
   A_r(:)=0.
   A_y(:)=0.
   B_xx(:,:)=0.
   B_xr(:,:)=0.
   B_rr(:,:)=0.
   B_yr(:,:)=0.
   B_yy(:,:)=0.
   !
   Do j_ant=1,N_antennas
      ant_x= ant_pos(1,j_ant) - core_pos(1)
      ant_y= ant_pos(2,j_ant) - core_pos(2)
      dist = sqrt(ant_x**2 + ant_y**2)
      If( dist .gt. D_max) cycle
      SinTh=ant_y/dist
      CosTh=ant_x/dist
      !theta = ATAN2(ant_y, ant_x)  ! sin(theta)=ant_y/dist ; cos(theta)=ant_x/dist
      j_d = dist/D_radial
      del_j = dist/D_radial - j_d
      Numb(j_d)  = Numb(j_d) + 1  ! should be at least 3 for a possibly sensible result
      Numb(j_d+1)= Numb(j_d+1) + 1  ! should be at least 3 for a possibly sensible result
      Rx = PulseTime_x(j_ant)/Time_errx(j_ant)
      Ry = PulseTime_y(j_ant)/Time_erry(j_ant)
      A_x(j_d)  =A_x(j_d)   + (1.-del_j)*Rx
      A_x(j_d+1)=A_x(j_d+1) +    del_j  *Rx
      A_y(j_d)  =A_y(j_d)   + (1.-del_j)*Ry
      A_y(j_d+1)=A_y(j_d+1) +    del_j  *Ry
      A_r(j_d)  =A_r(j_d)   + (1.-del_j)*(Rx * CosTh + Ry * Sinth)
      A_r(j_d+1)=A_r(j_d+1) +    del_j  *(Rx * CosTh + Ry * Sinth)
      B_xx(j_d,j_d)  =B_xx(j_d,j_d) + (1.-del_j)*(1.-del_j)/Time_errx(j_ant)
      B_xx(j_d,j_d+1)=B_xx(j_d,j_d+1) + (1.-del_j)*  del_j /Time_errx(j_ant)
      B_xx(j_d+1,j_d+1)=B_xx(j_d+1,j_d+1) + del_j *  del_j /Time_errx(j_ant)
      B_xr(j_d,j_d) =B_xr(j_d,j_d) + (1.-del_j)*(1.-del_j)*CosTh/Time_errx(j_ant)
      B_xr(j_d,j_d+1) =B_xr(j_d,j_d+1) + (1.-del_j)*del_j *CosTh/Time_errx(j_ant)
      B_xr(j_d+1,j_d+1) =B_xr(j_d+1,j_d+1) + del_j *del_j *CosTh/Time_errx(j_ant)
      !B_rx(j_d,j_d)  =B_rx(j_d,j_d)   + (1.-del_j)*(1.-del_j)*CosTh/Time_errx(j_ant)
      !B_rx(j_d,j_d+1)=B_rx(j_d,j_d+1) + (1.-del_j)*  del_j   *CosTh/Time_errx(j_ant)
      B_rr(j_d,j_d) =B_rr(j_d,j_d) + (1.-del_j)*(1.-del_j)*(CosTh*CosTh/Time_errx(j_ant)+SinTh*SinTh/Time_erry(j_ant))
      B_rr(j_d,j_d+1) =B_rr(j_d,j_d+1) + (1.-del_j)*del_j *(CosTh*CosTh/Time_errx(j_ant)+SinTh*SinTh/Time_erry(j_ant))
      B_rr(j_d+1,j_d+1) =B_rr(j_d+1,j_d+1) + del_j *del_j *(CosTh*CosTh/Time_errx(j_ant)+SinTh*SinTh/Time_erry(j_ant))
      B_yr(j_d,j_d) =B_yr(j_d,j_d) + (1.-del_j)*(1.-del_j)*SinTh/Time_erry(j_ant)
      B_yr(j_d,j_d+1) =B_yr(j_d,j_d+1) + (1.-del_j)*del_j *SinTh/Time_erry(j_ant)
      B_yr(j_d+1,j_d+1) =B_yr(j_d+1,j_d+1) + del_j *del_j *SinTh/Time_erry(j_ant)
      B_yy(j_d,j_d) =B_yy(j_d,j_d) + (1.-del_j)*(1.-del_j)/Time_erry(j_ant)
      B_yy(j_d,j_d+1) =B_yy(j_d,j_d+1) + (1.-del_j)*del_j /Time_erry(j_ant)
      B_yy(j_d+1,j_d+1) =B_yy(j_d+1,j_d+1) + del_j *del_j /Time_erry(j_ant)
      !B_xr(j_d+1,j_d)=B_xr(j_d,j_d+1)
      !B_yr(j_d+1,j_d)=B_yr(j_d,j_d+1)
   EndDo ! i_ant=1,N_antennas
   !  Take out entries with too few (<3) data (=antennas)
   !A_r(:)=0.
   !B_xr(:,:)=0.
   !B_yr(:,:)=0.
   i=0
   Njd(:)=0
   Numb_lim=3
   Do j_d=0,N_radial
      If(Numb(j_d) .lt.Numb_lim) cycle
      i=i+1
      Njd(i)=j_d
   EndDo ! j_d=1,N_dist
   N_points=i
   If(NFit.le.0) &
      write(2,*) 'N_points:', N_points,N_radial, D_max, D_radial
   !Flush(unit=2)
   !-----------------
   !  Construct matrices for matrix inversion A = B t
   !
   N = 3*N_points
   Allocate( Res(1:N,1), B(1:N,1:N), IPIV(1:N) )
   Res(:,1)=0.
   B(:,:)=0.
   Do i=1,N_points
      Res(i,1)=A_x(Njd(i))
      Res(i+N_points,1)=A_r(Njd(i))
      Res(i+2*N_points,1)=A_y(Njd(i))
      B(i,i)=B_xx(Njd(i), Njd(i))
      B(i,i+N_points)=B_xr(Njd(i), Njd(i))
      B(i+N_points,i+N_points)=B_rr(Njd(i), Njd(i))
      B(i+N_points,i+2*N_points)=B_yr(Njd(i), Njd(i))
      B(i+2*N_points,i+2*N_points)=B_yy(Njd(i), Njd(i))
      If(i.lt.N_points) Then
         B(i,i+1)=B_xx(Njd(i), Njd(i+1))
         B(i,i+1+N_points)=B_xr(Njd(i), Njd(i+1))
         B(i+N_points,i+1+N_points)=B_rr(Njd(i), Njd(i+1))
         B(i+N_points,i+1+2*N_points)=B_yr(Njd(i), Njd(i+1))
         B(i+2*N_points,i+1+2*N_points)=B_yy(Njd(i), Njd(i+1))
      EndIf
      If(i.gt.1) Then
         B(i,i-1+N_points)=B_xr(Njd(i-1), Njd(i))
         B(i+N_points,i-1+2*N_points)=B_yr(Njd(i-1), Njd(i))
      EndIf
   EndDo !i=1,N_points
   !Do i=1,N
   !   B(i,i)=B(i,i)+0.1
   !EndDo
   !Do i=1,N
   !   Write(2,*) Res(i,1),i, B(1:i,i)
   !EndDo
   !flush(unit=2)
   !
   !
   UPLO =  'U'  ! 'L'  ! 'U'  (Upper or lower half of A is used (L implies first index is .ge. second)
   NRHS = 1  !  The number of columns in matrix B
   !         ZHESV computes the solution to a complex system of linear equations
   !           B * X = Res,
   !         where B is an N-by-N Hermitian matrix and X and Res are N-by-NRHS matrices.
   ! Real(dp) :: Work(  !  WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   ! Integer :: LWork(N*N)  ! If LWORK = -1, then a workspace query is assumed
   ! INFO = 0: successful exit
   !       < 0: if INFO = -i, the i-th argument had an illegal value
   !       > 0: if INFO = i, D(i,i) is exactly zero.  The factorization has been completed, but the block diagonal
   !            matrix D is exactly singular, so the solution could not be computed.
   !On Linux:
   !LAPACKlib="/usr/lib/x86_64-linux-gnu/lapack/liblapack.a"
   !BLASlib="/usr/lib/x86_64-linux-gnu/blas/libblas.a"
   !export FFTLIB="-lm /home/olaf/NumLib/bin/libfftpack5.1d.a ${LAPACKlib} ${BLASlib}"  # FFT library, double precision
   !my Windows: lib is in the same directory as the FFT library.
   LWORK=-1
   Call dsysv (UPLO, N, NRHS, B, N, IPIV, Res, N, Dm, LWORK, INFO)
   LWORK=nint(Dm(1))
   Allocate (ReWork(LWORK))
   Call dsysv (UPLO, N, NRHS, B, N, IPIV, Res, N, ReWORK, LWORK, INFO)
   DeAllocate ( ReWork )
   !
   !write(2,*) 'Convergence criterium:',INFO
   !If(INFO.gt.0) Then
   !   Do i=1,N
   !      Write(2,*) i, IPIV(i), B(i,1:i)
   !   EndDo
   !EndIf
   DeAllocate( B,IPIV )
   !
   If(NFit.le.0) &
      write(2,"(A,4x,A,5x,A,9x,A,11x,A,5x,A)") ' distance', 'N','t_x', 't_r', 't_y'
   tx_fit(:)=0.
   ty_fit(:)=0.
   tr_fit(:)=0.
   j_d_prev=0
      tx_fit(j_d_prev) =Res(1,1)
      tr_fit(j_d_prev) =Res(N_points+1,1)
      ty_fit(j_d_prev) =Res(2*N_points+1,1)
   !Do j_d=0,N_radial
   Do i=1,N_points
      If(Njd(i).le. j_d_prev) cycle
      tx_fit(Njd(i)) =Res(i,1)
      tr_fit(Njd(i)) =Res(N_points+i,1)
      ty_fit(Njd(i)) =Res(2*N_points+i,1)
      Do j_d=j_d_prev+1,Njd(i)
         D=Njd(i)-j_d_prev
         tx_fit(j_d) =tx_fit(j_d_prev)*(Njd(i)-j_d)/D + Res(i,1)*(j_d-j_d_prev)/D  ! interpolate linearly
         tr_fit(j_d) =tr_fit(j_d_prev)*(Njd(i)-j_d)/D + Res(N_points+i,1)*(j_d-j_d_prev)/D  ! interpolate linearly
         ty_fit(j_d) =ty_fit(j_d_prev)*(Njd(i)-j_d)/D + Res(2*N_points+i,1)*(j_d-j_d_prev)/D  ! interpolate linearly
         j_d_prev=Njd(i)
      EndDo
      If(NFit.le.0) &
         write(2,"(I3,F7.1,I4,3G13.3)") i, Njd(i)*D_radial, Numb(Njd(i)), Res(i,1), Res(N_points+i,1), Res(2*N_points+i,1)
   EndDo ! i=1,N_points
   !  Extrapolate to have non-sero values
   Do j_d=Njd(N_points)+1,N_radial
      tx_fit(j_d) = 2*tx_fit(j_d-1) - tx_fit(j_d-2) ! interpolate linear extrapolation
      tr_fit(j_d) = 2*tr_fit(j_d-1) - tr_fit(j_d-2) ! interpolate linear extrapolation
      ty_fit(j_d) = 2*ty_fit(j_d-1) - ty_fit(j_d-2) ! interpolate linear extrapolation
   EndDo
   DeAllocate( Res )
   !
   chisq=0.
   NDoF=-N
   Do i =0,N_radial-2
      !write(2,*) 'distance bin #',i, Numb(i)
      Do j_ant=1,N_antennas
         ant_x= ant_pos(1,j_ant) - core_pos(1)
         ant_y= ant_pos(2,j_ant) - core_pos(2)
         dist = sqrt(ant_x**2 + ant_y**2)
         If( dist .gt. D_max) cycle
         SinTh=ant_y/dist
         CosTh=ant_x/dist
         theta = ATAN2(ant_y, ant_x)*180./pi  ! sin(theta)=ant_y/dist ; cos(theta)=ant_x/dist
         j_d = dist/D_radial
         If(i.eq. j_d) Then
            del_j = dist/D_radial - j_d
            !write(2,*) j_d,tr_fit(j_d),(1.-del_j), tr_fit(j_d+1),del_j
            !Flush(unit=2)
            t_r=tr_fit(j_d)*(1.-del_j) + tr_fit(j_d+1)*del_j
            t_x=tx_fit(j_d)*(1.-del_j) + tx_fit(j_d+1)*del_j + CosTh * t_r
            t_y=ty_fit(j_d)*(1.-del_j) + ty_fit(j_d+1)*del_j + SinTh * t_r
            !write(2,"(I3,2G13.3,';',3G13.3,';',3G13.3,';',13G13.3)") j_ant, dist, theta, &
            !   t_x, PulseTime_x(j_ant), Time_errx(j_ant), t_y, PulseTime_y(j_ant), Time_erry(j_ant) &
            !   ,t_r, del_j, j_d, tx_fit(j_d),(1.-del_j),tx_fit(j_d+1),del_j
            If(Numb(1) .ge.Numb_lim) Then
               chisq=chisq &
                  + (t_x- PulseTime_x(j_ant))**2/Time_errx(j_ant) + (t_y- PulseTime_y(j_ant))**2/Time_erry(j_ant)
               NDoF=NDoF+1
            EndIf
         EndIf
      EndDo
   EndDo !i =0,N_radial
   write(2,"(A,I4,A,1pG11.3,A,I4,A,1pG11.3,A,2(1pG11.3),I5)") 'NFit', NFit, &
      'chisq:', chisq,'/', NDoF , '= per deg freedom=', chisq/NDoF, ', core@',core_pos(1:2)
   !
   Do j_ant=1,N_antennas
      ant_x= ant_pos(1,j_ant) - core_pos(1)
      ant_y= ant_pos(2,j_ant) - core_pos(2)
      dist = sqrt(ant_x**2 + ant_y**2)
      !If( dist .gt. D_max) cycle
      SinTh=ant_y/dist
      CosTh=ant_x/dist
      !theta = ATAN2(ant_y, ant_x)*180./pi  ! sin(theta)=ant_y/dist ; cos(theta)=ant_x/dist
      j_d = dist/D_radial
      If(j_d .ge.  N_radial) j_d=N_radial-1   ! interpolate linear extrapolation) Then
      del_j = dist/D_radial - j_d
      !write(2,*) j_d,tr_fit(j_d),(1.-del_j), tr_fit(j_d+1),del_j
      !Flush(unit=2)
      t_r=tr_fit(j_d)*(1.-del_j) + tr_fit(j_d+1)*del_j
      t_x=tx_fit(j_d)*(1.-del_j) + tx_fit(j_d+1)*del_j + CosTh * t_r
      t_y=ty_fit(j_d)*(1.-del_j) + ty_fit(j_d+1)*del_j + SinTh * t_r
      !write(2,"(I3,2G13.3,';',3G13.3,';',3G13.3,';',13G13.3)") j_ant, dist, theta, &
      !   t_x, PulseTime_x(j_ant), Time_errx(j_ant), t_y, PulseTime_y(j_ant), Time_erry(j_ant) &
      !   ,t_r, del_j, j_d, tx_fit(j_d),(1.-del_j),tx_fit(j_d+1),del_j
      R(2*j_ant-1)=(t_x- PulseTime_x(j_ant))/sqrt( Time_errx(j_ant) )
      R(2*j_ant)  =(t_y- PulseTime_y(j_ant))/sqrt( Time_erry(j_ant) )
   EndDo
   !
   ChiSq=SUM( R(1:2*N_antennas)*R(1:2*N_antennas))
   NDoF=2*N_antennas - 2 - 3*N_points
   write(2,"(A,1pG11.3,A,I4,A,1pG11.3,A,2(1pG11.3),I5)") 'For position fitting, chisq:', chisq,'/', NDoF , &
      '= per deg freedom=', chisq/NDoF, ', core@',core_pos(1:2)
   !
End Subroutine FitPulseTimes

subroutine ufparm ( meqn, nvar, x )
!*****************************************************************************80
!! UFPARM is a user-supplied external routine.
!
!  Discussion:
!    The name of the routine, the argument list, and even whether
!       it is a function or subroutine, are left to the user.
!    NL2SOL simply passes the external reference from the calling
!       program through to the residual and jacobian routines.
!    If the user has no need for this facility, then a dummy
!       routine like this one may be used.
!
!  Modified:
!    07 February 2003
!
!  Parameters:
!    Input, integer ( kind = 4 ) MEQN, the number of functions.
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!    Input, real ( kind = 8 ) X(NVAR), the current value of the variables.
!
  implicit none
  integer ( kind = 4 ) meqn
  integer ( kind = 4 ) nvar
  real ( kind = 8 ) x(nvar)
  return
end subroutine ufparm
!----------------------------------------------
Subroutine FindPulseMax(N_samples, Shape_x, Time_max, Time_err, Val_max)
   Implicit none
   Integer, intent(in) :: N_samples
   double precision, intent(in) :: Shape_x(1:N_samples)
   double precision, intent(out) :: Time_max, Time_err, Val_max
   double precision :: SampleTime(1:N_samples)
   double precision :: Shape_x_pp(1:N_samples)
   Integer :: i_loc, i, safety
   double precision :: yp, tA, dt, B
   !
   Do i_loc=1,N_samples
      SampleTime(i_loc)=i_loc  ! *SamplingTime-tTrace_Offset_dwn  ! units [m]
   EndDo
   Call spline_cubic_set( N_samples, SampleTime, Shape_x, Shape_x_pp )
   i_loc=MaxLoc(Shape_x(:),1 )
   ! require that this does not lie near the beginning or end of the trace (defined as N_samples/10.)
   Safety=2+N_samples/10
   If(i_loc.le.Safety .or. i_loc.ge. N_samples-safety) Then
      Time_max=N_samples/2.
      Time_err=N_samples*2.
      Val_max=Shape_x(1)
   Else
      Yp=Shape_x(i_loc+1)-Shape_x(i_loc) - Shape_x_pp(i_loc+1)/6. -Shape_x_pp(i_loc)/3.    !CCorr_der(i_loc,j_corr,i_Peak)
      if(Yp.lt.0.) then
          i_loc=i_loc-1  !  shift to the point left of the real maximum
          Yp=Shape_x(i_loc+1)-Shape_x(i_loc) - Shape_x_pp(i_loc+1)/6. -Shape_x_pp(i_loc)/3. !CCorr_der(i_loc,j_corr,i_Peak)
      endif
      tA=Shape_x_pp(i_loc+1) - Shape_x_pp(i_loc)   ! 2 * A
      B=Shape_x_pp(i_loc)
      dt= ( -B - sqrt(B*B - 2.*tA * Yp) )/tA
      Time_max=i_loc+dt  ! position of the maximum
   EndIf
   Call spline_cubic_val( N_samples, SampleTime, Shape_x, Shape_x_pp, Time_max, Val_max) !, ypval, yppval )
   !write(2,*) '!FindPulseMax',i_loc, N_samples, Shape_x_pp(i_loc+1), Shape_x_pp(i_loc), Val_max
   !flush(unit=2)
   If( (Shape_x_pp(i_loc+1) + Shape_x_pp(i_loc)) .ne. 0.) then
      Time_err=-Val_max/(Shape_x_pp(i_loc+1) + Shape_x_pp(i_loc))  ! [/sample^2]
   Else
      Time_err=100.
   EndIf
   Return
End Subroutine FindPulseMax




!=====================================
!============================================
Subroutine ReadThetData(Base, DistMax)
!  - Read time traces for individual antennas where antenna function is already corrected for.
!  - Fold with impulse response to find pulse-time. Note: identical to just find position of max in time trace.
!  - Fit pulse times to reconstruct the core position.
   use constants, only : dp, pi, c_l
   Use ShowerData, only : N_antennas, N_timesamples, N_radial, core_pos, d_ant
   Use ShowerData, only : ant_pos, ant_toffset, ant_xtrace, ant_ytrace, ant_file, Sample_Offset
   Use ShowerData, only : SamplingTime, nu_min, nu_max
   Use ShowerData, only : Time_maxx, Time_errx, Val_maxx, Time_maxy, Time_erry, Val_maxy
   Use ShowerData, only : PulseTime_x, PulseTime_y, tx_fit, ty_fit, tr_fit
   Use ShowerData, only :  GroundLevel_sh, ZenithAngle_deg, ZenithAngle_rad, Azimuth_deg, EstCore_N, EstCore_E
   !use Atmosphere, only : AtmosphereInit
	implicit none
   !Integer, intent(out) :: N_antennas
   Character(len=*), intent(in) :: Base
   Real(dp), intent(inout) :: DistMax
   Integer :: nxx, i, i_t, N_samples, N, i_ant  ! , i_min, i_max, i_ant
   Character(len=10) :: Label
   !Character(len=20) :: Labela, labelb
   Real(dp) :: t_offset,  D_max, Ant_d, Zen_sh, Azi_sh ! , D_time
   Real(dp) :: a,b,c,d , ant_x, ant_y ! scratch variables
   Real(dp), allocatable :: Shape_x(:), Shape_x_pp(:), Shape_y(:), Shape_y_pp(:)
   Character(len=200) :: LineIn
   Real(dp) :: DistMin
   Integer :: j_ant, j_iterations, N_iterations
   !
!    write(4,"(3I5 )") N_ant, StParRange, GroundLevel, Zen_sh, Azi_sh, Core_N, Core_E, tTrace_dim_dwn
   write(2,*) 'Use regular radial antenna spacing.'
   N_antennas=75
   !
   !-----------------------------------------------
   !
   nxx=0
   ! To allow for antennas in a ring we distinguish
   !  j_ant running over all antenna files, and
   !  i_ant=1,N the accepted antennas.
   ! inner radius of ring is set to DistMin
   DistMin=DistMax/2.
   DistMin=0.
   write(2,*) 'Minimum distance of antennas to core=',DistMin,', max distance=',DistMax
   !
   N=0
   Do j_ant=1,N_antennas
      write(label,"(I2.2)") j_ant !Set value of label equal to j_ant
      OPEN(UNIT=4,STATUS='unknown',FILE=trim(Base)//'th_000'//trim(label)//'.csv',iostat=nxx) ! E13.4
      !write(2,*) 'reading file: "', trim(Base)//'th_000'//trim(label)//'.csv'
      If(j_ant.eq.1) Then
         If(nxx.ne.0) Then
            write(2,*) 'problems reading file:', trim(Base)//'th_00001.csv'
            STOP TRIM(Base)//'th_00001.csv is problematic'
         EndIf
         Read(4,*) label, Ant_d, N_timesamples, SamplingTime, t_offset, nu_min, nu_max, GroundLevel_sh, Zen_sh, Azi_sh
         !write(2,*) label, Ant_d, N_timesamples, SamplingTime, t_offset, nu_min, nu_max, GroundLevel_sh, Zen_sh, Azi_sh
         write(2,"(A,F7.2,A, I5,A,F6.3,A,F7.1,A)") 'Antenna separation=', Ant_d, '[m], trace length=', &
            N_timesamples,'[samples], sampling time=', SamplingTime,'[m], offset=', t_offset,'[m]'
         write(2,"(A,F7.2,A,F7.2,A, A,F6.3,A,F8.3,A,F8.3)") 'Frequency bandwidth=', nu_min,' --', nu_max,'[MHz],' &
            ,' ground level @', GroundLevel_sh,'[m], Zenith angle=', Zen_sh,'[deg], azimuth=', Azi_sh
         !write(2,*) '*** check dt_PSF,del_t', N_timesamples, SamplingTime
         ZenithAngle_deg = Zen_sh
         Azimuth_deg = Azi_sh
         ZenithAngle_rad=ZenithAngle_deg*pi/180.
         Allocate( ant_file(N_antennas), ant_pos(1:2,N_antennas), ant_toffset(N_antennas))
         !Allocate( Shape_x(N_antennas), Shape_x_pp(N_antennas))
         !Allocate( Shape_y(N_antennas), Shape_y_pp(N_antennas))
         Allocate( ant_xtrace(1:N_timesamples,N_antennas), ant_ytrace(1:N_timesamples,N_antennas) )
         Allocate( PulseTime_x(N_antennas), PulseTime_y(N_antennas) )
         !Allocate( Time_maxx(N_antennas), Time_errx(N_antennas), Val_maxx(N_antennas) )
         !Allocate( Time_maxy(N_antennas), Time_erry(N_antennas), Val_maxy(N_antennas) )
         Allocate( d_ant(N_antennas+1) )
         d_ant(:)=0.
      Else
         If(nxx.ne.0) exit
         Read(4,*) label, Ant_d
      EndIf
      If(Ant_d .lt. DistMin) Then
         Close(Unit=4)
         cycle
      ElseIf(Ant_d .le. DistMax) then
         N=N+1
      Else
         Close(Unit=4)
         exit
      EndIf
      d_ant(N)=Ant_d
      ant_pos(1,N) = ant_d
      ant_pos(2,N) = 0.
      ant_toffset(N)=-t_offset !  [m]
      PulseTime_x(N) = 0.
      PulseTime_y(N) = 0.
      Do i_t=1,N_timesamples
         Read(4,*) i, ant_xtrace(i_t,N), ant_ytrace(i_t,N)
         !write(2,*) i, ant_xtrace(i_t,i_ant), ant_ytrace(i_t,N)
         !Flush(unit=2)
      End Do ! i_t=1,N_t
      A=MaxVAL( ant_xtrace(:,N)*ant_xtrace(:,N),1 )
      write(2,*) '!E_max=', N, Ant_d, sqrt(A)
      Close(Unit=4)
   End Do ! i_d=1,N_d
   Sample_Offset= t_offset
!   ant_toffset(:) = ant_toffset(:) -0.25 +1.5* SamplingTime    ! [m]  Determined by fixing the peak in the PSF 1+1.5* SamplingTime=3
   ant_toffset(:) = ant_toffset(:) -0.15 +1.5* SamplingTime    ! [m]  Determined by fixing the peak in the PSF 1+1.5* SamplingTime=3
   write(2,*) 'Timing off-set correction=', +0.5 -1.0* SamplingTime,'[m]', +0.5 +1.0* SamplingTime,-0.25 +1.5* SamplingTime
   N_antennas=N
   d_ant(N_antennas+1)=d_ant(N_antennas)
   core_pos(1)= 0.
   core_pos(2)= 0.
   !write(2,*) '************ core position adjusted to:', core_pos
   DistMax = MaxVal( d_ant(1:N), 1 )
   !
   !stop 'ReadData'
   Return
End Subroutine ReadThetData
!=================================
Subroutine PrepareKernel(Label2, Label3, PSF_max, N_PSF, zeta_K, rs_0, phi_0)
!
   use constants, only : dp,pi
   Use ShowerData, only : ZenithAngle_deg
   implicit none
   Character(len=10), intent(out) :: Label2
   Character(len=13), intent(out) :: Label3(1:PSF_max)
   Integer, intent(in) :: PSF_max
   Integer, intent(out) :: N_PSF
   real(dp), intent(inout) :: zeta_K(1:PSF_max)
    real(dp), intent(in) :: rs_0, phi_0
!
   real(dp) :: zeta_f
    INTEGER :: i,k ,nxx
!    CHARACTER*12 REAL_C (3)
!    character*80 line
!    real(dp) :: r,z,T_o,t_r,nu
!    real(dp) :: zeta
!	 Character*80 :: lname
!	 Character*20 :: label
!
   nxx=0
   write(label2,"(I3.3)") NINT(ZenithAngle_deg*10.)
   If( (zeta_K(1)*1000. -10.) .lt. 10.) Then
      zeta_K(1)=20.
   EndIf
   !
   Do k=1,PSF_max
      zeta_f = zeta_K(k)*1000.
      If(zeta_f.lt.20) exit
      N_PSF=k
      Write(Label3(k),"(I3.3,i3.3,I2.2)") NINT(zeta_f/100.),NINT(rs_0),NINT(phi_0*10/pi) !Generates the names for the output files
   EndDo ! while (nxx.eq.0)
End Subroutine PrepareKernel
!=================================

Subroutine CorePosDistances()
    Use ShowerData, only : N_antennas, N_timesamples, N_radial, core_pos, d_ant, theta_ant, Weights_ant
    Use ShowerData, only :  GroundLevel_sh, ZenithAngle_deg, ZenithAngle_rad, Azimuth_deg, EstCore_N, EstCore_E
    use constants, only : dp, pi, c_l
    integer Tempdistance, Tempcount
    OPEN(UNIT=4,STATUS='unknown',FILE="AntennaPosInput.in") ! E13.4
    Read(4,*) N_antennas
    Allocate( d_ant(N_antennas) )
    Allocate( theta_ant(N_antennas) ) ! Degrees antanna is from positive x-axis
    Allocate (Weights_ant(N_antennas) )
    d_ant(:)=0.
    theta_ant(:) = 0.
    do j_ant = 1, N_antennas
        Read(4,*) d_ant(j_ant),Weights_ant(j_ant)
    end do
    Close(unit = 4)
    Return
End subroutine CorePosDistances
