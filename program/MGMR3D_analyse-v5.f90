module T_Analyze
    use RFootPars, only : SamplingTime_dwn      ! [m]
    use constants, only : dp, pi
    integer,save :: i_nu_ini,i_nu_max, tTrace_dim_dwn, nuTrace_dim_dwn
    real(dp),save :: tTrace_Offset_dwn, maxnu_dwn, nuTrace_step_dwn
    logical, save :: first=.true.
contains
!---------------------------------------------------------
    Subroutine Analyze(Ft)
!    use BigArrays, only : Ex_nui,Ey_nui,Er_nui,filt,test
!    use BigArrays, only : Ex_to,Ey_to,Er_to
!    use BigArrays, only : i_nu_ini,i_nu_max, ObsDist_dim, tTrace_dim_o,tTrace_step,padding
    use eventdata, only : N_ant
    use RFootPars, only : N_FitPar
    use FFT, only : FFTransform_su,DAssignFFT
    implicit none
    logical, Intent(in) :: Ft ! determines if this is run while fitting or after matching to data
    character*80 :: line=''
    character*5 :: comnd
    character*30 :: OFile=''
    real(dp) :: z=0.
    integer :: eof,i,nxx
    logical :: done=.false.
    !
    !write(2,*) 'calling Analyze:',Ft
    !flush(unit=2)
    call DwnSample_E() ! convert to time-spectra
    !
    !write(2,*) 'calling Analyze: Dwnsample done',Ft
    !flush(unit=2)
    call FFTransform_su(tTrace_dim_dwn)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !write(2,*) 'calling Analyze:',Ft
    !flush(unit=2)
    !write(2,"(A,I10,o16)") 'FFTransform_su(tTrace_dim_dwn)',tTrace_dim_dwn,tTrace_dim_dwn
    If(N_ant .gt. 0) then
        call Observab_List
    endif
!
    if(N_FitPar.le.0) Then
      Do
       read(*,"(A80)",iostat=eof) line
       !write(*,*) '>',trim(line),'<'
       if(eof.lt.0) exit
       if(trim(line).eq.'') then
           if(done) exit
           cycle
       endif
       read(line,*,iostat=nxx) comnd, z, OFile
       if(nxx.ne.0) Ofile=''
       write(*,*) 'command=',comnd,' args=',z, 'OutFilePrefix=',trim(OFile),'<'
       done=.true.
       if(trim(comnd).eq. 'grid') then
         call Observab_Grid(z, OFile)
       else if(trim(comnd).eq. 'theta') then
         call Observab_theta(z, OFile)
       else if(trim(comnd).eq. 'dist') then
         call Observab_dist(z, OFile)
       else
         exit
       endif
      enddo
    EndIf
9   continue
!
    call DAssignFFT()
    first=.false.
!
    return
    end
!---------------------------------------------------------
  subroutine DwnSample_E()
    use BigArrays, only : Ex_nu_dwn,Ey_nu_dwn,Er_nu_dwn,filt
    use BigArrays, only : Ex_to,Ey_to,Er_to, t_to
    use BigArrays, only : ObsDist_dim, tTrace_dim_o, tTrace_step
    use BigArrays, only : CEx,CEy,CEr, Ex_nu, Ey_nu, Er_nu
    use constants, only : c_l
    use RFootPars, only : nu_min, nu_max, padding, test
    use FFT, only : FFTransform_su,DAssignFFT, DownSamlple
    implicit none
!    character*80 line
!    real(dp) :: z
    integer :: i,idi
    real(dp), allocatable :: E(:)
    integer :: nu_dim, FF_dim
    !
    FF_dim=tTrace_dim_o+2*padding +1
    nu_dim=FF_dim/2 ; FF_dim=nu_dim*2
    if(SamplingTime_dwn.lt.tTrace_step) then
        SamplingTime_dwn=tTrace_step
        write(2,*) 'SamplingTime_dwn set equal to tTrace_step=',tTrace_step
    endif
    tTrace_dim_dwn=FF_dim*tTrace_step/SamplingTime_dwn -1
    nuTrace_dim_dwn=tTrace_dim_dwn/2  ; tTrace_dim_dwn=2*nuTrace_dim_dwn
    nuTrace_step_dwn=c_l/(tTrace_dim_dwn*SamplingTime_dwn)   ; maxnu_dwn=c_l/(2.*SamplingTime_dwn) ! in [GHz]
    tTrace_Offset_dwn=padding*tTrace_step
    if(first) write(2,203) padding, tTrace_dim_o, tTrace_step/c_l, tTrace_dim_o*tTrace_step/c_l,&
        SamplingTime_dwn/c_l, maxnu_dwn, nuTrace_step_dwn*1000.
203 format('Mesh in observer-time, number before start=',I5,', total #=',I4,&
        ', stepsize=',F5.3,'[ns], total length of calculated signal=',F6.1,'[ns]'/&
        'Down-sampling to',F5.3,'[ns], Frequency range upto ',F5.1,'[GHz], in steps of ',f5.1,'[MHz]')
    i_nu_ini= nu_min/(1000.*nuTrace_step_dwn)  ; i_nu_max=nu_max/(1000.*nuTrace_step_dwn)
    if(i_nu_ini .le. 0) i_nu_ini=1
    if(i_nu_max .ge. nuTrace_dim_dwn) i_nu_max=nuTrace_dim_dwn
    if(first) write(2,204) i_nu_ini*nuTrace_step_dwn*1000.,i_nu_max*nuTrace_step_dwn*1000.
204 format('Frequency filter from ',F5.1,'[MHz] till ',F6.0,'[MHz]')
!
    allocate(filt(i_nu_ini:i_nu_max))
    allocate(Ex_nu_dwn(i_nu_ini:i_nu_max,ObsDist_dim), Ey_nu_dwn(i_nu_ini:i_nu_max,ObsDist_dim), &
        Er_nu_dwn(i_nu_ini:i_nu_max,ObsDist_dim))
    allocate(Ex_nu(i_nu_ini:i_nu_max), Ey_nu(i_nu_ini:i_nu_max), Er_nu(i_nu_ini:i_nu_max))
    allocate(CEx(1:tTrace_dim_dwn), CEy(1:tTrace_dim_dwn), CEr(1:tTrace_dim_dwn))
    allocate(E(1:FF_dim))
    call FFTransform_su(FF_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !write(2,"(A,I10,o16)") 'FFTransform_su(FF_dim)',FF_dim,FF_dim
    !       For down-sampling
    filt(:)=0.      ! Frequency filter
    do i = i_nu_ini, i_nu_max
      filt(i)=1.
    end do
    Do idi =1,ObsDist_dim
        Call DownSamlple(Ex_to(1,idi), E,padding,tTrace_dim_o, FF_dim, filt, Ex_nu_dwn(i_nu_ini,idi), i_nu_ini,i_nu_max)
        Call DownSamlple(Ey_to(1,idi), E,padding,tTrace_dim_o, FF_dim, filt, Ey_nu_dwn(i_nu_ini,idi), i_nu_ini,i_nu_max)
        Call DownSamlple(Er_to(1,idi), E,padding,tTrace_dim_o, FF_dim, filt, Er_nu_dwn(i_nu_ini,idi), i_nu_ini,i_nu_max)
    enddo
    !
    call DAssignFFT()
    deallocate(E)
    deallocate(t_to, Ex_to, Ey_to, Er_to)
  end subroutine DwnSample_E
!-------------------------------------------------------------------------
   Subroutine Observab_dist(Antd, Ofile)  ! not unfolding antenna function
   use BigArrays, only : ObsDist_dim, ObsDist_step
   use BigArrays, only : ObsPlsTime
   use BigArrays, only : Ex_nu_dwn,Ey_nu_dwn,Er_nu_dwn
   use BigArrays, only : CEx,CEy,CEr, Ex_nu, Ey_nu, Er_nu
   use RFootPars, only : test
   use RFootPars, only : X_max, Zen_sh
   use constants, only : c_l  ! ,ci
   implicit none
   character(len=30), intent(inout) :: OFile
   real(dp), intent(in) :: Antd
   complex(dp) :: Cx
   real(dp) :: nu !,c
   character*80 :: line
   integer :: i,idi,ith,Nth
   real(dp) :: ddd,t_o,theta,StI,StQ,StU,StV
   real(dp) :: dphi, dph_y, th_y, t0_y, dph_r, th_r, t0_r, th0_r, dth_r, th0_x, dth_x
   !
   if(trim(OFile).eq.'') OFile='plot/dist_'
   call Intpl_Enu(Antd,idi)  ! Generates the E-field time traces, not unfolding antenna function
   if(idi.ge. ObsDist_dim) return
   write(line,'(I3.3)') idint(Antd)
   write(2,*) 'Calculate observables at distance=',Antd,'[m]',', file=',trim(OFile)//'tt-'//trim(line)//'.csv'
   OPEN(UNIT=4,STATUS='unknown',FILE=trim(OFile)//'tt-'//trim(line)//'.csv')
   write(4,"('!',2x,'t[ns]',7x,'Re(E_x)',9x,'Re(E_y)',9x,'Re(E_r)',9x,'Re(E_x E_r*)',4x,'Im(E_x E_r*)')")
   ! c=1000d0 *c_l
   DO i=1,tTrace_dim_dwn
       t_o=(i*SamplingTime_dwn - tTrace_Offset_dwn)/c_l
       Cx=CEx(i)*conjg(CEr(i))
       write(4,"(f9.1,5(' , 'E13.4))") T_o,real(CEx(i)),real(CEy(i)),real(CEr(i)),real(Cx),Imag(Cx) ! ,real(Cex(i)),imag(CEx(i)),real(Cey(i)),imag(CEy(i))
   enddo
   close(unit=4)
   OPEN(UNIT=4,STATUS='unknown',FILE=trim(OFile)//'nu-'//trim(line)//'.csv')
   write(2,*) 'file written: "', trim(OFile)//'nu-'//trim(line)//'.csv', '"'
   !idi=IDINT(Antd/ObsDist_Step)
   !if(idi.lt.1) idi=1
   !ddd=Antd/ObsDist_Step-idi
   !if(idi.ge. ObsDist_dim) return
   !T_o=(ddd*ObsPlsTime(idi+1)+(1.d0-ddd)*ObsPlsTime(idi))/SamplingTime_dwn
   !dphi       =-0.5 - t_o/nuTrace_dim_dwn
   !write(4,*) 'time-shift phase:', t_o,  t_o/nuTrace_dim_dwn,dphi
   t_o = maxloc(ABS(CEx),1)
   dphi       = - t_o/nuTrace_dim_dwn
   t0_y = maxloc(ABS(CEy),1)
   dph_y       = - t0_y/nuTrace_dim_dwn
   t0_r = maxloc(ABS(CEr),1)
   dph_r       = - t0_r/nuTrace_dim_dwn
   !dph_r       = dphi
   !write(4,*) 'time-shift phase:', t_o, dphi
   write(2,*) 'time-shift phase:', t_o, t0_y, t0_r, dphi, i_nu_ini*dphi
   Write(4,*) X_max, Zen_sh, Antd, t_o, t0_y, t0_r, ' !'
   write(4,"('!  nu[GHz]',5x,'Ex(ampl)',8x,'Ey(ampl)',8x,'Er(ampl)' )")
   th0_x=1.
   dth_x=0.
   th0_r=1.
   dth_r=0.
   DO i=i_nu_ini,i_nu_max
       nu=i*nuTrace_step_dwn
       theta=atan2(imagpart(ex_nu(i)),realpart(ex_nu(i)))/pi + i*dphi +1.
       theta=modulo(theta,2.d0) ! needed because 2*pi phase is lost when taking the atan
       If((theta-th0_x) .gt. 1.5) dth_x=dth_x+2.d0
       If((theta-th0_x) .lt. -1.5) dth_x=dth_x-2.d0
       th0_x=theta
       th_y=atan2(imagpart(ey_nu(i)),realpart(ey_nu(i)))/pi + i*dph_y +3.
       th_r=atan2(imagpart(er_nu(i)),realpart(er_nu(i)))/pi + i*dph_r +3.
       th_r=modulo(th_r,2.d0)
       If((th_r-th0_r) .gt. 1.5) dth_r=dth_r+2.d0
       If((th_r-th0_r) .lt. -1.5) dth_r=dth_r-2.d0
       th0_r=th_r
       write(4,"(f9.6,3(' , 'E13.4),3F8.2)") nu, &
         abs(ex_nu(i))/nuTrace_step_dwn,abs(ey_nu(i))/nuTrace_step_dwn,abs(Er_nu(i))/nuTrace_step_dwn &
            ,theta-dth_x, th_y, th_r-dth_r
   enddo
   close(unit=4)
!
   Nth=8
   Do ith =1,Nth
     theta=2.*pi*ith/Nth
     call GetStokes(Antd,theta,StI,StQ,StU,StV,PrntStks=.true.)
   enddo  ! ith loop over angles
    return
    end
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    Subroutine Observab_theta(thetaD, OFile)   ! not unfolding antenna function
    use BigArrays, only : ObsDist_dim, ObsDist_step
    use BigArrays, only : CEx, CEy, CEr
    use BigArrays, only : Ex_nu_dwn,Ey_nu_dwn,Er_nu_dwn
    use constants, only : dp, pi !,ci
      use FFT, only : FFTransform_B
    implicit none
    real(dp), intent(in) :: thetaD
    character(len=30), intent(inout) :: OFile
!    complex(dp) :: E_nu_int1(i_nu_ini:i_nu_max),E_nu_int2(i_nu_ini:i_nu_max)
!    complex(dp) :: Cx,Cy
    character*80 :: line
    integer :: i,idi
    real(dp) :: theta,Antd,StI,StQ,StU,StV
!
    if(trim(OFile).eq.'') OFile='plot/th_'
      write(2,*) 'Calculate observables at theta=',thetaD,'degree'
      write(line,'(I3.3)') idint(thetaD)
      OPEN(UNIT=4,STATUS='unknown',FILE=trim(OFile)//trim(line)//'.csv')
      write(2,*) 'file written: "',trim(OFile)//trim(line)//'.csv','"'
      write(4,"('!  d[m]',9x,'theta',10x,'I',13x,'Q/I',13x,'U/I',13x,'V/I' )")
      theta=thetaD*pi/180.
      Do idi =1,ObsDist_dim
        Antd=idi*ObsDist_Step
        call FFTransform_B(CEX, tTrace_dim_dwn, Ex_nu_dwn(i_nu_ini,idi), nuTrace_dim_dwn,i_nu_ini,i_nu_max)
        call FFTransform_B(CEy, tTrace_dim_dwn, Ey_nu_dwn(i_nu_ini,idi), nuTrace_dim_dwn,i_nu_ini,i_nu_max)
        call FFTransform_B(CEr, tTrace_dim_dwn, Er_nu_dwn(i_nu_ini,idi), nuTrace_dim_dwn,i_nu_ini,i_nu_max)
!        if(idi.eq.5) then
!            write(2,*) 'antenna distance=',Antd
!            Do i=1,tTrace_dim,10
!            write(2,*) i*tTrace_step,'[m], ',cex(i),abs(cex(i))
!            enddo
!        endif
        call GetStokes(Antd,theta,StI,StQ,StU,StV,PrntStks=.true.)
        write(4,"(f9.4,' , ',f9.3,4(' , 'E13.4))") Antd,theta*180/pi,StI,StQ/StI,StU/StI,StV/StI  ! ,S/StI ,sqrt(StQ*StQ+Stu*Stu+Stv*StV)
      enddo  ! idi loop over distances to the core
      close(unit=4)
    return
    end
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    Subroutine Observab_Grid(d_Grid_max, OFile)
    use eventdata, only : FileGrid, Norm_I, N_grid, d_grid, StI_max
    use BigArrays, only : ObsDist_dim, ObsDist_step
    use BigArrays, only : CEx,CEy,CEr, Ex_nu, Ey_nu, Er_nu
    use BigArrays, only : Ex_nu_dwn,Ey_nu_dwn,Er_nu_dwn
    use RFootPars, only : FShift_x,FShift_y
    use eventdata, only : Core_N, Core_E
    use constants, only : dp, pi !,ci
    !  use FFT, only : FFTransform_B
    implicit none
    real(dp), intent(in) :: d_Grid_max
    character(len=30), intent(inout) :: OFile
!    complex(dp) :: E_nu_int1(i_nu_ini:i_nu_max),E_nu_int2(i_nu_ini:i_nu_max)
!    complex(dp) :: Cx,Cy
    real*8 :: theta, ddd  !d_grid,
    integer :: i_grid, iant, Nth, ith  ! N_grid,
    character*8 :: Fname
    integer :: i,idi
    real(dp) :: Antd,StI,StQ,StU,StV,p_ang
    logical :: polar
    integer :: ix_grid,iy_grid
    real(dp) :: x_grid,y_grid  !, StI_max, Nrm
    complex(dp) :: ph_shft_2,d2,ph_shft_a,da
    real(dp) :: T_diff,T_ave, ProjectFact, T_pulse, FShiftx,FShifty, CoreN, CoreE
!
  Nth=8       !  take this as an input parameter
  d_grid=25   !  take this as an input parameter
  polar=.true.
  if(d_Grid_max .lt. 0.) polar=.false.
  if(abs(d_Grid_max) .gt. 1.) d_grid=abs(d_Grid_max)
  write(2,*) 'Calculate observables on a grid distance=',d_grid,' [m], polar=',polar ! ,'Norm_I:', Norm_I
  if(trim(OFile).ne.'') FileGrid=OFile
!
  !OPEN(UNIT=4,STATUS='unknown',FILE=trim(OFile)//'power.csv')
  OPEN(UNIT=4,STATUS='unknown',FILE=trim(FileGrid)//'_Stokes.csv')
  !write(4,"('!  d[m]',9x,'theta',10x,'I',15x,'Q',15x,'U',15x,'V' )")
  if(polar) then
    write(4,"('!  d[m]',6x,'theta',10x,'I',15x,'Q',15x,'U',15x,'V' )")
!      iant = 1
      N_grid=ObsDist_Step*ObsDist_dim/d_grid -1
      Do i_grid =1,N_grid
          Antd=i_grid*d_grid
          call Intpl_Enu(Antd,idi) ! not unfolding antenna function
          if(idi.ge. ObsDist_dim) cycle
          Do ith =1,Nth
            theta=2.*pi*ith/Nth
            call GetStokes(Antd,theta,StI,StQ,StU,StV,PrntFile=FileGrid)
            !write(4,*) Antd,' , ',theta*180./pi,' , ',StI!sum(abs(CEx)*abs(CEx))+sum(abs(CEy)*abs(CEy))
            write(4,"(f9.4,' , ',f9.3,4(' , 'E13.4))") Antd,theta*180./pi,StI,StQ,StU,StV
         enddo
!          if(iant.gt.NAnt_max) exit
      enddo
  else
      write(4,"('!  x[m]',7x,'y[m]',10x,'I',15x,'Q',15x,'U',15x,'V' )")
      N_grid=(3./4.)*ObsDist_Step*ObsDist_dim/d_grid -1
      StI_max=0.
      OPEN(UNIT=41,STATUS='unknown',FILE=trim(FileGrid)//'_StI.z')       ! make .z files for GLE-contour plots
      write(41,"('! nx ',I3,' ny ',I3,' xmin ',F8.2,' xmax ',F7.2,' ymin ',F8.2,' ymax ',F7.2 )") &
         2*N_grid+1,2*N_grid+1, -N_grid*d_grid,N_grid*d_grid, -N_grid*d_grid,N_grid*d_grid
      OPEN(UNIT=42,STATUS='unknown',FILE=trim(FileGrid)//'_StQ.z')       ! make .z files for GLE-contour plots
      write(42,"('! nx ',I3,' ny ',I3,' xmin ',F8.2,' xmax ',F7.2,' ymin ',F8.2,' ymax ',F7.2 )") &
         2*N_grid+1,2*N_grid+1, -N_grid*d_grid,N_grid*d_grid, -N_grid*d_grid,N_grid*d_grid
      OPEN(UNIT=43,STATUS='unknown',FILE=trim(FileGrid)//'_StU.z')       ! make .z files for GLE-contour plots
      write(43,"('! nx ',I3,' ny ',I3,' xmin ',F8.2,' xmax ',F7.2,' ymin ',F8.2,' ymax ',F7.2 )") &
         2*N_grid+1,2*N_grid+1, -N_grid*d_grid,N_grid*d_grid, -N_grid*d_grid,N_grid*d_grid
      OPEN(UNIT=44,STATUS='unknown',FILE=trim(FileGrid)//'_StV.z')       ! make .z files for GLE-contour plots
      write(44,"('! nx ',I3,' ny ',I3,' xmin ',F8.2,' xmax ',F7.2,' ymin ',F8.2,' ymax ',F7.2 )") &
         2*N_grid+1,2*N_grid+1, -N_grid*d_grid,N_grid*d_grid, -N_grid*d_grid,N_grid*d_grid
      OPEN(UNIT=45,STATUS='unknown',FILE=trim(FileGrid)//'_Ang.z')       ! make .z files for GLE-contour plots
      write(45,"('! nx ',I3,' ny ',I3,' xmin ',F8.2,' xmax ',F7.2,' ymin ',F8.2,' ymax ',F7.2 )") &
         2*N_grid+1,2*N_grid+1, -N_grid*d_grid,N_grid*d_grid, -N_grid*d_grid,N_grid*d_grid
!      iant = 1
      Do iy_grid =-N_grid,N_grid
         y_grid=iy_grid*d_grid
         Do ix_grid =-N_grid,N_grid
            x_grid=ix_grid*d_grid

      Call ObservablesPrep(y_grid, x_grid, Antd, theta, ProjectFact, T_pulse)
      !   Shift problems because in FitResults not the true antenna positions are written
      !      but (Antd, theta) while here (y_grid, x_grid)

            if(ProjectFact.le. 0.0) Then
               StI=0.d0
               StQ=0.d0
               StU=0.d0
               StV=0.d0
               p_ang=0.
            Else
               Call GetStokes(Antd,theta,StI,StQ,StU,StV)
               StQ=StQ/StI
               StU=StU/StI
               StV=StV/StI
               p_ang=atan2(stQ,StU)/(pi)  ! twice the polarization angle
               StI=StI/ProjectFact
            EndIf
            !write(4,*) x_grid,' , ',y_grid,' , ',StI!sum(abs(CEx)*abs(CEx))+sum(abs(CEy)*abs(CEy))
            write(4,"(f8.3,' , ',f8.3,4(' , 'E13.4),f8.3)") x_grid,y_grid,StI,StQ,StU,StV,p_ang, ProjectFact
            write(41,"(E13.4,1x)",advance='no') StI
            If(StI.gt.StI_max) StI_max=StI
            write(42,"(F6.3,1x)",advance='no') StQ
            write(43,"(F6.3,1x)",advance='no') StU
            write(44,"(F6.3,1x)",advance='no') StV
            write(45,"(F6.3,1x)",advance='no') p_ang
         enddo
         write(41,"(1x)")
         write(42,"(1x)")
         write(43,"(1x)")
         write(44,"(1x)")
         write(45,"(1x)")
!          if(iant.gt.NAnt_max) exit
      enddo
      close(unit=41)
      close(unit=42)
      close(unit=43)
      close(unit=44)
      close(unit=45)
      OPEN(UNIT=41,STATUS='unknown',FILE=trim(FileGrid)//'_StI.grd')       ! make .z files for GLE-contour plots
      write(41,*) 2*N_grid+1, -N_grid*d_grid, N_grid*d_grid, StI_max, d_grid, Norm_I, 0.0
      close(unit=41)
   endif
   close(unit=4)
   return
   end Subroutine Observab_Grid
!-------------------------------------------------------------------------
    Subroutine Observab_List
    ! calculate observables as used in the fitting routine; written to file "FitResult" in Subroutine "CompareRadioFoot"
    ! Voltages=.false. : calculate electric fields in the shower plane for fake antennas in the shower plane
    ! Voltages=.true. : Calculate voltages as measured by LOFAR antennas on the ground.
    !          LOFAR antenna positions should be transported to the shower plane for obtaining analytic result.
    !          The analytic result should be changed by some travel distance ratio to convert in to the antenna E-field.
    !          The antenna Efield needs to be pulled through an angle dependent antenna function.
    use constants, only : dp, pi !,ci
    use eventdata, only : distance_antenna, phi_antenna
    use eventdata, only : N_ant,StI_a,StQ_a,StU_a,StV_a
    use eventdata, only : St_I, St_Q, St_U, St_V, sigma_I, sigma_Q, sigma_U, sigma_V  ! just for printing
    use eventdata, only : Energy_sh
    use eventdata, only : N_scnt, x_scnt, y_scnt, S_scnt
    use BigArrays, only : ObsDist_dim, ObsDist_step, PenDepth
    use BigArrays, only : CEx,CEy,CEr, Ex_nu, Ey_nu, Er_nu
    use RFootPars, only : FShift_x,FShift_y, D_IMax, StParRange
    use RFootPars, only : X_0,lamx,X_max
    use RFootPars, only : N_FitPar
    use RFootPars, only : GroundLevel, Zen_sh, Azi_sh, Zen_B, Azi_B
    use eventdata, only : FileFitResult    ! 'plot/FitResult'
    use RFootPars, only : R_0,L_0, RL_param
    use eventdata, only : Voltages, Core_N, Core_E, Ant_N, Ant_E
    use eventdata, only : Eoff,Noff,vBE,vBN,vBU,vvBE,vvBN,vvBU, vBxvvB, RelMx_N, RelMx_E, RelMx_U
    ! use CrossProd,only : calc_alpha_vB
    !  use FFT, only : FFTransform_B
    implicit none
    character(len=40) :: OFile=''
    complex(dp) :: Cx,Cy
    real(dp) :: NPart, X_rh, W_tc
    integer :: i,idi,i_ant,Sample_Peak,m(1)
    real(dp) :: theta,dist,ant_x,ant_y, T_pulse, Pow(1:tTrace_dim_dwn), Ant_North, Ant_East
    Real(dp) :: ProjectFact, RelAnt_N, RelAnt_E, HorDis_AntMx, ZenAng, AziAng, Ang_Ux
   !Real(dp) :: BE, BN, BU, SN, SE, SU
   !Real(dp) :: rad, sZS, sZB, c, sVB
    external W_tc
    !
    if(N_FitPar.le.0) then
        If(FileFitResult .ne. 'plot/FitResult') then
            Ofile=FileFitResult
            !write(2,*) 'Ofile=',Ofile, N_FitPar
        endif
    endif
    if(first) write(2,"('window size of ',I3,' used for calculating observables')") StParRange
    Do i_ant=1,N_ant
      If(Voltages) Then ! Unfold antenna function, only when real antennas are specified
         Ant_North=Ant_N(i_ant)
         Ant_East=Ant_E(i_ant)
      Else
         Ant_North=distance_antenna(i_ant)*sin(phi_antenna(i_ant))
         Ant_East=distance_antenna(i_ant)*cos(phi_antenna(i_ant))
      EndIf
      Call ObservablesPrep(Ant_North, Ant_East, dist, theta, ProjectFact, T_pulse)
      If(ProjectFact.gt.0.) then !within calculational grid
         If(StParRange .gt.0) then
            DO i=1,tTrace_dim_dwn
                Cx=CEx(i) + cos(theta)*CEr(i)
                Cy=CEy(i) + sin(theta)*CEr(i)
                Pow(i)=DBLE(Cx*conjg(Cx) + Cy*conjg(Cy))
            enddo
            m = maxloc(Pow)
            Sample_Peak=m(1)
            call GetStokes(dist,theta,StI_a(i_ant),StQ_a(i_ant),StU_a(i_ant),StV_a(i_ant),PrntFile=OFile,SP=Sample_Peak)
         else
            call GetStokes(dist,theta,StI_a(i_ant),StQ_a(i_ant),StU_a(i_ant),StV_a(i_ant),PrntFile=OFile)
         endif
         StI_a(i_ant)=StI_a(i_ant)/ProjectFact
         StQ_a(i_ant)=StQ_a(i_ant)/ProjectFact
         StU_a(i_ant)=StU_a(i_ant)/ProjectFact
         StV_a(i_ant)=StV_a(i_ant)/ProjectFact
         If(first) write(2,"('i_ant=',i3,', d=',f8.1,'[m], phi=',f8.1,'[deg], I,Q/I,V/I,U/I=',e13.4,3F9.4,&
                                    '; Data:',e13.4,3F9.4)") &
           i_ant,dist,theta*180/pi,StI_a(i_ant), &
           StQ_a(i_ant)/StI_a(i_ant),StU_a(i_ant)/StI_a(i_ant),StV_a(i_ant)/StI_a(i_ant), & !,sqrt(StQ*StQ+Stu*Stu+Stv*StV)
           St_I(i_ant), St_Q(i_ant)/St_I(i_ant), St_U(i_ant)/St_I(i_ant), St_V(i_ant)/St_I(i_ant)
      Else
         StI_a(i_ant)=0.
         StQ_a(i_ant)=0.
         StU_a(i_ant)=0.
         StV_a(i_ant)=0.
         If(first) write(2,"('i_ant=',i3,', d=',f8.1,'[m], phi=',f8.1,'[deg] outside grid')") &
           i_ant,dist,theta*180/pi
      EndIf!(ProjectFact.gt.0.)
      !
      !if(i_ant.eq.98)
   enddo
    If(N_scnt .le. 0) return
    X_rh=PenDepth(0)
    write(2,*) 'Energy=',Energy_sh
    If(RL_param) then
      NPart=Energy_sh*( 1 - R_0*(X_max-X_rh)/L_0)**(1/(R_0*R_0)) * exp((X_max-X_rh)/(L_0*R_0))
    Else
      NPart=Energy_sh*( (X_rh-X_0)/(X_max-X_0))**((X_max-X_0)/lamx) * exp((X_max-X_rh)/lamx)
    EndIf
    Do i=1,N_scnt   ! particles on the ground in scintillator counters
        ant_x=FShift_x + x_scnt(i)
        ant_y=FShift_y + y_scnt(i)
        dist = sqrt(ant_x**2 + ant_y**2)
        S_scnt(i)= W_tc(dist,ant_x)* NPart
    enddo
    return
    end
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    Subroutine ObservablesPrep(Ant_North, Ant_East, dist, theta, ProjectFact, T_pulse)
    ! When 'Voltages'=.false. then Ant_East=x-position & Ant_North=y-position,
    ! calculate observables as used in the fitting routine; written to file "FitResult" in Subroutine "CompareRadioFoot"
    ! Voltages=.false. : Calculate electric fields in the shower plane for fake antennas in the shower plane
    ! Voltages=.true. : Calculate voltages as measured by LOFAR antennas on the ground.
    !          LOFAR antenna positions should be transported to the shower plane for obtaining analytic result.
    !          The analytic result should be changed by some travel distance ratio to convert in to the antenna E-field.
    !          The antenna Efield needs to be pulled through an angle dependent antenna function.
    use constants, only : dp, pi !,ci
    use BigArrays, only : ObsDist_dim !, ObsDist_step, PenDepth
    use RFootPars, only : FShift_x,FShift_y, D_IMax !, StParRange
    use RFootPars, only : Zen_sh, Azi_sh, Zen_B, Azi_B
    use eventdata, only : Voltages, Core_N, Core_E, Eoff,Noff, RelMx_N, RelMx_E, RelMx_U
    use eventdata, only : vBE,vBN,vBU,vvBE,vvBN,vvBU, vBxvvB, SN, SE, SU, sZS, Ang_Ux
    use CrossProd,only : calc_alpha_vB
    implicit none
    Real(dp), intent(out) :: ProjectFact, T_pulse, dist, theta
    Real(dp), intent(in) :: Ant_North, Ant_East
    integer :: idi
    real(dp) :: RelAnt_N, RelAnt_E, HorDis_AntMx, ZenAng, AziAng, ant_x,ant_y !, Pow(1:tTrace_dim_dwn)

   If(Voltages) Then ! Calculate voltages as measured by LOFAR antennas on the ground
      Noff= (FShift_x*vvBE -FShift_y*vBE)/vBxvvB  ! FShift_is the shift of the core in the antenna plane
      Eoff=-(FShift_x*vvBN -FShift_y*vBN)/vBxvvB
      RelMx_N=D_IMax*SN  ! Position shower current Max in Earth coordinates, placing Radio Shower core at (0,0,0)
      RelMx_E=D_IMax*SE
      RelMx_U=D_IMax*SU
      ! The above quantities really need to be calculated for a new fit run, not really each time this routine is called.
      RelAnt_N=Ant_North-Core_N-Noff  !=(-RelAnt_x*vvBE +RelAnt_y*vBE)/c
      RelAnt_E=Ant_East-Core_E-Eoff ! =(RelAnt_x*vvBN - RelAnt_y*vBN)/c
      !RelAnt_U=0.
      Ant_x= (RelAnt_N*vBN+ RelAnt_E*vBE)  ! equivalent to   dist*cos(theta); projected position on sh plane along v
      Ant_y= (RelAnt_N*vvBN+ RelAnt_E*vvBE)! equivalent to   dist*sin(theta)
      HorDis_AntMx=(RelAnt_N-RelMx_N)**2 + (RelAnt_E-RelMx_E)**2
      ProjectFact=(HorDis_AntMx + RelMx_U**2)/(D_IMax*D_IMax) ! Opening angle difference due to antennas on horizontal surface
      HorDis_AntMx=SQRT(HorDis_AntMx)
      ZenAng=ATAN2(HorDis_AntMx,RelMx_U)
      AziAng=atan2( RelAnt_N-RelMx_N , RelAnt_E-RelMx_E )  ! \phi=0 = East
      dist = sqrt(ant_x**2 + ant_y**2)
      theta = ATAN2(ant_y, ant_x)
      !If(first) write(2,*) i_ant, 'RelAnt_E,N=',RelAnt_E,RelAnt_N,', Sh plane Ant_x,y=',Ant_x,Ant_y, dist, theta*180./pi
      !write(2,*) 'horDist IMx=',HorDis_AntMx, ProjectFact
      call Intpl_Enu(dist,idi,T_pulse, ZenAng, AziAng, theta, Ang_Ux)  ! frequency to time spectrum unfolding the antenna function
   Else  ! Calculate electric fields in the shower plane
      !theta=phi_antenna(i_ant)
      !dist = distance_antenna(i_ant)
      ant_x=FShift_x + Ant_East ! dist*cos(theta)
      ant_y=FShift_y + Ant_North !dist*sin(theta)
      ProjectFact=1.
      dist = sqrt(ant_x**2 + ant_y**2)
      theta = ATAN2(ant_y, ant_x)
      call Intpl_Enu(dist,idi,T_pulse)  ! converts frequency spectra to time spectrum for this antenna, no antenna function
   EndIf

   if(idi.ge. ObsDist_dim) ProjectFact=-1.  ! outside calculation grid

   return
   end Subroutine ObservablesPrep
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine GetStokes(dist,theta,StI,StQ,StU,StV,PrntFile,PrntStks, SP)
    use BigArrays, only : CEx,CEy,CEr
    use constants, only : dp, pi, c_l  ! ,ci
    use eventdata, only : Fitting
    use RFootPars, only : StParRange
    use eventdata, only : Voltages
    implicit none
    real(dp), intent(in) :: dist,theta
    !complex(dp), intent(in) :: CEx(1:tTrace_dim_dwn),CEy(1:tTrace_dim_dwn),CEr(1tTrace_dim_dwn)
    real(dp), intent(out) :: StI,StQ,StU,StV
    Character(len=30),optional, intent(in) :: PrntFile
    integer,optional, intent(in) :: SP
    logical,optional, intent(in) :: PrntStks
    complex(dp) :: Cx,Cy
    real(dp) :: omega,Px,Py,tx,ty,TPx,TPy,dt,c,norm=17000.d0
    integer :: i,Unt, i_min, i_max
    logical :: Prnt
    character(len=7) :: Fname
    !
    c=1000.d0*c_l
    if(present(SP))then
        i_min=SP-StParRange/2
        if(i_min.lt.1) then
            write(2,*) 'peak to close to begin of trace for Stokes evaluation',i_min
            i_min=1
        endif
        i_max=i_min+StParRange-1
        if(i_max.gt.tTrace_dim_dwn) then
            i_min=tTrace_dim_dwn
            write(2,*) 'peak to close to end of trace for Stokes evaluation, range shortened to',i_max-i_min+1
        endif
    else
        i_min=1 ; i_max=tTrace_dim_dwn
    endif
    Prnt=.false.
    if(present(PrntFile))then
        if(PrntFile.ne.'') then
            Prnt=.true.
            !write(2,*) 'PrntFile=',PrntFile
        endif
    endif
    if(Prnt)then
        Unt=11
        !write(2,*) 'PrntUnt=',PrntUnt
        tx=theta
        if(theta.lt.0.0) tx=theta+2.d0*pi
        write(Fname,"(i3.3,'-',i3.3)") nint(dist),nint(tx*180/pi)
        !write(2,*) dist,nint(tx*180/pi),', file=',trim(PrntFile)//'ttrace-'//trim(Fname)//'.csv'
        OPEN(unit=Unt,STATUS='unknown', FILE=trim(PrntFile)//'ttrace-'//trim(Fname)//'.csv')
        write(Unt,"('!',2x,'t[mus]',7x,'Re(E_x)',9x,'Im(E_x)',9x,'Re(E_y)',9x,'Im(E_y)')")
    endif
     !   get Stokes
     StI=0.    ; StQ=0.  ; StU=0.  ; StV=0.
     Tx=0.   ; Ty=0. ; TPx=0.    ;TPy=0.
     DO i=i_min,i_max
       Cx=CEx(i) + cos(theta)*CEr(i)
       Cy=CEy(i) + sin(theta)*CEr(i)
       Px=DBLE(Cx*conjg(Cx)) ;   Py=DBLE(Cy*conjg(Cy))
       StI=StI+ Px+Py
       StQ=StQ+ Px-Py
       StU=StU + DREAL(Cx*conjg(Cy) + Cy*conjg(Cx))
       StV=StV + DIMAG(Cx*conjg(Cy) - Cy*conjg(Cx))
       Tx=Tx+i*Px    ; Ty=Ty+i*Py
       TPx=TPx+Px    ; TPy=TPy+Py
       if(prnt) write(Unt,"(f9.6,4(' , 'E13.4))") (i*SamplingTime_dwn-tTrace_Offset_dwn)/c, &
         norm*DReal(Cx),norm*DImag(Cx),norm*DReal(Cy),norm*DImag(Cy)
     enddo
     if(prnt) close(unit=Unt)
     StI=StI*SamplingTime_dwn ;   StQ= StQ*SamplingTime_dwn
     StU=StU*SamplingTime_dwn ;   StV= StV*SamplingTime_dwn
     If(Voltages) Then
      StQ=1. ; StU=0. ; StV=0.
     EndIf
     omega=pi*(i_nu_ini+i_nu_max)* nuTrace_step_dwn
     dt=(Ty/TPy - Tx/TPx)*SamplingTime_dwn/c_l
     if(present(PrntStks))then
         if(.not.PrntStks) return
         if(nint(theta*180/pi).eq. 90) then
          If(.not. Fitting) &
             write(2,"(' d=',f8.1,'[m], phi=',f8.1,'[deg], I,Q/I,V/I,U/I=',e13.4,3F9.4,', delay=',F7.2,'[ns]',F9.4)") &
            dist,theta*180/pi,StI,StQ/StI,StU/StI,StV/StI,atan(StV/StU)/omega,dt !,sqrt(StQ*StQ+Stu*Stu+Stv*StV)
         else
          If(.not. Fitting) &
             write(2,"(' d=',f8.1,'[m], phi=',f8.1,'[deg], I,Q/I,V/I,U/I=',e13.4,3F9.4)") &
            dist,theta*180/pi,StI,StQ/StI,StU/StI,StV/StI !,sqrt(StQ*StQ+Stu*Stu+Stv*StV)
         endif
     endif
    end subroutine GetStokes
!-------------------------------------
    Subroutine Intpl_Enu(Antd,idi,T_pulse, ZenAng, AziAng, theta, phi_pt)
    ! converts frequency spectra to time spectrum for this antenna
    ! When Zenith angle is specified, the antenna function is unfolded, i.e. option 'Voltages' was used in the input
    use BigArrays, only : ObsPlsTime
    use BigArrays, only : ObsDist_dim, ObsDist_step
    use BigArrays, only : CEx,CEy,CEr, Ex_nu, Ey_nu, Er_nu
    use BigArrays, only : Ex_nu_dwn,Ey_nu_dwn,Er_nu_dwn
    use RFootPars, only : nu_min, nu_max  ! in MHz
    use constants, only : dp,pi,ci
   use AntFunCconst, only : Freq_min, Freq_max,J_0p,J_0t,J_1p,J_1t
   use FFT, only : FFTransform_B
   implicit none
   real(dp), intent(in) :: Antd
   integer, intent(out) :: idi
   real(dp),optional, intent(in) :: ZenAng, AziAng, theta, phi_pt
   real(dp),optional, intent(out) :: T_pulse
   real*8 :: ddd
   integer :: i, i_freq
   complex(dp) :: ph_shft_2,d2,ph_shft_a,da
   real(dp) :: T_diff,T_ave, dnu, freq, dfreq
   !real(dp) :: Eoff,Noff,vBE,vBN,vBU,vvBE,vvBN,vvBU
   !Real(dp) :: BE, BN, BU, SN, SE, SU
   !Real(dp) :: rad, sZS, sZB, c, sVB, FShift_x,FShift_y
   Complex(dp) :: Ex, Ey, Er, E_p, E_t
   Real(dp) :: sTh, cTh, sPhi, cPhi
   !
   idi=IDINT(Antd/ObsDist_Step)
   if(idi.lt.1) idi=1
   ddd=Antd/ObsDist_Step-idi
   if(idi.ge. ObsDist_dim) return
   T_diff=(ObsPlsTime(idi+1)-ObsPlsTime(idi))/SamplingTime_dwn
   T_ave= ddd*T_diff
   if(present(T_pulse)) T_pulse=T_ave + ObsPlsTime(idi)/SamplingTime_dwn
   ph_shft_2=exp(-ci*pi*t_diff*i_nu_ini/nuTrace_dim_dwn)
   d2       =exp(-ci*pi*t_diff/nuTrace_dim_dwn)
   ph_shft_a=exp(ci*pi*t_ave *i_nu_ini/nuTrace_dim_dwn)
   da       =exp(ci*pi*t_ave /nuTrace_dim_dwn)
   If(present(ZenAng) .and. present(AziAng)) Then ! Unfold antenna function, only when Incoming angles are specified
      Call AntFun_Inv(ZenAng*180./pi, AziAng*180./pi) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
      dnu=(nu_max-nu_min)/(i_nu_max-i_nu_ini)   ! [MHz] Jones matrix is stored on 1MHz grid
      sTh =sin(theta)
      cTh =cos(theta)
      sPhi=sin(phi_pt)
      cPhi=cos(phi_pt)
      Er_nu(:)=0
      Do i=i_nu_ini,i_nu_max
         freq=i*dnu
         i_freq=Int(freq)
         dfreq=freq-i_freq
        Ex =ph_shft_a*( (1.-ddd)*Ex_nu_dwn(i,idi)+ddd*Ex_nu_dwn(i,idi+1)*ph_shft_2 )
        Ey =ph_shft_a*( (1.-ddd)*Ey_nu_dwn(i,idi)+ddd*Ey_nu_dwn(i,idi+1)*ph_shft_2 )
        Er =ph_shft_a*( (1.-ddd)*Er_nu_dwn(i,idi)+ddd*Er_nu_dwn(i,idi+1)*ph_shft_2 )
          Ex=Ex + cTh*Er
          Ey=Ey + sTh*Er
        E_p=cPhi*Ex + sPhi*Ey
        E_t=sPhi*Ex - cPhi*Ey
        Ex_nu(i)=((1.-dfreq)*J_1p(i_freq) + dfreq*J_1p(i_freq+1)) *E_p + &
                  ((1.-dfreq)*J_1t(i_freq) + dfreq*J_1t(i_freq+1)) *E_t
        Ey_nu(i)=((1.-dfreq)*J_0p(i_freq) + dfreq*J_0p(i_freq+1)) *E_p + &
                  ((1.-dfreq)*J_0t(i_freq) + dfreq*J_0t(i_freq+1)) *E_t
        ph_shft_2 = ph_shft_2*d2
        ph_shft_a = ph_shft_a*da
      enddo
      call FFTransform_B(CEX, tTrace_dim_dwn, Ex_nu, nuTrace_dim_dwn,i_nu_ini,i_nu_max)
      call FFTransform_B(CEy, tTrace_dim_dwn, Ey_nu, nuTrace_dim_dwn,i_nu_ini,i_nu_max)
      !call FFTransform_B(CEr, tTrace_dim_dwn, Er_nu, nuTrace_dim_dwn,i_nu_ini,i_nu_max)
      CEr(:)=0.
   Else
      Do i=i_nu_ini,i_nu_max
        Ex_nu(i)=ph_shft_a*( (1.-ddd)*Ex_nu_dwn(i,idi)+ddd*Ex_nu_dwn(i,idi+1)*ph_shft_2 )
        Ey_nu(i)=ph_shft_a*( (1.-ddd)*Ey_nu_dwn(i,idi)+ddd*Ey_nu_dwn(i,idi+1)*ph_shft_2 )
        Er_nu(i)=ph_shft_a*( (1.-ddd)*Er_nu_dwn(i,idi)+ddd*Er_nu_dwn(i,idi+1)*ph_shft_2 )
        ph_shft_2 = ph_shft_2*d2
        ph_shft_a = ph_shft_a*da
      enddo
      call FFTransform_B(CEX, tTrace_dim_dwn, Ex_nu, nuTrace_dim_dwn,i_nu_ini,i_nu_max)
      call FFTransform_B(CEy, tTrace_dim_dwn, Ey_nu, nuTrace_dim_dwn,i_nu_ini,i_nu_max)
      call FFTransform_B(CEr, tTrace_dim_dwn, Er_nu, nuTrace_dim_dwn,i_nu_ini,i_nu_max)
   EndIf
   Return
   end subroutine Intpl_Enu
end module T_Analyze
