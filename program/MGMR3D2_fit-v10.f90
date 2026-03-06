! Main program ----------------
!------------------------------
    Include 'Constants-v10.f90'
    Include 'MGMR3D2_Atmosphere.f90'
    Include 'MGMR3D2_PancakeFunct.f90'
    Include 'MGMR3D2_LateralDistr.f90'
    Include 'MGMR3D2_BA-v10.f90'
    Include 'MGMR3D2_RFootPars-v10.f90'  ! contains CrossProd; July 2020: fit parameters changed
    Include 'MGMR3D_FFT.f90'
    !Include 'C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\Imaging\LMA\LMA2019\Program/AntFuncCnst.f90' !Take this out when using the LINUX makefile
    !Include 'C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\Imaging\LMA\LMA2019\Program/AntFunct.f90' !Take this out when using the LINUX makefile
    Include '../../LOFLI/FORTRANsrc/AntFuncCnst.f90' !Take this out when using the LINUX makefile
    Include '../../LOFLI/FORTRANsrc/AntFunct.f90' !Take this out when using the LINUX makefile
    Include 'MGMR3D2_analyse-v10.f90'
    Include 'MakeGLEplots.f90'
    Include 'MGMR3D2-v10.f90'
    Include 'MGMR3D2_subr-v10.f90'
    Include 'MGMR3D2_Source-v10.f90'
    Include 'MGMR3D2_SourceKernel-v10.f90'
    include 'MGMR3D_spline.f90'
    Include 'MGMR3D2_SetParams-v10.f90'
    Include 'MGMR3D_Fit_RadioFoot-v4.f90'
    Include 'MGMR3D2_shower-v10.f90'
    !Include 'C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\Imaging\LMA\LMA2019\Program/nl2sol.f90' !Take this out when using the LINUX makefile
    Include '../../LOFLI/FORTRANsrc/nl2sol.f90' !Take this out when using the LINUX makefile
!------------------------------
	program MGMR3D2_fit
! subversion sept 18 concerns the following bugs:
!- Fitting of intensity now works even if there is no E-field specified
!- The default charge excess has been set to 0.22. This should correct for a cos(27deg) that was probably forgotten in generating fig 2 of the manual.
!- The implementation of the (vxB)xB current has been changed and is now hopefully more realistic
! Oct 2019: L,R parametrization for shower profile included
! Oct 2019: separate parameters Xb_0 and Xc_0 introduced to parametrize charge excess and drift vel.
    use eventdata, only : Fitting,ReadInput, ShPlane, NoisePower
    use eventdata, only : FileFitResult,FileShCurrent, PlotBase
    use eventdata, only : distance_antenna, phi_antenna, N_ant_max, N_ant, ZenithAngle_shower
    use eventdata, only : N_ant,St_I, St_Q, St_U, St_V, sigma_I, sigma_Q, sigma_U, sigma_V
    use eventdata, only : N_scnt, x_scnt, y_scnt, LORA, sigma_LORA, N_scnt_max
    use RFootPars, only : FitParam,N_FitPar_max,N_FitPar,AParMnm, Intensity_Weight, StParRange
    use RFootPars, only : step,stpv, line,N_frc,N_step_max,N_line_max, N_FitPar_bas, Fit_StI
    use RFootPars, only : Zen_sh, Azi_sh
    use eventdata, only : Voltages, Core_N, Core_E, Ant_N, Ant_E,Eoff,Noff, RelMx_N, RelMx_E, RelMx_U
    use eventdata, only : vBE,vBN,vBU,vvBE,vvBN,vvBU, vBxvvB, SN, SE, SU, sZS
    use RFootPars, only : release
    use constants, only : pi,dp
    use BigArrays, only : ObsDist_dim, ObsDist_Step
    use CrossProd,only : calc_alpha_vB
    use T_Analyze, only : Analyze  ! subroutine
    Use Systemcommands, only : MakeGLEplots
	 implicit none
    INTEGER DATE_T(8),i,WLength
    CHARACTER*12 REAL_C(3)
	character*80 :: lname
	character*40 :: fname, fnamex
    Character(len=10) ::  lnamex
    Character(len=250) ::  ldata
    character(len=3) :: txt2,txt3,txt4
    !character(len=1) :: txt1
    real(dp) :: ECore,NCore !,Eoff,Noff,vBE,vBN,vBU,vvBE,vvBN,vvBU
    character(Len=6)  :: date,time,dateId(100),timeId(100)
    integer ::     eventid, Oid,n, nxx
    integer :: N_inner,LBA_inner(20)
    real(dp) ::  el_arrival,x_antenna, y_antenna, Ea,Na, E_antenna, N_antenna
    Logical :: Ft=.false.
!    real(dp) :: ErrV,ErrU, a, av,ave,wa,wav,wave, phi, phi_cut, err_cut,chi2
!
!    OPEN(UNIT=2,STATUS='unknown',FILE='FitFoot.out')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T)
!    WRITE(2,230) DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
!230 FORMAT(3X,'FitFoot, run on ',I2,'/',I2,'/',I4,' , started at ',&
!          I2,':',I2,':',I2,'.',I3,1X,25(1H-))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    STOP
!    OPEN(UNIT=6,STATUS='unknown',FILE='fitt.out')
!    eventid, az_arrival, el_arrival, distance_antenna, phi_antenna, I, Q, U, V, sigma_I, sigma_Q, sigma_U, sigma_V
! where the az and el are in degrees with LOFAR / pipeline convention, and phi_antenna is in radians.
 !   read(4,*,iostat=nxx) lname
!    write(2,*) eventid,az_arrival
!
    release='v0.1, October 2023'
    write(*,*) 'starting'
    Fitting=.true.
    call SetParams
    ReadInput=.false.
    !
    If(step .or. stpv) then
        WLength=N_step_max
        N=N_FitPar_bas + 3*N_step_max
    elseif(line) then
        WLength=N_line_max+1
        N=N_FitPar_bas + 3*(N_line_max+1)
    else
        WLength=1
        N=N_FitPar_bas
    endif
    !
   CALL get_environment_variable("LIBRARY", ldata)
   !WRITE (*,*) "LIBRARY=",TRIM(ldata)
   !WRITE (2,*) "LIBRARY=",TRIM(ldata)
   !
   write(2,"('possible fit parameters:',/,10(A7,', '))") (AParMnm(i),i=1,N)
   Call GetNonZeroLine(lname)
   write(*,*) lname
   read(lname,*,iostat=nxx) FitParam
   if(FitParam(1).le.0) then
     N_FitPar=0
   else
     Do i=2,N_FitPar_max
         if(FitParam(i).lt.FitParam(i-1)) then
             N_FitPar=i-1
             exit
         endif
         !oid=mod(FitParam(i)-N_FitPar_bas+WLength-1,WLength) N_step_max
         oid=mod(FitParam(i)-N_FitPar_bas,WLength)
         !write(2,*) i,FitParam(i),N_FitPar_bas,WLength,oid, N_frc
         if(oid.gt.N_frc) then
             ! write(2,*) i,FitParam(i),N_FitPar_bas,WLength,N_frc
             write(2,*) FitParam(i),'this parameter out-of-range for force fitting',oid
             write(*,*) FitParam(i),'this parameter out-of-range for force fitting'
             stop
         endif
     enddo
     write(2,"('used fitting parameters:',/,10(A7,', '))") (AParMnm(FitParam(i)),i=1,N_FitPar)
   endif
   !
   !
   N_scnt=0
   Call GetNonZeroLine(lname)
   write(*,*) lname
   read(lname,*,iostat=nxx) fname,fnamex
   !write(*,*) fname,nxx,fnamex
   if(nxx.eq.0) then
     FileFitResult=trim(fnamex)
     !FileShCurrent=trim(fnamex)//'ShCurr'
     write(2,*) 'Pulse time traces written to file=',trim(FileFitResult)//'ttrace-ddd-ttt.csv'
   Else
     write(2,*) 'No Pulse time traces written to file=',trim(FileFitResult)//'ttrace-ddd-ttt.csv'
   endif
   write(2,*) 'Shower current structure written to file=',trim(FileShCurrent)//'.dat'
   !
   write(2,*) 'data read from file=',trim(fname)
   OPEN(UNIT=4,STATUS='old', FILE=trim(fname),err=9)
   write(2,*) 'data comparison written to file=',trim(FileFitResult)//'.dat'
   read(4,*,iostat=nxx) fnamex
   if(nxx.ne.0) Then
      write(2,*) 'Reading error on data file'
      Stop 'read error'
   EndIf
   If(fnamex.eq.'!Voltages:') then
      Voltages=.true.
      Fit_StI=.true.
      Call AntFieParGen()
      read(4,*) fnamex,oid,fname,Azi_sh,lnamex,Zen_sh
      write(2,*) trim(fnamex),oid,trim(fname),Azi_sh,trim(lnamex),Zen_sh
      call calc_alpha_vB  ! (vBE,vBN,vBU,vvBE,vvBN,vvBU)
      read(4,"(a80)") lname
      read(lname,*) lnamex,Core_E,txt2,Core_N
      write(2,*) 'Core_east=',Core_E,', Core_north=',Core_N
      N_ant=0
      Do ! i=1,10
         read(4,"(A250)",iostat=nxx) ldata
         if(nxx.ne.0) exit
         If(ldata(1:1).eq.'!') cycle
         N_ant=N_ant+1
         read(ldata,*,iostat=nxx) Ant_E(N_ant), Ant_N(N_ant), St_I(N_ant), St_Q(N_ant), St_U(N_ant), St_V(N_ant), &
          sigma_I(N_ant), sigma_Q(N_ant), sigma_U(N_ant), sigma_V(N_ant)
         If(St_I(N_ant).lt.0.) then
            write(2,*) 'data for antenna deleted (#, E, N, I)', N_ant, Ant_E(N_ant), Ant_N(N_ant), St_I(N_ant)
            N_ant=N_ant-1
            cycle
         EndIf
         if(N_ant.ge.(N_ant_max)) then
             write(2,*) 'antenna_number overflow',oid
             exit
         else
             if(Intensity_Weight) then
               If(NoisePower.gt. 0.) then
                  N=ABS(StParRange)
                  sigma_I(N_ant)= sqrt((St_I(N_ant)+NoisePower)*NoisePower*4/N)
                  sigma_Q(N_ant)= sqrt((St_I(N_ant)+NoisePower)*NoisePower*4/N)
                  sigma_U(N_ant)= sqrt((St_I(N_ant)+NoisePower)*NoisePower*4/N)
                  sigma_V(N_ant)= sqrt((St_I(N_ant)+NoisePower)*NoisePower*4/N)
                  If(N_ant.eq.1) &
                     write(2,"('noise power of ',g10.3,' used for Stokes error estimates, N=',i3,', I,sigma(1)=',2g10.3, &
                        '; @ (N,E)=',2g10.3,'[m]')") NoisePower, N, St_I(1), sigma_I(N_ant), Ant_N(N_ant), Ant_E(N_ant)
               Else
                  If(N_ant.eq.1) &
                     write(2,"('10% error bars of Stokes with noise power of ',g10.3,', I,sigma(1)=',2g10.3, &
                        '; @ (N,E)=',2g10.3,'[m]')") NoisePower, St_I(1), sigma_I(N_ant), Ant_N(N_ant), Ant_E(N_ant)
                 sigma_I(N_ant)= (0.1*St_I(N_ant)-NoisePower)  ! NoisePower is negative
                 sigma_Q(N_ant)= (0.1*St_I(N_ant)-NoisePower)  ! NoisePower is negative
                 sigma_U(N_ant)= (0.1*St_I(N_ant)-NoisePower)  ! NoisePower is negative
                 sigma_V(N_ant)= (0.1*St_I(N_ant)-NoisePower)  ! NoisePower is negative
               EndIf
             endif
         endif
         !write(2,*) N_ant,Ant_E(N_ant), Ant_N(N_ant), St_I(N_ant), sigma_I(N_ant)
      EndDo
      !Zen_sh=90.-el_arrival
      !call calc_alpha_vB
!    ZenithAngle_shower=0.2
    else
        backspace(unit=4)
        read(4,*) fnamex,oid,fname,Azi_sh,lnamex,Zen_sh
        write(2,*) trim(fnamex),oid,trim(fname),Azi_sh,trim(lnamex),Zen_sh
        call calc_alpha_vB ! (vBE,vBN,vBU,vvBE,vvBN,vvBU)
        read(4,"(a80)") lname
        read(lname,"(1x,a7)") fnamex
        ShPlane=.true.
        if(TRIM(fnamex) .eq.'core_E=') ShPlane=.false.  ! This is east in reality
        if(TRIM(fnamex) .eq.'core_x=') then  ! This is east in reality
            ShPlane=.false.
            write(2,*) 'An obsolete option, core_x is read as core_E'
        EndIf
        If(.not. ShPlane) Then
            !write(2,*) '>',lname,'<'
            read(lname,*) lnamex,ECore,txt2,NCore,txt3,Eoff,txt4,Noff
            !write(2,*) lnamex,ECore,txt2,NCore,txt3,Eoff,txt4,Noff
            read(4,"(a80)") lname
        endif
        N_ant=1
        Do ! i=1,10
            !write(2,*) N_ant
            if(ShPlane) then
                read(4,*,iostat=nxx) x_antenna, y_antenna, St_I(N_ant), St_Q(N_ant), St_U(N_ant), St_V(N_ant), &
                   sigma_I(N_ant), sigma_Q(N_ant), sigma_U(N_ant), sigma_V(N_ant)
                if(nxx.ne.0) exit
            else
                read(4,"(A250)",iostat=nxx) ldata
                if(nxx.ne.0) exit
                read(ldata,*,iostat=nxx) E_antenna, N_antenna, St_I(N_ant), St_Q(N_ant), St_U(N_ant), St_V(N_ant), &
                   sigma_I(N_ant), sigma_Q(N_ant), sigma_U(N_ant), sigma_V(N_ant)
                if(nxx.ne.0) then
                    read(ldata,*,iostat=nxx) txt2
                    !write(2,*) N_ant,txt2
                    if(txt2.eq.'!x_') then
                        exit
                    else
                        write(2,*) ldata
                        cycle
                    endif
                endif
                Ea=E_antenna-ECore-Eoff
                Na=N_antenna-NCore-Noff
                x_antenna=Ea*vBE + Na*vBN  !
                y_antenna=Ea*vvBE + Na*vvBN !
            endif
            !write(2,*) N_ant, ShPlane, x_antenna, y_antenna
            if(N_ant.ge.(N_ant_max-1)) then
                write(2,*) 'antenna_number overflow',oid
                exit
            else
                if(Intensity_Weight) then
                  If(NoisePower.gt. 0.) then
                     N=ABS(StParRange)
                     If(N_ant.eq.1) &
                        write(2,"('noise power of ',g10.3,' used for Stokes error estimates, N=',i3,', I(1)=',g10.3)") &
                        NoisePower, N, St_I(1)
                    sigma_I(N_ant)= sqrt((St_I(N_ant)+NoisePower)*NoisePower*4/N)
                    sigma_Q(N_ant)= sqrt((St_I(N_ant)+NoisePower)*NoisePower*4/N)
                    sigma_U(N_ant)= sqrt((St_I(N_ant)+NoisePower)*NoisePower*4/N)
                    sigma_V(N_ant)= sqrt((St_I(N_ant)+NoisePower)*NoisePower*4/N)
                  Else
                     If(N_ant.eq.1) &
                        write(2,"('10% error bars of Stokes with noise power of ',g10.3,', I(1)=',g10.3)") NoisePower, St_I(1)
                    sigma_I(N_ant)= (0.1*St_I(N_ant)-NoisePower)  ! NoisePower is negative
                    sigma_Q(N_ant)= (0.1*St_I(N_ant)-NoisePower)  ! NoisePower is negative
                    sigma_U(N_ant)= (0.1*St_I(N_ant)-NoisePower)  ! NoisePower is negative
                    sigma_V(N_ant)= (0.1*St_I(N_ant)-NoisePower)  ! NoisePower is negative
                  EndIf
                endif
                distance_antenna(N_ant)=sqrt(x_antenna**2 + y_antenna**2)
                phi_antenna(N_ant)=ATAN2(y_antenna, x_antenna)
!                if(sigma_I(N_ant) .lt. St_I(N_ant)/100.) then
!                    sigma_I(N_ant)=St_I(N_ant)/10.
!                    sigma_Q(N_ant)=St_Q(N_ant)/2.
!                    sigma_U(N_ant)=St_U(N_ant)/2.
!                    sigma_V(N_ant)=St_V(N_ant)/2.
!                endif
!                if(St_Q(N_ant)/St_I(N_ant).gt.-.2) write(2,*) N_ant, x_antenna, y_antenna
!                if(St_V(N_ant)/St_I(N_ant).gt. 0.4) write(2,*) N_ant, x_antenna, y_antenna
                !write(2,*) N_ant, distance_antenna(N_ant),phi_antenna(N_ant)
                if(distance_antenna(N_ant) .lt. (ObsDist_dim-20)* ObsDist_Step) then
                    N_ant=N_ant+1
                else
                    write(2,*) 'Antenna lies outside grid',distance_antenna(N_ant),'>', ObsDist_dim* ObsDist_Step
                endif
            endif
        enddo
        N_ant=N_ant-1
        !
        write(2,*) 'reading antennas is finished'
        N_scnt=0
        Do ! i=1,10
            N_scnt=N_scnt +1
            read(4,*,iostat=nxx) x_scnt(N_scnt), y_scnt(N_scnt), LORA(N_scnt), sigma_LORA(N_scnt)
            if(nxx.ne.0) exit
            if(N_scnt.ge.(N_scnt_max-1)) then
                write(2,*) 'antenna_number overflow',oid
                exit
            endif
        enddo
        write(2,*) 'reading scintillators is finished'
        N_scnt=N_scnt-1
    endif
    N_scnt=0        ! do not fit particles on the ground
    close(unit=4)
    !
    write(2,*) 'event number=',oid,' , # antennas=',N_ant,' , # scintillators=',N_scnt,', Zenith angle=',ZenithAngle_shower
    If(.not. Fitting) goto 9
!
    call Fit_RadioFoot
    !write(2,*) 'call Fit_RadioFoot done!', Fitting
    !flush(unit=2)
    !Ft=.not.Fitting
    !Call  Analyze(Ft)
    goto 8
    !
9   Fitting=.false.
    write(2,*) 'Input-data file could not be found!'
    call calc_alpha_vB !(vBE,vBN,vBU,vvBE,vvBN,vvBU)
    call MGMR3D
    !Call  Analyze(Fitting) does not work as all intermediate results are cleaned after leaving MGMR
8 continue
    Call MakeGLEplots(Submit=.true.)
    end program MGMR3D2_fit
!
!------------------------------
!    include 'C:/OlafsUtil/LSQ/nl2sol.f90'
!    include '../NumLib/LSQ/nl2sol.f90'
!------------------------------
!------------------------------
