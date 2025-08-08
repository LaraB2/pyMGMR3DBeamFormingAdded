    subroutine SetParams
!    all length and times are measured in meters where  1 ns=0.3 meter=1GHz^{-1}
    use constants, only : pi,dp,c_l
    use RFootPars, only : release
    use RFootPars, only : R_0,L_0, R_02,L_02, RL_param, AParMnmBas
    use RFootPars, only : step,N_frc,h_frc,Force,alpha_frc,N_step_max, N_FitPar_bas
    use RFootPars, only : line,  h_frcL,ForceL,alpha_frcL,D_ESmooth,N_line_max,AParMnm, AlternativeSmooth
    use RFootPars, only : stpv !,  h_frc_tr,Frc_tr,alpha_frc_tr
    use RFootPars, only : rh0,X_0,lamx,X_max,lam_tc,lam_100,XDepAlpha,J0t,J0Q, alpha_vB, u0, F_over_beta, a_ChX
    use RFootPars, only : X_02, lamx2, X_max2, Energy_sh2
    use RFootPars, only : NTo, MoliereRadius, IntegrateCurrent, PancakeIncField
    use RFootPars, only : RnrmA, RnrmB, FShift_x,FShift_y, test, SelectFh, Intensity_Weight, vDrift2, F_lim
    use RFootPars, only : nu_min, nu_max, padding, SamplingTime_dwn, StParRange
    use RFootPars, only : GroundLevel, Zen_sh, Azi_sh, Zen_B, Azi_B, Fit_StI
    use BigArrays, only : IndRefOrho, TopAtmExpon
    use BigArrays, only : ObsDist_dim, ObsDist_Step, AtmHei_step, AtmHei_dim, LamInt_dim
    use BigArrays, only : nuTrace_step, nuTrace_dim, tTrace_step, tTrace_dim_b, tTrace_dim_o, tTrace_Offset
    use eventdata, only : Voltages, Energy_sh, NoisePower, OutFileLabel
    implicit none
!    logical first
    INTEGER DATE_T(8),i
    CHARACTER*12 REAL_C (3)
!
    character*4 :: E_field_Param
    character*80 lineTXT
    character(len=80) :: OutFile
    Character(len=7) ::  txt1,txt2,txt3
    integer :: k,nxx,indx(1:N_step_max)
    character*1 :: Digit
    real(dp) :: SamplingTime,h(0:N_step_max),F(0:N_step_max),a(0:N_step_max),dfx
    real(dp) :: sin_alpha
    NAMELIST /ShPars/ OutFileLabel,test,AtmHei_dim,AtmHei_step, &
        SelectFh, lam_tc ,lam_100, XDepAlpha, IntegrateCurrent, PancakeIncField,  &
        ObsDist_dim, ObsDist_Step,  tTrace_step, lamx, u0, a_ChX, J0Q, padding, D_ESmooth, AlternativeSmooth, u0, &
        F_lim, nu_min,nu_max, SamplingTime, StParRange, Voltages, rh0,MoliereRadius, J0t, GroundLevel, X_0, X_max, &
        RnrmA, RnrmB, Zen_sh, Azi_sh, Zen_B, Azi_B, Intensity_Weight, NoisePower, Energy_sh, RL_param, R_0,L_0, Fit_StI, &
        X_02, lamx2, X_max2, R_02,L_02, Energy_sh2
    !
    Test=.false.
    AtmHei_dim=2000d0 ; AtmHei_step=10.d0 ! [m]
    Energy_sh=1.d9  ! energy in [GeV]
    Energy_sh2=-1.
    ObsDist_dim=42 ; ObsDist_Step=10.
    !    tTrace_Offset=60 ; NTo=150 ; tTrace_step=0.005d0
    tTrace_Offset=10
    tTrace_step=0.02d0 ! checked stable against decreasing tTrace_step, better not be much larger
    tTrace_dim_b=NTo+tTrace_Offset
    tTrace_dim_o=tTrace_dim_b*tTrace_dim_b/4
    nu_min=30.  ; nu_max= 80. ; SamplingTime=5 ! [ns]
    StParRange=11  ! Stokes parameter sampling
    Zen_sh=0d0 ; Azi_sh=0d0; X_0=50.0d0 ; X_max=700   ! in [g/cm^2]
    J0t=12.d0 ; Zen_B=90d0 ; Azi_B=-90d0  ! J0t=B[mu T]*0.3
    !  J0T = 14.77, zen_B=22.19 , azi_B=-90.   ! direction magnetic field at LOFAR (49.5 mu T)
    !  J0T = 16.8 , zen_B=27    , azi_B=-90.   ! direction magnetic field at GRAND (56 mu T)
    F_lim=1. ! limiting force in units of 100 keV/m applied in fitting
    !
    !pi=2._dp *asin(1._dp)  ! This is the real constant pi, now set in constants
    OutFile='MGMR3D_fit-'
!    first=.true.
    read(*,NML = ShPars)
    OPEN(UNIT=2,STATUS='unknown',FILE=TRIM(OutFile)//TRIM(OutFileLabel)//'.out')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    write(2,"(3x,5(1H-),1x,'MGMR3D_fit release of ',A22,25(1H-))") release
    CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T)
    WRITE(2,230) DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
230 FORMAT(3X,5(1H-),1x,'run on ',I2,'/',I2,'/',I4,' , started at ',&
          I2,':',I2,':',I2,'.',I3,1X,25(1H-))
    ! write(2,'("Radio Emission for (Analytically) Parametrized EAS in the presence of Atmospheric Electric Fields")')
    ! write(2,'(5x,"MGMR3D")')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    write(2,200) AtmHei_dim,AtmHei_step,AtmHei_dim*AtmHei_step/1000.
200 format('Mesh in height, number=',I4,', stepsize=',f4.1,'[m], max height=',F4.1,'[km]')
!    h0=8000 ! meter

!    read(*,*) rh0,MoliereRadius, lam_tc ,lam_100, XDepAlpha, IntegrateCurrent
!
!    read(*,*) ObsDist_dim, ObsDist_Step
    if(ObsDist_dim.eq.1) ObsDist_Step=2.*ObsDist_Step
    if(ObsDist_dim.gt.200) ObsDist_dim=200      ! Just for safety that not accidentaly much too much space is used
    write(2,202) ObsDist_dim, ObsDist_Step
202 format('Mesh in distance from core, number=',I3,', stepsize=',F5.1,'[m]')
    if(test) write(2,201) LamInt_dim
201 format('Grid for \lambda integration, Max number steps=',I4)
   If(RL_param) Then
      !If(L_0.gt.0.1  .and. R_0 .lt. 3) then
      AParMnmBas(2)='    R_0'
      AParMnmBas(3)='    L_0'
         lamx=R_0
         X_0=L_0
      !EndIf
      Write(2,*) '(L,R) parametrization used for long. sh. profile with L=', X_0,', R=',lamx
     ! Write(*,*) '(L,R) parametrization used for long. sh. profile with L=', -X_0,', R=',lamx
   Else
      Write(2,*) 'Gaisser-Hillas parametrization used for long. sh. profile with X_0=', X_0,', lambda=',lamx
   Endif
!
    SamplingTime_dwn=SamplingTime*c_l  ! units [m]
!    nuTrace_dim=tTrace_dim_o/2
!    nuTrace_step=c_l/(tTrace_dim_o*tTrace_step)   ; maxnu=c_l/(2.*tTrace_step) ! in [GHz]
!    write(2,203) tTrace_Offset, NTo, tTrace_step, Nto*tTrace_step, maxnu, nuTrace_step*1000.
!203 format('Mesh in observer-time, number before start=',I3,', total #=',I4,&
!        ', stepsize=',F5.3,'[m/c], total length of calculated signal=',F6.1,'[m/c]'/&
!        'Frequency range upto ',F5.1,'[GHz], in steps of ',f5.1,'[MHz]')
!    i_nu_ini= nu_min/(1000.*nuTrace_step)  ; i_nu_max=nu_max/(1000.*nuTrace_step)
!    if(i_nu_max .ge. nuTrace_dim) i_nu_max=nuTrace_dim
!    write(2,204) nu_min,i_nu_max*nuTrace_step*1000.
!204 format('Frequency filter from ',F5.1,'[MHz] till ',F6.0,'[MHz]')
!
    Call GetNonZeroLine(lineTXT)
    ! write(2,*) lineTXT
    read(lineTXT,*,iostat=nxx) FShift_x, FShift_y, alpha_vB
    If(Azi_B .eq. 0d0) then ! B would be pointing to the east
        if(nxx.eq.0) then
            sin_alpha = sin(alpha_vB*pi/180.d0)
        else
            sin_alpha = 1.d0
        endif
        write(2,205) alpha_vB, J0t,sin_alpha,J0t*sin_alpha
205 format('alpha_vB=',F5.1,'deg, net Lorentz force = (',F5.2,' [keV/m]) x (sin(alpha)= ',F5.2,' ) =',F6.2,' [keV/m]')
    endif
    ! write(2,*) nxx, FShift_x, FShift_y, alpha_vB
!
    write(2,NML = ShPars)
!
    Call GetNonZeroLine(lineTXT)
    read(lineTXT,*,iostat=nxx) E_field_Param !, dfx
    if(nxx.ne.0) dfx=0.d0
    ! write(*,*) E_field_Param, dfx
    !if(dfx.gt.0.01) D_ESmooth=dfx
    !    write(2,"('E-fiels smoothing parameter=',f6.2)") D_ESmooth
    !    read(*,*) E_field_Param
    If(E_field_Param .eq. 'step' .or. E_field_Param .eq. 'stpv') then
        write(2,"(' E-fields entered in steped-layer mode, smoothing parameter=',f6.2)") D_ESmooth
        Do i=1,N_step_max
            write(Digit,"(I1)") i
            AParMnm(N_FitPar_bas+i)='h_frc-'//digit
            AParMnm(N_FitPar_bas+i+N_step_max)='Force-'//digit
            AParMnm(N_FitPar_bas+i+N_step_max*2)='a_frc-'//digit
        enddo
        N_frc=N_step_max
        Do i=1,N_step_max
            read(*,"(A80)") lineTXT
            !write(*,*) 'input:',lineTXT,':end'
            read(lineTXT,*,iostat=nxx) h(i),F(i),a(i)
            if(nxx.ne.0 .or. h(i).le.0.) then
                N_frc=i-1
                exit
            endif
        enddo
        Do i=1,N_frc
            h_frc(i)=h(i)*1000.
            Force(i)=F(i)
            alpha_frc(i)=a(i)*pi/180.
        enddo
        If(E_field_Param .eq. 'stpv') then  ! true values are entered, not the change
            stpv=.true.
            write(2,*) 'Derived step E-fields from true fields:'
            call stpv2step(h_frc,Force,alpha_frc,h,F,a)
            Do i=1,N_frc
                write(2,*) h(i),F(i),a(i)*180./pi
            enddo
        else
            step=.true.
            ! Calculate the true values of the electric fields
            write(2,*) 'Derived true E-fields from E-field steps:'
            call step2stpv(h_frc,Force,alpha_frc,h,F,a)
            Do i=1,N_frc
                write(2,*) h(i),F(i),a(i)*180./pi
            enddo
        endif

!        write(2,"(1x,I1,10F7.0)" ) N_frc, h_frc
    elseIf(E_field_Param .eq. 'line') then
        write(2,*) 'E-fields entered in line mode'
        line=.true.
        Do i=0,N_line_max
            write(Digit,"(I1)") i
            AParMnm(N_FitPar_bas+i+1)='Force-'//digit
            AParMnm(N_FitPar_bas+i+2+N_line_max)='a_frc-'//digit
            AParMnm(N_FitPar_bas+i+3+2*N_line_max)='h_frcL'//digit
        enddo
        N_frc=0
        read(*,*) h_frcL(0),ForceL(0),alpha_frcL(0)
        Do i=1,N_line_max
            read(*,"(A80)") lineTXT
!            write(*,*) 'input:',lineTXT,':end'
            read(lineTXT,*,iostat=nxx) h_frcL(i),ForceL(i),alpha_frcL(i)
            if(nxx.ne.0) exit
!            write(*,*) i,h_frcL(i),ForceL(i),alpha_frcL(i)
            if(h_frcL(i).le.h_frcL(i-1)) exit
            N_frc=i
        enddo
        If(N_frc.eq.0) write(2,*) 'need at least 2 input lines in line-mode'
        Do i=0,N_frc
            h_frcL(i)=h_frcL(i)*1000.
            alpha_frcL(i)=alpha_frcL(i)*pi/180.
        enddo
!        write(2,"(1x,I1,10F7.0)" ) N_frc, h_frcL
    else
        write(2,*) 'No E-field specified'
        N_frc=0
    endif
!
    end subroutine SetParams
!
    subroutine order_height(h_frc,N_step_max,N_frc,indx)
        !use RFootPars, only : N_frc,h_frc,N_step_max ! does not work, maybe because h_frc is not 'save' only in common and thus treated like a local variable here
        use constants, only : dp
        implicit none
        real(dp), intent(inout) :: h_frc(1:N_step_max)
        integer, intent(in) :: N_frc,N_step_max
        integer, intent(out) :: indx(1:N_step_max)
        integer :: i,j,k
        real(dp) :: h0
        h0=0.
        Do j=1,N_frc
            h_frc(j)=abs(h_frc(j))
            if(h_frc(j).gt.h0 ) then
                h0=h_frc(j)
                k=j
            endif
        enddo
        indx(1)=k
        Do i=2,N_frc
            h0=0.
            Do j=1,N_frc
                if(h_frc(j).gt.h0 .and. h_frc(j).lt.h_frc(indx(i-1))) then
                    h0=h_frc(j)
                    k=j
                endif
            enddo
            indx(i)=k
            ! Kndex(k)=i
        enddo   !   Now heights are ordered, heighest first, and all are positive.
    end subroutine order_height
!
    subroutine GetNonZeroLine(lineTXT)
    implicit none
    character*80, intent(inout) :: lineTXT
    integer :: eof
        lineTXT=''
        Do while (trim(lineTXT).eq.'')
            read(*,"(A80)",iostat=eof) lineTXT
            if(eof.lt.0) exit
        enddo
    end subroutine GetNonZeroLine
!
    subroutine stpv2step(h_frc,Force,alpha_frc,h,F,a)
    use constants, only : pi,dp
    use RFootPars, only : step,stpv,N_frc,N_step_max ! note that these are true values now
    !use RFootPars, only : h_frc,Force,alpha_frc ! does not work because these are local in the calling routine
    implicit none
    real(dp), intent(out) :: h(0:N_step_max),F(0:N_step_max),a(0:N_step_max)
    real(dp), intent(inout) :: h_frc(1:N_step_max)
    real(dp), intent(in) :: Force(1:N_step_max)
    real(dp), intent(in) :: alpha_frc(1:N_step_max)
    real(dp) :: dfx,dfy,F_s(1:N_step_max),a_s(1:N_step_max)
    integer :: i,indx(1:N_step_max)
        ! convert true values that are entered to the change
        if(step) write(2,*) 'use the stpv option when calling this routine'
        call order_height(h_frc,N_step_max,N_frc,indx)
        Do i=1,N_frc
            h(i)=h_frc(indx(i))
            F(i)=Force(indx(i))
            a(i)=alpha_frc(indx(i))
        enddo
        F(0)=0.d0
        a(0)=0.d0
        Do i=1,N_frc
            dfx=F(i)*cos(a(i)) - F(i-1)*cos(a(i-1))
            dfy=F(i)*sin(a(i)) - F(i-1)*sin(a(i-1))
            F_s(i)=sqrt(dfx*dfx+dfy*dfy)
            a_s(i)=atan2(dfy,dfx)
        enddo
        Do i=1,N_frc
            F(i) = F_s(i)
            a(i) = a_s(i)
        enddo
    end subroutine stpv2step
!
    subroutine step2stpv(h_frc,Force,alpha_frc,h,F,a)
    use constants, only : pi,dp
    use RFootPars, only : step,stpv,N_frc,N_step_max
    implicit none
    real(dp), intent(out) :: h(0:N_step_max),F(0:N_step_max),a(0:N_step_max)
    real(dp), intent(inout) :: h_frc(1:N_step_max)
    real(dp), intent(in) :: Force(1:N_step_max)
    real(dp), intent(in) :: alpha_frc(1:N_step_max)
    real(dp) :: dfx,dfy
    integer :: i,indx(1:N_step_max)
        ! Calculate the true values of the electric fields from the changes
        if(stpv) write(2,*) 'use the step option when calling this routine'
        call order_height(h_frc,N_step_max,N_frc,indx)
        F(0)=0.d0
        a(0)=0.d0
        Do i=1,N_frc
            dfx=F(i-1)*cos(a(i-1)) + Force(indx(i))*cos(alpha_frc(indx(i)))
            dfy=F(i-1)*sin(a(i-1)) + Force(indx(i))*sin(alpha_frc(indx(i)))
            F(i)=sqrt(dfx*dfx+dfy*dfy)
            a(i)= atan2(dfy,dfx)
            h(i)=h_frc(indx(i))
        enddo
    end subroutine step2stpv
