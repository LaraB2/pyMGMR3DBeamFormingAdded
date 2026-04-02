    subroutine Fit_RadioFoot
    use eventdata, only : N_ant  ! distance_antenna, phi_antenna,
    use eventdata, only : N_scnt
    use eventdata, only : St_I, St_Q, St_U, St_V, sigma_I, sigma_Q, sigma_U, sigma_V
!    use RFootPars, only : rh0,X_0,lamx,X_max,lamy,X_may,lam_tc,lam_100,XDepAlpha,J0x,J0y,J0Q
    !use RFootPars, only : NTo, MoliereRadius
    use RFootPars, only : FitParam,N_FitPar_max,N_FitPar,APar,AParMnm, Fit_StI
    implicit none
!    integer, intent(in) :: N_ant
  integer ( kind = 4 ) :: meqn  ! Number of data points
  integer ( kind = 4 ), parameter :: v_dim = 20000 ! dim=93 + n*p + 3*n + p*(3*p+33)/2 ,n-meqn, p=nvar
  integer ( kind = 4 ) :: nvar    ! number of parameters
  integer ( kind = 4 ) iv(60+N_FitPar)
  external CompareRadioFoot  ! subroutine that compares FitFunction to Data
  external ufparm  ! Dummy external routine
  integer ( kind = 4 ) uiparm(1),i    ! Not really used
  real ( kind = 8 ) urparm(1)    ! Not really used
  real ( kind = 8 ) v(v_dim) ! dim=93 + n*p + 3*n + p*(3*p+33)/2 ,n-meqn, p=nvar
  real ( kind = 8 ) x(N_FitPar) ! parameters that are optimized
  integer ::  error
!
!  Set the initial solution estimate.
!
    Do i=1,N_FitPar
        X(i)=APar(FitParam(i))
    enddo
    nvar=N_FitPar
    Meqn=N_ant*4 + N_scnt
    If(Fit_StI) Then
      N_scnt=0
      Meqn=N_ant
    EndIf
    i=93 + Meqn*nvar + 3*Meqn + nvar*(3*nvar+33)/2
    if(v_dim.lt.i) then
        write(2,*) 'v_dim too small, should be',i
        write(*,*) 'Fit abandones, v_dim too small, should be',i
        return
    endif

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
!    iv(18)=1
! iv(nfcov).... iv(40) is the number of calls made on calcr when
!             trying to compute covariance matrices.
! iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
!             calcj) so far done (including those used for computing
!             the covariance).
! iv(ngcov).... iv(41) is the number of calls made on calcj when
!             trying to compute covariance matrices.
! iv(niter).... iv(31) is the number of iterations performed.
!
    if(N_FitPar .gt. 0) then
      If(meqn .gt. nvar) then   ! otherwise the system is underdetermined
         write(2,*) 'nl2sol diagnostics, iv 1,15,16',IV(1),IV(15),IV(16)
          !V(44)=1.d+0
          !V(40)=1.d+0
         write(2,*) 'nl2sol diagnostics, v 32,36,40,44', V(32), V(36), V(40), V(44)
        call nl2sno ( meqn, nvar, x, CompareRadioFoot, iv, v, uiparm, urparm, ufparm, error )
        !
        write(2,"('Result, chi^2/ndf=',F9.2)") 2*v(10)/(meqn-nvar)
        write(2,"('# fie calls=',i3,' and # of iterations=',i3)") iv(6),iv(31)
        Do i=1,N_FitPar
            APar(FitParam(i))=X(i)
            write(2,"(1x,A7,'=',ES10.3)") AParMnm(FitParam(i)),x(i)
        enddo
      endif
      write(2,*) 'Best fit recalculation'
      write(*,*) 'Best fit recalculation:',x(1:nvar)
      N_FitPar=0
    endif
    call CompareRadioFoot ( meqn, nvar, x, N_FitPar, v, uiparm, urparm, ufparm )
  return
  end subroutine Fit_RadioFoot

subroutine CompareRadioFoot ( meqn, nvar, x, nf, r, uiparm, urparm, ufparm )
    use constants, only : dp, pi
    use eventdata, only : distance_antenna, phi_antenna, F_max
    use eventdata, only : FileGrid, Norm_I, N_grid, d_grid, StI_max, FileFitResult
    use eventdata, only : N_ant,StI_a,StQ_a,StU_a,StV_a  ! result analytic calculation
    use eventdata, only : St_I, St_Q, St_U, St_V, sigma_I, sigma_Q, sigma_U, sigma_V
    use eventdata, only : N_scnt, x_scnt, y_scnt, LORA, sigma_LORA, S_scnt
    use eventdata, only : Voltages, Core_N, Core_E, Ant_N, Ant_E
    use eventdata, only : Eoff,Noff,vBE,vBN,vBU,vvBE,vvBN,vvBU
!    use RFootPars, only : rh0,X_0,lamx,X_max,lamy,X_may,lam_tc,lam_100,XDepAlpha,J0t,J0y,J0Q
    use RFootPars, only : RnrmA, RnrmB, X_max, F_lim, PenFacHeight  ! NTo, MoliereRadius,
    use RFootPars, only : FitParam,N_FitPar_max,N_FitPar,APar, Fit_StI
    use RFootPars, only : FShift_x,FShift_y !, D_IMax
!*****************************************************************************80
!
!  Discussion:
!    Given the value of the vector X, this routine computes the value of F(X),
!   F(X)=VoU - sin(\eta)/[cos(\eta)+ 2 cos(\phi) E_c/E_t]
!   with \eta=X(1) and E_c/E_t = X(2)
!    the vector function whose norm we are trying to minimize.
!
!  Modified:
!    24 July 2016
!
!  Author:
!    Olaf Scholten
!
!  Parameters:
!    Input, integer ( kind = 4 ) MEQN, the number of functions.
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!    Input, real ( kind = 8 ) X(NVAR), the current value of the variables.
!    Input, integer ( kind = 4 ) NF, the number of times the residual routine
!       has been called so far.
!    Output, real ( kind = 8 ) R(MEQN), the residual vector, that is, the
!       value of the functions for the given input value of the variables.
!    Input, integer ( kind = 4 ) UIPARM(*), a user array.
!    Input, real ( kind = 8 ) URPARM(*), a user array.
!    Input, external UFPARM, an external reference to a user subroutine
!       or function.
!
   implicit none
   integer ( kind = 4 ) meqn
   integer ( kind = 4 ) nvar
   integer ( kind = 4 ) nf
   real ( kind = 8 ) r(meqn),chisq_I,chisq_Q,chisq_U,chisq_V,chisq, W, chisq_scnt
   external ufparm
   integer ( kind = 4 ), intent(in) :: uiparm(1)
   real ( kind = 8 ), intent(in) :: urparm(1)
   real ( kind = 8 ) x(nvar)
   real ( kind = 8 ) :: norm,A,B, norm_scnt,dist,theta,ant_x,ant_y
   real(dp), parameter :: LORA2Part=40., LOFAR2MGMR=1.d8
   integer :: i, i_ant
   logical :: fitratio=.true.
   common /W_TrCurr/ Norm_tc,Moli_tc,s_tc,dalpha,Nrm_alpha(0:41)
   real*8 :: Norm_tc, Moli_tc, s_tc, dalpha, Nrm_alpha, p_ang
   Real(dp) :: RelAnt_N, RelAnt_E, Antx,Anty
   !
!    N_ant=meqn/4
    write(2,*)
    if(N_FitPar.gt.0) then
        Do i=1,N_FitPar
            APar(FitParam(i))=X(i)
        enddo
        write(2,200) NF, X
        write(*,200) NF, X
200     format('Fitting-try',I3,':',15(G12.6,','))
    endif
    call MGMR3D
    chisq_I=0. ; chisq_Q=0. ; chisq_U=0. ; chisq_V=0. ; chisq=0.
    norm=1.
    W=1.d0
    if(F_max .gt. F_lim) then
        W=1.d0 + 4.*(F_max/F_lim-1.d0)**2
        write(2,"('Penalty factor for large force=',F5.1,', F_max=',F5.2)") W,F_max
    endif
    W=W* PenFacHeight
    !Write(2,*) 'W,fitratio=',W,fitratio
    if(fitratio) then
      A=0.  ; B=0.
      Do i=1,N_ant
!          A=A+ StI_a(i)/sigma_I(i)
!          B=B+ St_I(i)/sigma_I(i)  ! sets Sum[R(i)]=0
          A=A+ (StI_a(i)/sigma_I(i))**2
          B=B+ St_I(i)*StI_a(i)/(sigma_I(i)**2)  ! sets d/dN[Sum(R(i)**2)]=0
      enddo
      norm=B/A
      !Write(2,*) 'B,A=',B,A
      If(Fit_StI) Then
         Do i=1,N_ant
             R(i)=(St_I(i)- norm*StI_a(i))/sigma_I(i)
             If(R(i) .lt. -20.) Then
                  Write(2,*) 'Data point changed, i=',i, St_I(i), sigma_I(i), R(i)
                  sigma_I(i)=10.*ABS(R(i))
             Endif
             !If(ABS(R(i)).gt. 3.) R(i)=SIGN(3.d0,R(i))
             R(i)=R(i)*W
             chisq_I=chisq_I + R(i)**2
         enddo
      Else
         Do i=1,N_ant
             R(4*i-3)=(St_I(i)- norm*StI_a(i))/sigma_I(i)
             R(4*i-2)=(St_Q(i)- StQ_a(i)*St_I(i)/StI_a(i) )/sigma_Q(i)
             R(4*i-1)=(St_U(i)- StU_a(i)*St_I(i)/StI_a(i) )/sigma_U(i)
             R(4*i  )=(St_V(i)- StV_a(i)*St_I(i)/StI_a(i) )/sigma_V(i)
             !If(ABS(R(4*i-3)).gt. 3.) R(4*i-3)=SIGN(3.d0,R(4*i-3))
             R(4*i-3)=R(4*i-3)*W
             chisq_I=chisq_I + R(4*i-3)**2
             chisq_Q=chisq_Q + R(4*i-2)**2
             chisq_U=chisq_U + R(4*i-1)**2
             chisq_V=chisq_V + R(4*i)**2
         enddo
      Endif
    else
      A=0.  ; B=0.
      If(Fit_StI) Then
         Do i=1,N_ant
             A=A+ (StI_a(i)/sigma_I(i))**2
             B=B+ St_I(i)*StI_a(i)/(sigma_I(i)**2)
         enddo
      Else
         Do i=1,N_ant
             A=A+ (StI_a(i)/sigma_I(i))**2+ (StQ_a(i)/sigma_Q(i))**2+ (StU_a(i)/sigma_U(i))**2+ (StV_a(i)/sigma_V(i))**2
             B=B+ St_I(i)*StI_a(i)/(sigma_I(i)**2) + &
                  St_Q(i)*StQ_a(i)/(sigma_Q(i)**2) + St_U(i)*StU_a(i)/(sigma_U(i)**2) + St_V(i)*StV_a(i)/(sigma_V(i)**2)
         enddo
      Endif
      norm=B/A
      If(Fit_StI) Then
         Do i=1,N_ant
             R(i)=W*(St_I(i)- norm*StI_a(i))/sigma_I(i)
             chisq_I=chisq_I + R(i)**2
         enddo
      Else
         Do i=1,N_ant
             R(4*i-3)=W*(St_I(i)- norm*StI_a(i))/sigma_I(i)
             R(4*i-2)=W*(St_Q(i)- norm*StQ_a(i))/sigma_Q(i)
             R(4*i-1)=W*(St_U(i)- norm*StU_a(i))/sigma_U(i)
             R(4*i  )=W*(St_V(i)- norm*StV_a(i))/sigma_V(i)
             chisq_I=chisq_I + R(4*i-3)**2
             chisq_Q=chisq_Q + R(4*i-2)**2
             chisq_U=chisq_U + R(4*i-1)**2
             chisq_V=chisq_V + R(4*i)**2
         enddo
      Endif
    endif
    Norm_I=norm
    !write(2,*) 'Norm_I:', Norm_I
    !Flush(Unit=2)
    !
    chisq_scnt=0.
    If(N_scnt.gt.0) then
        A=0. ; B=0.
        Do i=1,N_scnt
          ! A=A+ StI_a(i)/sigma_I(i)
          ! B=B+ St_I(i)/sigma_I(i)  ! sets Sum[R(i)]=0
          A=A+ (S_scnt(i)/sigma_LORA(i))**2
          B=B+ LORA(i)*S_scnt(i)/(sigma_LORA(i)**2)  ! sets d/dN[Sum(R(i)**2)]=0
        enddo
        norm_scnt=B/A
        Do i=1,N_scnt
            R(4*N_ant+i)=(LORA(i)- norm_scnt*S_scnt(i))/sigma_LORA(i)
            chisq_scnt=chisq_scnt + R(4*N_ant+i)**2
        enddo
        write(2,*) 'ScintillatorFit: ReducedChisq=',chisq_scnt/N_scnt, '; norm=',norm_scnt*LORA2Part
        OPEN(UNIT=4,STATUS='unknown',FILE=trim(FileFitResult)//'_scnt.dat' ) !'plot/FitResult.dat')
        write(4,"('!',14x,'x [m],',20x,'y [m]',18x,'p, p_calc, sigma_p')")
        Do i=1,N_scnt
          write(4,*) x_scnt(i), y_scnt(i),LORA(i), norm_scnt*S_scnt(i),sigma_LORA(i)
        enddo
        close(unit=4)
    endif
    !
    chisq=chisq_I + chisq_Q+chisq_U+chisq_V + chisq_scnt
    !write(2,*) 'RadioFit-ReducedChisq',' all=',(chisq)/meqn,'; I=',(chisq_I)/N_ant,'; Q=',(chisq_Q)/N_ant,&
    !     '; U=',(chisq_U)/N_ant,'; V=',(chisq_V)/N_ant, '; norm=',norm*LOFAR2MGMR
    write(2,"('RadRedChisq; all=',g11.5,'; I=',g11.5,'; Q=',g11.5,'; U=',g11.5,'; V=',g11.5, '; norm=',g11.5, &
      '; ReNorm=',g11.5)") &
        chisq/meqn, chisq_I/N_ant, chisq_Q/N_ant, chisq_U/N_ant,chisq_V/N_ant,norm*LOFAR2MGMR, &
        norm*LOFAR2MGMR/(Moli_tc*Moli_tc*(RnrmA-RnrmB* X_max/300))
    write(*,*) 'ReducedChisq=',(chisq)/meqn
    OPEN(UNIT=4,STATUS='unknown',FILE=trim(FileFitResult)//'.dat' ) !'plot/FitResult.dat')
    write(4,"('!CoreDist[m], phi[rad]',3x,'I',11x,'I_calc',8x,'sigma_I',7x,'Q',9x,'Q_calc',4x,'sigma_Q',3x,&
        'U',9x,'U_calc',4x,'sigma_U',3x,'V',9x,'V_calc',4x,'sigma_V',2x,'Ant_x',2x,'Ant_y',2x,'2ph_pol/pi')")
    OPEN(UNIT=14,STATUS='unknown',FILE=trim(FileFitResult)//'-n.dat' ) !'plot/FitResult.dat')
    write(14,"('!CoreDist[m], phi[rad]',3x,'I',11x,'I_calc',8x,'sigma_I',7x,'Q',9x,'Q_calc',4x,'sigma_Q',3x, &
        'U',9x,'U_calc',4x,'sigma_U',3x,'V',9x,'V_calc',4x,'sigma_V',2x,'Ant_x',2x,'Ant_y',2x,'2ph_pol/pi')")
    !write(2,*) 'Fit_RadioFoot:CompareRadioFoot',Core_N,Noff,Core_E,Eoff
    Do i_ant=1,N_ant
      If(Voltages) Then ! Unfold antenna function, only when real antennas are specified
         Antx=Ant_E(i_ant)
         Anty=Ant_N(i_ant)
         RelAnt_N=Anty -Core_N-Noff  !=(-RelAnt_x*vvBE +RelAnt_y*vBE)/c
         RelAnt_E=Antx -Core_E-Eoff ! =(RelAnt_x*vvBN - RelAnt_y*vBN)/c
         !RelAnt_U=0.
         Ant_x= (RelAnt_N*vBN+ RelAnt_E*vBE)  ! equivalent to   dist*cos(theta); projected position on sh plane along v
         Ant_y= (RelAnt_N*vvBN+ RelAnt_E*vvBE)! equivalent to   dist*sin(theta)
      Else
         theta=phi_antenna(i_ant)
         dist = distance_antenna(i_ant)
         Antx=dist*cos(theta)
         Anty=dist*sin(theta)
         ant_x=FShift_x + Antx
         ant_y=FShift_y + Anty
      EndIf
      !write(2,*) 'Fit_RadioFoot:CompareRadioFoot',theta,dist,FShift_x,FShift_y,ant_x,ant_y
      dist = sqrt(ant_x**2 + ant_y**2)
      theta = ATAN2(ant_y, ant_x)
      If(Fit_StI) Then
         If(Voltages) Then
            write(4,"(f9.3,f10.6,3e14.5,9f10.5,1x, 3F7.1)") dist, theta,St_I(i_ant), norm*StI_a(i_ant),sigma_I(i_ant) &
               , 0.0 , 0.0, 0.1 &
               , 0.0 , 0.0, 0.1 &
               , 0.0 , 0.0, 0.1, Antx,Anty, 0.0
            write(14,"(f9.3,f10.6,3e14.5,9g12.4,1x, 3F7.1)") dist, theta,St_I(i_ant), norm*StI_a(i_ant),sigma_I(i_ant) &
               , 0.0 , 0.0, 0.1*sigma_I(i_ant) &
               , 0.0 , 0.0, 0.1*sigma_I(i_ant) &
               , 0.0 , 0.0, 0.1*sigma_I(i_ant), Antx,Anty, 0.0
         Else
            p_ang=atan2(stQ_a(i_ant),StU_a(i_ant))/pi  ! twice the polarization angle
            write(4,"(f9.3,f10.6,3e14.5,9f10.5,1x, 2F7.1, F7.2)") dist, theta,St_I(i_ant), norm*StI_a(i_ant),sigma_I(i_ant) &
               , StQ_a(i_ant)/StI_a(i_ant) ,StQ_a(i_ant)/StI_a(i_ant), 0.1 &
               , StU_a(i_ant)/StI_a(i_ant),StU_a(i_ant)/StI_a(i_ant), 0.1 &
               , StV_a(i_ant)/StI_a(i_ant) ,StV_a(i_ant)/StI_a(i_ant), 0.1, Antx,Anty, p_ang
            write(14,"(f9.3,f10.6,3e14.5,9g12.4,1x, 2F7.1, F7.2)") dist, theta,St_I(i_ant), norm*StI_a(i_ant),sigma_I(i_ant) &
               , norm*StQ_a(i_ant) ,norm*StQ_a(i_ant), 0.1*sigma_I(i_ant) &
               , norm*StU_a(i_ant) ,norm*StU_a(i_ant), 0.1*sigma_I(i_ant) &
               , norm*StV_a(i_ant) ,norm*StV_a(i_ant), 0.1*sigma_I(i_ant), Antx,Anty, p_ang
         EndIf
      Else
         p_ang=atan2(stQ_a(i_ant),StU_a(i_ant))/pi  ! twice the polarization angle
         write(4,"(f9.3,f10.6,3e14.5,9f10.5,1x, 2F7.1, F7.2)") dist, theta,St_I(i_ant), norm*StI_a(i_ant),sigma_I(i_ant) &
            ,St_Q(i_ant)/St_I(i_ant),StQ_a(i_ant)/StI_a(i_ant),sigma_Q(i_ant)/St_I(i_ant) &
            ,St_U(i_ant)/St_I(i_ant),StU_a(i_ant)/StI_a(i_ant),sigma_U(i_ant)/St_I(i_ant) &
            ,St_V(i_ant)/St_I(i_ant),StV_a(i_ant)/StI_a(i_ant),sigma_V(i_ant)/St_I(i_ant), Antx,Anty, p_ang
         write(14,"(f9.3,f10.6,3e14.5,9g12.4,1x, 2F7.1, F7.2)") dist, theta,St_I(i_ant), norm*StI_a(i_ant),sigma_I(i_ant) &
            ,St_Q(i_ant),norm*StQ_a(i_ant),sigma_Q(i_ant) &
            ,St_U(i_ant),norm*StU_a(i_ant),sigma_U(i_ant) &
            ,St_V(i_ant),norm*StV_a(i_ant),sigma_V(i_ant), Antx,Anty, p_ang
      EndIf
    enddo
    close(unit=4)
    close(unit=14)
      OPEN(UNIT=41,STATUS='unknown',FILE=trim(FileGrid)//'_StI.grd')       ! needed to properly normalize and make .z files for GLE-contour plots
      StI_max=StI_max*Norm_I
      Do i_ant=1,N_ant
         If(St_I(i_ant).gt.StI_max) StI_max=St_I(i_ant)
      Enddo
      StI_max=StI_max/Norm_I
      ! get position shower core on grid (Antx,y where Ant_x,y=0,0)
      If(Voltages) Then ! Unfold antenna function, only when real antennas are specified
         Anty = Core_N+Noff  !=(-RelAnt_x*vvBE +RelAnt_y*vBE)/c
         Antx = Core_E+Eoff ! =(RelAnt_x*vvBN - RelAnt_y*vBN)/c
      Else
         Antx=-FShift_x
         Anty=-FShift_y
      EndIf
      write(41,*) 2*N_grid+1, -N_grid*d_grid, N_grid*d_grid, StI_max, d_grid, Norm_I, Antx,Anty, 0.0
      close(unit=41)
    return
end subroutine CompareRadioFoot

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
