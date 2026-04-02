Module ShowerStruct
   use constants, only : dp
!   Integer, parameter :: AtmHei_dim=2000
!   real(dp), parameter :: AtmHei_step=10.d0 ! [m]
!   Real(dp) :: Xi(0:AtmHei_dim), PenDepth(0:AtmHei_dim+1)
   Integer, parameter :: ParDim=3
   !Real(dp) :: X_max(1:ParDim), J_max(1:ParDim)
   !Real(dp) :: J_total, J_ratio(1:ParDim)
   !Real(dp) :: R_0(1:ParDim), L_0(1:ParDim)
   !Real(dp) :: X_0(1:ParDim), Lamx(1:ParDim)
   Logical :: RL_param
   Real(dp), allocatable :: GrdWght(:,:)  ! GrdWght(0:N_tr,0:N_z)
   Real(dp), allocatable :: RIntAmpl(:)  ! GrdWght(0:N_z)
   Integer :: Nt_grd, Nz_grd
   Real(dp) :: delt_grd, delz_grd, zlow_grd
!============================================
!Contains
!---------------------------
End Module ShowerStruct
! ---------------------------------------------
!---------------------------
Subroutine ShowCurr2(NF, zlow_grd, delt_grd, N_tr, LongProfile, J_max, X_max, R_0, L_0, RL_param )
!  Angles are taken into account when calculating PenDepth(i_atm)
!   Real(dp) :: Xi(0:AtmHei_dim), PenDepth(0:AtmHei_dim)
   use constants, only : dp
   Use Atmosphere, only : AtmHei_dim, AtmHei_step
   Use Atmosphere, only : PenDepth
   Use ShowerStruct, only : ParDim
   !Use ShowerStruct, only : X_max, EnergyCR, R_0, L_0, X_0, Lamx
   Implicit none
   Real(dp), intent(out) :: LongProfile(0:N_tr), zlow_grd, delt_grd
   Integer, intent(in) :: NF, N_tr
   Real(dp), intent(in) :: X_max(1:ParDim), J_max(1:ParDim)
   Real(dp), intent(in) :: R_0(1:ParDim), L_0(1:ParDim)  ! /equivalent to/    Real(dp) :: X_0, Lamx
   Logical, intent(in) :: RL_param
   Real(dp), parameter :: R_0min=0.1, L_0min=20
   !Profile parameters on input
   Real(dp) :: R(1:3), L(1:3)  ! /equivalent to/    Real(dp) :: X_0, Lamx
   Integer :: i_tr, i_atm, i_z, i, i_sub
   Real(dp) :: NPart, X_rh, zeta
   !
   R(:)=R_0(:)
   L(:)=L_0(:)
   Do i_sub=1,ParDim
      If(J_max(i_sub).eq.0) cycle  ! do not process when particle number (~Energy) =0
      If(RL_param) then
         If(R_0(i_sub).lt.R_0min) Then
            R(i_sub) = R_0min + (R_0(i_sub)-R_0min)*(R_0(i_sub)-R_0min)
            If(L_0(i_sub).le.L_0min) Then
               L(i_sub)=L_0min + (L_0(i_sub)-L_0min)*(L_0(i_sub)-L_0min)
            EndIf
            If(NF.le.0) write(2,*) '*** EnergyCR, X_max, L_0, R_0:', J_max(i_sub), X_max(i_sub), L(i_sub), R(i_sub), &
                     ', (R,L)_input=', R_0(i_sub), L_0(i_sub)
         Else
            If(NF.le.0) write(2,*) 'EnergyCR, X_max, L_0, R_0:', J_max(i_sub), X_max(i_sub), L_0(i_sub), R_0(i_sub)
         EndIf
      Else
         If(NF.le.0) write(2,*) 'EnergyCR, X_max, X_0, Lamx:', J_max(i_sub), X_max(i_sub), R(i_sub), L(i_sub)
      EndIf
   EndDo
   ! Calculate longitudinal shower profile
   LongProfile(:) = 0.
   Do i_tr=0,N_tr
      zeta= i_tr*delt_grd + zlow_grd   !distance above ground along shower axis
      i_atm= (i_tr*delt_grd + zlow_grd)/AtmHei_step
      X_rh = PenDepth(i_atm)
      Do i_sub=1,ParDim
         If(J_max(i_sub).eq.0) cycle  ! do not process when particle number (~Energy) =0
         If(RL_param) then
            If(( 1 - R(i_sub)*(X_max(i_sub)-X_rh)/L(i_sub)).le.0) then
               NPart=0
            ElseIf( LOG10( 1 - R(i_sub)*(X_max(i_sub)-X_rh)/L(i_sub))/(R(i_sub)*R(i_sub)) .lt. -50. ) Then
               NPart=0
            ElseIf( (X_max(i_sub)-X_rh)/(L(i_sub)*R(i_sub)) .lt. -50. ) Then
               NPart=0
            Else
              NPart=J_max(i_sub)*( 1 - R(i_sub)*(X_max(i_sub)-X_rh)/L(i_sub))**(1/(R(i_sub)*R(i_sub))) &
                     * exp((X_max(i_sub)-X_rh)/(L(i_sub)*R(i_sub)))
            endif
         Else ! Non RL form where [R <=> X_0] & [L <=> lamx]
            If(X_rh.le.R(i_sub)) Then
               NPart=0
            else
               NPart=J_max(i_sub)*( (X_rh-R(i_sub))/(X_max(i_sub)-R(i_sub)))**((X_max(i_sub)-R(i_sub))/L(i_sub)) &
                     * exp((X_max(i_sub)-X_rh)/L(i_sub))
            Endif
         EndIf
         LongProfile(i_tr) = LongProfile(i_tr) + NPart
      EndDo
   EndDo
   !
   Return
End Subroutine ShowCurr2
!+++++++++++++++++++++++++++++++++++++++
Subroutine Fit_ShSruct(FitOption, R_Int_data, R_GrdWght, N_tr, N_z, del_tr, del_z, z_low, J_max, X_max, R_0, L_0, FitParams, &
         RL_Switch, qual, error, RInt_I, LongProfile)
! Fit X_max, EnergyCR, [R_0, L_0  --or-- X_0, lamx]
   use constants, only : dp
   Use ShowerStruct, only : RL_param, ParDim
   !Use ShowerStruct, only : X_max, EnergyCR, R_0, L_0, X_0, Lamx
   Use ShowerStruct, only : RIntAmpl, GrdWght, Nt_grd, Nz_grd, delt_grd, delz_grd, zlow_grd
   implicit none
!    integer, intent(in) :: N_ant
   !Integer, parameter :: ParDim=3
   Integer, intent(in) :: FitOption
   !  =1 fit PSF
   !  =2 fit with weights
   !  =3 fit with full kernel
   Real(dp), intent(in) :: R_Int_data(0:N_z)
   Real(dp), intent(in) :: R_GrdWght(0:N_tr,0:N_z)
   Integer, intent(in) :: N_tr, N_z
   Real(dp), intent(in) :: del_tr, del_z, z_low
   Real(dp), intent(inout) :: X_max(1:ParDim), J_max(1:ParDim), R_0(1:ParDim), L_0(1:ParDim)
   logical, intent(in) :: RL_Switch
   Real(dp), intent(out) :: qual
   integer, intent(out) ::  error
   Real(dp), intent(out) :: RInt_I(0:N_z), LongProfile(0:N_tr)
   Integer, intent(in) :: FitParams(12)
   !Real(dp) :: X_max(1:ParDim), J_max(1:ParDim)
   Real(dp) :: J_total, J_ratio(1:ParDim)
   !Real(dp) :: R_0(1:ParDim), L_0(1:ParDim)
   !Integer, parameter :: ParDim=3
   Integer :: N_FitPar,i,j,k,m
   integer ( kind = 4 ) :: meqn  ! Number of data points
   integer ( kind = 4 ) :: v_dim  ! dim=93 + n*p + 3*n + p*(3*p+33)/2 ,n-meqn, p=nvar
   integer ( kind = 4 ) :: nvar    ! number of parameters
   integer ( kind = 4 ) iv(60+12)
   external CompareR_Intf  ! subroutine that compares FitFunction to Data
   external ufparm  ! Dummy external routine
   integer ( kind = 4 ) uiparm(0:12)    ! FitOption, FitParams
   real ( kind = 8 ), allocatable :: urparm(:)    ! Not really used
   real ( kind = 8 ), allocatable :: v(:) ! dim=93 + n*p + 3*n + p*(3*p+33)/2 ,n-meqn, p=nvar
   real ( kind = 8 ) x(7) ! parameters that are optimized
!
!  Set the initial solution estimate.
!
   Nz_grd=N_z
   uiparm(0)=FitOption  !  Determines the kernel, PSF/full  and x/r
   !
   Allocate( RIntAmpl(0:Nz_grd) )
   RIntAmpl(0:Nz_grd) = R_Int_data(0:N_z)
   Nt_grd=N_tr
!   If( Allocated(GrdWght) ) Deallocate( GrdWght )  ! GrdWght may be needed for a separate calculation of profile
   Allocate( GrdWght(0:Nt_grd,0:Nz_grd) )
   GrdWght(0:Nt_grd,0:Nz_grd)=R_GrdWght(0:N_tr,0:N_z)
   delt_grd=del_tr ; delz_grd= del_z ; zlow_grd= z_low
   !
   !write(2,"(A,100G13.4)") '!R_Int_data(i_z):', R_Int_data(0:Nz_grd)
   !write(2,"(A,100G13.4)") '!delt_grd:', delt_grd, delz_grd, zlow_grd,Nt_grd,Nz_grd, GrdWght(Nt_grd/2,Nz_grd/2)
   !
   Meqn=Nz_grd+1
   Allocate( urparm(1:Meqn+2+Nt_grd+4*ParDim) )
   !write(2,*) '!urparam', Meqn+2+Nt_grd+4*ParDim, Meqn,Nt_grd,ParDim
   J_total=J_max(1)
   J_ratio(1)=1
   J_ratio(2:ParDim) = J_max(2:ParDim)/J_max(1)
   m=2+meqn+Nt_grd
   urparm(m) = J_total
   urparm(m+1:m+ParDim) = J_ratio(1:ParDim)
   urparm(m+ParDim+1:m+2*ParDim) = X_max(1:ParDim)
   urparm(m+2*ParDim+1:m+3*ParDim) = R_0(1)  !:ParDim)
   urparm(m+3*ParDim+1:m+4*ParDim) = L_0(1)  !:ParDim)
   RL_param = RL_Switch
   N_FitPar = 0
   m=0
   Do i=1,12
      j=FitParams(i)
      !write(2,*) '!FitParams1', i,j,m
      !If(j.le.m) exit  !  do not supply the same parameter twice and exit for negative values
      If(j.le.0) exit  !   exit for negative values
      m=0
      Do k=1,i-1
         If( FitParams(k) .eq. j) m=-1
      EndDo
      If(m.ne. 0) exit
      k=(j-1)/ParDim
      m=j-ParDim*k
      uiparm(i)=FitParams(i)
      !write(2,*) '!FitParams2', i,j,k,m
      Select Case(k)
         Case(0)
            X(i)=J_ratio(m)
            Write(2,"(A,I2,A,I1,A,1pG13.4)") ' -- Fitting parameter ',j,', J_ratio(',m,') =',X(i)
            If(m.eq.1) Then
               write(2,*) '*** This parameter should not be fitted!!'
               stop 'Not a fit parameter, J_ratio(1)'
            EndIf
         Case(1)
            X(i)=X_max(m)
            Write(2,"(A,I2,A,I1,A,1pG13.4)") ' -- Fitting parameter ',j,', X_max(',m,') =',X(i)
         Case(2)
            X(i)=R_0(m)
            Write(2,"(A,I2,A,I1,A,1pG13.4)") ' -- Fitting parameter ',j,', X0 or R(',m,') =',X(i)
            If(m.gt.1) Then
               write(2,*) '*** This parameter should not be fitted!!'
               stop 'Not a fit parameter, R(2:3)'
            EndIf
         Case(3)
            X(i)=L_0(m)
            Write(2,"(A,I2,A,I1,A,1pG13.4)") ' -- Fitting parameter ',j,', Lam_x or L(',m,') =',X(i)
            If(m.gt.1) Then
               write(2,*) '*** This parameter should not be fitted!!'
               stop 'Not a fit parameter, L(2:3)'
            EndIf
         Case Default
            !If(j.gt.ParDim*4) exit
            Exit
      End Select
      N_FitPar=i
      m=j
   EndDo
   !
   nvar=N_FitPar
   v_dim=93 + Meqn*nvar + 3*Meqn + nvar*(3*nvar+33)/2
   !write(2,*) 'Fitting; N_z, N_tr, Meqn, N_FitPar, v_dim are:',N_z, N_tr, Meqn, N_FitPar, v_dim
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
         !write(2,*) 'nl2sol diagnostics, iv 1,15,16',IV(1),IV(15),IV(16)
          !V(44)=1.d+0
          !V(40)=1.d+0
         !write(2,*) 'nl2sol diagnostics, v 32,36,40,44', V(32), V(36), V(40), V(44)
        call nl2sno ( meqn, nvar, x, CompareR_Intf, iv, v, uiparm, urparm, ufparm, error )
        !
        write(2,"('Result, chi^2/ndf=',g15.5)") 2*v(10)/(meqn-nvar)
        write(2,"('# fie calls=',i3,' and # of iterations=',i3)") iv(6),iv(31)
      endif
      write(2,*) 'Best fit recalculation'
      write(*,"(10G13.5)") 'Best fit:',x(1:nvar)
      m=2+meqn+Nt_grd
      J_total=urparm(m)
      J_ratio(1:ParDim) = urparm(m+1:m+ParDim)
      X_max(1:ParDim)   = urparm(m+ParDim+1:m+2*ParDim)
      R_0(1:ParDim)     = urparm(m+2*ParDim+1)  !:m+3*ParDim)
      L_0(1:ParDim)     = urparm(m+3*ParDim+1)  !:m+4*ParDim)
      J_max(1:ParDim)   = J_total * J_ratio(1:ParDim)
    endif
    i=-1    ! Renenerate best fitting result
    call CompareR_Intf ( meqn, nvar, x, i, v, uiparm, urparm, ufparm )
    RInt_I(0:Nz_grd) = urparm(1:meqn)
    LongProfile(0:Nt_grd) =urparm(1+meqn:1+meqn+Nt_grd)
    qual=(meqn/(meqn-nvar))*100.*sqrt( sum(V(1:Meqn)*V(1:Meqn))/sum(RInt_I(0:Nz_grd)*RInt_I(0:Nz_grd)) )
   !write(2,"(A,100G13.4)") '!LongProfile(i_z):', LongProfile(0:Nt_grd)
   write(2,*) 'chi^2/df=', qual/(meqn-nvar)
   DeAllocate( RIntAmpl, GrdWght, urparm, v )
   return
end subroutine Fit_ShSruct
!
subroutine CompareR_Intf( meqn, nvar, x, nf, r, uiparm, urparm, ufparm )
    use constants, only : dp, pi
   Use ShowerStruct, only : RIntAmpl, Nt_grd, Nz_grd, delt_grd, delz_grd, zlow_grd, GrdWght
   !Use ShowerStruct, only : X_max, J_max, R_0, L_0 !, X_0, Lamx
   Use ShowerStruct, only : RL_param, ParDim !, ShowCurr
   Use PSF, only : hcnt_l, hcnt_u, Grd_Kx , Grd_Kr
   !Use PSF, only : GrdPSFWghtx, GrdPwrWghtx, GrdPSFWghtr, GrdPwrWghtr
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
   integer ( kind = 4 ), intent(in) :: meqn
   integer ( kind = 4 ), intent(in) :: nvar
   real ( kind = 8 ), intent(in) :: X(*) !  x(nvar)
   integer ( kind = 4 ), intent(in) :: nf
   real ( kind = 8 ), intent(out) :: r(meqn)
   integer ( kind = 4 ), intent(in) :: uiparm(0:12)
   real ( kind = 8 ), intent(inout) :: urparm(*)
   Real(dp) :: X_max(1:ParDim), J_max(1:ParDim)
   Real(dp) :: J_total, J_ratio(1:ParDim)
   Real(dp) :: R_0(1:ParDim), L_0(1:ParDim)
   external ufparm
   !
   real(dp) :: W, ChiSq
   integer :: i,j,k,m, i_z, i_h
   logical :: fitratio=.true.
   !common /W_TrCurr/ Norm_tc,Moli_tc,s_tc,dalpha,Nrm_alpha(0:41)
   !real*8 :: Norm_tc, Moli_tc, s_tc, dalpha, Nrm_alpha, p_ang
   !Real(dp) :: RelAnt_N, RelAnt_E, Antx,Anty

!    N_ant=meqn/4
!   write(2,*) 'Next Fitting round'
   m=2+meqn+Nt_grd
   Do i=1,nvar
      j=uiparm(i)
      urparm(m+j)= X(i)
      !k=(j-1)/ParDim
      !Select Case(k)
      !   Case(0)
      !      J_ratio(j-ParDim*k) = X(i)
      !   Case(1)
      !      X_max(j-ParDim*k) = X(i)
      !   Case(2)
      !      R_0(1:ParDim) = X(i)  ! fit probably not stabe when all are fitted independently
      !   Case(3)
      !      L_0(1:ParDim) = X(i)
      !   Case Default
      !      write(2,*) '*** should never reach here!!!!'
      !End Select
   EndDo
   J_total=urparm(m)
   J_ratio(1:ParDim) = urparm(m+1:m+ParDim)
   X_max(1:ParDim)   = urparm(m+ParDim+1:m+2*ParDim)
   R_0(1:ParDim)     = urparm(m+2*ParDim+1)  !:m+3*ParDim)
   L_0(1:ParDim)     = urparm(m+3*ParDim+1)  !:m+4*ParDim)
   J_max(1:ParDim)   = J_total * J_ratio(1:ParDim)
   !write(2,*) 'params',X(1:Nvar)
   !write(2,*) J_max(1:ParDim)
   !write(2,*) X_max(1:ParDim)
   !write(2,*) R_0(1:ParDim)
   !write(2,*) L_0(1:ParDim)
   Flush(Unit=2)
!   EnergyCR=X(1)
!   X_max = X(2)
!   If(RL_param) Then
!      R_0 = X(3)
!      L_0 = X(4)
!   Else
!      X_0 = X(3)
!      Lamx= X(4)
!   EndIf
    !
    !Call FoldWghtCurr(Nt_grd, Nz_grd, RInt_I, LongProfile)
   IF(meqn .ne. (Nz_grd+1) ) write(2,*) '****(meqn .ne. (Nz_grd+1) ):',meqn, (Nz_grd+1)
   !Call ShowCurr(NF, Nt_grd, urparm(1+meqn) )
   i_z=nf
   Call ShowCurr2(i_z, zlow_grd, delt_grd, Nt_grd, urparm(1+meqn), J_max, X_max, R_0, L_0, RL_param )
   !write(2,"(A,100G13.4)") '!profile(t)',urparm(1+meqn:1+meqn+Nt_grd)
   ! RInt_I(0:Nz_grd) = urparm(1:meqn)
   ! LongProfile(0:Nz_grd) =urparm(1+meqn:2*meqn)
   Select Case (uiparm(0))
      CASE (3)
         !write(2,*) '!hcnt_l:hcnt_u,Nt_grd,Nz_grd',hcnt_l,hcnt_u,Nt_grd,Nz_grd
         Do i_z=0,Nz_grd
            urparm(1+i_z)=0.
            Do i_h=hcnt_l,hcnt_u
               W=SUM(Grd_Kx(i_h, 0:Nt_grd, i_z)*urparm(1+meqn:1+meqn+Nt_grd) )
               urparm(1+i_z)= urparm(1+i_z) + W*W
            EndDo
            urparm(1+i_z)=sqrt(urparm(1+i_z))*delt_grd
            !R(i_z+1)=RIntAmpl(i_z)-urparm(1+i_z)  !  urparm(1:2*(Nz_grd+1))
         EndDo
         !write(2,"(A,100G13.4)") '!R_calc(z)',urparm(1:1+Nz_grd)
      CASE (13)
         !write(2,*) '!hcnt_l:hcnt_u,Nt_grd,Nz_grd',hcnt_l,hcnt_u,Nt_grd,Nz_grd
         Do i_z=0,Nz_grd
            urparm(1+i_z)=0.
            Do i_h=hcnt_l,hcnt_u
               W=SUM(Grd_Kr(i_h, 0:Nt_grd, i_z)*urparm(1+meqn:1+meqn+Nt_grd) )
               urparm(1+i_z)= urparm(1+i_z) + W*W
            EndDo
            urparm(1+i_z)=sqrt(urparm(1+i_z))*delt_grd
            !R(i_z+1)=RIntAmpl(i_z)-urparm(1+i_z)  !  urparm(1:2*(Nz_grd+1))
         EndDo
         !stop
      CASE DEFAULT
         Do i_z=0,Nz_grd
            urparm(1+i_z)=SUM( GrdWght(0:Nt_grd,i_z)*urparm(1+meqn:1+meqn+Nt_grd)*delt_grd )   ! GrdWght(0:Nt_grd,0:Nz_grd)
            !write(2,*) '!CASE DEFAULT:', i_z, GrdWght(Nt_grd/2,i_z), &
            !   urparm(1+meqn+Nt_grd/2), delt_grd
            !R(i_z+1)=RIntAmpl(i_z)-urparm(1+i_z)  !  urparm(1:2*(Nz_grd+1))
         EndDo
   END SELECT    !
   !write(2,"(A,100G13.4)") '!R_calc(z)',urparm(1:1+Nz_grd)
   ChiSq= SUM( RIntAmpl(0:Nz_grd)-urparm(1:Nz_grd+1) )
   !W= sqrt(SUM(urparm(1:Nz_grd+1)*urparm(1:Nz_grd+1)))
   W= SUM( urparm(1:Nz_grd+1) )
   If(W.eq.0.) W=1.
   !write(2,*) 'urparm(m):', urparm(m), ChiSq, W, Nz_grd, urparm((Nz_grd+1)/2) !, 1. + ChiSq/W
   !Flush(unit=2)
   W= 1. + ChiSq/W
   urparm(m)=W*urparm(m)  ! energy
   urparm(1:1+meqn+Nt_grd) = urparm(1:1+meqn+Nt_grd) *W  ! renormalize R and driving current
   R(1:Nz_grd+1)= RIntAmpl(0:Nz_grd) - urparm(1:Nz_grd+1)  !  urparm(1:2*(Nz_grd+1))
    !write(2,*) 'SUM(R(1:Nz_grd+1)):', SUM(R(1:Nz_grd+1))
   !
   ChiSq= SUM( R(1:meqn)*R(1:meqn) )
   W =SUM( RIntAmpl(0:Nz_grd)*RIntAmpl(0:Nz_grd) )
   write(2,200) NF, ChiSq, ChiSq/W, X(1:nvar)
200     format('Fitting-try=',I3,', Chi^2=',g12.5,', Chi^2/powr=',g12.5,'; Fit Parameters=',12(1pG12.6,','))
   Flush(Unit=2)
   !
!   Flush(unit=2)
    !write(2,"(A,100G13.4)") '! W, R(i_z):',W, R(1:1+Nz_grd)
    !write(2,"(A,100G13.4)") '        (i_z):', (i_z,i_z=0,Nz_grd)
    !write(2,"(A,100G13.4)") 'RIntAmpl(i_z):', RIntAmpl(:)
    !write(2,"(A,100G13.4)") '  RInt_I(i_z):', urparm(1:1+Nz_grd)
    !write(2,"(A,100G13.4)") 'LongProf(i_z):', urparm(meqn+1:meqn+1+Nz_grd)
    !stop 'CompareR_Intf'
    return
end subroutine CompareR_Intf
