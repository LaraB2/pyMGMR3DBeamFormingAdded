Module CurrentFolding !------------------------
contains
Subroutine UnFoldJ(Base, GrdWght, pow, N_tr,N_zf, del_tr, del_z, z_low, Frac, J_fit)
   ! J(t_r) = (WW)^-1 (PRW)
   ! WW  = \sum_\zeta  W(t_r, \zeta)  *  W(t'_r, \zeta)
   ! PRW = \sum_\zeta  W(t_r, \zeta)  *  PSF_pow(\zeta)
   ! the t_r and zeta_f scales have the same off-set z_low which is implicitly understood
   use constants, only : dp
   Implicit none
   Real(dp), intent(in) :: GrdWght(0:N_tr,0:N_zf)
   Real(dp), intent(in) :: pow(0:N_zf)  ! really its sqrt
   Integer, intent(in) :: N_tr,N_zf
   Real(dp), intent(in) :: del_tr, del_z, z_low
   Real(dp), intent(in) :: Frac
   character(Len=*), intent(in)  :: Base
   Real(dp), intent(out) :: J_fit(0:N_tr)
   !
   Real(dp), allocatable :: WW(:,:)
   Real(dp), allocatable :: PRW(:)
   !Real(dp), allocatable :: GrdWght_r(:,:)
   Real(dp), allocatable :: WF(:,:)
   !Integer :: i_r, i_r2, i_tr, itr_low, itr_hi, i_z, NDc, i_c, j_c
   Integer :: i_r, i_tr, i_z, NDc, i_c, j_c
   Real(dp) :: t_r, chisq, totpow
   Real(dp) :: Residual(0:N_zf)
   !===============================================
   !         ZHESV computes the solution to a complex system of linear equations
   !           A * X = B,
   !         where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS matrices.
   ! https://netlib.org/lapack/explore-html/d3/d9d/group__complex16_h_esolve_gac2def3254215d3a7c56bc162483193d4.html#gac2def3254215d3a7c56bc162483193d4
   !=======================================================
   Character(len=1) :: UPLO
   Integer :: i, NRHS, Neq, LWork, INFO
   Integer, allocatable :: IPIV(:)
   Real(dp), allocatable :: ReB(:,:)
   Real(dp), allocatable :: ReWork(:)
   Real(dp) :: D, Dm(1), W
   Integer, parameter :: NDc_max=40
   Integer(dp) :: Dc(1:NDc_max+1) !, Diz(1:NDc_max+1)
   !
   ! Get grid points Di
   i_c=1
   Dc(i_c)=0  ! grid follows t_r (=D_current) Real distance is Di(i)*del_tr
   !Diz(i_c)=0  ! grid follows t_r (=D_current) Real distance is Di(i)*del_tr
   i_r=0
   W=GrdWght(Dc(i_c),i_r)
   !write(2,*) '!?',del_tr, del_z, Frac,  1./Frac  ! procedure depends strongly on del_ value
   Do i_tr=0,N_tr  ! loop over D_current
      t_r=i_tr*del_tr
      i_z=NINT(t_r/del_z)  ! del_z=del_tr is assumed ?
     ! Write(2,"(A,10g12.3)") '!?',i_tr,i_r,GrdWght(i_tr,i_z),GrdWght(Dc(i_c),i_z), GrdWght(i_tr,i_r), Dc(i_c) &
     !    , W, GrdWght(i_tr,i_z)/W, GrdWght(Dc(i_c),i_z)/W, GrdWght(i_tr,i_r)/W
      !If(GrdWght(Dc(i_c),i_z)/GrdWght(i_tr,i_z) .lt. Frac) Then ! new gridpoint for current folding
      If((ABS(GrdWght(i_tr,i_r)/W) .lt. Frac) .or. (ABS(GrdWght(i_tr,i_r)/W) .gt. 1./Frac )) Then ! new gridpoint for current folding
     !    Write(2,"(A,10g12.3)") '!selected'
         i_c=i_c+1
         Dc(i_c)=i_tr
         i_r=i_z
         W=GrdWght(Dc(i_c),i_z)
         !Diz(i_c)=i_z
         If(i_c.eq.NDc_max) exit
      !Else ! just continue
      !   Write(2,"(A,10g12.3)") '!not selected'
      EndIf
   EndDo
   NDc=i_c+1
   Dc(NDc)=N_tr  ! add last point at top of atmosphere
   write(2,"(A,F5.2,I4,A, 50I4)") 'For frac=',Frac, NDc,' current  gridpoints at',Dc(1:NDc)
   !--------------------------
   ! WF = GrdWght(t_r,:)  * F_i(t_r)  where F_i is a triangle function =1 @ t_r=i*del_r; =0 @ t_r=(i \pm 1)*del_r
   !   identical to GrdWght_r, just for checking
   !N_r=NDc
   Allocate( WF(1:NDc,0:N_zf) )
   WF(:,:)=0.
   write(2,*) 'NDc:', NDc
   flush(unit=2)
   i_c=2
   Do i_tr=0,N_tr
      If(i_tr.gt.Dc(i_c)) Then
         i_c=i_c+1
      EndIf
      d=(i_tr-Dc(i_c-1))*1./(Dc(i_c)-Dc(i_c-1))  ! between 0 and 1
      WF(i_c-1,:)=WF(i_c-1,:) + (1.-d)*GrdWght(i_tr,:)*del_tr
      If(i_c.lt.NDc) Then
         WF(i_c,:)=WF(i_c,:) + d*GrdWght(i_tr,:)*del_tr
      EndIf
   EndDo
   !
   Allocate( PRW(1:NDc), WW(1:NDc,1:NDc) )
   Do i_c=1,NDc ! Loop over  grid points for current determination
      Do j_c=1,NDc  !
         WW(i_c,j_c)= del_z*SUM( WF(i_c,:)* WF(j_c,:) ) ! Integrate over zeta_f
      EndDo
      PRW(i_c)= del_z*SUM( WF(i_c,:)*pow(:) ) ! Integrate over zeta_f
   EndDo !
   !
   Allocate( ReB(1:NDc,1), IPIV(1:NDc) )
   UPLO =  'L'  ! 'U'  (Upper or lower half of A is used (L implies first index is .ge. second)
   NRHS = 1  !  The number of columns in matrix B
   Neq = NDc-1    ! take out last current grid point
   ReB(:,1)=PRW(:)
   !
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
   Call dsysv (UPLO, Neq, NRHS, WW, NDc, IPIV, ReB, NDc, Dm, LWORK, INFO)
   !                 ^ dimension for inversion (= number of eqs to solve)
   !                                ^ first dimension in declaration of WW
   !                                                ^ first dimension in declaration of ReB
   LWORK=nint(Dm(1))
   Allocate (ReWork(LWORK))
   Call dsysv (UPLO, Neq, NRHS, WW, NDc, IPIV, ReB, NDc, ReWORK, LWORK, INFO)
   DeAllocate ( ReWork )
   !
   write(2,*) 'Convergence criterium:',INFO
   If(INFO.gt.0) Then
      Do i=1,NDc
         Write(2,*) i, IPIV(i), WW(i,1:i)
      EndDo
   EndIf
   DeAllocate( PRW, WW )
   !
   ReB(Neq+1:NDc,1)=0.
   write(2,*) 'J_fit(Neq):',Neq, ReB(1:Neq,1)
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(Base)//'-PWL-currentPoints.dat')
   Do i_c=1,Neq
      write(4,"(I3,5G13.5)") i_c,(Dc(i_c)*del_tr + z_low)/1000., ReB(i_c,1)
   Enddo
   Close(unit=4)
   !
   ! Construct current profile
   i_c=2
   J_fit(:)=0.
   Do i_tr=0,N_tr
      If(i_tr.gt.Dc(i_c)) Then
         i_c=i_c+1
      EndIf
      d=(i_tr-Dc(i_c-1))*1./(Dc(i_c)-Dc(i_c-1))  ! between 0 and 1
      J_fit(i_tr) = J_fit(i_tr) + (1.-d)*ReB(i_c-1,1)
      If(i_c.lt.NDc) Then
         J_fit(i_tr) = J_fit(i_tr) + d*ReB(i_c,1)
      EndIf
   EndDo
   !
   Chisq=0.
   TotPow=0.
   Do i_z=0,N_zf
      Residual(i_z)=pow(i_z) - SUM(WF(1:Neq,i_z)*ReB(1:Neq,1))
      !write(2,"(A,I4, 5G13.4)") 'Residual:', i_z, Residual(i_z), pow(i_z), 100.*Residual(i_z)/pow(i_z)
      chisq=Chisq+ Residual(i_z)*Residual(i_z)
      TotPow= TotPow + pow(i_z)*pow(i_z)
      write(2,"(20G13.4)") 'Folding Fit', i_z, pow(i_z), ', % off:',100.* Residual(i_z)/pow(i_z)
   EndDo
   DeAllocate( WF )
   DeAllocate( ReB, IPIV )
   write(2,*) 'ChiSq:', chisq, TotPow, 100.*sqrt( chisq/TotPow),'% in amplitude'
   !
   Return
End Subroutine UnFoldJ
!=================== Obsolete:
End Module CurrentFolding !------------------------
