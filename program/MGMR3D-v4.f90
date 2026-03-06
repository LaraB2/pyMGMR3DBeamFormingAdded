!------------------------------------------------
    Subroutine MGMR3D
!    subroutine ttret
!    all length and times are measured in meters where  1 ns=0.3 meter=1GHz^{-1}
      use constants
      use BigArrays
    use RFootPars, only : NTo, MoliereRadius,test
    use eventdata, only : Fitting,ReadInput
    use T_Analyze, only : Analyze
    implicit none
!    logical first
!    integer, parameter :: NAnt_max=1  ! 900
    INTEGER DATE_T(8),i,k ,nxx,idi, CD_i ! ,inu,ith,Nth, Ntr
    CHARACTER*12 REAL_C (3)
    character*80 line
    real(dp) :: r,z,T_o,t_r,nu
    real(dp) :: Exk,AxDk,Eyk,AyDk,Ark,Erhk ! ,StI,StQ,StU,StV,,X_rh,  a,b,c,dfhl,dfha,dr,fh,NN,rh,theta
    real(dp) :: ha,hb,sqha,sqhb,ddd,s,t_shft,lam,dlamb
    common /W_TrCurr/ Norm_tc,Moli_tc,s_tc,dalpha,Nrm_alpha(0:41)
      real(dp) :: Norm_tc,Moli_tc,s_tc,dalpha,Nrm_alpha ! ,alpha
    common / time / zeta_c,hcto_c
      integer :: N_branch
      real(dp) :: T_oi,zeta_c,hcto_c,zeta_d,hcto_d ! ,Antd
!
    If(ReadInput) call SetParams
!    first=.true.
!    read(*,*) out_file  ,plt_file
!    Test=.true.
!    Test=.false.
    Moli_tc=MoliereRadius
!
!
!===========================================================================================
    CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T) !-----------------------------
    write(*,232) 'Initialize',(DATE_T(i),i=5,8)               !-----------------------------
    write(2,232) 'Initialize',(DATE_T(i),i=5,8)               !-----------------------------
    Flush(unit=2)
!===========================================================================================
!232 FORMAT(3X,A10,'@',I2,':',I2,':',I2,'.',I3,1X,25(1H-))
!
    call AssignDim()   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !   Initialize Basis time-trace grid for observer time after 'Cherenkov time'=t_ch(idi) for line emissions
    k=2
    t_tb(1)=-0.5*tTrace_step   !  half stepsize moved to avoid touching the Cherenkov time.
!    DO i=2,k
!        t_tb(i)=(i-1.5)*tTrace_step   !  half stepsize moved to avoid touching the Cherenkov time.
!    enddo
    DO i=k,tTrace_dim_b
        t_tb(i)=t_tb(i-1)+tTrace_step*(1.+(i-k)/2.)  !  half stepsize moved to avoid touching the Cherenkov time.
    enddo
    !   Initialize Observer time-trace grid after 'shower-ground'=0 for EAS emission
    !   needs to be regular to allow for fourier transform
    DO i=1,tTrace_dim_o
        t_to(i)=i*tTrace_step   !  starts after shouwer touching ground.
    enddo
    write(2,"('last true basis-grid time=',f8.3,', last observer time=',f8.3,'[ns]')") t_tb(NTo)/c_l,t_to(tTrace_dim_o)/c_l
    !
    Line2Obsrv(1)=0.5*ObsDist_Step  ! not used (yet) -------------------------------------------------
    k=1+150/ObsDist_Step
    if( (k+2) .ge. ObsDist_dim) then
        write(2,*) '***ERROR*** ObsDist_dim should be larger than',k+2
        stop 'ObsDist_dim'
    endif
    DO i=2,k ! regular grid in the shower plane upto about 150 m to cover Cherencov,
        Line2Obsrv(i)=Line2Obsrv(i-1) + ObsDist_Step
    enddo
    Do i=k+1,ObsDist_dim
        Line2Obsrv(i)=Line2Obsrv(i-1) + ObsDist_Step*(1.+(i-k-1.)/2.)
    enddo
    ! ------------------------------------------------------------------------------------------------------
!
    call Initialize_shower
!
    k=1     ! First index to be used in rs integral in "LateralInt"
    Line2Core(1)=CoreDist_A(0)/2.
    Line2Core(2)=CoreDist_A(0)
    Do i=1,CoreDist_Dim
        Line2Core(2*i+1)=(CoreDist_A(i-1) + CoreDist_A(i))/2.
        Line2Core(2*i+2)=CoreDist_A(i)
    enddo
    Line2Core(k-1)=-Line2Core(k)
    Line2Core(2*CoreDist_Dim+3)=3.*CoreDist_A(CoreDist_Dim)/2. - CoreDist_A(CoreDist_Dim-1)/2.
!
!===========================================================================================
    CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T) !-----------------------------
    if(test) write(*,232) 'pancakeInt',(DATE_T(i),i=5,8)               !-----------------------------
    if(test) write(2,232) 'pancakeInt',(DATE_T(i),i=5,8)               !-----------------------------
!===========================================================================================
!
    allocate(Ex_tb(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim), Ex_spld(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim))
    allocate(AxD_tb(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim), AxD_spld(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim))
    allocate(Ey_tb(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim), Ey_spld(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim))
    allocate(AyD_tb(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim), AyD_spld(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim))
    allocate(Ar_tb(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim), Ar_spld(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim))
    allocate(Erh_tb(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim), Erh_spld(1:tTrace_dim_b,ObsDist_dim,0:CoreDist_Dim))
    Ex_tb(:,:,:)=0.d0
    AxD_tb(:,:,:)=0.d0
    Ey_tb(:,:,:)=0.d0
    AyD_tb(:,:,:)=0.d0
    Ar_tb(:,:,:)=0.d0
    Erh_tb(:,:,:)=0.d0
!  t-t_retarded data for plotting only
   if(ObsDist_dim.gt.1 ) Then
      k=min(ObsDist_dim,45)
      OPEN(UNIT=4,STATUS='unknown',FILE='plot/t_tret.dat')
      write(4,"(1x,50F7.1)") AtmHei_step,((idi-0.5)*ObsDist_Step, idi=1,k),0.0
      Do i = 0, AtmHei_dim
         z=i*AtmHei_step !     t_r=-z
         !\lambda=hcto - \sqrt{(R-N\zeta)^2+(N^2-1) d^2} \;. \eqlab{Lam_zeta}
         !lambda(i)=hcto_c + sqrt( (sqrt(ObsDist*ObsDist+z*z)-(1.+xi(i))*z)**2+((1.+xi(i))**2-1.)*ObsDist*ObsDist)
         Do idi=1,k
            ObsDist=(idi-0.5)*ObsDist_Step !        ro=(id-0.5)*ObsDist_Step
            T_obs(idi) = (-z+sqrt(ObsDist*ObsDist+z*z)*(1.+xi(i)))/0.3  ! 0.3 to convert from [m] to [ns]
         EndDo
         write(4,*) z/1000.,(T_obs(idi), idi=1,k)
      EndDo !i = 0, AtmHei_dim
      Close(unit=4)
   EndIf
!
    Do idi=1,ObsDist_dim
        ! write(2,*) 'idi set to',idi
        ObsDist=(idi-0.5)*ObsDist_Step !        ro=(id-0.5)*ObsDist_Step
        s=ObsDist*(1.+xi(0))
        !   Find Branch points (if any)
        nxx=-1
        do i = 0, AtmHei_dim
            z=i*AtmHei_step
            t_r=-z
            T_obs(i)=t_r+sqrt(ObsDist*ObsDist+z*z)*(1.+xi(i))
            if(T_obs(i).lt.s) then  ! find branchpoint corresponding to minimal T_obs
              s=T_obs(i)
              r=T_obs(i)
              k=i
            else if(T_obs(i).gt.s) then  ! we are on the side where T_obs increases again
              if(T_obs(i).gt.r) then   ! find branchpoint corresponding to maximal T_obs
                r=T_obs(i)
                nxx=i
              endif
            endif
        enddo
        !write(2,*) 'MGMR3D: ObsDist=',ObsDist
        call FindZetaC(k,zeta_c,hcto_c)
        !   If(.not. Fitting) write(2,211) ObsDist,zeta_c/1000.,-hcto_c
        !211 format(/,'distance to core=',f5.1,'[m]',/,' Cherenkov height=',F7.3,'[km], first signal received at t=',f8.4,'[m/c]')
        if(nxx.gt.1 .and. (Nxx .lt. AtmHei_dim)) then
          N_branch=4
          call FindZetaC(nxx,zeta_d,hcto_d)
          !    If(.not. Fitting) write(2,212) zeta_d/1000.,-hcto_d
          !212 format(' second Cherenkov height=',F7.3,'[km], signal received before t=',f8.4,'[m/c]')
        endif
        z_ch(idi)=zeta_c  ; t_ch(idi)=-hcto_c
        !write(2,*) 'Cherenkov',idi,zeta_c,t_ch(idi),k,T_obs(AtmHei_dim)
        if(k.ge.AtmHei_dim) then   ! No cherenkov solution found below max height
            t_ch(idi)=s
        endif
        !write(2,*) idi,k,r,s,nxx,t_ch(idi)
        !
        if(ObsDist_dim.eq.1) OPEN(UNIT=4,STATUS='unknown',FILE='plot/t_tret.dat')
        do i = 0, AtmHei_dim
            z=i*AtmHei_step
            t_r=-z
            !\lambda=hcto - \sqrt{(R-N\zeta)^2+(N^2-1) d^2} \;. \eqlab{Lam_zeta}
            lambda(i)=hcto_c + sqrt( (sqrt(ObsDist*ObsDist+z*z)-(1.+xi(i))*z)**2+((1.+xi(i))**2-1.)*ObsDist*ObsDist)
            !        if(i.lt.nt) z=z+dz  ! distance dz behind front
            if(ObsDist_dim.eq.1) write(4,*) i,T_obs(i),lambda(i) ! ,t_r+sqrt(ObsDist*ObsDist+r*r)*(1.+xi(i)) !&
            !                ,t_r+sqrt(d*d+z*z)*(1.+xi(i+1))
        enddo
        if(ObsDist_dim.eq.1) close(unit=4)
        !
        Do CD_i=0,CoreDist_Dim  ! Distance to the shower core, called rs in radial integral
        !   In performing this calculation the same same distance to observer (geometry) is kept
        !   Keeping the same geometry thus keeps the same Cherenkov conditions.
        !   What changes is the pancake thikness parameter that depends on distance to core.
            !   get integration grid in lambda -----------------------------------------------
            lam=PancakeThi(CD_i)
            do  ! construct conditionally a grid for the lambda (square of pancake thickness) integrals
                dlamb=lam/20.
                dL_I(0)=0.
                zeta_I(1,0)=zeta_c  ; zeta_I(2,0)=zeta_c
                if(ObsDist_dim.eq.1) write(2,*) 'Integration grid:, i, lambda, weight, zeta1, zeta2'
                !    L_I(LamInt_dim)=NT*dz
                LamGrid_Nrt=sqrt(1.*lam)/dlamb
                LamGrid_max=LamInt_dim
                ha=0. ; sqha=0.       ! starting value of grid in lambda
                Do i=1,LamInt_dim   ! Check integrate over pancake thickness
                      hb=i*i*dlamb*dlamb ! next grid-point in lambda, quadratic grid
                      ddd=hb-ha
                      if(ddd.gt.lam/8.) then  ! zero crossing of derivative pancake function is at lam/2.
                        ddd=lam/8.
                        hb=ha+ddd
                      endif
                      sqhb=sqrt(hb)
                      if(i.le.LamGrid_Nrt) then ! d sqrt(lambda)
                        L_I(i)=((sqhb+sqha)/2.)**2      ! central value in sqrt(lambda)
                        JL_I(i)=sqhb+sqha ! sqrt(lambI)  !Jacobian
                        dL_I(i)=sqhb-sqha     ! difference of the sqrt(lambda) over interval
                        !        End_I(i)=-sqha ! +sqrt(h_c)  =length last interval/JL_I(i)
                      else      ! linear integral
                        L_I(i)=(hb+ha)/2.
                        JL_I(i)= 1.
                        dL_I(i)=hb-ha
                        !        End_I(i)=-ha  ! +h_c  =length last intervall
                      endif
                      ha=hb  ; sqha=sqhb
                      call getZeta(L_I(i),zeta_I(1,i-1),zeta_I(1,i))
                      !      if(ObsDist_dim.eq.1) write(2,*) i,l_i(i),JL_I(i)*dL_I(i),zeta_I(1,i),zeta_I(2,i)
                      if((zeta_I(1,i).lt.0.) .and. (zeta_I(2,i) .gt. AtmHei_dim*AtmHei_step)) then
                        LamGrid_max=i  ; exit
                      endif
                      !      if(hb .gt. (-t_ch(idi)+NTo*dto)) exit  !  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& check why this does not stop
                enddo
                !stop '1' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !    if(idi.eq.22) write(2,*) '  i',i
                if(i.gt.10) exit ! Guarantee a minimal number of integration points
                lam=lam/2.
            enddo
            !stop '2' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(CD_i.le.1) then
              !     If(.not. Fitting) write(2,*) 'lam-max:',i,PancakeThi(CD_i),lam
              if(i.ge.LamInt_dim .and. (.not. Fitting) ) then
                !   write(2,*) 'Lambda grid exceeded, d=',ObsDist,hb,t_ch(idi)+NTo*tTrace_step,zeta_I(1,LamInt_dim),zeta_I(2,LamInt_dim)
                !write(*,*) 'extent Lam-grid@ zeta=',zeta_I(1,LamInt_dim),'t=',(t_ch(idi)+L_I(LamInt_dim))/tTrace_step
              endif
            endif
            !    if(Ndi.eq.1) close(unit=4)
            !
            !write(2,*) 'MGMR3D: @',CD_i,', first (h-h_c) integation steps:', L_I(1), L_I(2),zeta_I(1,1),zeta_I(1,2)
            !
            !   Calculate point-like (in radial extent) fields --------------------------------
            DO i=2,NTo
              T_o=t_ch(idi)+t_tb(i)
              ! Note that around d=140 m an inflection occurs and an extra branchpoint should be considered.
              !    if(idi.eq.22 ) write(2,*) '  T_o',i
              call CurrField(Exk,AxDk,Eyk,AyDk,Ark,Erhk, CD_i,T_o)
              !      if(ObsDist_dim.eq.1) write(2,*) T_o,'Exk=',Exk
              !   if(idi.eq.22 ) write(2,*) '  T_o',T_o,exk,axdk
              Ex_tb(i,idi,CD_i) =Exk
              AxD_tb(i,idi,CD_i)=AxDk
              Ey_tb(i,idi,CD_i) =Eyk
              AyD_tb(i,idi,CD_i)=AyDk
              Ar_tb(i,idi,CD_i) =Ark
              Erh_tb(i,idi,CD_i)=Erhk
              !      Ex(i)=sin(2.*pi*T_o/1.)  !lambda=6 m should be 50 MHz, 1 m = 300 MHz
              !      Ey(i)=0. ; if(i.eq.(Nto+NToi/2)) Ey(i)=100.
              !  if(int(T_o/0.1)==65 .and. Obsdist.gt.400. .and. CD_i==0) write(2,*) 'a',Ar_tb(i,idi,CD_i),Ark
            enddo
            !if(idi==41 .and. CD_i==0) write(2,*) 'b',Ar_tb(5,idi,CD_i),t_ch(idi)+t_tb(5)
            !
            !     make a smooth pulse-tail to zero
            k=tTrace_Offset/2
            DO i=1,k
              Ex_tb(NTo+i,idi,CD_i)  =Exk*(k-i)/k
              AxD_tb(NTo+i,idi,CD_i) =AxDk*(k-i)/k  ! make a smooth pulse-tail to zero
              Ey_tb(NTo+i,idi,CD_i)  =Eyk*(k-i)/k    ! make a smooth pulse-tail to zero
              AyD_tb(NTo+i,idi,CD_i) =AyDk*(k-i)/k  ! make a smooth pulse-tail to zero
              Ar_tb(NTo+i,idi,CD_i)  =Ark*(k-i)/k    ! make a smooth pulse-tail to zero
              Erh_tb(NTo+i,idi,CD_i) =Erhk*(k-i)/k    ! make a smooth pulse-tail to zero
            enddo
            !301 Format(f5.1,4f9.1,i2)
            !302 Format(f6.2,6e13.5)
            !
            call spline_cubic_set( tTrace_dim_b, t_tb(1), Ex_tb(1,idi,CD_i), Ex_spld(1,idi,CD_i) )
            call spline_cubic_set( tTrace_dim_b, t_tb(1), AxD_tb(1,idi,CD_i), AxD_spld(1,idi,CD_i) )
            call spline_cubic_set( tTrace_dim_b, t_tb(1), Ey_tb(1,idi,CD_i), Ey_spld(1,idi,CD_i) )
            call spline_cubic_set( tTrace_dim_b, t_tb(1), AyD_tb(1,idi,CD_i), AyD_spld(1,idi,CD_i) )
            call spline_cubic_set( tTrace_dim_b, t_tb(1), Ar_tb(1,idi,CD_i), Ar_spld(1,idi,CD_i) )
            call spline_cubic_set( tTrace_dim_b, t_tb(1), Erh_tb(1,idi,CD_i), Erh_spld(1,idi,CD_i) )
            !t_shft=0.5*tTrace_step
            !Do i=k,tTrace_dim-1       ! tTrace_Offset,NTo+tTrace_Offset
            !    call spline_cubic_val( tTrace_dim-tTrace_Offset, t_tb(k), Ex_tb(k,idi,CD_i), Ex_spld, t_tb(i)+t_shft, Ex_tb_h(i))
            !enddo
            if(test) then
              write(line,"(i2.2,'-',i2.2)") idi,CD_i
              OPEN(UNIT=4,STATUS='unknown',FILE='plot/Ex_dd'//trim(line)//'.dat')
              ! write(2,*) 't_shft',t_shft,trim(line)
              write(4,"('!',9x,'t[m/c]',20x,'E_x',23x,'t_tb',23x,'A_r',23x,'A_yD',23x,'A_r',23x,'E_rh' )")
              DO i=1,tTrace_dim_b-1
                  T_o=(t_ch(idi) +t_tb(i))/c_l  ! in [n sec]
                  write(4,*) T_o,Ex_tb(i,idi,CD_i),t_tb(i)/c_l,Ar_tb(i,idi,CD_i), Ex_spld(i,idi,CD_i),Ar_spld(i,idi,CD_i)
                  !write(4,*) T_o,Ex(i),AxD(i),Ey(i),AyD(i),Ar(i),Erh(i) ! ,real(Cex(i)),imag(CEx(i)),real(Cey(i)),imag(CEy(i))
              enddo
              close(unit=4)
            endif
            !
        enddo ! CD_i
        !
        ! stop
    enddo   ! idi loop over distances to the core
    if(test) then
        write(2,*)'t_ch(1:ObsDist_dim) [ns], ObsDist_dim=', ObsDist_dim,', step=',ObsDist_Step,'[m]'
        write(2,"(10f7.3/,10f7.3)") t_ch/c_l
        write(2,*)'z_ch(1:ObsDist_dim) [m]'
        write(2,"(10f7.0/,10f7.0)") z_ch
    endif
    !stop
    !
    if(ObsDist_dim.eq.1) stop
    !
    deallocate(L_I, JL_I, dL_I) !; deallocate(End_I) ;
    deallocate(zeta_I, Line2Obsrv)
    deallocate(T_obs, lambda)
    deallocate(xi, dxi, ddxi)
    deallocate(Ix, Iy, IQ)
    deallocate(alpha_tr, Ix_int, Iy_Int)
    !
!===========================================================================================
    CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T) !-----------------------------
    if(test) write(*,232) 'radial-int',(DATE_T(i),i=5,8)               !-----------------------------
    if(test) write(2,232) 'radial-int',(DATE_T(i),i=5,8)               !-----------------------------
!===========================================================================================
!
    allocate(Ex_to(1:tTrace_dim_o,ObsDist_dim) , Ey_to(1:tTrace_dim_o,ObsDist_dim) )
    allocate(Er_to(1:tTrace_dim_o,ObsDist_dim), ObsPlsTime(1:ObsDist_dim) )
    !
    Ex_to(:,:)=0.d0 ;  Ey_to(:,:)=0.d0 ;  Er_to(:,:)=0.d0
    !
    ObsPlsTime(1)=t_ch(1)
    Do idi =1,ObsDist_dim
        ObsDist=idi*ObsDist_Step   ! Note: this is inbetween the grid points used before
        !if(.not. fitting) write(*,*) 'lateral:',idi,ObsDist
        call LateralInt(idi)
        !
        if(idi.gt.1) then
            call Get_ObsPlsDelay(idi)
        endif
        if(test) then
            write(line,"(i2.2)") idi
            OPEN(UNIT=4,STATUS='unknown',FILE='plot/Ex_int-d'//trim(line)//'.dat')
            write(4,"('!',9x,'t[ns]',20x,'E_x',23x,'E_y',23x,'E_r',2x,'@ ObsDist=',F9.3 )") ObsDist
            DO i=1,tTrace_dim_o-1
!              T_o=t_to(i)/c_l  ! in [n sec]
              write(4,*) t_to(i)/c_l,Ex_to(i,idi),Ey_to(i,idi),Er_to(i,idi)
            enddo
            close(unit=4)
        endif
        !
    enddo
    if(test) then
        write(2,*) 'Observer pulse timings(1:ObsDist_dim) [ns]:'
        write(2,"(20f6.1)") ObsPlsTime/c_l
    endif
    deallocate(Z_ch, T_ch) !;
    deallocate(t_tb, Ex_tb, AxD_tb, Ey_tb, AyD_tb, Ar_tb, Erh_tb)
    deallocate(Ex_spld, AxD_spld, Ey_spld, AyD_spld, Ar_spld, Erh_spld)
!
!===========================================================================================
    CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T)  !-----------------------------
    if(.not. fitting) write(*,232) 'Analyze',(DATE_T(i),i=5,8)  !-----------------------------
    if(test) write(2,232) 'Analyze',(DATE_T(i),i=5,8)          !-----------------------------
!===========================================================================================
!
    Call  Analyze(fitting)
!
    CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T)
    !write(*,232) 'endgame',(DATE_T(i),i=5,8)
    !write(2,232) 'endgame',(DATE_T(i),i=5,8)
232 FORMAT(3X,'start ',A10,' @ ',I2,':',I2,':',I2,'.',I3,1X,25(1H-))
    Flush(unit=2)
!
9    continue
    call DAssignArr()
    deallocate(PenDepth)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T)
    WRITE(2,231) DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
231 FORMAT(3X,'MGMR3D, run on ',I2,'/',I2,'/',I4,' , stopped at ',&
          I2,':',I2,':',I2,'.',I3,1X,25(1H-))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    end subroutine MGMR3D
!    end subroutine ttret

!------------------------------
    Include 'MGMR3D_subr.f90'
    include 'MGMR3D_spline.f90'
!------------------------------
!------------------------------
