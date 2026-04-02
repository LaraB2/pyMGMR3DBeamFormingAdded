!
Subroutine SourceMoments(CD_i, zeta, h, fh,dfh, SEx, SAxD, SEy, SAyD, SAr, SErh)
!  Rdlambda_hdr=dlambda_hdr/lambda
!  Int_dh=Int_dh + dz* dhdz * dfh ! should be =0 {= [f(h=0)=0] - [f(h=\infty)=0] } when integrated over the complete height,
!     any finite value is due to integration-step-size problems which is why it is more stable to include the subtraction
!     at large observer distances where 'cald' is almost constant over the integration range,
!     the integral over dfh*Jx should yield zero, but this is from a cancellation of relatively large numbers.
!   Calculate the source @ height=\zeta as function of tmh=(t-h) for a fixed observerdistance to this shower ray
!   A backtrace integration is performed over observer position
!   CD_i (CoreDist_A(CD_i)) determines the position of the pencil shower w.r.t. the core of the real shower which determines the pancake thickness.
!   zeta= height in the atmosphere of the back-trace voxel
!   Shower time varies with h, distance behind the shower front.
!
   use BigArrays, only : Ix,Iy,IQ, Ix_int,Iy_Int
   use Atmosphere, only : AtmHei_step, AtmHei_dim
   use Pancakefunction, only : CoreDist_Dim, CoreDist_A, alpha_tr, PancakeThi
   use constants, only : dp
   Use Pancakefunction, only : fieh
   implicit none
   !    real*8 :: Ex(0:CoreDist_Dim),AxD(0:CoreDist_Dim),Ey(0:CoreDist_Dim),AyD(0:CoreDist_Dim),Ar(0:CoreDist_Dim),nt_r,t_o
   real(dp), intent(in) :: zeta, h
   real(dp), intent(out) :: fh, dfh
   real(dp), intent(out) :: SEx, SAxD, SEy, SAyD, SAr, SErh
   integer, intent(in) :: CD_i
   real(dp) :: nt_r, di,Jx,Jx_int,Jy,Jy_int,dJx,dJy,JQ
   !h_c,h,di,R,NN,dr,calD,Jx,Jx_int,Jy,Jy_int,dJx,dJy,JQ, t_r
   !real(dp) :: Intx,Inty,Intr,B,Lbd,weigh ! ,dh,max_h,fh_ce
   real(dp) :: dfhl,dfha,alpha,  dfhdr, dlambda_hdr, Rdlambda_hdr, lambda
   integer :: i !, idi, i_t, i_h, N_h
   !
   lambda=PancakeThi(CD_i)
   if(CD_i.eq.0) then  ! calculate d(lambda)/dr
     dlambda_hdr=0.5*(PancakeThi(1)-PancakeThi(0))/(CoreDist_A(1)-CoreDist_A(0))
   elseif(CD_i.eq.CoreDist_Dim) then
     dlambda_hdr=(PancakeThi(CD_i)-PancakeThi(CD_i-1))/(CoreDist_A(CD_i)-CoreDist_A(CD_i-1))
   else
     dlambda_hdr=(PancakeThi(CD_i+1)-PancakeThi(CD_i-1))/(CoreDist_A(CD_i+1)-CoreDist_A(CD_i-1))
   endif
   Rdlambda_hdr=dlambda_hdr/lambda
   !
   SEx = 0.
   SAxD= 0.
   SEy = 0.
   SAyD= 0.
   SAr = 0.
   SErh=  0.
   If(h.lt.0.) Return
   nt_r=zeta-h     ! - t_r (=negative t_r) to have a positive number
   i=(nt_r/AtmHei_step)
   if(i.lt.1) Return
   di=nt_r/AtmHei_step-i
   !
   JQ= di*IQ(i+1) + (1.-di)*IQ(i)
   Jx= di*Ix(i+1) + (1.-di)*Ix(i)
   Jx_int= di*Ix_int(i+1) + (1.-di)*Ix_int(i)  ! contribution from moving dipole
   dJx= (Ix(i+1) - Ix(i))/AtmHei_step  !  seems to give a vanishing O(10^{-4}) contribution to E
   Jy= di*Iy(i+1) + (1.-di)*Iy(i)
   !write(2,*) 'jy',jy,i,di
   Jy_int= di*Iy_int(i+1) + (1.-di)*Iy_int(i)
   dJy= (Iy(i+1) - Iy(i))/AtmHei_step  !  seems to give a vanishing contribution to E
   alpha=di*alpha_tr(i+1) + (1.-di)*alpha_tr(i) ! the parameter that determines the width of the pancake function
   !
   call fieh(h,CD_i,alpha,fh,dfhl,dfha)  ! fh is being defined
   !    write(2,*) 'fieh2',h,CD_i,alpha,fh,dfhl,dfha
   dfh=dfhl + dfha *(alpha_tr(i+1)-alpha_tr(i))/AtmHei_step
   dfhdr = - (fh + dfhl*h ) * Rdlambda_hdr
   !
   ! Need to check the sign in the second (small) term, d/dh(J_x), even changing signs hardly reflect in results !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SEx = (dfh*Jx - fh *dJx) ! (dfh*Jx/cald - fh *dJx/cald)
   SAxD= Jx_int*fh   ! proportional to integrated current, thus like the moving dipole
   SEy = (dfh*Jy - fh *dJy)
   SAyD= Jy_int*fh
   SAr =  -(fh *JQ)     ! still needs to be times dw/dr
   SErh=  -(dfhdr *JQ)    ! still needs to be times w(r)
   !
   SAxD=Jx
   !
   Return
End Subroutine SourceMoments
!------------------------------

Subroutine PencilSource_zeta(CD_i, zeta_0)
!   Calculate the source @ height=\zeta as function of tmh=(t-h) for a fixed observerdistance to this shower ray
!   A backtrace integration is performed over observer position
!   CD_i (CoreDist_A(CD_i)) determines the position of the pencil shower w.r.t. the core of the real shower which determines the pancake thickness.
!   zeta= height in the atmosphere of the back-trace voxel
!   Shower time varies with h, distance behind the shower front.
!
    use BigArrays, only :  ObsDist_Step, ObsDist_dim
   use Atmosphere, only : AtmHei_step, AtmHei_dim, xi,dxi
   use Pancakefunction, only : CoreDist_Dim, CoreDist_A, alpha_tr, PancakeThi
    use BigArrays, only : t_tb, Ex_tb,AxD_tb,Ey_tb,AyD_tb,Ar_tb, Erh_tb, t_ch, Z_ch
    use BigArrays, only : Ix,Iy,IQ
    use RFootPars, only : NTo
    use constants, only : dp
    use eventdata, only : PlotBase
   Use Pancakefunction, only : fieh
    implicit none
!    real*8 :: Ex(0:CoreDist_Dim),AxD(0:CoreDist_Dim),Ey(0:CoreDist_Dim),AyD(0:CoreDist_Dim),Ar(0:CoreDist_Dim),nt_r,t_o
    real(dp), intent(in) :: zeta_0
    integer, intent(in) :: CD_i
    real(dp) :: h_c, h_0,di, R_0,NN_0, dr, t_r  ! , zeta, h,Jx,Jx_int,Jy,Jy_int,dJx,dJy,JQ
    real(dp) :: fh,dfh,dfhl,dfha,alpha, dlambda_hdr, Rdlambda_hdr, dfhdr, Int_dh, weight, norm, normc
    real(dp) :: Ex,AxD,Ey,AyD,Ar,Erh, Rx, Sx_h, Sx, Sxc, zeta_s, d_z, h_s, NN_s, R_s, ObsDist
    Real(dp) :: S1x_h, S1x, S1xc
    Real(dp) :: SEx, SAxD, SEy, SAyD, SAr, SErh
    integer :: i, idi, i_t, i_h, N_h, i_z, N_z, idi_N
    real(dp) :: hto, T_o, dhdz, lambda, dto, d_h,nrm, OD_max
    real(dp) :: DistMax=400.d0
    Real(dp) :: Kx, Kxc, Kx_h, h_mx, zeta_mx, h_mn, zeta_mn, h_mn2, zeta_mn2, dr_s, ddr_s, CalD_s, dCalD_s
    Real(dp) :: Kx_stk0, Kx_stk1, Kx_stk2
    real(dp) :: Zeta_Ch
    Integer :: iz_lo, iz_hi, i_sign
    character(len=10) :: txt
      Real(dp) :: tr_s, WeiFie, SxW, WeiFiec, SxWc
!
   Ex=0.   ; Ey=0. ; Erh=0.
   !zeta=zeta_0
   lambda=PancakeThi(CD_i)
   if(CD_i.eq.0) then  ! calculate d(lambda)/dr
     dlambda_hdr=0.5*(PancakeThi(1)-PancakeThi(0))/(CoreDist_A(1)-CoreDist_A(0))
   elseif(CD_i.eq.CoreDist_Dim) then
     dlambda_hdr=(PancakeThi(CD_i)-PancakeThi(CD_i-1))/(CoreDist_A(CD_i)-CoreDist_A(CD_i-1))
   else
     dlambda_hdr=(PancakeThi(CD_i+1)-PancakeThi(CD_i-1))/(CoreDist_A(CD_i+1)-CoreDist_A(CD_i-1))
   endif
   Rdlambda_hdr=dlambda_hdr/lambda
   !
   write(2,*) 'Source_zeta:', zeta_0, ', CD_i, lambda:', CD_i, lambda
   write(txt,"(I3.3,i2.2)") NINT(zeta_0/100.),CD_i
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'zeta_'//trim(txt)//'.csv')
   write(4,*) ' Source_zeta:', zeta_0, ', CD_i, lambda:', CD_i, lambda,' !'
   write(4,*) '!     h     Rx(backtrace antennas),  Sxc(MGMR, cut antenna domain),  Sxc(MGMR, full antenna domain), Source !'
   OPEN(UNIT=14,STATUS='unknown',FILE=TRIM(PlotBase)//'zetaW_'//trim(txt)//'.csv')
   write(14,*) ' Source_zeta:', zeta_0, ', CD_i, lambda:', CD_i, lambda,' !'
   write(14,*) '!     -tr_s     weifiec,  weifie !'
   !N_h=5
   N_h=25
   d_h=5.*lambda/N_h
   N_z=N_h*5
   d_z=d_h*20.
   If((ObsDist_dim-1)*ObsDist_Step .lt. DistMax) DistMax=(ObsDist_dim-1)*ObsDist_Step
   !
   nrm=1.e6
   norm=0.
   normc=0.
   idi_N = ObsDist_dim ! Cherenkov heights and distances are not defined for larger distances
   Do idi=1, idi_N
      ObsDist=(idi-0.5)*ObsDist_Step !        ro=(id-0.5)*ObsDist_Step
      R_0=sqrt(ObsDist*ObsDist +zeta_0*zeta_0)
      weight=2.*ObsDist*ObsDist_Step   !*(DistMax-ObsDist)
      norm=norm+weight/(R_0*R_0)
      If( (ObsDist.lt.DistMax) ) normc=normc+weight/(R_0*R_0)
   EndDo ! idi=1, idi_N
   norm=norm*nrm
   normc=normc*nrm
   !
   Do i_h=-N_h/3, N_h
      h_0=i_h*d_h
      !h=h_0
      !write(2,*) 'Source_set; z,h=',zeta_0,h_0
      !Flush(unit=2)
      t_r=h_0-zeta_0 ! should be negative, all quantities at the source are given
      i=(zeta_0/AtmHei_step)
      if(i.lt.1) return
      di=zeta_0/AtmHei_step-i
      NN_0=1.+ di*xi(i+1) + (1.-di)*xi(i)
      dr=   di*dxi(i+1) + (1.-di)*dxi(i)
      !
      Rx=0.
      OD_max=-1.
      Sx=0.
      Sxc=0.
      S1x=0.
      S1xc=0.
      Kx=0.
      Kxc=0.
      !write(2,*) 'dto:', 'dto, ObsDist, Rx/norm, Ex/(NN_0*R_0), Ex, cald '
      ! Determine weighting function for the current(tr_z)
      If(i_h.eq.N_h/3) Then
      i_t=i
      SxW=0.
      SxWc=0.
      Do i_z=1,AtmHei_dim,1
         tr_s=-i_z*AtmHei_step
         WeiFie=0.
         WeiFiec=0.
         Do idi=1, idi_N
            ObsDist=(idi-0.5)*ObsDist_Step !        ro=(id-0.5)*ObsDist_Step
            R_0=sqrt(ObsDist*ObsDist +zeta_0*zeta_0)
            t_o=NN_0*R_0+t_r !observer is summed over
            !
            Call Findhs(t_o, ObsDist, tr_s, h_s, NN_s, dr_s)
            If(h_s.le.0) cycle
            !
            zeta_s=h_s-tr_s
            alpha=alpha_tr(i_z) ! the parameter that determines the width of the pancake function
            call fieh(h_s,CD_i,alpha,fh,dfhl,dfha)  ! fh is being defined
            dfh=dfhl + dfha *(alpha_tr(i_z+1)-alpha_tr(i_z))/AtmHei_step  ! depends on retarded time
            R_s=sqrt(zeta_s*zeta_s +ObsDist*ObsDist)
            CalD_s=NN_s*( R_s - NN_s*zeta_s - dr_s*R_s*R_s)
            WeiFie=WeiFie+ 2.*ObsDist*ObsDist_Step*dfh/((NN_s*R_s-CalD_s)*NN_0*R_0)
            If(ObsDist.lt.DistMax) Then
               WeiFiec=WeiFiec+ 2.*ObsDist*ObsDist_Step*dfh/((NN_s*R_s-CalD_s)*NN_0*R_0)
            EndIf
            !write(2,*) h_s, fh, dfh, zeta_s
         EndDo ! idi=1, idi_N
         !write(2,*) 'WeiFie:',-tr_s,WeiFie
         write(14,*) -tr_s, WeiFiec, WeiFie
         SxW=SxW+ WeiFie*Ix(i_z)*AtmHei_step
         SxWc=SxWc+ WeiFiec*Ix(i_z)*AtmHei_step
      EndDo !i_z=1, AtmHei_dim
      write(2,*) 'SxW',SxWc/normc,SxW/norm
      EndIf
      Do idi=1, idi_N
         ObsDist=(idi-0.5)*ObsDist_Step !        ro=(id-0.5)*ObsDist_Step
         R_0=sqrt(ObsDist*ObsDist +zeta_0*zeta_0)
         t_o=NN_0*R_0+t_r !observer is summed over
         hto=h_0-t_o  ! =-(NN_0*R-zeta_0)
         dto=T_o-t_ch(idi)  ! =+t_tb(i_t)
         Zeta_Ch=z_ch(idi)
         If(dto.lt.0.) cycle
         Do i_t=2,nto
            If((dto-t_tb(i_t)).lt. 0.) Then ! Find point where (t_tb(i_t) about = dto)
               exit
            EndIf
         Enddo
         weight=2.*ObsDist*ObsDist_Step   !*(DistMax-ObsDist)
         !
         Sx_h = 0.
         S1x_h = 0.
         Do i_z=-N_z,N_z
            zeta_s=zeta_0+i_z*d_z
            i=(zeta_s/AtmHei_step)
            if((i.lt.1) .or. (i.ge.AtmHei_dim) ) cycle
            di=zeta_s/AtmHei_step-i
            NN_s=1.+ di*xi(i+1) + (1.-di)*xi(i)
            R_s=sqrt(zeta_s*zeta_s +ObsDist*ObsDist)
            h_s=h_0 + zeta_s - zeta_0 + NN_0* R_0 - NN_s * R_s
            If(h_s.le.0.) cycle
            Call SourceMoments(CD_i, zeta_s, h_s, fh,dfh, SEx, SAxD, SEy, SAyD, SAr, SErh)
            Sx_h = Sx_h + d_z*SEx/(NN_s*R_s)
            S1x_h = S1x_h + d_z*SAxD*dfh/(NN_s*R_s)  !  SAxD = Jx with the present hack
         EndDo ! i_z=-N_h,N_h
         !write(2,*) 'Sx_h', zeta_0 - N_z*d_z, zeta_0 + N_z*d_z, Sx_h, S1x_h
         !
         ! here was the now parked piece of code.
         Sx_h=Sx_h*weight/(NN_0 * R_0 )
         S1x_h=S1x_h*weight/(NN_0 * R_0 )
         Sx=Sx+Sx_h
         S1x=S1x+S1x_h
         !
         If((i_t.lt.nto) .and. (ObsDist.lt.DistMax) ) Then
            OD_max=ObsDist
            Ex= ((dto-t_tb(i_t-1))* Ex_tb(i_t,idi,CD_i) - (dto-t_tb(i_t))*Ex_tb(i_t-1,idi,CD_i) )/(t_tb(i_t)-t_tb(i_t-1)) ! this pulse interpolation is primitive
            Rx=Rx + Ex*weight/(NN_0*R_0)  ! radial increas i.s.o. fall-off for backtracking
            Sxc=Sxc+Sx_h
            S1xc=S1xc+S1x_h
            Kxc=Kxc+Kx_h
         EndIf
         !write(2,*) 'dto:', dto, ObsDist, Rx/norm, Ex/(NN_0*R_0), Ex, cald !, R_0 , NN_0*zeta_0 , dr*R_0*R_0 ! Ex,Ex_tb(i_t,idi,CD_i),Ex_tb(i_t-1,idi,CD_i), cald
         !
      EndDo ! idi=1,ObsDist_dim
      !stop 'Source_zeta'
      Call SourceMoments(CD_i, zeta_0, h_0, fh,dfh, SEx, SAxD, SEy, SAyD, SAr, SErh)
      If(norm.gt.0.) Then
         Rx =Rx/normc              ! 'experimental' backtrace source (finite frequency, limited antenna range
         Sxc= Sxc/normc          ! Theoretical Backtrace source, integrated over the 'experimental' domain of antenna distances
         Sx = Sx/norm            ! Theoretical Backtrace source, Integrated over a large antenna distances
         S1xc= S1xc/normc          ! Theoretical Backtrace source, integrated over the 'experimental' domain of antenna distances
         S1x = S1x/norm            ! Theoretical Backtrace source, Integrated over a large antenna distances
         Kxc= Kxc/normc          ! Theoretical Backtrace source, integrated over the 'experimental' domain of antenna distances
         Kx = Kx/norm            ! Theoretical Backtrace source, Integrated over a large antenna distances
         SEx=SEx * lambda*10./nrm  !  (R*R)  ! true source at this height
      EndIf
      !
      !write(2,*) 'Source_zeta; @h=',h_0, Rx,';', Sxc, Kxc, Sx, Kx,';', SEx, OD_max ! , dfh*Jx*nrm, fh *dJx*nrm
      write(2,*) 'Source_zeta; @h=',h_0, Rx,';', Sxc,  Sx, ';', SEx, OD_max, S1xc, S1x ! , dfh*Jx*nrm, fh *dJx*nrm
      write(4,*) h_0,', ', Rx,', ', Sxc,', ', S1xc,', ', SEx,', ', Kxc, ', '
      !stop 'Source_zeta1'
   EndDo ! i_h=1,N_h
   Close(Unit=4)
   Close(Unit=14)
   !stop 'Source_zeta2'
!    EndDo
    return
    !
         h_mx=0.
         h_mn=AtmHei_step*AtmHei_dim
         h_mn2=h_mn
         !d_z=AtmHei_step
         iz_lo=-zeta_Ch/d_z-1
         iz_hi=(h_mn-zeta_Ch)/d_z+1
         i_sign=-1
         Kx_h=0.
         Do i_z=iz_lo, iz_hi  ! park this piece of code, this way does not reall (or rather -really not-) work.
            zeta_s=zeta_Ch-0.5*d_z+i_z*d_z  ! symmetric arounf Cherenkov
            if(i_z.gt.0) i_sign=+1 ! to take care of sign flip of CalD at the cherenkov point
            i=(zeta_s/AtmHei_step)
            if((i.lt.1) .or. (i.ge.AtmHei_dim) ) cycle
            di=zeta_s/AtmHei_step-i
            NN_s=1.+ di*xi(i+1) + (1.-di)*xi(i)
            dr_s=   di*dxi(i+1) + (1.-di)*dxi(i)
            ddr_s=   (dxi(i+1) - dxi(i))/AtmHei_step
            R_s=sqrt(zeta_s*zeta_s +ObsDist*ObsDist)
            h_s=h_0 + zeta_s - zeta_0 + NN_0* R_0 - NN_s * R_s
            If(h_s.le.0.) cycle
            Call SourceMoments(CD_i, zeta_s, h_s, fh,dfh, SEx, SAxD, SEy, SAyD, SAr, SErh)
            CalD_s=NN_s*( R_s - NN_s*zeta_s - dr_s*R_s*R_s)
            dCalD_s=dr_s*(R_s - 3*NN_s*zeta_s - dr_s*R_s*R_s)+NN_s*(zeta_s/R_s-NN_s-ddr_s*R_s*R_s)
            Kx_h=Kx_h + d_z*i_sign*(fh*SAxD)*dCalD_s/(CalD_s*CalD_s)
            If(h_s.gt.h_mx) Then
               h_mx=h_s
               zeta_mx=zeta_s
               Kx_stk1=(fh*SAxD)/CalD_s
               h_mn2=h_s
               zeta_mn2=zeta_s
               Kx_stk2=(fh*SAxD)/CalD_s
            EndIf
            If(h_s.lt.h_mn2 .and. zeta_s.gt.zeta_mx) Then
               h_mn2=h_s
               zeta_mn2=zeta_s
               Kx_stk2=(fh*SAxD)/CalD_s
            EndIf
            If(h_s.lt.h_mn .and. zeta_s.le.zeta_Ch) Then
               h_mn=h_s
               zeta_mn=zeta_s
               Kx_stk0=fh/CalD_s
            EndIf
            !write(2,*) 'zeta,h,fh,calD',zeta_s,h_s,';', (fh*SAxD),CalD_s, dCalD_s
         EndDo ! i_z=-N_h,N_h
         !write(2,*) 'Kx_mx',h_mx, zeta_mx, Kx_stk1,';', Kx_h, Sx_h
         !write(2,*) zeta_Ch-0.5*d_z+iz_lo*d_z, zeta_Ch-0.5*d_z+iz_hi*d_z , Kx_h , zeta_Ch
         !write(2,*) 'Kx_mn', h_mn, zeta_mn, Kx_stk0,';', h_mn2, zeta_mn2, Kx_stk2,';', Kx_h, Sx_h,zeta_0-N_z*d_z, zeta_0+N_z*d_z
         Kx_h=Kx_h*weight/(NN_0 * R_0 )
         Kx=Kx+Kx_h
!
End Subroutine PencilSource_zeta
!------------------------------
!------------------------------
Subroutine WideSource_zeta(zeta_0, rs_0)
!   Calculate the source @ height=\zeta as function of tmh=(t-h) for a fixed observerdistance to this shower.
!   A back-beaming integration is performed over observer position using the electric fields at these antennas
!     or the source structure for this shower.
!   CD_i (CoreDist_A(CD_i)) determines the position of the pencil shower w.r.t. the core of the real shower which determines the pancake thickness.
!   zeta= height in the atmosphere of the back-trace voxel.
!   Shower time varies with h, distance behind the shower front.
!
    use Atmosphere, only : xi,dxi, AtmHei_dim, AtmHei_step
   use Pancakefunction, only : PancakeThi, CoreDist_Dim
    use BigArrays, only : t_ch, t_to, Ex_to,Ey_to,Er_to
    ! Ex_to(i=1,tTrace_dim_o;idi=1,ObsDist_dim) Field as function of time for different observer positions
    !    where t_to(i)=i*tTrace_step   !  starts after shouwer touching ground.
   !
   use BigArrays, only : Line2Core
    use BigArrays, only : ObsDist_dim, tTrace_dim_b, tTrace_dim_o, tTrace_step, ObsDist_Step
    use constants, only : dp, pi
   use eventdata, only : PlotBase
   Use Pancakefunction, only : Find_CD_Interpol
   use LateralDistribution, only : W_tc
    implicit none
!    real*8 :: Ex(0:CoreDist_Dim),AxD(0:CoreDist_Dim),Ey(0:CoreDist_Dim),AyD(0:CoreDist_Dim),Ar(0:CoreDist_Dim),nt_r,t_o
    real(dp), intent(in) :: zeta_0, rs_0  ! zeta_0 is the distance from ground-core, rs_0 is distance from shower axis of focus point
    real(dp) :: h_0,R_0, phi_0, NN_0,dr_0, d_phi0 ! all _0 mark the quantities at the focus point of retraced rays
    real(dp) :: rc_a, t_a ! all _0 mark the quantities at the focus point of retraced rays
    real(dp) ::  zeta_s, h_s, NN_s
   Real(dp) :: rc_r, dc_r, xc_r, rr_a, dr_a, xr_a, phi_r, d_pr, c_pr, lim_ca, R_r, h_r, CD_d
   Integer :: i_pr, ic_r, NC_r, ir_a
    real(dp) :: di, fh,dfh0,dfh1, weight, norm, normc, nrm, wei, Factr
    real(dp) :: Rx_h, Rx, Sx_h, Sx, Sxc, d_z, Ex,Ey, Erh, S1x, S1xc, S1x_h, dSx_h, dS1x_h
    Real(dp) :: SEx, S1Ex, SEx0, SEx1, SAxD0, SAxD1, SEy, SAyD, SAr, SErh
    integer :: i, idi, i_t, i_h, N_h, i_z, idi_N, i_phi, N_phir, N_phi0, Nz_lo,Nz_hi, CD_i
    real(dp) :: dhdz, lambda, d_h, OD_max, d_t, dwr, drs, rcr_max
    !integer, save :: il_start = 1, CD_i_prev=-1
    real(dp) :: DistMax=500.d0
    character(len=10) :: txt
    real(dp)  RW_z, RE_z

!
! lateral integral
! create grid at a fixed d0 as distance between core and observer
!
   call Find_CD_Interpol(rs_0,i,di)
   lambda=(1.-di)*PancakeThi(i)+ di*PancakeThi(i+1)
   Ex=0.   ; Ey=0. ; Erh=0.
   !
   lambda=5.
   write(2,*) 'wideSource_zeta:', zeta_0, rs_0, lambda
   !
   write(txt,"(I3.3,i2.2)") NINT(zeta_0/100.)
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(PlotBase)//'Wzeta_'//trim(txt)//'.csv')
   write(4,*) ' WideSource_zeta:', zeta_0, ' !'
   write(4,*) '!     h     Rx(backtrace antennas),  Sxc(MGMR, cut antenna domain),  Sxc(MGMR, full antenna domain), Source !'
   N_h=20
   !d_h=5.*lambda/N_h
   d_h=10./N_h
   d_z=d_h*20.
   d_z=10.
   !d_h=d_h/3.
   !N_h=N_h*3
   N_phir=9 ! same as for E-field calc
   N_phi0=1
   d_phi0=pi/N_phi0
   If((ObsDist_dim-1)*ObsDist_Step .lt. DistMax) DistMax=(ObsDist_dim-1)*ObsDist_Step
   !
   nrm=1.e6
   norm=0.
   normc=0.
   idi_N = ObsDist_dim-1
   NC_r=CoreDist_Dim*2+2
   rcr_max= Line2Core(NC_r)
   Do idi=1, idi_N
      rc_a=idi*ObsDist_Step !
      R_0=sqrt(rc_a*rc_a +zeta_0*zeta_0)
      weight=2.*rc_a*ObsDist_Step   !*(DistMax-ObsDist)
      norm=norm+weight/R_0 ! /(R_0*R_0)
      If(rc_a.lt.DistMax) normc=normc+weight/R_0 ! /(R_0*R_0)
   EndDo ! idi=1, idi_N
   norm=norm*nrm
   normc=normc*nrm
   RW_z=0.
   RE_z=0.
   !
   Do i_h=-N_h,N_h
   !Do i_h=1,1
      h_0=i_h*d_h
      !tr_0=h_0-zeta_0 ! should be negative, all quantities at the source are given
      !write(2,*) 'Source_zeta; z,h=',zeta_0,h
      !Flush(unit=2)
      i=(zeta_0/AtmHei_step)
      if(i.lt.1) return
      di=zeta_0/AtmHei_step-i
      NN_0=1.+ di*xi(i+1) + (1.-di)*xi(i)
      dr_0=   di*dxi(i+1) + (1.-di)*dxi(i)
      !
      Rx=0.
      OD_max=-1.
      Sx=0.
      Sxc=0.
      S1x=0.
      S1xc=0.
      !write(2,*) 'dto:', 'dto, ObsDist, Rx/norm, Ex/(NN*R), Ex, cald '
      Do idi=1, idi_N
      !Do idi=36, 36 ! 50
         rc_a=idi*ObsDist_Step ! Distance of antenna to shower axis
         !phi_0=0. ; i_phi=0  ; d_phi0=pi
         weight=2.*rc_a*ObsDist_Step*d_phi0   !*(DistMax-ObsDist)
         Do i_phi=1,N_phi0  ! angle between ObsDist and rs_0
            phi_0=(i_phi-0.5)*d_phi0
            If(N_phi0.eq.1) phi_0=0.
            R_0=sqrt(zeta_0*zeta_0+ rc_a*rc_a + rs_0*rs_0 - 2*rs_0*rc_a*cos(phi_0)) ! Distance between antenna and source in shower
            t_a=NN_0*R_0 + h_0-zeta_0 ! time in antenna is derived quantity
         !If(i_h.eq.0) write(2,*) 't_shift:', NN_0*R_0-zeta_0, zeta_0, rc_a, NN_0, R_0
            Factr=1.d0/(NN_0*R_0)
            If(t_a.lt.0.) cycle
            !
            i_t=t_a/tTrace_step
            If(i_t.ge.tTrace_dim_o) cycle
            Rx_h = 0.
            If(rc_a.lt.DistMax) Then
               d_t=t_a/tTrace_step - i_t  ! between 0 and 1
               OD_max=rc_a
               Ex= ((1.-d_t)* Ex_to(i_t,idi) + d_t*Ex_to(i_t+1,idi) )   ! this pulse interpolation is primitive
               Rx_h= Ex*weight*Factr  ! radial fall-off for backtracking
            !   write(2,*) '!Rx_h,Ex:', Rx_h, i_t,idi,d_t,Ex, Ex_to(i_t,idi), Ex_to(i_t+1,idi), tTrace_step
            EndIf
            !write(2,*) '!Rx_h:', Rx_h/normc, i_t,idi,d_t,Ex, t_a, NN_0*R_0, rc_a, rs_0, zeta_0,i_phi, tTrace_step
            Rx=Rx+Rx_h
            !write(2,*) 'dto:', dto, ObsDist, Rx/norm, Ex/(NN*R), Ex, cald !, R , NN*zeta , dr*R*R ! Ex,Ex_tb(i_t,idi,CD_i),Ex_tb(i_t-1,idi,CD_i), cald
            !
            Sx_h = 0.
            S1x_h = 0.
            Nz_lo=t_a/(xi(1)*d_z)-30
            If(Nz_lo.lt.1) Nz_lo=1
            Nz_hi=t_a/(xi(AtmHei_dim)*d_z)+1
            If(Nz_hi.gt. int(AtmHei_dim*AtmHei_step/d_z)) Nz_hi= int(AtmHei_dim*AtmHei_step/d_z)
            !write(2,*) 'Nz_lo, Nz_hi:',Nz_lo, Nz_hi, d_z, t_a-xi(1)*Nz_lo*d_z, Nz_lo*d_z, Nz_hi*d_z, t_a
            !write(2,*) 'rc_a',rc_a
            Do i_z=1,Nz_hi ! Integrate over source as required for obtaining the E-field, r_s(i_rs), phi_s
            !Do i_z=216,218 ! Integrate over source as required for obtaining the E-field, r_s(i_rs), phi_s
               zeta_s=i_z*d_z
!If(zeta_s.lt.6000) cycle   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!If(zeta_s.gt.6500) cycle   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               i=(zeta_s/AtmHei_step)
               if((i.lt.1) .or. (i.ge.AtmHei_dim) ) cycle
               di=zeta_s/AtmHei_step-i
               NN_s=1.+ di*xi(i+1) + (1.-di)*xi(i)
               ! Perform radial integral over different rays in shower, use notation where
               !  rc_a: Distance antenna to shower axis=core
               !  rc_r: Distance ray to core ; ic_r, dc_r, xc_r (pojection along rc_a direction)
               !  rr_a: Distance antenna to shower ray ; ir_a, dr_a, xr_a
               !  phi_r: angle between rc_a and rc_r or rc_a and rr_a depending on ray integration ; i_pr, d_pr
               !  R_r: Distance from antenna to position on ray at height zeta_s
               !  lim_ca=rc_a/2. half way point
               !
               lim_ca=rc_a/2.
               d_pr=pi/N_phir
               Do i_pr=1,N_phir       ! angle with antenna direction
                  phi_r=(i_pr-0.5)*d_pr
                  c_pr=cos(phi_r)
                  dSx_h=0.
                  dS1x_h=0.
                  Do ic_r=1,Nc_r     ! distance ray to shower axis, use predefined grid
                     rc_r=Line2Core(ic_r)
                     dc_r=(Line2Core(ic_r+1)-Line2Core(ic_r-1))/2.
                     xc_r=rc_r*c_pr
                     If(xc_r.gt. lim_ca) exit ! getting closer to antenna than core
                     If(Line2Core(ic_r+1)*c_pr .gt. lim_ca) then ! next step would move out of range, take care of proper area (approx)
                        dc_r=lim_ca/c_pr-(rc_r-dc_r/2.)  ! remaining distance along the center line betwee observer and core
                     endif
                     R_r=sqrt(zeta_s*zeta_s + rc_a*rc_a + rc_r*rc_r - 2*xc_r*rc_a)
                     h_r=t_a + zeta_s - NN_s * R_r
                     If(h_r.le.0.) cycle
                     call Find_CD_Interpol(rc_r,CD_i,CD_d)
                     Call SourceMoments(CD_i, zeta_s, h_r, fh,dfh0, SEx0, SAxD0, SEy, SAyD, SAr, SErh)
                     Call SourceMoments(CD_i+1, zeta_s, h_r, fh,dfh1, SEx1, SAxD1, SEy, SAyD, SAr, SErh)
                     SEx=(1.-CD_d)*SEx0 + CD_d*SEx1
                     S1Ex=(1.-CD_d)*dfh0*SAxD0 + CD_d*dfh1*SAxD1
                     wei=W_tc(rc_r,dwr)
                     dSx_h = dSx_h+ weight*Factr*d_z*SEx*wei*dc_r*d_pr/(NN_s*R_r)
                     dS1x_h = dS1x_h+ weight*Factr*d_z*S1Ex*wei*dc_r*d_pr/(NN_s*R_r)
                     !write(2,*) '!dSx_h,dic_r:', dSx_h, SEx, wei, h_r, rc_r,ic_r, NN_s*R_r,zeta_s, rc_a !,'fh', fh,dfh
                  EndDo ! ic_r=0,CoreDist_Dim-1CD_i
                  Sx_h = Sx_h + 2*dSx_h
                  S1x_h = S1x_h + 2*dS1x_h
!If(abs(dSx_h/normc).gt. 1.0)  write(2,*) '!dSx_h,ic_r:',i_z, i_pr, dSx_h/normc, zeta_s,h_r
                  !
                  dr_a=ObsDist_Step ! to be able to take care of rs=0
                  dSx_h=0.
                  dS1x_h=0.
                  Do ir_a=1, ObsDist_dim      ! distance ray to antenna, use regular grid
                     rr_a=(ir_a-0.5)*ObsDist_Step
                     xr_a=rr_a*c_pr
                     If(xr_a.gt. lim_ca) exit ! getting closer to core than antenna
                     If((ir_a+0.5)*ObsDist_Step*c_pr .gt. lim_ca) then ! next step would move out of range, take care of proper area (approx)
                        dr_a=lim_ca/c_pr-(ir_a-1)*ObsDist_Step  ! remaining distance along the center line betwee observer and core
                     endif
                     R_r=sqrt(zeta_s*zeta_s + rr_a*rr_a)
                     h_r=t_a + zeta_s - NN_s * R_r
                     If(h_r.le.0.) cycle
                     rc_r=sqrt(rc_a*rc_a + rr_a*rr_a - 2*xr_a*rc_a)  !
                     If(rc_r .gt. rcr_max) cycle
                     call Find_CD_Interpol(rc_r,CD_i,CD_d)
                     Call SourceMoments(CD_i, zeta_s, h_r, fh,dfh0, SEx0, SAxD0, SEy, SAyD, SAr, SErh)
                     Call SourceMoments(CD_i+1, zeta_s, h_r, fh,dfh1, SEx1, SAxD1, SEy, SAyD, SAr, SErh)
                     SEx=(1.-CD_d)*SEx0 + CD_d*SEx1
                     S1Ex=(1.-CD_d)*dfh0*SAxD0 + CD_d*dfh1*SAxD1
                     wei=W_tc(rc_r,dwr)*rr_a/rc_r
                     dSx_h = dSx_h+ weight*Factr*d_z*SEx*wei*dr_a*d_pr/(NN_s*R_r)
                     dS1x_h = dS1x_h+ weight*Factr*d_z*S1Ex*wei*dr_a*d_pr/(NN_s*R_r)
                     !write(2,*) '!dSx_h,ir_a:',i_z,ir_a,i_pr, dSx_h,zeta_s, h_r, SEx, SEx0, SEx1,CD_d
                  EndDo ! ic_r=1,CoreDist_Dim      ! distance to shower axis, use predefined grid
                  Sx_h = Sx_h + 2*dSx_h
                  S1x_h = S1x_h + 2*dS1x_h
                  !If(abs(dSx_h/normc).gt. 1.0) write(2,*) '!dSx_h,ir_a:',i_z,i_pr, dSx_h/normc,zeta_s,h_r
               EndDo ! i_pr=1,      ! angle with antenna direction
               !write(2,*) '!dSx_h,ir_a:',i_z,ir_a,i_pr, dSx_h,zeta_s, h_r, SEx, SEx0, SEx1,CD_d
            !If(abs(Sx_h/normc).gt. 1.0) write(2,*) '!Sx_h,rc_a:',i_z, Sx_h/normc,zeta_s, rc_a
            EndDo ! i_z=-N_h,N_h
            !write(2,*) '!Rx_h, Sx_h:', rc_a, i_phi,t_a,';', Rx_h/normc, Sx_h/normc,', normc=',normc
            Sx=Sx+Sx_h
            S1x=S1x+S1x_h
            !stop 'WSource_zeta0'
            If (rc_a.lt.DistMax)  Sxc=Sxc+Sx_h
            If (rc_a.lt.DistMax)  S1xc=S1xc+S1x_h
            !
         EndDo ! i_phi=1,N_phi0
         !
      EndDo ! idi=1,ObsDist_dim


      If(norm.gt.0.) Then
         Rx =Rx/normc              ! back-beamed source from calculated E fields (almost inite frequency, limited antenna range)
         Sxc= Sxc/normc          ! Theoretical Backtrace=beamed source, integrated over the 'experimental' domain of antenna distances, using currents
         ! should expect (Rx = Sxc)  within the intrinsic time resolution on 0.02 in the code. True for pos h, but need a fine Phi_r resolution in integral
         !  For neg h there appears to be a constant difference.
         !  This constant offset is due to a too sharp zeta-cutoff in "RayFields_t_obs", now corrected.
         S1xc= S1xc/normc          ! Theoretical Backtrace source, integrated over the 'experimental' domain of antenna distances, using only (df/dh*I_x, ignoring the second f*dI_x/dh)
         Sx = Sx/norm            ! Theoretical Backtrace source, Integrated over a large antenna distances
         S1x = S1x/norm            ! Theoretical Backtrace source, Integrated over a large antenna distances
         !SEx=SEx * lambda*10./nrm  !  (R*R)  ! true source at this height
      EndIf
      RE_z=RE_z+  Rx*Rx*d_h
      RW_z=RW_z+  Sxc*Sxc*d_h
      !
      SEx=0.
      write(2,*) 'WSource_zeta; @h=',h_0, Rx, Sxc, S1xc, OD_max ! , dfh*Jx*nrm, fh *dJx*nrm
      flush(Unit=2)
      write(4,*) h_0,', ', Rx,', ', Sxc,', ', Sx,', ', S1xc
      !If(i_h .eq. -N_h/2) stop 'WideSource_zeta-i_h'
   EndDo ! i_h=1,N_h
   Close(Unit=4)
   write(2,"(A,2(1pG13.4),A,2(1pG13.4))") 'sq-integrated E gives R^2, R:', RE_z, sqrt(RE_z), ', dt=',d_h
   write(2,"(A,2(1pG13.4),A,2(1pG13.4))") 'sq-integrated Weight gives R^2, R:', RW_z, sqrt(RW_z), ', dt=',d_h
   !stop 'WSource_zeta2'
!    EndDo
    return
End Subroutine WideSource_zeta
!------------------------------
