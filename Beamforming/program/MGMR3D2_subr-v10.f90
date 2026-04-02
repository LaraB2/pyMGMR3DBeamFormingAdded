!------------------------------
Subroutine LateralInt(idi)
   use BigArrays, only : t_tb, Ex_tb,AxD_tb,Ey_tb,AyD_tb,Ar_tb, Erh_tb,t_ch
   use BigArrays, only : Ex_spld,AxD_spld,Ey_spld,AyD_spld,Ar_spld, Erh_spld
   use BigArrays, only : t_to, Ex_to,Ey_to,Er_to
   use BigArrays, only : ObsDist_dim, tTrace_dim_b, tTrace_dim_o, ObsDist_Step
   use BigArrays, only : Line2Core
   use Pancakefunction, only : CoreDist_Dim
   use constants, only : ci,pi,dp
   Use Pancakefunction, only : Find_CD_Interpol
   use LateralDistribution, only : W_tc
   implicit none
   integer, intent(in) :: idi
   integer :: id,CD_i
   Integer :: i_pr, ic_r, NC_r, ir_a, N_phi
   Real(dp) :: rc_a, rc_r, dc_r, xc_r, rr_a, dr_a, xr_a, phi_r, d_pr, c_pr, lim_ca
   real(dp)  :: gth, Wei,weid,rcr_max,did
   real(dp)  :: Int,d0,int2,intx,wrs,wro,dwr,CD_d,tshft ! ,CoreDist,r_s
   real(dp) :: DEx(tTrace_dim_o),DAxD(tTrace_dim_o),DEy(tTrace_dim_o),DAyD(tTrace_dim_o), &
         DAr(tTrace_dim_o),DErh(tTrace_dim_o)
   !
   ! lateral integral
   ! create grid at a fixed d0 as distance between core and observer
   rc_a= idi*ObsDist_Step
   lim_ca=rc_a/2.
   int=0.  ; intx=0. ;  int2=0.
   NC_r=CoreDist_Dim*2+2
   rcr_max= Line2Core(NC_r)
   !   Nth=15  ; dth=pi/Nth
   N_phi=9
   d_pr=pi/N_phi
   gth=2.
   Do i_pr=1,N_phi       ! angle with antenna direction
      phi_r=(i_pr-0.5)*d_pr
      c_pr=cos(phi_r)
      Do ic_r=1,Nc_r     ! distance ray to shower axis, use predefined grid
         rc_r=Line2Core(ic_r)
         dc_r=(Line2Core(ic_r+1)-Line2Core(ic_r-1))/2.
         xc_r=rc_r*c_pr
         If(xc_r.gt. lim_ca) exit ! getting closer to antenna than core
         If(Line2Core(ic_r+1)*c_pr .gt. lim_ca) then ! next step would move out of range, take care of proper area (approx)
            dc_r=lim_ca/c_pr-(rc_r-dc_r/2.)  ! remaining distance along the center line betwee observer and core
         endif
         rr_a=sqrt(rc_a*rc_a + rc_r*rc_r - 2*xc_r*rc_a)
         wei=dc_r*d_pr*W_tc(rc_r,dwr)
         weid=dc_r*d_pr*dwr
         !
         id=IDINT(rr_a/ObsDist_Step + 0.5)  ;  did=rr_a/ObsDist_Step+0.5-id
         if(id.ge.ObsDist_dim) exit
         if(id.lt. 1) then
          id=1  ; did=0.
         endif
         call Find_CD_Interpol(rc_r,CD_i,CD_d)
         !
         tshft=((1.-did)*t_ch(id)+did*t_ch(id+1))   ! need to apply time-shift
         call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Ex_tb(1,id,CD_i),   Ex_spld(1,id,CD_i),    Ex_tb(1,id,CD_i+1),   Ex_spld(1,id,CD_i+1), &
            Ex_tb(1,id+1,CD_i), Ex_spld(1,id+1,CD_i),  Ex_tb(1,id+1,CD_i+1), Ex_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DEx )
         call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            AxD_tb(1,id,CD_i),   AxD_spld(1,id,CD_i),    AxD_tb(1,id,CD_i+1),   AxD_spld(1,id,CD_i+1), &
            AxD_tb(1,id+1,CD_i), AxD_spld(1,id+1,CD_i),  AxD_tb(1,id+1,CD_i+1), AxD_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DAxD )
         Ex_to(:,idi)=Ex_to(:,idi) + gth*wei *DEx(:) - gth*DAxD(:)*weid
         call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Ey_tb(1,id,CD_i),   Ey_spld(1,id,CD_i),    Ey_tb(1,id,CD_i+1),   Ey_spld(1,id,CD_i+1), &
            Ey_tb(1,id+1,CD_i), Ey_spld(1,id+1,CD_i),  Ey_tb(1,id+1,CD_i+1), Ey_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DEy )
         call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            AyD_tb(1,id,CD_i),   AyD_spld(1,id,CD_i),    AyD_tb(1,id,CD_i+1),   AyD_spld(1,id,CD_i+1), &
            AyD_tb(1,id+1,CD_i), AyD_spld(1,id+1,CD_i),  AyD_tb(1,id+1,CD_i+1), AyD_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DAyD )
         Ey_to(:,idi)=Ey_to(:,idi) + gth*wei *DEy(:) - gth*DAyD(:)*weid
         call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Ar_tb(1,id,CD_i),   Ar_spld(1,id,CD_i),    Ar_tb(1,id,CD_i+1),   Ar_spld(1,id,CD_i+1), &
            Ar_tb(1,id+1,CD_i), Ar_spld(1,id+1,CD_i),  Ar_tb(1,id+1,CD_i+1), Ar_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DAr )
         call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Erh_tb(1,id,CD_i),   Erh_spld(1,id,CD_i),    Erh_tb(1,id,CD_i+1),   Erh_spld(1,id,CD_i+1), &
            Erh_tb(1,id+1,CD_i), Erh_spld(1,id+1,CD_i),  Erh_tb(1,id+1,CD_i+1), Erh_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DErh )
         Er_to(:,idi)=Er_to(:,idi) + gth*(DAr(:)*weid + DErh(:)*wei)*c_pr
      enddo     ! ic_r=1,Nc_r     ! distance ray to shower axis, use predefined grid
      !
      dr_a=ObsDist_Step ! to be able to take care of rs=0
      Do ir_a=1, ObsDist_dim      ! distance ray to antenna, use regular grid
         rr_a=(ir_a-0.5)*ObsDist_Step
         xr_a=rr_a*c_pr
         If(xr_a.gt. lim_ca) exit ! getting closer to core than antenna
         If((ir_a+0.5)*ObsDist_Step*c_pr .gt. lim_ca) then ! next step would move out of range, take care of proper area (approx)
            dr_a=lim_ca/c_pr-(ir_a-1)*ObsDist_Step  ! remaining distance along the center line betwee observer and core
         endif
         !R_r=sqrt(zeta_s*zeta_s + rc_a*rc_a + rr_a*rr_a - 2*xr_a*rc_a)
         !R_r=sqrt(zeta_s*zeta_s + rr_a*rr_a)
         !h_r=t_a + zeta_s - NN_s * R_r
         !If(h_r.le.0.) cycle
         rc_r=sqrt(rc_a*rc_a + rr_a*rr_a - 2*xr_a*rc_a)  !
         If(rc_r .gt. rcr_max) cycle
        wei=dr_a*d_pr*W_tc(rc_r,dwr)*rr_a/rc_r  !
        weid=dr_a*d_pr*dwr*rr_a/rc_r          !
         call Find_CD_Interpol(rc_r,CD_i,CD_d)
!
        id=ir_a
        tshft=t_ch(id)
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Ex_tb(1,id,CD_i), Ex_spld(1,id,CD_i), Ex_tb(1,id,CD_i+1), Ex_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DEx )
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            AxD_tb(1,id,CD_i),   AxD_spld(1,id,CD_i),    AxD_tb(1,id,CD_i+1),   AxD_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DAxD )
        Ex_to(:,idi)=Ex_to(:,idi) + gth*wei *DEx(:) - gth*DAxD(:)*weid
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Ey_tb(1,id,CD_i),   Ey_spld(1,id,CD_i),    Ey_tb(1,id,CD_i+1),   Ey_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DEy )
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            AyD_tb(1,id,CD_i),   AyD_spld(1,id,CD_i),    AyD_tb(1,id,CD_i+1),   AyD_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DAyD )
        Ey_to(:,idi)=Ey_to(:,idi) + gth*wei *DEy(:) - gth*DAyD(:)*weid
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Ar_tb(1,id,CD_i),   Ar_spld(1,id,CD_i),    Ar_tb(1,id,CD_i+1),   Ar_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DAr )
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Erh_tb(1,id,CD_i),   Erh_spld(1,id,CD_i),    Erh_tb(1,id,CD_i+1),   Erh_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DErh )
        Er_to(:,idi)=Er_to(:,idi) + gth*(DAr(:)*weid + DErh(:)*wei)*(rc_a-xr_a)/rc_r
      enddo     ! id
!      write(2,*) irs,ith,id,did,wei,Ex_nui((inui+inum)/2,idi)
    enddo       ! ith
!    write(2,*) 'Int W(r)=',Ex_nui(inui,idi),Ex_nui(inum,idi)
!    write(2,*) 'd0=',d0,int,int2,int+int2
    ! stop
    return
    end subroutine LateralInt
!-------------------
!------------------------------
    subroutine LateralInt_x(idi)  ! obsolete
    use BigArrays, only : t_tb, Ex_tb,AxD_tb,Ey_tb,AyD_tb,Ar_tb, Erh_tb,t_ch
    use BigArrays, only : Ex_spld,AxD_spld,Ey_spld,AyD_spld,Ar_spld, Erh_spld
    use BigArrays, only : t_to, Ex_to,Ey_to,Er_to
    use BigArrays, only : ObsDist_dim, tTrace_dim_b, tTrace_dim_o, ObsDist_Step
    use BigArrays, only : Line2Core
   use Pancakefunction, only : CoreDist_Dim
   use constants, only : ci,pi,dp
   Use Pancakefunction, only : Find_CD_Interpol
   use LateralDistribution, only : W_tc
    implicit none
    integer, intent(in) :: idi
    integer :: Nrs,irs,ith,nth,id,CD_i
    real(dp)  :: dth,theta,drs,rs,y_s,x_s,x_o,ro,Wei,weid,rs_max,gth,did
    real(dp)  :: Int,d0,int2,intx,wrs,wro,dwr,CD_d,tshft ! ,CoreDist,r_s
    real(dp) :: DEx(tTrace_dim_o),DAxD(tTrace_dim_o),DEy(tTrace_dim_o),DAyD(tTrace_dim_o), &
        DAr(tTrace_dim_o),DErh(tTrace_dim_o)
!
! lateral integral
! create grid at a fixed d0 as distance between core and observer
    d0= idi*ObsDist_Step
    int=0.  ; intx=0. ;  int2=0.
    nrs=CoreDist_Dim*2+2
    rs_max= Line2Core(nrs)
 !   Nth=15  ; dth=pi/Nth
    Nth=7  ; dth=pi/Nth
    gth=1.
    Do ith=0,Nth    ! for theta=0 point lies inbetween sh-core and observer
      theta=ith*dth
      do irs = 1, nrs   ! distance from shower core; CoreDist_A(0,CoreDist_Dim)
        rs=Line2Core(irs)
        drs=(Line2Core(irs+1)-Line2Core(irs-1))/2.
        !rs=(irs-0.5)*drs
        x_s=rs*cos(theta)
        if(x_s.gt. d0/2.) exit
        wrs=drs
        !wrs=0.     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        If((rs+drs)*cos(theta) .gt. d0/2.) then
            wrs=d0/(2.*cos(theta))-(rs-drs/2.)  ! remaining distance to the center line betwee observer and core
        endif
        wei=wrs*dth*W_tc(rs,dwr)
        weid=wrs*dth*dwr
        y_s=rs*sin(theta)
        ro=sqrt(y_s*y_s+(x_s-d0)*(x_s-d0))
        !
        id=IDINT(ro/ObsDist_Step + 0.5)  ;  did=ro/ObsDist_Step-(id-0.5)
        if(id.ge.ObsDist_dim) exit
        if(id.lt. 1) then
          id=1  ; did=0.
        endif
        call Find_CD_Interpol(rs,CD_i,CD_d)
        !
        tshft=((1.-did)*t_ch(id)+did*t_ch(id+1))   ! need to apply time-shift
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Ex_tb(1,id,CD_i),   Ex_spld(1,id,CD_i),    Ex_tb(1,id,CD_i+1),   Ex_spld(1,id,CD_i+1), &
            Ex_tb(1,id+1,CD_i), Ex_spld(1,id+1,CD_i),  Ex_tb(1,id+1,CD_i+1), Ex_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DEx )
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            AxD_tb(1,id,CD_i),   AxD_spld(1,id,CD_i),    AxD_tb(1,id,CD_i+1),   AxD_spld(1,id,CD_i+1), &
            AxD_tb(1,id+1,CD_i), AxD_spld(1,id+1,CD_i),  AxD_tb(1,id+1,CD_i+1), AxD_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DAxD )
        Ex_to(:,idi)=Ex_to(:,idi) + gth*wei *DEx(:) - gth*DAxD(:)*weid
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Ey_tb(1,id,CD_i),   Ey_spld(1,id,CD_i),    Ey_tb(1,id,CD_i+1),   Ey_spld(1,id,CD_i+1), &
            Ey_tb(1,id+1,CD_i), Ey_spld(1,id+1,CD_i),  Ey_tb(1,id+1,CD_i+1), Ey_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DEy )
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            AyD_tb(1,id,CD_i),   AyD_spld(1,id,CD_i),    AyD_tb(1,id,CD_i+1),   AyD_spld(1,id,CD_i+1), &
            AyD_tb(1,id+1,CD_i), AyD_spld(1,id+1,CD_i),  AyD_tb(1,id+1,CD_i+1), AyD_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DAyD )
        Ey_to(:,idi)=Ey_to(:,idi) + gth*wei *DEy(:) - gth*DAyD(:)*weid
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Ar_tb(1,id,CD_i),   Ar_spld(1,id,CD_i),    Ar_tb(1,id,CD_i+1),   Ar_spld(1,id,CD_i+1), &
            Ar_tb(1,id+1,CD_i), Ar_spld(1,id+1,CD_i),  Ar_tb(1,id+1,CD_i+1), Ar_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DAr )
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Erh_tb(1,id,CD_i),   Erh_spld(1,id,CD_i),    Erh_tb(1,id,CD_i+1),   Erh_spld(1,id,CD_i+1), &
            Erh_tb(1,id+1,CD_i), Erh_spld(1,id+1,CD_i),  Erh_tb(1,id+1,CD_i+1), Erh_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DErh )
        Er_to(:,idi)=Er_to(:,idi) + gth*(DAr(:)*weid + DErh(:)*wei)*cos(theta)
      enddo     ! irs
!
      Do id = 1, ObsDist_dim   ! distance from observer in direction of shower core for theta=0.
        ro=(id-0.5)*ObsDist_Step
        x_o=ro*cos(theta)
        if(x_o.gt. d0/2.) exit
        wro=ObsDist_Step
        !wro=0. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        If((ro+ObsDist_Step)*cos(theta) .gt. d0/2.) then
            wro=d0/(2.*cos(theta))-(ro-ObsDist_Step/2.)  ! remaining distance to the center line betwee observer and core
        endif
        x_s=d0-x_o  ;  y_s=ro*sin(theta)  ; rs=sqrt(x_s*x_s+y_s*y_s)
        if(rs.gt. rs_max) cycle
        wei=wro*dth*W_tc(rs,dwr)*ro/rs  ! replaced wrs by wro
        weid=wro*dth*dwr*ro/rs          ! replaced wrs by wro
!
        call Find_CD_Interpol(rs,CD_i,CD_d)
        tshft=t_ch(id)
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Ex_tb(1,id,CD_i), Ex_spld(1,id,CD_i), Ex_tb(1,id,CD_i+1), Ex_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DEx )
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            AxD_tb(1,id,CD_i),   AxD_spld(1,id,CD_i),    AxD_tb(1,id,CD_i+1),   AxD_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DAxD )
        Ex_to(:,idi)=Ex_to(:,idi) + gth*wei *DEx(:) - gth*DAxD(:)*weid
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Ey_tb(1,id,CD_i),   Ey_spld(1,id,CD_i),    Ey_tb(1,id,CD_i+1),   Ey_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DEy )
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            AyD_tb(1,id,CD_i),   AyD_spld(1,id,CD_i),    AyD_tb(1,id,CD_i+1),   AyD_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DAyD )
        Ey_to(:,idi)=Ey_to(:,idi) + gth*wei *DEy(:) - gth*DAyD(:)*weid
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Ar_tb(1,id,CD_i),   Ar_spld(1,id,CD_i),    Ar_tb(1,id,CD_i+1),   Ar_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DAr )
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Erh_tb(1,id,CD_i),   Erh_spld(1,id,CD_i),    Erh_tb(1,id,CD_i+1),   Erh_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DErh )
        Er_to(:,idi)=Er_to(:,idi) + gth*(DAr(:)*weid + DErh(:)*wei)*(x_s/rs)
      enddo     ! id
!      write(2,*) irs,ith,id,did,wei,Ex_nui((inui+inum)/2,idi)
      gth=2. ; if(ith.eq.(Nth-1)) gth=1.
    enddo       ! ith
!    write(2,*) 'Int W(r)=',Ex_nui(inui,idi),Ex_nui(inum,idi)
!    write(2,*) 'd0=',d0,int,int2,int+int2
    ! stop
    return
    end subroutine LateralInt_x
!-------------------
subroutine Get_ObsPlsDelay(idi)
    use BigArrays, only : t_to, Ex_to,Ey_to,Er_to
    use BigArrays, only : tTrace_dim_o, tTrace_step, ObsPlsTime
    use constants, only : dp
    use RFootPars, only : test
    implicit none
    integer, intent(in) :: idi
    integer :: i,j,k,nth,id,inu,CD_i,m(1)
    real(dp) :: Overlap(100), Norma, Normb
    !
    k=ObsPlsTime(idi-1)/tTrace_step
    If((k+300).ge.tTrace_dim_o) then
        write(2,*) 'ERROR, time-trace too short to contain full peak at distance#',idi
        stop 'ERROR, time-trace too short'
    endif
    Do i=1,100
        Overlap(i)=0.0
        Norma=0.
        Normb=0.
        Do j=1,200
            Overlap(i)=Overlap(i) + Ex_to(k+j,idi-1)*Ex_to(k+j+i,idi) + Ey_to(k+j,idi-1)*Ey_to(k+j+i,idi)
            Norma = Norma + Ex_to(k+j,idi-1)*Ex_to(k+j,idi-1) + Ey_to(k+j,idi-1)*Ey_to(k+j,idi-1)
            Normb = Normb + Ex_to(k+j+i,idi)*Ex_to(k+j+i,idi) + Ey_to(k+j+i,idi)*Ey_to(k+j+i,idi)
        enddo
        Overlap(i)=Overlap(i)/sqrt(Norma*Normb)
    enddo
    m = maxloc(Overlap)
    !if(test) then
    !    write(2,*) idi,m
    !    write(2,"(20f6.3)") Overlap
    !endif
    ObsPlsTime(idi)=ObsPlsTime(idi-1) + m(1)*tTrace_step
end subroutine Get_ObsPlsDelay
!------------------------------
!------------------------------
subroutine getZeta(lamb,zetai,zeta)
   ! lamb is a measure for the distance behind the shower front (difference with Cherenkov height)
   ! zeta is the emission height in the atmosphere (including the distance behind the shower front)
   use Atmosphere, only : AtmHei_dim, AtmHei_step
   implicit none
   real*8 :: zetaa,zetab,lamb,da,db,zetai(2),zeta(2),zet
   integer :: i,k, branch
   common / time / zeta_c,hcto_c
      real*8 :: zeta_c,hcto_c
   !
   !       lower branch starting from the Cherenkov point
   branch=-1
   !write(2,*) 'getZeta-lower branch:',k,zet,lamb
   call GetZet_Lam(zet,lamb,branch)
   zeta(1)=zet
   !
   !       upper branch
   branch=+1
   !write(2,*) 'getZeta-upper branch:',k,zet,lamb
   call GetZet_Lam(zet,lamb,branch)
   zeta(2)=zet
   !
   return
end subroutine getZeta
!------------------------------
subroutine GetZet_Lam(zet,lam,branch)
    ! Solves for \zeta:  n\sqrt{(-\beta t +h)^2 + (1-\beta^2 n^2)d^2} =n(R-n\zeta)
    ! \zeta = (h + n \sqrt{h^2+(1-n^2)d^2})/(1-n^2)
    ! given (-t+h)= hcto_c - lam where lam is input parameter
    ! input: value of zet is used as start value for the search between values zeta and zetb
    ! input: da and db are the differences between lama, lamb values for zeta and zetb with lam
    use Atmosphere, only : xi,dxi, AtmHei_dim, AtmHei_step
    use BigArrays, only : ObsDist
    implicit none
    real*8, intent(inout) :: zet
    real*8, intent(in) :: lam
    integer, intent(in) :: branch  ! should be +1 or -1
    real*8 :: dd,NN,di  ! ,dr
    real*8, parameter :: Converg=1.e-6
    integer, parameter :: MaxIter=4
    integer :: i,k, brnch
    common / time / zeta_c,hcto_c
      real*8 :: zeta_c,hcto_c
    Real*8 zetx, hxto, root
    Real*8 dr,R,calD,hmto,h0,lamx,dsqldz
!
    hxto=hcto_c - lam
    zetx=zet
    If(branch.gt.0) then
      brnch=+1
    Else
      brnch=-1
    EndIf
    !write(2,*) 'GetZet_Lam: calling=', zet,lam,brnch
    Do k=1,MaxIter  ! find value of `zet' that zeros `DD' for a given value of `lam', all have units of length
      i=zetx/AtmHei_step
      if(i.lt.1) then
         i=1
         di=0.
      else if(i.lt.AtmHei_dim) then
        di=zetx/AtmHei_step-i
      else
        i=AtmHei_dim-1
        di=1.
      endif
      NN=1.+ di*xi(i+1) + (1.-di)*xi(i)
      R=sqrt(zetx*zetx+ObsDist*ObsDist)
      !hxto=zetx-NN*R
      !dhdz=-(NN*R-zetx)*(R-NN*zetx - dr*R*R)/(R*hxto)
      ! branch*dhdz/sqrt(hcto_c-hxto)
      root=(hxto*hxto + (1-nn*nn)*ObsDist*ObsDist)
      If(root.lt.0) root=0.
      zet= (hxto - brnch*nn*sqrt(root) )/(1-nn*nn)
      !write(2,*) 'GetZet_Lam:', k,i,nn, zet, zetx, root
      !
      !checking retarded distance
      !
      If(abs(zet-zetx).lt.converg) then
         zetx=zet  ! done for later convergence testing
         exit
      Else If(abs(zet-zetx) .lt. 5*AtmHei_step) then  ! factor is set relatively arbitrarily
         zetx=(zetx+zet)/2
      Else
         zetx=zet
      End If
    enddo
    i=zetx/AtmHei_step
    if(i.lt.(AtmHei_dim-1) .and. i.gt.1) then
      di=zetx/AtmHei_step-i
      NN=1.+ di*xi(i+1) + (1.-di)*xi(i)
      R=sqrt(zetx*zetx+ObsDist*ObsDist)
      hxto=zetx- NN*R
      lamx=hcto_c - hxto
      DD=lam -lamx
      If(ABS(dd).gt. 1.e-6) Then  ! refine using derivative
         !write(2,*) 'GetZet_Lam: lam=', lam,'zeta=', zet, zetx, DD, '=lam-lamx'
         dr=   di*dxi(i+1)+(1.-di)*dxi(i)  ! change in refractivity
         dsqldz=-(NN*R-zetx)*(R-NN*zetx - dr*R*R)/(2.*R*hxto*sqrt(lamx))  ! d\sqrt(h_c-h)/d\zeta
         zet=zetx + (sqrt(lamx)-sqrt(lam))/dsqldz   ! refine using derivative in sqrt(lam)
         !write(2,*) sqrt(lamx), sqrt(lam), dsqldz, (sqrt(lamx)-sqrt(lam))/dsqldz, zetx, zeta_c
         i=zet/AtmHei_step
         if(i.lt.(AtmHei_dim-1) .and. i.gt.1) then
            di=zet/AtmHei_step-i
            NN=1.+ di*xi(i+1) + (1.-di)*xi(i)
            R=sqrt(zet*zet+ObsDist*ObsDist)
            lamx=hcto_c - (zet- NN*R)
            If(ABS(lam -lamx).gt. abs(DD)) Then
               !write(2,*) 'GetZet_Lam, itr:lam,x,diff=',lam, lamx, lam -lamx, 'old:',DD &
               !   ,'; zeta,x,c=', zet, zetx, zeta_c,  dsqldz, (AtmHei_dim-1.)*AtmHei_step
               Zet=zetx ! new value worse, do not accept derivative change
            EndIf
         Else  ! new value out of bound, do not accept derivative change
            Zet=zetx ! new value out of bound, do not accept derivative change
         EndIf
      EndIf
    EndIf
    !
end subroutine GetZet_Lam
!------------------------------
!------------------------------
subroutine AnalyzeGeometry(step)
   ! checks the expression for dh/d\zeta
   use Atmosphere, only : xi,dxi, AtmHei_dim, AtmHei_step
   use BigArrays, only : ObsDist
   implicit none
   real*8, intent(in) :: step
   !real*8, intent(in) :: lam
   real*8 :: NN,di  ! ,dr
   integer :: i,k, branch
   common / time / zeta_c,hcto_c
   real*8 :: zeta_c,hcto_c
   Real*8 hxto,zetx
   Real*8 dr,R, dhdz,root, calD
   !
   Do k=-10,9
      zetx=zeta_c+(k+0.5)*step
      branch=-1
      If(zetx .lt. zeta_c) branch=+1
      i=zetx/AtmHei_step
      if(i.lt.1) then
         cycle
      else if(i.lt.AtmHei_dim) then
        di=zetx/AtmHei_step-i
      else
        cycle
      endif
      NN=1.+ di*xi(i+1) + (1.-di)*xi(i)
      dr=   di*dxi(i+1)+(1.-di)*dxi(i)  ! change in refractivity
      R=sqrt(ObsDist*ObsDist+zetx*zetx)
      hxto=zetx-NN*R
      dhdz=-(NN*R-zetx)*(R-NN*zetx - dr*R*R)/(R*hxto)
      cald=NN*(R - NN*zetx - dr*R*R) ! retarded distance, including the needed factor n=NN
      !
      write(2,*) 'Investigate \zeta dependence:',k,zetx,hcto_c-hxto &
         ,dhdz,branch*dhdz/sqrt(hcto_c-hxto), dhdz/cald, cald
   End Do
   !
   !
end subroutine AnalyzeGeometry
!------------------------------
subroutine dhdzeta(lam,zet,dhdz)
   ! checks the expression for dh/d\zeta
   use Atmosphere, only : xi,dxi, AtmHei_dim, AtmHei_step
   use BigArrays, only : ObsDist
   implicit none
   real*8, intent(in) :: zet
   real*8, intent(in) :: lam
   real*8 :: NN,di  ! ,dr
   integer :: i
   common / time / zeta_c,hcto_c
   real*8 :: zeta_c,hcto_c
   Real*8 hxto
   Real*8 dr,R, dhdz
!
   dhdz=0.
   hxto=hcto_c - lam
   i=zet/AtmHei_step
   if((i.lt.1) .or. (i.ge.AtmHei_dim) ) return
   di=zet/AtmHei_step-i
   NN=1.+ di*xi(i+1) + (1.-di)*xi(i)
   dr=   di*dxi(i+1)+(1.-di)*dxi(i)  ! change in refractivity
   !
   R=sqrt(zet*zet+ObsDist*ObsDist)
   dhdz=-(NN*R-zet)*(R-NN*zet - dr*R*R)/(R*hxto)
   !
   !
end subroutine dhdzeta
!------------------------------
!------------------------------
  subroutine FindZetaC(k,zeta,hcto)
!
!     In general: zeta= negative retarded time
!     Zeta_C= the zeta (=height from which signal is emitted) for which the retarded distance CalD equals zero
!     Cald is equal to the the denominator of eq 34, and may be different from eq 35 by a factor n.
!     The programmed form for CalD rather follows eq 36 where aslo the change in index of refraction with height is accounted for
!  This subroutine is called at line 124 of      Subroutine MGMR3D
!     where k is index in atmospheric height for which
!            T_obs(i)=-z+sqrt(ObsDist*ObsDist+z*z)*(1.+xi(i))  ! z=i*AtmHei_step  (i=atm height)
!     is minimal (i=k), i.e. close to the cherenkov distance
!
!
   use Atmosphere, only : AtmHei_dim, AtmHei_step
    implicit none
    real*8 :: zeta1,zeta2,calDa,calDb,calD,zeta,hcto ! ,N_s
    integer :: i,k
    real*8, parameter :: Converg=1.e-10
!
    if(k.ge.AtmHei_dim) then
! at large height $xi*zeta=constant$ and $xi'=-xi*zeta/zeta^2$
      zeta1=(AtmHei_dim-1)*AtmHei_step
      zeta2=10*AtmHei_dim*AtmHei_step+5.
    else
      zeta1=(k-1)*AtmHei_step
      zeta2=(k+1)*AtmHei_step
    endif
    call GetCalD(zeta1,calDa,hcto)
    call GetCalD(zeta2,calDb,hcto)
    if(k.ge.AtmHei_dim) then
      if(calDa*calDb .gt. 0.) then
        zeta=zeta2
        return
      endif
    endif
    Do i=1,10
      zeta=(zeta1*calDb-zeta2*calDa)/(calDb-calDa)
      call GetCalD(zeta,calD,hcto)
     if(calDa*calD.le. 0.) then
        zeta2=zeta
        calDb=calD
      else if(calDb*calD.le. 0.) then
        zeta1=zeta
        calDa=calD
      else
        write(2,*) 'FindZetaC:',zeta1,calDa,zeta2,calDb,zeta,hcto,calD
        exit
      endif
    enddo
    if(k.ge.AtmHei_dim) write(2,*) 'z',zeta,cald
    return
  end subroutine FindZetaC
!------------------------------
    subroutine GetCalD(zeta,calD,hcto)
    use Atmosphere, only : xi,dxi, PenDepth, IndRefOrho, AtmHei_step, AtmHei_dim, TopAtmExpon
    use BigArrays, only : ObsDist
    use constants, only : dp
    implicit none
    real(dp) :: zeta,calD,hcto,NN,R2,di,dr,h0  ! ,hcto2
    integer :: i
    i=zeta/AtmHei_step
    di=zeta/AtmHei_step-i
    if(i.ge.AtmHei_dim) then
      if(i.ge. AtmHei_dim*2) then
        hcto=-IndRefOrho*PenDepth(0)  ! (IndRefOrho == Index of refraction over density)(PenDepth == penetration depth in [g/cm^2])
        calD=ObsDist*ObsDist/(2.*zeta)
!  N.B. at d=140 m there is an inflection point (at 12+ km height) in t-ret v.s. t-obs which causes a 3rd branch.
        return
      endif
!      NN=1.+rh0 * (1./zeta)* h0*(1.- exp(-zeta/h0)) ! mean index of refraction along path
      h0= TopAtmExpon
      NN= 1. + IndRefOrho*(PenDepth(0)-PenDepth(AtmHei_dim)*exp(-(zeta-AtmHei_dim*AtmHei_step)/h0))/zeta
      dr=IndRefOrho * (-PenDepth(0)+PenDepth(AtmHei_dim)*(zeta/h0+1)*exp(-(zeta-AtmHei_dim*AtmHei_step)/h0))/(zeta*zeta)   !  ((h0+zeta)*exp(-zeta/h0) - h0)/(zeta*zeta)
    else
      NN=1.+ di*xi(i+1)+(1.-di)*xi(i)   ! Mean index of refraction from height zeta to ground
      dr=   di*dxi(i+1)+(1.-di)*dxi(i)  ! change in refractivity
    endif
    R2=ObsDist*ObsDist+zeta*zeta
    calD=sqrt(R2) - NN*zeta- dr*R2 ! corresponds to eq 35 of the notes where calD=D/n and a term is added to account for varying inderx of refraction
    hcto=-sqrt(dr*R2*dr*R2 + (NN*NN-1.)*ObsDist*ObsDist)  !=-\sqrt{(R-N\zeta)^2+(N^2-1) d^2}
!    hcto2=-sqrt( (sqrt(R2)-NN*zeta)**2+(NN**2-1.)*d*d)
!    write(2,*) 'hcto-hcto2=',hcto-hcto2
    !write(2,*) 'GetCalD:',calD,zeta, hcto, 'h0to - hcto =', -ObsDist*sqrt(NN*NN-1.) - hcto
    end subroutine GetCalD
!------------------------------
!------------------------------
!------------------------------
!------------------------------
Subroutine RayFields_t_obs(idi,CD_i)
!   Calculate the field (i.e. integral over \zeta for all observer times t_o) for a single ray = pencil shower
! idi (=ObsDist) determines the distance to the pencil shower
! CD_i (CoreDist_A(CD_i)) determines the position of this pencil shower w.r.t. the core of the real shower
    use BigArrays, only : Ex_tb,AxD_tb,Ey_tb,AyD_tb,Ar_tb, Erh_tb, t_tb, t_ch
    use BigArrays, only : Ix,Iy,IQ, Ix_int,Iy_Int, ObsDist !, ObsDist_dim, ObsDist_Step
    use RFootPars, only : NTo
    use Atmosphere, only : xi,dxi, AtmHei_dim, AtmHei_step
   use Pancakefunction, only : CoreDist_Dim, CoreDist_A, alpha_tr, PancakeThi
    use constants, only : dp
   Use Pancakefunction, only : fieh
    implicit none
!    real*8 :: Ex(0:CoreDist_Dim),AxD(0:CoreDist_Dim),Ey(0:CoreDist_Dim),AyD(0:CoreDist_Dim),Ar(0:CoreDist_Dim),nt_r,t_o
    integer, intent(in) :: idi, CD_i
    real(dp) :: Ex,AxD,Ey,AyD,Ar,Erh
    real(dp) :: t_o
    real(dp) :: h_c,h,zeta,di,R,NN,dr,calD,Jx,Jx_int,Jy,Jy_int,dJx,dJy,JQ, nt_r
    real(dp) :: Intx,Inty,Intr,B,Lbd,weigh ! ,dh,max_h,fh_ce
    real(dp) :: fh,dfh,dfhl,dfha,alpha,lam, dlambda_hdr, dfhdr, Int_dh
    integer :: i, i_z, iz1, iz2, branch, i_t
    real(dp) :: hto, zet1, zet2, dhdz, dz
    !integer, save :: il_start = 1, CD_i_prev=-1
    common / time / zeta_c,hcto_c
    real(dp) :: zeta_c,hcto_c,lambda
    logical :: subtr
!
   DO i_t=2,NTo
      T_o=t_ch(idi)+t_tb(i_t)
      h_c=hcto_c+T_o
      Ex=0.   ; Ey=0. ; Ar=0. ; AxD=0.    ; AyD=0.    ; Erh=0.
      lambda=PancakeThi(CD_i)
      if(CD_i.eq.0) then  ! calculate d(lambda)/dr
        dlambda_hdr=0.5*(PancakeThi(1)-PancakeThi(0))/(CoreDist_A(1)-CoreDist_A(0))
      elseif(CD_i.eq.CoreDist_Dim) then
        dlambda_hdr=(PancakeThi(CD_i)-PancakeThi(CD_i-1))/(CoreDist_A(CD_i)-CoreDist_A(CD_i-1))
      else
        dlambda_hdr=(PancakeThi(CD_i+1)-PancakeThi(CD_i-1))/(CoreDist_A(CD_i+1)-CoreDist_A(CD_i-1))
      endif
      !
      !write(2,*) 'RayFields_t_obs: i_t=', i_t
      !flush(unit=2)
      lam=hcto_c + T_o  !  corresponds to h=0
      branch=-1
      zet1=zeta_c
      call GetZet_Lam(zet1,lam,branch)
      !If(zet1.lt.0.) zet1=0.
      branch=+1
      zet2=zeta_c
      call GetZet_Lam(zet2,lam,branch)
      !if(zet2.gt. AtmHei_step*AtmHei_dim) zet2=AtmHei_step*AtmHei_dim


      Intx=0. ;   Inty=0. ; Intr=0. ; B=0.
      dz=2.*AtmHei_step     ! will be halved later
      !dz=AtmHei_step     ! results in no observable difference
      i_z=(zet2-zet1)/dz
      If(i_z.lt.10) i_z=10   ! minimum number of integration steps
      dz=(zet2-zet1)/(2.*i_z)   ! should minimize effects of missing (over estimating) integral due to contributions at the end tips
      iz1=MAX( (zet1-zeta_c)/dz-50, -zeta_c/dz-0.5 )
      iz2=MIN( (zet2-zeta_c)/dz+50, (AtmHei_step*AtmHei_dim-zeta_c)/dz+0.5 )
      !iz1=(zet1-zeta_c)/dz-150
      !iz2=(zet2-zeta_c)/dz+150 ! no observable difference
      !write(2,*) 'CurrField_zeta:',iz1,iz2, zet1, zet2, zeta_c, dz
      !flush(unit=2)
      !If(cd_i.le.2 .and. idi.lt.3 .and. i_t.lt.10) write(2,*) 'RayFields_t_obs: i_t=', i_t, t_tb(i_t), &
      !      iz1*dz+zeta_c, iz2*dz+zeta_c, idi, cd_i
      !flush(unit=2)
      Do i_z=iz1, iz2
         zeta=zeta_c +(i_z+0.5)*dz
         i=(zeta/AtmHei_step)
         di=zeta/AtmHei_step-i
!If(zeta.lt.6000) cycle !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!If(zeta.gt.6500) cycle !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         If(i.lt.1) Then
            cycle
         Else If(i.ge.AtmHei_dim) then
            cycle
         endif
         R=sqrt(ObsDist*ObsDist +zeta*zeta)
         NN=1.+ di*xi(i+1) + (1.-di)*xi(i)
         dr=   di*dxi(i+1) + (1.-di)*dxi(i)
         hto=  zeta-NN*R
         h= hto + T_o
         If(h.lt.0.) cycle  ! Even infront of the shower front
         If(h.gt.10.*lambda) cycle  !  Too far from shower front
         dhdz=(R-NN*zeta - dr*R*R)/R
         !cald=NN*(R - NN*zeta - dr*R*R) ! retarded distance, including the needed factor n=NN
         weigh=dz/(NN*R) ! =dz* dhdz/cald
         !
         nt_r=zeta-h     ! - t_r (=negative t_r) to have a positive number
         i=(nt_r/AtmHei_step)
         if(i.lt.1) cycle
         di=nt_r/AtmHei_step-i
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
         dfhdr = - (fh + dfhl*h ) *dlambda_hdr/lambda
         Int_dh=Int_dh + dz* dhdz * dfh ! should be =0 when integrated over the complete height, any finite value is due to
         !                               integration-step-size problems which is why it is more stable to include the subtraction
         !  at large observer distances where 'cald' is almost constant over the integration range,
         !    the integral over dfh*Jx should yield zero, but this is from a cancellation of relatively large numbers.
         !
         !    fh=0. ! Just for checking the effects of this term, effect of this is invusible.
         ! Need to check the sign in the second (small) term, d/dh!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         Ex =Ex  + weigh * (dfh*Jx -fh *dJx) ! (dfh*Jx/cald - fh *dJx/cald)
         AxD=AxD + weigh*Jx_int*fh   ! proportional to integrated current, thus like the moving dipole
         Ey =Ey  + (dfh*Jy - fh *dJy) * weigh
         AyD=AyD + weigh*Jy_int*fh
         Ar =Ar  -(fh *JQ) * weigh    ! still needs to be times dw/dr
         Erh =Erh  -(dfhdr *JQ) * weigh    ! still needs to be times w(r)
      EndDo !  i_z=-iz1, iz2
      Ex_tb(i_t,idi,CD_i) =Ex
      AxD_tb(i_t,idi,CD_i)=AxD
      Ey_tb(i_t,idi,CD_i) =Ey
      AyD_tb(i_t,idi,CD_i)=AyD
      Ar_tb(i_t,idi,CD_i) =Ar
      Erh_tb(i_t,idi,CD_i)=Erh
   enddo  ! i=2,NTo
   return
End Subroutine RayFields_t_obs
!------------------------------
