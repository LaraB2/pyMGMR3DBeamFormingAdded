!------------------------------
Module Pancakefunction !------------------------
   use constants, only : dp
   implicit none
   real(dp),save :: dalpha=1./10.   ! Normalization of pancake function
   real(dp),save :: Nrm_alpha(0:41)
   Integer, save :: SelectFh=4
   integer, parameter :: CoreDist_Dim=20  ! number of line-showers per shower
   real(dp), save :: PancakeThi(0:CoreDist_Dim)
   real(dp), save :: CoreDist_A(0:CoreDist_Dim)    ! Distances of pencil showers to core
   real*8, allocatable, save :: alpha_tr(:)
contains
!------------------------------
Subroutine PancakeInit(XDepAlpha, PenDepth, AtmHei_dim, ObsDist_dim, ObsDist_Step, lam_100, lam_tc)
!  Initialize for pancake function
   use constants, only : dp
    implicit none
   Real(dp), intent(in) :: XDepAlpha
   Real(dp), intent(in) :: PenDepth(0:AtmHei_dim+1)
   Integer, intent(in) :: AtmHei_dim, ObsDist_dim
   Real(dp), intent(in) :: ObsDist_Step, lam_100, lam_tc
   Integer :: i, j, CD_i
   Real(dp) :: one = 1.0
   Real(dp) :: DerLambda=0.
   Real(dp) :: DerAlpha=0.
   Real(dp) :: Alpha, z, a, b, c, ddd, scrtch(0:40)
   !
   !  Initialize the normalization of the pancake
   write(2,*) '!PancakeInit:', XDepAlpha, AtmHei_dim, ObsDist_dim, ObsDist_Step, lam_100, lam_tc
   Nrm_alpha(:)=1.d0 !  just for set-up, used in fieh
   l_alpha: Do i=0,40
      alpha=1.d0+i * dalpha
      z=0.
      l_idi: do j=1,10000
         a=j*alpha/40.d0    ! distance
         !call fieh(a,1,alpha,ddd,b,c)  !    subroutine fieh(h,CoreDist_i,alpha,fh,dfhl,dfha)
         call fiehLam(a,one, DerLambda,alpha, DerAlpha,ddd,b,c)  !    subroutine fiehLamFx(h,Lambda, DerLambda, Alpha, DerAlpha, fh,dfh,dfhdr)
         z=z+ddd*alpha/40.d0
      enddo l_idi
      scrtch(i)=1.d0/z
   enddo l_alpha
   Nrm_alpha(0:40)=scrtch(0:40)
   Nrm_alpha(41)=Nrm_alpha(40)
   !
   !   Initialize Radial dependence of pancake thickness
    b=(ObsDist_dim)*ObsDist_Step/100.
    c=lam_100*b/lam_tc
    a=c**(1./(CoreDist_Dim))
    b=(ObsDist_dim)*ObsDist_Step/c    ! =100. * lam_tc/lam_100  !
    CoreDist_A(0)=b    ; PancakeThi(0)=lam_tc
    Do CD_i=1,CoreDist_Dim  ! calculate pancake thickness-parameter as function of distance to shower-core
      PancakeThi(CD_i)=PancakeThi(CD_i-1)*a  ! Later this will be multiplied by alpha
      CoreDist_A(CD_i)=PancakeThi(CD_i)*b/lam_tc
    enddo
    if(CoreDist_A(CoreDist_Dim).lt.ObsDist_dim*ObsDist_Step) CoreDist_A(CoreDist_Dim)=ObsDist_dim*ObsDist_Step
    !
    !  Depth-dependent factor to increase pancakethickness
    allocate(alpha_tr(0:AtmHei_dim+1))
    ! 0.0048 is the value for lofar, but should be angle dependent
    ! may be as large as 0.3 for lightning situations
    alpha_tr(0:AtmHei_dim)=1.d0 + XDepAlpha*sqrt(PenDepth(0:AtmHei_dim)/600.) + 0.0048
    !
   Return
End Subroutine PancakeInit
!------------------------------
Subroutine fiehcore(h,Lambda, Alpha, fh,dfhl,dfha)

!   Calculate the value of the pancake function, fh, and
!       the derivaties v.s. lambda (=integration variable in line-emission i.e. NOT 'lam'), dfhl, and
!       the derivaties v.s. alpha, dfha.
!   h=distance behind front
   use constants, only : dp
    implicit none
    real(dp), intent(in) :: h,Lambda, Alpha
    real(dp), intent(out) :: fh,dfhl,dfha
    real(dp) :: x,lam,di,nrm,a,b,dfdx,d
    real(dp), parameter :: c=1.d0, xScale=10.d0 , xScale3=30.d0
    integer :: i
    lam=alpha*Lambda
    x=h/lam
    if(x.gt.20.) then
      fh=0.0d0 ; dfhl=0. ; dfha=0.
      return
    endif
    i=IDINT((alpha-1.)/dalpha)
    di=alpha-1.-i*dalpha
    nrm=((1.-di)*Nrm_alpha(i) + di*Nrm_alpha(i+1))/lam
    !
    select case (SelectFh)
    !---------------------------------------------------
    case (1)
    !------  fh=x^alpha e^(-2x)  ----------------------
        fh =Nrm*x**(alpha-1) *exp(-2.*x) ! =fh/x
        dfhl=fh*(alpha-2*x)/lam   !  =df/dh at fixed alpha
        fh=fh*x
        dfha=(Nrm_alpha(i+1)-Nrm_alpha(i))/dalpha  ! =df/da at fixed h
        dfha=(dfha/Nrm + log(x) - (1-2*x/alpha))*fh
        ! PT=lam/alpha
        ! dfhdPT=fh/PT  - dfdx * x/PT =  alpha/lam *fh (1. - (alpha-2 x))
        ! from line integration:       dfhdr = - (1. + alpha - 2. * h/(alpha*lam) ) * fh*dlambda_hdr/(lam=pancakethicjness)
    !---------------------------------------------------
    case (3)
    !------  fh=x^2 /(exp(sqrt(x))+c*alpha)  -------------
        x=x*xScale3       ! x=10 h / lam
        a=sqrt(x)
        b=exp(a)
        d=c*alpha
        fh =Nrm*x/(b+d)
        dfdx=fh*(2.d0 - 0.5d0*a*b/(b+d))
        dfhl=dfdx*xScale3/lam     ! (dfh/dx) * (dx/dh=10/lam)
        fh=fh*x
        dfha=(Nrm_alpha(i+1)-Nrm_alpha(i))/dalpha  ! =df/da at fixed h
        dfha=(dfha/Nrm - c/(b+d) )*fh - dfdx*x/alpha ! (dx/da=-x/alpha)
    !---------------------------------------------------
    case (4)
    !------  fh=x /(exp(sqrt(x))+1)  -------------
        x=x*xScale
        a=sqrt(x)
        b=exp(a)
        fh =Nrm /(b+c)
        dfdx=fh*(1.d0- 0.5d0*a*b/(b+c))
        dfhl=dfdx*xScale/lam
        fh=fh*x
        !dfha=(Nrm_alpha(i+1)-Nrm_alpha(i))/dalpha  ! =df/da at fixed h
        dfha=fh*xScale/alpha - dfdx*x*xScale/alpha  ! uncorrected again;  Correctly corrected ?????????????????????????????????????????????????????????????????????
    !---------------------------------------------------
    case default
    !------  fh=x^alpha /(exp(sqrt(x))+1)  -------------
        x=x*xScale
        a=sqrt(x)
        b=exp(a)
        fh =Nrm*x**(alpha-1) /(b+c)
        dfdx=fh*(alpha- 0.5d0*a*b/(b+c))
        dfhl=dfdx*xScale/lam
        fh=fh*x
        dfha=(Nrm_alpha(i+1)-Nrm_alpha(i))/dalpha  ! =df/da at fixed h
        dfha=(dfha/Nrm + log(x))*fh - dfdx*x*xScale/alpha
    !---------------------------------------------------
    !
    end select
    !
    return
end subroutine fiehcore
!------------------------------================================
Subroutine fiehLam(h,Lambda, DerLambda, Alpha, DerAlpha, fh,dfh,dfhdr)
!   Calculate the value of the pancake function, fh, and
!       the derivaties v.s. lambda (=integration variable in line-emission i.e. NOT 'lam'), dfhl, and
!       the derivaties v.s. alpha, dfha.
!   h=distance behind front
   use constants, only : dp
   implicit none
   real(dp), intent(in) :: h,Lambda, DerLambda, Alpha, DerAlpha
   real(dp), intent(out) :: fh,dfh,dfhdr
   real(dp) :: dfhl,dfha
   !
   Call fiehcore(h,Lambda, Alpha, fh,dfhl,dfha)
   !
   dfh=dfhl + dfha * DerAlpha
   dfhdr = - (fh + dfhl*h ) * DerLambda/Lambda
!            dfh=dfhl + dfha *(alpha_tr(i+1)-alpha_tr(i))/AtmHei_step
!            dfhdr = - (fh + dfhl*h ) *dlambda_hdr/lambda

    return
end subroutine fiehLam
!------------------------------================================
subroutine fieh(h,CoreDist_i,alpha,fh,dfhl,dfha)
   use constants, only : dp
   implicit none
   integer, intent(in) :: CoreDist_i
   real(dp), intent(in) :: h, Alpha
   real(dp), intent(out) :: fh, dfhl, dfha
   real(dp) :: Lambda
   !
   Lambda = PancakeThi(CoreDist_i)
   !
   Call fiehcore(h,Lambda, Alpha, fh,dfhl,dfha)
   !
   Return
   !
End subroutine fieh
!------------------------------================================
Subroutine GetLambda(rc_r, Lambda, DerLambda, CD_i, CD_d)
!
   use constants, only : dp
   implicit none
   real(dp), intent(in) :: rc_r
   real(dp), intent(out) :: Lambda, DerLambda, CD_d
   integer, intent(out) :: CD_i
   real(dp) :: nt_r, di,Jx,Jx_int,Jy,Jy_int,dJx,dJy,JQ
   real(dp) :: dlambda0, dlambda1
   integer :: i !, idi, i_t, i_h, N_h
   !
   call Find_CD_Interpol(rc_r,CD_i,CD_d)
   Lambda=(1.-CD_d)*PancakeThi(CD_i) + CD_d*PancakeThi(CD_i+1)
   if(CD_i.eq.0) then  ! calculate d(lambda)/dr
     DerLambda=0.5*(PancakeThi(1)-PancakeThi(0))/(CoreDist_A(1)-CoreDist_A(0))
   elseif(CD_i.ge.(CoreDist_Dim-1)) then
     DerLambda=(PancakeThi(CoreDist_Dim)-PancakeThi(CoreDist_Dim-1))/(CoreDist_A(CoreDist_Dim)-CoreDist_A(CoreDist_Dim-1))
   else
     dlambda0=(PancakeThi(CD_i+1)-PancakeThi(CD_i-1))/(CoreDist_A(CD_i+1)-CoreDist_A(CD_i-1))
     dlambda1=(PancakeThi(CD_i+2)-PancakeThi(CD_i))/(CoreDist_A(CD_i+2)-CoreDist_A(CD_i))
     DerLambda=(1.-CD_d)*dlambda0 + CD_d*dlambda1
   endif
   Return
End Subroutine GetLambda
!------------------------------
subroutine Find_CD_Interpol(CD,CD_i,CD_d)
   use constants, only : dp
   implicit none
    real(dp), intent(in) :: CD
    integer , intent(out) :: CD_i
    real(dp), intent(out) :: CD_d
    integer :: i
    CD_i=CoreDist_Dim-1
    Do i=1,CoreDist_Dim
      if(CD.lt.CoreDist_A(i)) then
        CD_i=i-1
        exit
      endif
    enddo
    CD_d=(CD-CoreDist_A(CD_i))/(CoreDist_A(CD_i+1)-CoreDist_A(CD_i))
    if(CD_d.lt.0.) CD_d=0.
    ! write(2,*) 'cd',cd,cd_i,cd_d,CoreDist_A(CD_i),CoreDist_A(CD_i+1)
end subroutine Find_CD_Interpol
!------------------------------
End Module Pancakefunction !------------------------
!---------------------------------------------------------
Subroutine Findhs(t_a,r_a,tr_s,h_s,NN_s,dr_s)
!  For a given t_a (antenna time), r_a (distance antenna to the ray considered), and t_r (retarded time in the shower)
!     find the emission height, or
!     NN_s*(h_s-tr_s) = sqrt((t_a-tr_s)*(t_a-tr_s)-NN_s*NN_s*r_a*r_a)
!     zeta_s = h_s-tr_s
!        or
!     NN_s^2 [zeta_s^2 + r_a^2] = (t_a-tr_s)^2
   use Atmosphere, only : AtmHei_step, AtmHei_dim, xi,dxi
   use constants, only : dp
   implicit none
   real(dp), intent(in) :: t_a,r_a,tr_s
   real(dp), intent(out) :: h_s,NN_s,dr_s
   Integer :: i, k
   real(dp) :: di,RT, zeta_s, h
   Real(dp), parameter :: Conv=1.d-6
   !
   h_s=0.
   h=-1.
   Do k=1,10
      zeta_s=-tr_s + h_s
      i=(zeta_s/AtmHei_step)
      if(i.lt.1) Then
         NN_s=1.+ xi(1)
         dr_s=   dxi(1)
      ElseIf(i.ge.AtmHei_dim)  Then
         NN_s=1.+ xi(AtmHei_dim)
         dr_s=   dxi(AtmHei_dim)
      Else
         di=zeta_s/AtmHei_step-i
         NN_s=1.+ di*xi(i+1) + (1.-di)*xi(i)
         dr_s=   di*dxi(i+1) + (1.-di)*dxi(i)
      EndIf
      RT=(t_a-tr_s)*(t_a-tr_s)-NN_s*NN_s*r_a*r_a
      If(RT.lt.0.) exit
      h=tr_s + sqrt(RT)/NN_s
      If(abs(h-h_s).lt. conv) exit
      h_s=h
      If(k.gt.5) write(2,*) 'iteration#',k, h, h_s, t_a,r_a,tr_s
   EndDo ! k=1,10
   !
   Return
End Subroutine Findhs
!------------------------------
