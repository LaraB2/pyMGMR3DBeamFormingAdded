      module FFT
    use constants, only : dp,pi,ci
  implicit none
  integer ( kind = 4 ), save :: lensav
  integer ( kind = 4 ), save :: lenwrk
    integer ( kind = 4 ), save :: N_Time
    integer ( kind = 4 ), save :: N_nu
!  Allocate the work arrays.
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave
      contains
!      ------------------------------------
    subroutine FFTransform_su(N_t)
! Initialize FFT transform and allocate workspace
! N_t= Length of input array in time domain.
  implicit none
  integer , intent (in) :: N_t
  integer ( kind = 4 ) ier
!  integer ( kind = 4 ) inR
!  integer ( kind = 4 ) lenR
  N_Time=N_T
  N_nu=N_Time/2
  lensav = N_t + int ( log ( real ( N_t, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4
  lenwrk = N_t
  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )
!    write(2,*) 'FFTransform_su',lenwrk,lensav,N_t
  call Rfft1i ( N_t, wsave, lensav, ier )
  return
  end subroutine FFTransform_su
!      ------------------------------------
    subroutine DAssignFFT()
       deallocate ( work )
       deallocate ( wsave )
    end subroutine DAssignFFT
!      ------------------------------------
  subroutine DownSamlple(E_t,E,padding,tTrace_dim,FF_dim, filt, E_nu_dwn, inui,inum)

!  subroutine FFTransform_FF(N_t,A,N_f,f,t_shft,C, wsave, lensav, work, lenwrk)
! Fourier transform to frequency, apply filter and timeshift
! A= input array in time domain, length=N_t
! F=input filter function (complex) in frequency
! C_nu=output array (complex) in frequency, length=N_f; after filter
  implicit none
  integer, intent(in) :: inui,inum,padding,tTrace_dim, FF_dim
  complex(dp), intent(out) :: E_nu_dwn(inui:inum)
  REAL(dp), intent(in) :: E_t(1:tTrace_dim)
  REAL(dp), intent(inout) :: E(1:FF_dim)
  complex(dp), intent(in) :: filt(inui:inum)
  integer :: i,nu_dim, inuip
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inR
  integer ( kind = 4 ) lenR
!
!  Compute the FFT coefficients.
!
  inR = 1
  lenR = FF_dim
  nu_dim=ff_dim/2
  E(1:padding)=0.0d0 ;  E(padding+1:tTrace_dim+padding)=E_t(1:tTrace_dim) ;  E(tTrace_dim+padding+1:FF_dim)=0.0d0
  call Rfft1f (FF_dim, inR, E, FF_dim, wsave, lensav, work, lenwrk, ier )
  If(inui.eq.0) then
   InuIp=InuI+1
   E_nu_dwn(0)=filt(0)* E(1)
  Else
   InuIp=InuI
  EndIf
  !E_nu_dwn(InuI)=0.
  DO I=inuip,inum
    E_nu_dwn(I)=filt(I)*CMPLX(E(2*I),E(2*I+1))
  ENDDO
  return
  end subroutine DownSamlple
!
!      ------------------------------------
  subroutine FFTransform_FF_x(A,F,t_shft,Cnu,inui,inum)
!  subroutine FFTransform_FF(N_t,A,N_f,f,t_shft,C, wsave, lensav, work, lenwrk)
! Fourier transform to frequency, apply filter and timeshift
! A= input array in time domain, length=N_t
! f=input filter function (complex) in frequency
! t_shft= shift of pulse in time in units of sample-time
! c=output array (complex) in frequency, length=N_f; after filter
  !1! use BigArrays, only : tTrace_dim_o, nuTrace_dim  ! seems not to be used in MGMR
  !  use constants, only : dp,pi,ci
  implicit none
  integer, intent(in) :: inui,inum
  complex(dp), intent(out) :: Cnu(inui:inum)
  !1! REAL(dp), intent(in) :: A(1:tTrace_dim_o),t_shft
  REAL(dp), intent(in) :: A(1:N_Time),t_shft
  integer :: i
  !1! REAL(dp) :: D(1:tTrace_dim_o)
  !1! complex(dp) :: C(0:nuTrace_dim),F(0:nuTrace_dim),ph_shft,ipi
  REAL(dp) :: D(1:N_Time)
  complex(dp) :: C(0:N_nu),F(0:N_nu),ph_shft,ipi
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inR
  integer ( kind = 4 ) lenR
  ipi=ci*pi  ! This is the real constant i*pi
!
!  Compute the FFT coefficients.
!
!    write(2,*) 'FFTransform_FF',lenwrk,lensav,N_t
  D(:)=A(:)
  inR = 1
!1!   lenR = tTrace_dim_o
!1! !  call Cfft1f ( N_t, inR, B, lenR, wsave, lensav, work, lenwrk, ier )
!1!   call Rfft1f ( tTrace_dim_o, inR, D, lenR, wsave, lensav, work, lenwrk, ier )
!1!   CALL R2C(D,C,tTrace_dim_o,nuTrace_dim)
  lenR = N_Time
!  call Cfft1f ( N_t, inR, B, lenR, wsave, lensav, work, lenwrk, ier )
  call Rfft1f ( N_Time, inR, D, lenR, wsave, lensav, work, lenwrk, ier )
  CALL R2C(D,C,N_Time,N_nu)
!
!    write(2,*) '0',abs(B(0)),ph(B(0)), f(0)
!    write(2,*) '1',abs(B(1)),ph(B(1)), f(1)
!    write(2,*) 'n/2',abs(B(N_f/2)),ph(B(N_f/2)) , f(n_f/2)
!    write(2,*) 'n-1',abs(B(N_f-1)),ph(B(N_f-1))
!    write(2,*) 'n',abs(B(N_f)),ph(B(N_f))
  DO I=inui,inum
!1!     ph_shft=exp(ipi*t_shft*i/nuTrace_dim)
    ph_shft=exp(ipi*t_shft*i/N_nu)
    Cnu(I)=C(I)*F(I) *ph_shft
  ENDDO
  return
  end subroutine FFTransform_FF_x
!
!      ------------------------------------
!      ------------------------------------
  Subroutine RFTransform_CF(A,Cnu)
  ! transform real time trace to complex frequency
   !  time and frequency dimension should have been set in a previous call to "RFTransform_su(N_T)"
  implicit none
  complex(dp), intent(out) :: Cnu(0:N_nu)
  REAL(dp), intent(in) :: A(1:N_Time)
  !
  integer :: i
  REAL(dp) :: D(1:N_Time)
  complex(dp) :: C(0:N_nu)
  integer ( kind = 4 ) :: ier, inr=1, lenR
!
!  Compute the FFT coefficients.
!
  !  write(*,*) 'FFTransform_FF',lenwrk,lensav
  D(:)=A(:) !  *Hann(:)
  call Rfft1f ( N_Time, inR, D, N_Time, wsave, lensav, work, lenwrk, ier )
  CALL R2C(D,Cnu,N_Time,N_nu)
!
  return
  End Subroutine RFTransform_CF
!
!      ------------------------------------
  Subroutine RFTransform_CF2RT(Cnu,RD)
! Fourier transform back to time domain
! Cnu=input array (complex) in frequency, length=N_f
! RD= Real output array in time domain, length=N_t
!  use constants, only : dp,ci
  implicit none
  complex(dp), intent(in) :: Cnu(0:N_nu)
  real(dp), intent(out) :: RD(1:N_Time)
  !
  integer ( kind = 4 ) :: ier, inr=1, lenR
  !
  lenR = N_Time
!
!  Compute inverse FFT of coefficients.
  CALL C2R(RD,Cnu,N_Time,N_nu)
  call Rfft1b ( N_Time, inR, RD, lenR, wsave, lensav, work, lenwrk, ier )
!
  return
  End Subroutine RFTransform_CF2RT
!
!      ------------------------------------
  subroutine FFTransform_B(CD,tTrace_dim,Cnu,nuTrace_dim,inui,inum)
! Fourier transform back to time domain
! C=output array (complex) in frequency, length=N_f; after filter
! CD= complex output array in time domain, length=N_t; after filter
  !use constants, only : dp,ci
  implicit none
  integer, intent(in) :: inui,inum,tTrace_dim, nuTrace_dim
  complex(dp), intent(in) :: Cnu(inui:inum)
  complex(dp), intent(out) :: CD(1:tTrace_dim)
!  complex *16 Ci
    INTEGER*4 i
    REAL*8 D(1:tTrace_dim),ID(1:tTrace_dim)
    complex*16 C(0:nuTrace_dim),IC(0:nuTrace_dim)
     integer ( kind = 4 ) ier
  integer ( kind = 4 ) inR
  integer ( kind = 4 ) lenR
!  ci=complex(0.,1.)  ! This is the real constant i*pi
  inR = 1
  lenR = tTrace_dim
  IC(:)=0.
  c(:)=0.
  DO I=inui,inum
      !  write(2,*) '!FFTransform_BI', i, Cnu(i)
      !flush(Unit=2)
    C(i)=Cnu(I)
    IC(I)= + ci*Cnu(i)
  ENDDO
!
!  Compute inverse FFT of coefficients.
!
  CALL C2R(D,C,tTrace_dim,nuTrace_dim)
  call Rfft1b ( tTrace_dim, inR, D, lenR, wsave, lensav, work, lenwrk, ier )
  CALL C2R(ID,IC,tTrace_dim,nuTrace_dim)
  call Rfft1b ( tTrace_dim, inR, ID, lenR, wsave, lensav, work, lenwrk, ier )
!
  DO I=1,tTrace_dim
    CD(i)=D(i) + ci*ID(i)
  ENDDO
  return
  end subroutine FFTransform_B
!
!------------------------------------------------
SUBROUTINE R2C(R,C,Nr,Nc)
! Nc=Nr/2 note that zero and max frequency components are always real
    IMPLICIT NONE
    INTEGER Nr,Nc,I
    REAL*8 R(1:Nr)
    COMPLEX*16 C(0:Nc)
  C(0)=R(1)
  DO I=1,Nc-1
    C(I)=CMPLX(R(2*I),R(2*I+1))
  ENDDO
  C(Nc)=R(2*Nc)
  RETURN
END
SUBROUTINE  C2R(R,C,Nr,Nc)
! Nc=Nr/2 note that zero and max frequency components are always real
    IMPLICIT NONE
    INTEGER Nr,Nc,I
    REAL*8 R(1:Nr)
    COMPLEX*16 C(0:Nc)
  R(1)=DREAL(C(0))
  DO I=1,Nc-1
    R(2*I)=DREAL(C(I))
    R(2*I+1)=DIMAG(C(I))
  ENDDO
  R(2*Nc)=DREAL(C(Nc))
  RETURN
END
!      ------------------------------------
	function ph(zzz)
   ! use constants, only : pi
	implicit none
	complex*16, intent(in) :: zzz
	real*8 :: re,im,ph
!	common / constants / CI,pi
!	  COMPLEX*16 CI
!	  REAL*8 pi
	re=REALpart(zzz)
	im=IMAGpart(zzz)
	if(im.eq.0) then
		ph=0.
	else
		ph=datan(im/re)
	endif
	if(re.lt.0.) ph=ph+pi
	if(ph.gt.pi) ph=ph-2*pi
	ph=ph/pi
	return
	end
!      ------------------------------------
      end module FFT
! ==========================================================
