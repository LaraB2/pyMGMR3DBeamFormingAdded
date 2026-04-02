Module PSF
   use constants, only : pi,dp, c_l
   Integer, parameter :: PSF_max=25
   Integer, parameter :: t_dim=512
   Integer, save :: N_psf
   Character(len=10), save :: Label2
   Character(len=13), save :: Label3(1:PSF_max)
   Character(len=20), save :: Base
   Integer, save :: hcnt_l, hcnt_u !, N_zf, N_tr
   Real(dp), save, allocatable :: Grd_Kx(:,:,:)
   Real(dp), save, allocatable :: GrdPSFWghtx(:,:)     ! GrdPSFWght(0:N_tr,0:N_zf) )  ! GrdPSFWght(i_tr,i_zf) ; put on rectangular grid where linear interpolation of current is assumed
   Real(dp), save, allocatable :: GrdPwrWghtx(:,:)     ! GrdPwrWght(0:N_tr,0:N_zf) )  ! GrdPwrWght(i_tr,i_zf) ; put on rectangular grid where linear interpolation of current is assumed
   Real(dp), save, allocatable :: PSF_tracex(:,:)
   Real(dp), save, allocatable :: Grd_Kr(:,:,:)  !  Obsolete, but referred to in shower structure
contains
!---------------------------
!---------------------------
!---------------------------
   !----------------------------------------------------
   Subroutine ReadKernel( N_PSF, dt_PSF, z_low, del_z, del_tr, N_zf, N_tr, AtmHei_dim, AtmHei_step &
   , Base, KernBase &
   , Grd_K &
   , GrdPSFWght &    ! GrdPSFWght(0:N_tr,0:N_zf) )  ! GrdPSFWght(i_tr,i_zf) ; put on rectangular grid where linear interpolation of current is assumed
   , GrdPwrWght &   ! GrdPwrWght(0:N_tr,0:N_zf) )  ! GrdPwrWght(i_tr,i_zf) ; put on rectangular grid where linear interpolation of current is assumed
   , PSF_trace )
   !----------------------------------------------------
   !  Read Kernel Inp_K(i_h,0:t_r,if_PSF) (depending on shower time and zeta_f=Zf_PSF(1):Zf_PSF(N_PSF))
   !  i_h=hcnt_l, hcnt_u   ;  h=(i_h-0.5)*(dt_PSF = del_t)
   !  i_tr=it_min,N_tr+it_min      ;  t_r=i_tr*dN_tr*AtmHei_step
   !  if_PSF=0,N_PSF       ;  zeta_f=Zf_PSF(if_PSF)
   !
   use KernelFilter, only : FreqFilterKernel
   Use ShowerData, only : SamplingTime, nu_min, nu_max
	implicit none
   Integer, intent(in) :: N_PSF
   Real(dp), intent(out) :: dt_PSF
   Real(dp), intent(out) :: z_low
   Real(dp), intent(in) :: del_z
   Real(dp), intent(in) :: del_tr
   Integer, intent(out) :: N_zf
   Integer, intent(out) :: N_tr
   Integer, intent(in) :: AtmHei_dim
   Real(dp), intent(in) :: AtmHei_step
   Character(len=25), intent(in) :: KernBase, Base
   Real(dp), allocatable, intent(out) :: Grd_K(:,:,:)
   Real(dp), allocatable, intent(out) :: GrdPSFWght(:,:)     ! GrdPSFWght(0:N_tr,0:N_zf) )  ! GrdPSFWght(i_tr,i_zf) ; put on rectangular grid where linear interpolation of current is assumed
   Real(dp), allocatable, intent(out) :: GrdPwrWght(:,:)     ! GrdPwrWght(0:N_tr,0:N_zf) )  ! GrdPwrWght(i_tr,i_zf) ; put on rectangular grid where linear interpolation of current is assumed
   Real(dp), allocatable, intent(out) :: PSF_trace(:,:)
   !Character(len=10), intent(in) :: Label2
   !Character(len=13), intent(in) :: Label3(1:N_PSF)
   Integer :: i, k, nxx, if_PSF, ih_psf, irt_psf, N_z, i_zf, i_tr, ih_1, ih_2, odN_tr, oih_1, oih_2, i_h
   Character(len=3) :: Label
   Integer :: dN_tr, it_min, padding, NuDim_new
   Real(dp) :: Norm, d_h, Zf_PSF(1:N_PSF), od_h
   Real(dp), allocatable ::  Inp_K(:,:,:), temp(:,:), K_filt(:,:), RMS_InpK(:,:), PSF_unfilt(:,:)
   Character(len=220) :: Line
   !
   Real(dp) :: X,Y, z , W ! scratch
   Real(dp) :: t_r, d_tr
   Real(dp) :: zeta_f, d_zeta, zeta_low, Zwin
   !
   Do if_PSF=1,N_PSF  ! label of PSF file
      OPEN(UNIT=4,STATUS='old',FILE=TRIM(KernBase)//'SrcKrnl-'//TRIM(label2)//TRIM(Label3(if_psf))//'.z',iostat=nxx) ! E13.4
      write(2,*) 'Reading file "',TRIM(KernBase)//'SrcKrnl-'//TRIM(label2)//TRIM(Label3(if_psf))//'.z','"'
      If(nxx.ne.0) Then
         write(2,*) 'problems reading file:', TRIM(KernBase)//'SrcKrnl-'//TRIM(label2)//TRIM(Label3(if_psf))//'.z'
         STOP 'SrcKrnl file is problematic'
      EndIf
      read(4,*) label
      ! assumed that for all files: dN_tr, hcnt_l, hcnt_u, dt_PSF  are identical
      read(4,"(A220)") Line
      If(line(1:5) .eq. '! Gen') Then
         Write(2,*) Line
         read(4,"(A220)") Line
      EndIf
      !Write(2,*) Line
      read(Line,*,iostat=nxx) label, N_z, dN_tr, it_min, ih_1, ih_2, Zf_PSF(if_PSF), d_h ! i_atm=i_tr*dN_tr+it_min
      If(ih_1.lt.1-ih_2) then
         write(2,*) '******** -ih_1 should be smaller then ih_2', ih_1, ih_2
         Stop 'ih_1, ih_2'
      EndIf
      !ih_1=1-ih_2  !  Needed for downsampling
      If(.not. Allocated(Inp_K)) Then
         Allocate( Inp_K(1-ih_2:ih_2, 0:AtmHei_dim/dN_tr, 1:N_PSF), RMS_InpK(0:AtmHei_dim/dN_tr,1:N_PSF) )
         odN_tr=dN_tr ; oih_1=ih_1 ; oih_2=ih_2 ; od_h=d_h
      Else
         If((odN_tr.ne.dN_tr) .or. (oih_1.ne.ih_1) .or. (oih_2.ne.ih_2) .or. (od_h.ne.d_h) ) Then
            Write(2,*) '*** Mismatch in generating parameters: (', odN_tr, oih_1, oih_2, od_h,') differs from (', &
               odN_tr, oih_1, oih_2, od_h,')'
            Stop 'File SrcKrnl mismatch'
         EndIf
      EndIf
      Inp_K(:, 0:it_min,if_PSF)=0.
      RMS_InpK(:,if_PSF)=0.
      Do i_zf=0, N_z  !   i_atm=(itf_0+it_min)*dN_tr
         Read(4,*) Inp_K(ih_1:ih_2, i_zf+it_min,if_PSF)
         RMS_InpK(i_zf+it_min,if_PSF)=sqrt( SUM(Inp_K(ih_1:ih_2, i_zf+it_min,if_PSF)*Inp_K(ih_1:ih_2, i_zf+it_min,if_PSF)) ) ! used for interpolation
      EndDo
      Close(Unit=4)
      !Write(2,"(30g13.4)") '!K(-4:5, 20, i)=',ih_1,ih_2, Inp_K(-4:5, 20, if_PSF)
      !Write(2,"(A,30g13.4)") '!RMS_InpK(20:30, i)=', N_z, RMS_InpK(20:30, if_PSF)
      !Write(2,"(A,30g13.4)") '!Inp_K(1,20:30, i)=', 1, Inp_K(1,20:30, if_PSF)
      !Write(2,"(A,30g13.4)") '!Inp_K(1,20:30, i)=', ih_1+k+1, Inp_K(ih_1+k+1,20:30, if_PSF)
      !Write(2,"(A,30g13.4)") '!Inp_K(1,20:30, i)=', ih_2-k-1, Inp_K(ih_2-k-1,20:30, if_PSF)
   EndDo
   !---------------------------------------------------------------
   !  Produce gridded & frequency-filtered & downsampled kernel  Grd_K(i_h, i_tr, i_zf) with
   !  i_h=hcnt_l, hcnt_u   ;  h=(i_h-0.5)*dt_PSF
   !  i_zf=0,N_zf          ;  zeta_f=i_zf*del_z + z_low  ! focal height for interferometry
   !  i_tr=0,N_tr          ;  t_r=i_tr*del_tr + Z_low    ! the contribution of the current at z=t_r to the interferometry trace
   !
   !  old sampletime =d_h, new sampletime = del_t
   !  Frequency filter between nu_min & nu_max [MHz]
   dt_PSF=SamplingTime ! =5*c_l  ! 5~ns converted to units [m]
   Padding=1+50*dt_PSF/d_h  ! allow for at least 50 new samples left and right of the region of interest.
   write(2,*) 'SampleTime_old, SampleTime_new:', d_h, dt_PSF, &
      '[m], nu_min, nu_max:', nu_min, nu_max,'[MHz]'
   hcnt_l = -40./dt_PSF+0.5  ! corresponds to h=-40 [m]
   hcnt_u = +50./dt_PSF+0.5  ! corresponds to h=+50 [m]; needed to have a minimal frequency resolution for filtering
   If(hcnt_u*dt_PSF .gt. (ih_2+Padding)*d_h) Then
      Padding= hcnt_u*dt_PSF/d_h  -ih_2 +1
      write(2,*) '** Padding increased to ',Padding
   EndIf
   !
   Z_low=Zf_PSF(1)  !
   N_zf=(Zf_PSF(N_PSF)-Z_low)/del_z  !  zeta_f=i_zf*del_z  + Z_low
   N_tr=(Zf_PSF(N_PSF)-Z_low)/del_tr  !  t_r = i_tr*del_tr + Z_low; new binning in t_r
   if_psf=1
   Allocate( Grd_K(hcnt_l:hcnt_u, 0:N_tr, 0:N_zf) )
   Allocate( temp(ih_1:ih_2, 0:N_tr) )
   Allocate( PSF_unfilt(ih_1:ih_2, 0:N_tr) )
   Allocate( PSF_trace(hcnt_l:hcnt_u, 0:N_zf) )
   write(2,*) '!ReadKernel, hcnt_l,hcnt_u,ih_1,ih_2:', hcnt_l,hcnt_u,ih_1,ih_2
   write(2,"(A,20G13.4)") '!in time [m], hcnt_l,hcnt_u,ih_1,ih_2:', hcnt_l*dt_PSF,hcnt_u*dt_PSF, ih_1*d_h,ih_2*d_h &
      , (ih_1-Padding)*d_h, (ih_2+Padding)*d_h
   Flush(unit=2)
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(Base)//'BeamKer-Dc.dat')
   k=25 !=5.5km
   k=NINT((5500.-z_low)/del_z)
   write(2,*) '******************** k=',k,';', TRIM(Base)//'BeamKer-Dc.dat'
   !k=50 !=10.5km
   !k=45 !=  km
   !k=3
   write(4,"(A,F7.0,i4,A,F6.3,A,F6.3,A)") '! un-filtered kernel at fixed D_c=',k*del_tr + z_low,k, &
         ', h=-0.5,+0.5*',d_h, '(*',dt_PSF,'), & h=0. & same for filtered '
   Do i_zf=0,N_zf
      zeta_f=i_zf*del_z + z_low
      If(if_psf .lt. N_PSF) Then
         If(zeta_f .gt. Zf_PSF(if_psf+1) ) Then
            if_psf=if_psf+1 ! Assume read-in grid is more sparse then requested grid
         EndIf
      EndIf
      Call Interpol(zeta_f,if_psf, zeta_f, AtmHei_dim, AtmHei_step, dN_tr, ih_1,ih_2, Zf_PSF, Inp_K, RMS_InpK, &
            PSF_unfilt(ih_1:ih_2,i_zf) )
      Do i_tr=0,N_tr  ! assume Inp_K is depends mainly or var=t_r/zeta_f and zeta_f for fixed h
         t_r=i_tr*del_tr + z_low
!Subroutine Interpol(zeta_f,if_psf, t_r, AtmHei_dim, AtmHei_step, dN_tr, ih_1, ih_2, Zf_PSF, Inp_K, RMS_InpK, temp )
         Call Interpol(zeta_f,if_psf, t_r, AtmHei_dim, AtmHei_step, dN_tr, ih_1, ih_2, Zf_PSF, Inp_K, RMS_InpK, &
             temp(ih_1, i_tr) )
      EndDo
      !-------------------------------------------------------------
      ! Perform frequency filtering of the Kernel
      !write(2,"(A,50G12.4)") '!-1temp(-2:+2,,)',N_tr-1, temp(-2:+2, N_tr-1)
      !write(2,"(A,50G12.4)") '!trtemp(-2:+2,,)',N_tr, temp(-2:+2, N_tr)
      !write(2,*) '!zeta:',i_zf,zeta_f, del_tr, z_low, N_tr*del_tr + z_low
      If(temp(0, N_tr) .gt. 0.5) Stop 'temp(Nz)'
      Call FreqFilterKernel(temp, ih_2, N_tr, d_h, padding, dt_PSF, nu_min, nu_max, &
         K_filt, NuDim_new)
      If(hcnt_l.lt. -NuDim_new) hcnt_l=-NuDim_new
      If(hcnt_u.gt.  NuDim_new) Then
         write(2,*) '**** Insufficient paddig, hcnt_u:', hcnt_u, NuDim_new
         stop  '**hcnt_u**'
         hcnt_u= NuDim_new
      EndIf
      Grd_K(hcnt_l:hcnt_u, 0:N_tr, i_zf)=K_filt(hcnt_l:hcnt_u, 0:N_tr)
      i_tr=k
      x=1.e5 ! (i_tr*del_tr + z_low)/zeta_f
      write(4,"(30g13.4)") i_zf, zeta_f/1000., x*temp(0, i_tr), x*temp(1, i_tr), x*(temp(0, i_tr)+ temp(1, i_tr))/2. &
         , x*K_filt(0, i_tr), x*K_filt(1, i_tr), x*(K_filt(0, i_tr)+ K_filt(1, i_tr))/2. &
         , (i_tr*del_tr + z_low - zeta_f)/1000.
      If(i_zf.eq.k) Then
         !Stop 'interp::'
         OPEN(UNIT=14,STATUS='unknown',FILE=TRIM(Base)//'BeamKer-Db.dat')
         write(14,"(A,F7.0,i4,A,F6.3,A,F6.3,A)") '! un-filtered kernel at fixed D_b=',i_zf*del_z + z_low,i_zf, &
               ', h=-0.5,+0.5*',d_h, '(*',dt_PSF,'), & h=0. & same for filtered '
         Do i_tr=0,N_tr
            write(14,"(30g13.4)") i_tr, (i_tr*del_tr + z_low)/1000. &
               , x*temp(0, i_tr), x*temp(1, i_tr), x*(temp(0, i_tr)+ temp(1, i_tr))/2. &
               , x*Grd_K(0, i_tr, i_zf), x*Grd_K(1, i_tr, i_zf), x*(Grd_K(0, i_tr, i_zf)+ Grd_K(1, i_tr, i_zf))/2. &
               , (i_tr*del_tr + z_low - zeta_f)/1000.
         EndDo
         Close(unit=14)
         OPEN(UNIT=14,STATUS='unknown',FILE=TRIM(Base)//'BeamKer-PSF.dat')
         write(14,"(A,F7.0,i4,A,F6.3,A,F6.3,A)") '! un-filtered kernel PSF(t_b) at fixed D_b=D_c=',i_zf*del_z + z_low,i_zf, &
               ',  & same for filtered '
         Do i_h=ih_1,ih_2
            If(i_h.le.hcnt_l .or. i_h.gt.hcnt_u) Then
               write(14,"(30g13.4)") i_h, (i_h-0.5)*d_h/c_l, x*temp(i_h, i_zf)  &
                  , (i_h-0.5)*dt_PSF/c_l, 0.0
            Else
               write(14,"(30g13.4)") i_h, (i_h-0.5)*d_h/c_l, x*temp(i_h, i_zf)  &
                  , (i_h-0.5)*dt_PSF/c_l, x*Grd_K(i_h, i_zf, i_zf)
            EndIf
         EndDo
         Close(unit=14)
      EndIf
      !
      ! ---------------------------------
   EndDo
   Close(unit=4)
   !
   !Write(2,"(A,30g13.4)") '!PSF_unfilt(-4:5, 10)=',  PSF_unfilt(-4:5, 10)
   !Write(2,"(A,30g13.4)") '!PSF_unfilt(-4:5, 20)=',  PSF_unfilt(-4:5, 20)
   !-------------------------------------------------------------
   ! Perform frequency filtering of the Pulse Spread Function
   DeAllocate(K_filt)
   !write(2,*) '!PSF-filt:', ih_2, N_zf, d_h, padding, dt_PSF, nu_min, nu_max, NuDim_new
   Call FreqFilterKernel(PSF_unfilt, ih_2, N_zf, d_h, padding, dt_PSF, nu_min, nu_max, &
         K_filt, NuDim_new)
   !Write(2,"(A,30g13.4)") '!K_filt(-4:5, 10)=',  K_filt(-4:5, 10)
   !Write(2,"(A,30g13.4)") '!K_filt(-4:5, 20)=',  K_filt(-4:5, 20)
   !  K_filt  is a tempory array for storing un-normalized traces for z_f=t_r
   Deallocate( Inp_K, temp, PSF_unfilt)
   !Short circuit the PSF:
   !-------------------------------------------------
   !  Get weights
   Allocate( GrdPwrWght(0:N_tr,0:N_zf), GrdPSFWght(0:N_tr,0:N_zf) )
   irt_psf=(10000-z_low)/del_z  ! take PSF at height of 10. km, If variable timing shifts
   If(irt_psf.gt.N_zf) irt_psf=N_zf
   write(2,*) 'hcnt_l:hcnt_u,N_tr,N_zf',hcnt_l,hcnt_u,N_tr,N_zf, irt_psf


   OPEN(UNIT=14,STATUS='unknown',FILE=TRIM(Base)//'BeamPSF-Dc.dat')
   !k=27 !=5.5km, set earlier
   write(14,"(A,F7.0,i4,A,F6.3,A,F6.3,A)") '! filtered K*PSF at fixed D_c=',k*del_tr + z_low,k, &
         '),  & same for filtered '


   Do i_zf=0,N_zf ! Loop over focal height on the shower axis
      !***************                                Short circuit the PSF to the one at 10 km:
      K_filt(hcnt_l:hcnt_u, i_zf)=K_filt(hcnt_l:hcnt_u, irt_psf)
      !*********************************************
      x=sqrt(SUM(K_filt(hcnt_l:hcnt_u, i_zf)*K_filt(hcnt_l:hcnt_u, i_zf)) )
      If(x.gt.0.)  PSF_trace(hcnt_l:hcnt_u, i_zf)=K_filt(hcnt_l:hcnt_u, i_zf)/x
      !Write(2,"(30g13.4)") '!PSF norm=',i_zf, x, K_filt(-4:5, i_zf)
      !
      Do i_tr=0,N_tr
         GrdPwrWght(i_tr,i_zf)=sqrt(SUM(Grd_K(hcnt_l:hcnt_u, i_tr, i_zf)*Grd_K(hcnt_l:hcnt_u, i_tr, i_zf)))
         GrdPSFWght(i_tr,i_zf)= SUM(Grd_K(hcnt_l:hcnt_u, i_tr, i_zf)*PSF_trace(hcnt_l:hcnt_u, i_zf))
      EndDo
      !write(2,*) '!GrdPSFWght(i_zf,i_zf)',GrdPSFWght(i_zf,i_zf),x
      i_tr=Maxloc(GrdPSFWght(0:N_tr,i_zf),1)-1
      !write(2,*) '!GrdPSFWght', x, i_zf, MaxVal(GrdPSFWght(0:N_tr,i_zf)), i_tr, GrdPwrWght(i_tr,i_zf)
      !If(GrdPwrWght(i_tr,i_zf).gt.2)
      ! write(2,"(50G12.4)") '!Grd_K(hcnt_l:0,,)', Grd_K(hcnt_l:0, i_tr, i_zf)
      !If(GrdPwrWght(i_tr,i_zf).gt.2)
      ! write(2,"(50G12.4)") '!Grd_K(0:hcnt_u,,)', Grd_K(0:hcnt_u, i_tr, i_zf)
      !write(2,"(50G12.4)") '!Grd_K(-2:+2,,)', Grd_K(-2:+2, i_tr, i_zf)
      !If(GrdPwrWght(i_tr,i_zf).gt.2) write(2,"(50G12.4)") '!Grd_K(-2:+2,,)', Grd_K(-2:+2, i_tr-1, i_zf)
      !  Just for plotting -------------
      zeta_f=i_zf*del_z + z_low
      i_tr=k
      write(14,"(30g13.4)") i_zf, zeta_f/1000., GrdPSFWght(i_tr,i_zf)/x &
         , (i_tr*del_tr + z_low - zeta_f)/1000.
      If(i_zf.eq.k) Then
         !Stop 'interp::'
         OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(Base)//'BeamPSF-Db.dat')
         write(4,"(A,F7.0,i4,A,F6.3,A,F6.3,A)") '! filtered K*PSF at fixed D_b=',i_zf*del_z + z_low,i_zf, &
               '),  & same for filtered '
         Do i_tr=0,N_tr
            write(4,"(30g13.4)") i_tr, (i_tr*del_tr + z_low)/1000., GrdPSFWght(i_tr,i_zf)/x &
               , (i_tr*del_tr + z_low - zeta_f)/1000.
         EndDo
      EndIf
      Do if_PSF=1,N_psf
         If( ((zeta_f-Zf_PSF(if_PSF)).ge.-del_z/2.) .and. ((zeta_f-Zf_PSF(if_PSF)).lt.del_z/2.) ) Then
            OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(Base)//'KPSF-'//TRIM(label2)//TRIM(Label3(if_psf))//'.dat') ! E13.4
            write(2,*) 'writing file "',TRIM(Base)//'KPSF-'//TRIM(label2)//TRIM(Label3(if_psf))//'.dat','"'
            write(4,*) '!', x, nudim_new, dt_PSF, hcnt_l, hcnt_u, zeta_f
            flush(unit=2)
            Do i_h=hcnt_l,hcnt_u !   h= (i_h-0.5)*SampleTime_new
               write(4,"(I5,3(1pG13.4))") i_h, PSF_trace(i_h, i_zf), (i_h-0.5)*SamplingTime
            EndDo
            Close(Unit=4)
            !
            OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(Base)//'KPwr-'//TRIM(label2)//TRIM(Label3(if_psf))//'.dat') ! E13.4
            write(2,*) 'writing file "',TRIM(Base)//'KPwr-'//TRIM(label2)//TRIM(Label3(if_psf))//'.dat','"'
            write(4,*) '!', x, N_tr, zeta_f
            flush(unit=2)
            Do i_tr=0,N_tr  !   i_atm=(i_tr+it_min)*dN_tr  ! dK_CohPwr not really used anywhere
               write(4,"(F8.1,3(1pG13.4))") i_tr*del_tr + Z_low, &
                  GrdPwrWght(i_tr,i_zf)/x, GrdPwrWght(i_tr,i_zf)/x, GrdPSFWght(i_tr,i_zf)/x
            EndDo
            Close(Unit=4)
            exit
         EndIf
      EndDo
      !  Just for plotting -------------
   EndDo ! i_zf=0,N_zf ! Loop over focal height on the shower axis
   Close(Unit=14)
   !
   Deallocate( K_filt)
   !-------------------------------------------------
   Return
End Subroutine ReadKernel
!---------------------------
Subroutine TestWeight(Base, GrdPSFWght, GrdPwrWght, del_tr, N_tr, del_z, N_z, z_low, Ix, AtmHei_step, AtmHei_dim, &
         Grd_K, hcnt_l, hcnt_u, ttrace_zf, ih_min, ih_max, PSF_trace )
   Use Systemcommands, only : MakeGLEplots
   Use ShowerData, only : SamplingTime
   use FFT, only : FFTransform_su, DAssignFFT, RFTransform_CF, RFTransform_CF2RT
   Implicit none
   Character(len=25), intent(in) :: Base
   Real(dp), intent(in) :: GrdPSFWght(0:N_tr,0:N_z), GrdPwrWght(0:N_tr,0:N_z), Ix(0:AtmHei_dim)
   Real(dp), intent(in) :: del_tr, del_z, AtmHei_step, z_low
   Integer, intent(in) :: N_tr, N_z, AtmHei_dim
   Real(dp), intent(in) :: Grd_K(hcnt_l:hcnt_u,0:N_tr,0:N_z)
   Integer, intent(in) :: hcnt_l, hcnt_u
   real(dp), intent(in) :: ttrace_zf( ih_min:ih_max,0:N_z)
   Integer, intent(in) :: ih_min, ih_max
   real(dp), intent(in) :: PSF_trace(hcnt_l:hcnt_u, 0:N_z)
   !
   Integer :: i_tr, i_atm, i_z, i, ih1, ih2, Nh, nuh
   Real(dp) :: PSF_I, Pwr_I, t_r, rat, Zeta_max, Zeta_step, R_h(hcnt_l:hcnt_u,0:N_z), A, B, C, D
   Character(len=30) :: PlotType, PlotName, PlotDataFile
   Complex(dp) :: Cnu_R(0:(hcnt_u-hcnt_l)/2), Cnu_E(0:(hcnt_u-hcnt_l)/2)
   !
   !   write(2,*) 'Reading file "',TRIM(Base)//'SrcKrnl-'//TRIM(label2)//TRIM(Label3(if_psf))//'.z','"'
   !write(2,*) del_tr, N_tr, del_z, N_z, z_low
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(Base)//'BeamPwr-tst.dat')
   write(2,*) 'writing file "',TRIM(Base)//'BeamPwr-tst.dat','" ',z_low, N_z
   flush(unit=2)
   write(4,"('!',A7,3g13.4,A4,20g13.4)") 'zeta[km]', 'W_PSF*I_x', 'W_Pwr*I_x', 'I_x' &
               , 'i_z', 'W_PSF(i_tr=', 'N_tr/2-5','....','N_tr/2+5', ',i_z)'
   Do i=0,N_tr
      t_r=i*del_tr + z_low
      i_z= (t_r-z_low)/del_z
      PSF_I=0.
      Pwr_I=0.
      R_h(hcnt_l:hcnt_u,i_z) = 0.
      Do i_tr=1,N_tr
         i_atm= (i_tr*del_tr + z_low)/AtmHei_step
         PSF_I=PSF_I + GrdPSFWght(i_tr,i_z)*Ix(i_atm)*del_tr
         Pwr_I=Pwr_I + GrdPwrWght(i_tr,i_z)*Ix(i_atm)*del_tr
         R_h(hcnt_l:hcnt_u,i_z) = R_h(hcnt_l:hcnt_u,i_z)+ Grd_K(hcnt_l:hcnt_u, i_tr, i_z)*Ix(i_atm)*del_tr
         !write(2,"(2F8.2,3g13.4,20g13.4)") t_r/1000., i_tr*del_tr + z_low, Ix(i_atm), i_tr &
         !      , Grd_K(0:2, i_tr, i_z)
      EndDo
      !write(2,"(A,F7.2,3g13.4,20g13.4)") '!R_h(0:2, i_z)',t_r/1000., R_h(0:2, i_z)
      !write(2,*) 'TestWeight, i_z, PSF*I, Pwr*I:', i_z, PSF_I, Pwr_I
      i_atm = t_r/AtmHei_step
      write(4,"(F7.2,3g13.4,I4,20g13.4)") t_r/1000., PSF_I, Pwr_I, Ix(i_atm) &
               , i_z, (GrdPSFWght(i_tr,i_z), i_tr=N_tr/2-5,N_tr/2+5)
   EndDo
   Close(unit=4)
   !
   !---------------------------
   ! Testing fourier spectrum; kernel*I =?= E_sum
   !
   !ih1=MAX(ih_min, hcnt_l)
   !ih2=MIN(ih_max, hcnt_u)
   !i_z=25
   !write(2,*) '!FFT test:',I_z, i_z*del_z+z_low,ih1,ih2,ih1*SamplingTime,ih2*SamplingTime
   !nuh=(ih2-ih1+1)/2
   !Nh=2*nuh
   !ih1=ih2-Nh+1
   !
   !Call FFTransform_su(Nh)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !Call RFTransform_CF(R_h(ih1:ih2,i_z), Cnu_R(0))  ! transform to frequency
   !Call RFTransform_CF(ttrace_zf( ih1:ih2,i_z), Cnu_E(0))  ! transform to frequency
   !Call DAssignFFT()      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !Do i=0,nuh
   !   write(2,"(I4,3g13.4,20g13.4)") i,500.*i/nuh,i/(nuh*SamplingTime),ABS(Cnu_R(i)),ABS(Cnu_E(i))
   !EndDo
   !---------------------------
   ! Testing kernel*I =?= E_sum
   !ih1=-20
   !ih2=+20
   !write(2,*) '!Grd_K*I test:',ih1,ih2,ih1*SamplingTime,ih2*SamplingTime
   !Do i_z=1,N_z
   !   A=sqrt( SUM(R_h(ih1:ih2,i_z) * R_h(ih1:ih2,i_z)) )
   !   B=sqrt( SUM(ttrace_zf( ih1:ih2,i_z) * ttrace_zf( ih1:ih2,i_z)) )
   !   If(A*B .eq. 0.) Then
   !      write(2,"(A,20(1pg13.4))") '!Grd_K*I AB=0:',I_z, i_z*del_z+z_low, A, B
   !      exit
   !   EndIf
   !   rat=SUM(R_h(ih1:ih2,i_z) * ttrace_zf( ih1:ih2,i_z)) /(A*B)
   !   C= SUM(R_h(ih1:ih2,i_z) * PSF_trace(ih1:ih2,i_z))
   !   D= SUM(ttrace_zf( ih1:ih2,i_z) * PSF_trace( ih1:ih2,i_z))
   !   write(2,"(A,20(1pg13.4))") '!Grd_K*I test:',I_z, i_z*del_z+z_low, A, B, A/B, rat, C, D, C/D
   !EndDo
   !
   ih1=MAX(ih_min, hcnt_l)
   ih2=MIN(ih_max, hcnt_u)
   write(2,*) '!Grd_K*I test:',ih1,ih2,ih1*SamplingTime,ih2*SamplingTime
   Do i_z=1,N_z
      A=sqrt( SUM(R_h(ih1:ih2,i_z) * R_h(ih1:ih2,i_z)) )
      B=sqrt( SUM(ttrace_zf( ih1:ih2,i_z) * ttrace_zf( ih1:ih2,i_z)) )
      If(A*B .eq. 0.) Then
         write(2,"(A,20(1pg13.4))") '!Grd_K*I AB=0:',I_z, i_z*del_z+z_low, A, B
         exit
      EndIf
      rat=SUM(R_h(ih1:ih2,i_z) * ttrace_zf( ih1:ih2,i_z)) /(A*B)
      C= SUM(R_h(ih1:ih2,i_z) * PSF_trace(ih1:ih2,i_z))
      D= SUM(ttrace_zf( ih1:ih2,i_z) * PSF_trace( ih1:ih2,i_z))
      write(2,"(A,20(1pg13.4))") '!Grd_K*I test:',I_z, i_z*del_z+z_low, A, B, A/B, rat, C, D, C/D
   EndDo
   ! ---------------------------------
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(Base)//'BeamPwr-tstK.dat')
   write(2,*) 'writing file "',TRIM(Base)//'BeamPwr-tstK.dat','" ',z_low, N_z
   Zeta_step = N_z * del_z / 10.
   write(4,"(A,17G13.4)") '!    i_h,  Kern*I_x(i_h,zeta_f[km]=',(i_z*Zeta_step, i_z=3,10)
   Do  i=ih1,ih2  ! same format & layout as in Beam_PSF-v20.f90, line 163
      write(4,"(22(1pG13.4))") (i-0.5)*SamplingTime, (R_h(i,NINT(i_z*Zeta_step/del_z)), i_z=1,10) & !, PSF_trace(i)
         , R_h(i,10), ttrace_zf( i,10), R_h(i,N_z-10), ttrace_zf( i,N_z-10)
   EndDo ! i_zf=0,N_zf
   Close(Unit=4)
   !
   Zeta_step = N_z * del_z / 10.
   !If( (17*1000.-z_low)/del_z .gt. N_z) Then
   !   write(2,*) 'Maximum height reaches only till ',N_z*del_z+z_low,'[km], not 17'
   !   write(*,*) 'Maximum height reaches only till ',N_z*del_z+z_low,'[km], not 17'
   !   stop
   !EndIf
   ! TestWeight:          70          70   1000.0000000000000        200.00000000000000        80.000000000000000
   OPEN(UNIT=4,STATUS='unknown',FILE=TRIM(Base)//'BeamPwr-w.dat')
   write(2,*) 'writing file "',TRIM(Base)//'BeamPwr-w.dat','" ',z_low, N_z
   flush(unit=2)
   write(4,*) z_low, Zeta_step, N_psf, ' "'//TRIM(label2)//'" !'
   write(4,*) (' "'//TRIM(Label3(i))//'" ', i=1,N_psf), ' "'//TRIM(label2)//'" !'
   write(4,"('!',(20(1pG13.5)))") 't_r', 'W_PSF(i_tr,','Zeta=z_low', '....', 'z_low +','10*z_step)'
   Do i_tr=0,N_tr
      t_r=i_tr*del_tr + z_low
      !write(4,*) t_r/1000., (GrdPSFWght(i_tr,NINT((i*1000.-z_low)/del_z))*1.e5, i=1,17)
      write(4,"(20(1pG13.5))") t_r/1000., (GrdPSFWght(i_tr,NINT(i*Zeta_step/del_z)), i=0,10)
   EndDo
   Close(unit=4)
   PlotType='Test_Weight'
   PlotName=trim(Base)//'Test_Weight'
   PlotDataFile=trim(Base)
   Call MakeGLEplots(PlotType, PlotName, PlotDataFile) !, Submit, CleanPlotFile
   Return
End Subroutine TestWeight
! ---------------------------------------------
End Module PSF
!==========================================================
Subroutine Interpol(zeta_f,if_psf, t_r, AtmHei_dim, AtmHei_step, dN_tr, ih_1, ih_2, Zf_PSF, Inp_K, RMS_InpK, temp )
   use constants, only : dp
   Implicit none
   Integer, intent(in) :: if_psf, ih_1, ih_2
   Real(dp), intent(in) :: zeta_f, t_r
   Integer, intent(in) :: dN_tr, AtmHei_dim
   Real(dp), intent(in) :: AtmHei_step
   Real(dp), intent(in) :: Zf_PSF(*), Inp_K(ih_1:ih_2, 0:AtmHei_dim/dN_tr,*), RMS_InpK(0:AtmHei_dim/dN_tr,*)
   Real(dp), intent(out) :: temp(ih_1:ih_2)
   Real(dp) :: norm
   Integer :: itf_0, itf_1, padding, NuDim_new
   Real(dp) :: d_zeta, tr_0, tr_1, dtf_0, dtf_1, del_tf, var
   !
   d_zeta=(zeta_f-Zf_PSF(if_psf))/(Zf_PSF(if_psf+1)-Zf_PSF(if_psf))
   del_tf=dN_tr*AtmHei_step  !  step size used on input file
   var=(t_r/zeta_f)
   tr_0=var*Zf_PSF(if_psf)
   tr_1=var*Zf_PSF(if_psf+1)
   !tr_0=(t_r-zeta_f)+Zf_PSF(if_psf)
   !tr_1=(t_r-zeta_f)+Zf_PSF(if_psf+1)
   If(tr_0.le.0.) tr_0=0
   If(tr_1.le.0.) tr_1=0
   itf_0=INT(tr_0/del_tf)  ! step size on input will be smaller than new grid spacing in t_t
   dtf_0=tr_0/del_tf- itf_0
   If(itf_0.ge.AtmHei_dim/dN_tr) Then
      itf_0 = AtmHei_dim/dN_tr-1
      dtf_0 = 1.
   EndIf
   itf_1=INT(tr_1/del_tf)
   dtf_1=tr_1/del_tf - itf_1
   If(itf_1.ge.AtmHei_dim/dN_tr) Then
      itf_1 = AtmHei_dim/dN_tr-1
      dtf_1 = 1.
   EndIf
   temp(ih_1:ih_2)=(1.-d_zeta)*(1.-dtf_0)*Inp_K(ih_1:ih_2, itf_0,  if_PSF) + &
                  (1.-d_zeta)*dtf_0  *   Inp_K(ih_1:ih_2, itf_0+1,if_PSF) + &
                  d_zeta  *   (1.-dtf_1)*Inp_K(ih_1:ih_2, itf_1,  if_PSF+1) + &
                  d_zeta  *    dtf_1  *  Inp_K(ih_1:ih_2, itf_1+1,if_PSF+1)
   norm=sqrt( SUM(temp(ih_1:ih_2)*temp(ih_1:ih_2)) )
   var=(1.-d_zeta)*(1.-dtf_0)*RMS_InpK(itf_0,if_PSF) + (1.-d_zeta)*dtf_0*RMS_InpK(itf_0+1,if_PSF) + &
      d_zeta*(1.-dtf_1)*RMS_InpK(itf_1,  if_PSF+1)     +  d_zeta * dtf_1 *RMS_InpK(itf_1+1,if_PSF+1)
   !write(2,"(30(1pG11.4))") '!Interpol',d_zeta, tr_0, tr_1, dtf_0, dtf_1, if_PSF, itf_0, itf_1, norm, var &
   !   , RMS_InpK(itf_0,if_PSF), RMS_InpK(itf_0+1,if_PSF), RMS_InpK(itf_1,  if_PSF+1),RMS_InpK(itf_1+1,if_PSF+1)
   If(norm.ne.0.) Then
      temp(ih_1:ih_2)=temp(ih_1:ih_2)*var/norm
      norm=var
   EndIf
   Return
End Subroutine Interpol
!----------------------------
Subroutine FoldCurrent1(A, Ix, AI, N_tr, del_tr, z_low)
   use Atmosphere, only : AtmHei_step, AtmHei_dim
   use constants, only : dp
   Implicit none
   Real(dp), Intent(in) :: A(0:N_tr), Ix(0:AtmHei_dim)
   Real(dp), Intent(out) :: AI
   Integer, intent(in) :: N_tr
   Real(dp), Intent(in) :: del_tr, z_low
   Integer :: i_tr, i_atm
   AI=0.
   Do i_tr=0,N_tr
      i_atm= (i_tr*del_tr + z_low)/AtmHei_step
      AI=AI + A(i_tr)*Ix(i_atm)*del_tr
   EndDo
   Return
End Subroutine FoldCurrent1
!----------------------------
Subroutine FoldCurrent2(A, ih1, ih2, Ix, AI, N_tr, del_tr, z_low)
   use Atmosphere, only : AtmHei_step, AtmHei_dim
   use constants, only : dp
   Implicit none
   Real(dp), Intent(in) :: A(ih1:ih2,0:N_tr), Ix(0:AtmHei_dim)
   Real(dp), Intent(out) :: AI(ih1:ih2)
   Integer, intent(in) :: N_tr, ih1, ih2
   Real(dp), Intent(in) :: del_tr, z_low
   Integer :: i_tr, i_atm
   AI(:)=0.
   Do i_tr=0,N_tr
      i_atm= (i_tr*del_tr + z_low)/AtmHei_step
      AI(ih1:ih2)=AI(ih1:ih2) + A(ih1:ih2,i_tr)*Ix(i_atm)*del_tr
   EndDo
   Return
End Subroutine FoldCurrent2
