!---------------------------------------------------------
      module KernelFilter
      contains
Subroutine FreqFilterKernel(K_in, NuDim_in, N_tr, SampleTime_old, padding, SampleTime_new, nu_min, nu_max, &
      K_h, NuDim_new)
   use constants, only : dp, pi, c_l
   use FFT, only : FFTransform_su,DAssignFFT, RFTransform_CF, RFTransform_CF2RT
   implicit none
   Integer, intent(in) :: N_tr
   Integer, intent(in) :: NuDim_in ! half the length of input trace [samples]
   REAL(dp), intent(in) :: K_in(1-NuDim_in:NuDim_in,0:N_tr)
   real(dp), intent(in) :: SampleTime_old ! step size before FreqFilt [m]
   integer, intent(in) :: padding ! number of zero samples added at begin and end input trace
   real(dp), intent(in) :: SampleTime_new  ! sampling time trace after FreqFilt [ns]; typically 5 for LOFAR
   real(dp), intent(in) :: nu_min, nu_max ! block filter frequencies [MHz]
   Real(dp), allocatable, intent(out) :: K_h(:,:)     ! K_h(1-nudim_new:nudim_new,0:N_tr)
   Integer, intent(out) :: NuDim_new ! half the length of output trace [samples]
   integer :: i_tr,i
   integer :: tDim_in
   real(dp), allocatable :: tTrace(:)  ! scratch array
   integer :: FF_dim, i_nu_ini, i_nu_max  !, tTrace_Offset_dwn
   real(dp) :: NuSample_new
   integer :: tdim_new ! length of input trace [samples]
   complex(dp), allocatable :: filt(:)
   complex(dp), allocatable :: Cnu(:)
   complex(dp), allocatable :: K_nu(:,:)
   logical, save :: First=.true.
    !
    tDim_in=NuDim_in*2
    FF_dim=(NuDim_in+padding)*2 ! dimension time trace before frequency filtering, padding 'padding' number of zeros at begin and end
    !
    if(SampleTime_new.lt.SampleTime_old) then ! new should be down sampled from input, not up
        write(2,*) 'SampleTime_new should not be smaller than SampleTime_old=',SampleTime_new ,SampleTime_old
        Stop 'FreqFilterKernel:SampleTime_new'
    endif
    tdim_new=FF_dim*SampleTime_old/SampleTime_new  ! cover the same length in time [m] for the two traces
    nudim_new=tdim_new/2  ; tdim_new=2*NuDim_new
    NuSample_new=c_l/(tdim_new*SampleTime_new)    ! in [GHz]
    !write(2,"(A,I4,A,F6.3,A,F7.2,A,F7.2)") 'padding=',padding, ', sampling time=', SampleTime_old, '[m], input trace length=', &
    !  tDim_in*SampleTime_old, '[m] after padding:', FF_dim*SampleTime_old
    i_nu_ini= NINT( nu_min/(1000.*NuSample_new) )  ; i_nu_max=NINT( nu_max/(1000.*NuSample_new) )
    if(i_nu_ini .lt. 0) i_nu_ini=0
    if(i_nu_max .ge. NuDim_new) i_nu_max=NuDim_new
    If(First) Then
       write(2,"(A,F5.3,A,F7.2,A,F7.2,A,F5.1,A,F7.1,A,F5.1,A)") 'Down-sampling to samples of ',SampleTime_new/c_l,&
         '[ns] and trace length of ',tdim_new*SampleTime_new,'[m]=',tdim_new*SampleTime_new/c_l, &
         '[ns], frequency filter from ',i_nu_ini*NuSample_new*1000.,' till ',i_nu_max*NuSample_new*1000., &
         '[MHz], steps of ',NuSample_new*1000.,'[MHz]'
      First=.false.
    EndIf
   !
   allocate(filt(0:nudim_new))
   allocate(tTrace(1:FF_dim))
   allocate(Cnu(0:NuDim_in+padding))
   allocate(K_nu(0:nudim_new,0:N_tr))
   If(.not. Allocated(K_h)) allocate(K_h(1-nudim_new:nudim_new,0:N_tr))
   !write(2,*) '!nudim_new:', nudim_new,nudim_in, ',  i_nu_ini:', i_nu_ini,i_nu_max
   !
   call FFTransform_su(FF_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   filt(:)=0.      ! Frequency filter
   filt(i_nu_ini: i_nu_max)=1.
   !
!   write(2,*) 'FreqFilterKernel:', FF_dim, N_tr, SampleTime_old, ' nu_max_old=',0.5*c_l/SampleTime_old
!   write(2,*) 'FreqFilterKernel:', tdim_new*SampleTime_new, FF_dim*SampleTime_old
!   write(2,*) (Real(NuSample_new*i_tr),i_tr =0,nudim_new)
   Do i_tr =0,N_tr
      Do i=0,padding-1
         tTrace(i+1)=i*K_in(1-NuDim_in,i_tr)/padding
         tTrace(FF_dim-i)=i*K_in(NuDim_in,i_tr)/padding
      EndDo
      !tTrace(1:padding)=0.0d0
      tTrace(padding+1:padding+tDim_in)=K_in(1-NuDim_in:NuDim_in,i_tr)
      !tTrace(tDim_in+padding+1:FF_dim)=0.0d0
   !   write(2,"(A,20(1pg13.4))") '!FreqFilterKernel1:',i_tr,tTrace(padding-1:padding),';',tTrace(padding+1:padding+3)
   !   write(2,"(A,20(1pg13.4))") '!FreqFilterKernel2:',tTrace(FF_dim-padding-3:FF_dim-padding-1),';' &
   !         ,tTrace(FF_dim-padding:FF_dim-padding+1)
      Call RFTransform_CF(tTrace,Cnu(0))  ! transform to frequency
      K_nu(0:nudim_new, i_tr)=filt(0:nudim_new)*Cnu(0:nudim_new) ! apply filter and reduce nu-trace length
!      If(i_tr.lt.5) Then
!         write(2,*) real(abs(Cnu(0:nudim_new)))
!      EndIf
   EndDo !  i_tr =1,ObsDist_dim
   Call DAssignFFT()      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
!   write(2,*) 'FreqFilterKernel:', tdim_new, NuSample_new, SampleTime_new, ' nu_max_new=',0.5*c_l/SampleTime_new
   Call FFTransform_su(tdim_new)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   Do i_tr =0,N_tr
      Call RFTransform_CF2RT(K_nu(0, i_tr),K_h(1-nudim_new, i_tr)) !transform to real-time trace
   EndDo !  i_tr =1,ObsDist_dim
   Call DAssignFFT()      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
   Deallocate(tTrace,filt,Cnu)
!   Flush(Unit=2)
!   stop
   Return
   End subroutine FreqFilterKernel
end module KernelFilter
!-------------------------------------------------------------------------
