Module Systemcommands
contains
Subroutine MakeGLEplots(PlotType, PlotName, PlotDataFile, Submit, CleanPlotFile, mode)
   character(LEN=*), intent(in), optional :: PlotType, PlotName, PlotDataFile
   character(LEN=*), intent(inout), optional :: CleanPlotFile
   Logical, intent(in), optional :: Submit
   Character(len=*), intent(in), optional :: mode  ! interesting only at the first call
   Character(len=40), save :: BatchFile=''
   character*100 :: shellin
   Character(len=26), save :: CallGLE='call GLE -d jpg -r 300 -o '  ! space at end is important
   Character(len=5), save :: ext='.jpg ' ! space at end is important
   Character(len=13) :: Sys
   Character(len=13), save :: GLEFolder, Wfolder
   Logical, save :: Windows=.false.
   Character(len=1), save :: FChar='/'
   Integer :: i, k, StrL
   !
   !inquire(unit=10, opened=itsopen)  ; if ( itsopen ) then okay  ; logical itsopen
   If(BatchFile.eq. '' ) Then  ! Open unit 10
      CALL GET_ENVIRONMENT_VARIABLE("OS", Sys) !
      !write(2,*) 'MakeGLEplots:Sys=', Sys
      If(Sys(1:7).eq.'Windows') Windows=.true.
      If( present(PlotName) ) Then
         StrL=LEN_TRIM(PlotName)
         k=1
         Do i=1,StrL
            If((PlotName(i:i) .eq. '/') .or. (PlotName(i:i) .eq. '\')) Then
               exit
            Else
               k=i
            EndIf
         EndDo
         If(k.gt.30) k=30.
         BatchFile='Afig'//TRIM(PlotName(1:k))
      Else
         BatchFile='Afig'
      EndIf
      If( present(mode) ) Then
         BatchFile=TRIM(BatchFile)//TRIM(mode)
      EndIf
      If(windows) then
         FChar='\'
         CallGLE='call GLE -d jpg -r 300 -o '  ! space at end is important
         ext='.jpg '
         GLEFolder='%GLEscr%\'
         Wfolder='%Workdir%/'
         BatchFile=TRIM(BatchFile)//'.bat'
         OPEN(UNIT=10,STATUS='REPLACE',ACTION='WRITE',FILE=TRIM(BatchFile) )
      Else
         CallGLE='gle -d pdf -o '  ! space at end is important
         ext='.pdf '
         GLEFolder='${M_GLEscr}/'
         Wfolder='${Workdir}/'
         BatchFile=TRIM(BatchFile)//'.sh'
         Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE=TRIM(BatchFile) )   ! It is assumed that this scipt is launched in the FlashFolder
         Write(10,"(7(A,/) )") &
            '#!/bin/bash', &   !  '#!/bin/bash -v', &
            '#', '#',  'source ${MGMR_Base}/MGMR_SC.sh'
      EndIf
   EndIf
   !
   !write(2,*) 'MakeGLEplots:BatchFile=', BatchFile
   If(present(PlotType)) Then  ! PlotName, PlotDataFile
      !write(2,*) 'MakeGLEplots:PlotType=', PlotType
      !write(2,*) 'MakeGLEplots:PlotName=', PlotName
      !write(2,*) 'MakeGLEplots:PlotDataFile=', PlotDataFile
      If(Windows) Then
            Write(10,"(4(A,1x))") TRIM(CallGLE), trim(PlotName)//ext,  TRIM(GLEFolder)//trim(PlotType)//".gle", &
                  TRIM(Wfolder)//trim(PlotDataFile)
      Else
            Write(10,"(A)") &
            "gle -d pdf -o "//trim(PlotName)//".pdf "//TRIM(GLEFolder)//trim(PlotType)//".gle "&
                  //TRIM(Wfolder)//TRIM(PlotDataFile)
      EndIf
   EndIf
   !
   If(present(CleanPlotFile)) Then  ! PlotName, PlotDataFile
      StrL=LEN_TRIM(CleanPlotFile)
      Do i=1,StrL
         If((CleanPlotFile(i:i) .eq. '/') .or. (CleanPlotFile(i:i) .eq. '\')) Then
            CleanPlotFile(i:i)=FChar
         EndIf
      EndDo
      If(Windows) Then
        write(10,"('del ',A)") TRIM(CleanPlotFile)
      Else
        write(10,"('rm ',A)") TRIM(CleanPlotFile)
      EndIf
   EndIf
   Flush(unit=10)
   !
   If(present(Submit)) Then
      If(Submit) Then
         If(windows) then
            Close(unit=10)
            shellin = TRIM(BatchFile)
            CALL system(shellin)
            !write(2,*) 'system command submitted'
            Write(2,*) 'GLE should be running now'
         Else
            Close(unit=10)
            shellin = 'chmod 755 '//TRIM(BatchFile)
            CALL system(shellin)
            shellin = 'nohup ./'//TRIM(BatchFile)//'  >'//TRIM(BatchFile)//'.log 2>&1  & '
            !shellin = 'source ./'//TRIM(BatchFile)  ! does not work
            !shellin = 'source '//TRIM(BatchFile)
            CALL system(shellin)
            write(2,*) 'shellin:',shellin
            Write(2,*) TRIM(BatchFile)//'  is submitted'
            Write(*,*) TRIM(BatchFile)//'  is submitted'
         EndIf
         BatchFile=''
      EndIf
   EndIf
   Return
End Subroutine MakeGLEplots
!======================
Subroutine CreateNewFolder(Folder, Message)
! Search for folder definition in FileName and create the folder when such is present
   Implicit none
   Character(len=*), intent(inout) :: Folder
   Character(len=*), intent(out) :: Message
!   Character(len=*), intent(in) :: FileName
   Integer :: i, StrL, nxx
   Character(len=10) :: Sys
   Logical :: Windows
   Character(len=1) :: FChar='/'
   CALL GET_ENVIRONMENT_VARIABLE("OS", Sys) !
   If(Sys(1:7).eq.'Windows') Windows=.true.
   If(Windows) FChar='\'
!   !
!   Call GETCWD(FlashFolder)  ! Get current working directory
   StrL=LEN_TRIM(Folder)
   If((Folder(StrL:StrL) .eq. '/') .or. (Folder(StrL:StrL) .eq. '\')) Then
      Folder(StrL:StrL)=' '
      StrL=StrL-1
   EndIf
   If(StrL.le.0) Then
      Write(2,*) '*** All files will be placed in top working directory ***'
      Return
   EndIf
   If(Windows) then
      Call EXECUTE_COMMAND_LINE("mkdir "//TRIM(Folder) , WAIT=.true., EXITSTAT=nxx, CMDSTAT=i)
      !Call EXECUTE_COMMAND_LINE("mkdir /? "//TRIM(Folder) , WAIT=.true., EXITSTAT=nxx, CMDSTAT=i)
   Else
      Call EXECUTE_COMMAND_LINE("mkdir -v  "//TRIM(Folder) , WAIT=.true., EXITSTAT=nxx, CMDSTAT=i)
   EndIf
      !Call EXECUTE_COMMAND_LINE("touch "//'files/'//TRIM(Simulation)//'_Structure.dat' , WAIT=.true., &
      !      EXITSTAT=nxx, CMDSTAT=i, CMDMSG=lname)
   If(nxx.eq.0 ) then
      Message= "Created folder: "//TRIM(Folder)
   Else
      Message= 'folder "'// TRIM(Folder)//'" probably existed already'
   EndIf
   Folder(StrL+1:StrL+1)='/'
   !flush(unit=2)
   Return
End Subroutine CreateNewFolder
End Module Systemcommands
