#!/bin/bash
# 
# 
source ${PWD}/../MGMR_SC.sh

#------- Run MGMR3D to generate a sample shower profile ----------------------------
export LIBRARY="${FFTlib}"
source ${M_ProgDir}/RunMGMR.sh MGMR3D2_fit-v10 MGMR3D2-SKA.in

#exit
#--------------- Run Beaming Kernel generator ----------------------------
export LIBRARY="${FFTlib} ${BLASlib}"
source ${M_ProgDir}/RunMGMR.sh BeamForm-GenerateKernel Kernel-SKA.in

#exit
#--------------- Run Beaming to extract source currents---------------------------
export LIBRARY="${FFTlib} ${BLASlib}"
source ${M_ProgDir}/RunMGMR.sh BeamForm-ExtractCurrent BeamCurr-SKA.in
exit
   --- ---------- -------- ------------- -   nohup ./NewBeamForm-SKA.sh  >NewBeamForm-SKA.log 2>&1  & 


# read -rsp $'Press enter to continue...\n'
# cd plot
#gle -d pdf -o ../FitStokes.pdf ${M_ProgDir}FitStokes.GLE ${RunFolder}/plot/FitResult
#gle -d jpg -r 200 -o ../FitStokes-map.jpg ${ProgDir}FitStokes-map.GLE ${RunFolder}/plot/
#gle -d pdf -o ../sh-current.pdf ${M_ProgDir}sh-current.GLE ${RunFolder}/plot/
#exit

