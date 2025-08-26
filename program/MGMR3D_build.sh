#!/bin/bash
# 
# 
#read -rsp $'Press enter to continue...\n'
# gfortran -c ~/NumLib/LSQ/nl2sol.f90 
# cp nl2sol.o ~/NumLib/LSQ/nl2sol.o

export MG_Base="~/pyMGMR/"
#ProgDir="~/../olaf/"
FFD="${MG_Base}/program"
cd "${FFD}"
#echo ${FFD}
#
make -f "${FFD}/MGMR3D_fit-makefile-v5.make"
# rm  .mod
#read -rsp $'Press enter to continue...\n'
