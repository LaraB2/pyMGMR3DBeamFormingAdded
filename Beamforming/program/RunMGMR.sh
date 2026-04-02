#!/bin/bash
#
#
Prog="$1"
echo "Program: $1"
echo "Input from: $2"
echo "program run is: ${Prog} with argument $2"
#
cd ${M_ProgDir}
gfortran  ${FCFLAGS} -fopenmp -o ${Prog} ${Prog}.f90 ${LIBRARY} 
rm *.mod
#
cd ${RunFolder}
export input=' <'$2
echo "executing command: '${M_ProgDir}/${Prog}  ${input} '"
${M_ProgDir}/${Prog}  <$2   # This is apparently very different from    ${LL_bin}/${Prog}  ${input}
#
