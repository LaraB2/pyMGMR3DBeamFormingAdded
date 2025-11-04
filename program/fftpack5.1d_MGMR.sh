#!/bin/bash
#
#set source = fftpack5.1d  in  ${FFTPackBase}
# Make sure the folder  ${FFTPackBase}/bin/  exists
#source  ${LL_Base}/ShortCuts.sh
   FFTPackBase="/Users/users/scholten/NumLib"  # subfolder /FFTPACK should contain the FFT library, double precision
   export LL_Base="/Users/users/scholten/LOFLI"
   export LL_bin=${LL_Base}/bin          # folder containing f90split.make
   export LL_src=${LL_Base}/FORTRANsrc   # folder containing f90split_OS.f90

cd ${LL_src}
make -f ${LL_scripts}/f90split.make

cd ${FFTPackBase}/FFTPACK/
mkdir temp
cd temp
${LL_bin}/f90split ../fftpack5.1d.f90
#
for FILE in `ls -1 *.f90`;
do
  gfortran -fallow-argument-mismatch -c -O3 $FILE
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done
rm *.f90
# May want to use the following gfortran compiler option:
#-std=std    and use the 'legacy' option
#    Specify the standard to which the program is expected to conform, which may be one 
# of ‘f95’, ‘f2003’, ‘f2008’, ‘f2018’, ‘f2023’, ‘gnu’, or ‘legacy’. 
# The default value for std is ‘gnu’, which specifies a superset of the latest Fortran standard that includes all of the extensions supported by GNU Fortran, although warnings are given for obsolete extensions not recommended for use in new code. The ‘legacy’ value is equivalent but without the warnings for obsolete extensions, and may be useful for old nonstandard programs. The ‘f95’, ‘f2003’, ‘f2008’, ‘f2018’, and ‘f2023’ values specify strict conformance to the Fortran 95, Fortran 2003, Fortran 2008, Fortran 2018 and Fortran 2023 standards, respectively; errors are given for all extensions beyond the relevant language standard, and warnings are given for the Fortran 77 features that are permitted but obsolescent in later standards. The deprecated option ‘-std=f2008ts’ acts as an alias for ‘-std=f2018’. It is only present for backwards compatibility with earlier gfortran versions and should not be used any more. ‘-std=f202y’ acts as an alias for ‘-std=f2023’ and enables proposed features for testing Fortran 202y. As the Fortran 202y standard develops, implementation might change or the experimental new features might be removed.

#
ar cr libfftpack5.1d.a *.o
rm *.o
rm *.out
#
mv libfftpack5.1d.a ${FFTPackBase}/bin
cd ..
rmdir temp
#
echo "Library installed as ${FFTPackBase}/bin/libfftpack5.1.a."
