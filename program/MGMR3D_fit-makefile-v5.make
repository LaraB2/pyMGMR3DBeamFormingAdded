FC = gfortran
FCFLAGS = -ggdb -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -finit-real=inf

LL_Base = ../LightningPros
Lib_Base = ../Library/bin
F_NLSOL = $(LL_Base)
F_AntFie = $(LL_Base)

EXECUTE = MGMR3D_fit-v5
LIBRARY = $(Lib_Base)/libfftpack5.1d.a

# Only link the files that are NOT included via 'include' statements
OBJS = constants.o Atmosphere.o nl2sol.o AntFuncCnst.o AntFunct.o MGMR3D_fit-v5.o

all: $(EXECUTE)

$(EXECUTE): $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXECUTE) $(OBJS) -lm $(LIBRARY)

constants.o: Constants.f90
	$(FC) $(FCFLAGS) -c Constants.f90 -o constants.o

Atmosphere.o: Atmosphere.f90 constants.o
	$(FC) $(FCFLAGS) -c Atmosphere.f90 -o Atmosphere.o

nl2sol.o: $(F_NLSOL)/nl2sol.f90
	$(FC) $(FCFLAGS) -c $(F_NLSOL)/nl2sol.f90 -o nl2sol.o

AntFuncCnst.o: $(F_AntFie)/AntFuncCnst.f90
	$(FC) $(FCFLAGS) -c $(F_AntFie)/AntFuncCnst.f90 -o AntFuncCnst.o

AntFunct.o: $(F_AntFie)/AntFunct.f90 AntFuncCnst.o
	$(FC) $(FCFLAGS) -c $(F_AntFie)/AntFunct.f90 -o AntFunct.o

# The main file compiles and brings in the included logic
MGMR3D_fit-v5.o: MGMR3D_fit-v5.f90 constants.o Atmosphere.o
	$(FC) $(FCFLAGS) -c MGMR3D_fit-v5.f90 -o MGMR3D_fit-v5.o

clean:
	rm -f *.mod *.o $(EXECUTE)