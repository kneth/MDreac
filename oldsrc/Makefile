FC	= gfortran -cpp -std=gnu
CC	= gcc
CPP	= /lib/cpp
OPTFLAG	= -O6 -fomit-frame-pointer -fno-strength-reduce -malign-double -funroll-loops
#OPTFLAG	= -g
#OPTFLAG	= -g -qstrict -O3 -qdpc -qhot -qipa
CPPFLAG	= -DPBC_LIST -UVERLET -L/usr/lib
CFLAG	= -I$(PVM_ROOT)/include $(ARCHCFLAGS)
INCLUDE	=
LIB	= -L$(PVM_ROOT)/lib/$(PVM_ARCH)
MLIB	= -L/home/kneth/Numeric/mass
MV	= mv
CP	= cp
RM	= rm -f

PROGRAM	= mdreac parmdreac geninit pargeninit append rateconst sf3d reacener diffcoeff

DIST	= /Users/kneth/bin

all: mdreac geninit

parmdreac: mdreac.f mdreac.inc mdreac.common randgen.f randgen.inc
	$(CPP) -P -DPVM_VER -UGEN_INIT $(CPPFLAG) mdreac.f mdreac.tmp.f
	f2c -w mdreac.tmp.f $(CFLAG)
	$(CC) $(OPTFLAG) $(CFLAG) -o parmdreac mdreac.tmp.c $(LIB) -lfpvm3 -lpvm3 -lf2c -lm
	$(RM) mdreac.tmp.f
	strip parmdreac
	$(MV) parmdreac $(DIST)


mdreac: mdreac.F mdreac.inc mdreac.common randgen.f randgen.inc
	$(FC) -UPVM_VER -UGEN_INIT $(CPPFLAG) $(OPTFLAG) -o mdreac mdreac.F


geninit: mdreac.F mdreac.inc mdreac.common randgen.f randgen.inc
	$(FC) -UPVM_VER -DGEN_INIT $(CPPFLAG) $(OPTFLAG) mdreac.F -o geninit
	strip geninit


pargeninit: mdreac.f mdreac.inc mdreac.common randgen.f randgen.inc
	$(CPP) -P -DPVM_VER -DGEN_INIT $(CPPFLAG) mdreac.f geninit.tmp.f
	$(FC) $(OPTFLAG) $(CFLAG) -o pargeninit geninit.tmp.f $(LIB) $(MLIB) -lpvm3 -lfpvm3
	$(RM) geninit.tmp.f
	strip pargeninit
	$(MV) pargeninit $(DIST)

clean:
	$(RM) core *~ *.o *.tmp.f
