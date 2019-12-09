###################################
#   C++ sample program
###################################
#
#CXX	= g++
CXX	= icpc 
#
CXXFLAGS  = -O3 -DNDEBUG -funroll-loops 
#CXXFLAGS  = -g 
#
# Linux
BLASLIB	= -lblas -latlas
LAPACKLIB = -llapack -llapack_atlas 
#
ifeq "$(shell uname)" "Darwin"
BLASLIB	= -lblas
LAPACKLIB = -llapack
endif
#
# FOR INTEL MKL ON SHPYNX
ifeq "$(shell hostname)" "sphynx.ccs.tulane.edu" 
MKLROOT=/usr/local/opt/intel/composerxe/mkl
endif
#
ifneq "$(MKLROOT)" ""
BLASLIB	=-L$(MKLROOT)/lib/intel64/ -DMKL -mkl
CXXFLAGS  += -DMKL
LAPACKLIB =
endif
#
#
EX73_OBJ=ex73.o stokeslet2d.o
EX74_OBJ=ex74.o stokeslet2dLapack.o stokeslet2d.o
EX75_OBJ=ex75.o stokeslet2dEigen3.o 
#
TARGET=ex73 ex74 ex75
#
all:	$(TARGET)	

ex73	:$(EX73_OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $(EX73_OBJ)

ex74	:$(EX74_OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $(EX74_OBJ) $(BLASLIB) $(LAPACKLIB)

ex75	:$(EX75_OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $(EX75_OBJ) 

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	rm -f *.o $(TARGET) 

