###################################
#   C++ sample program
###################################
#
CXX	= g++-8
#
CXXFLAGS  = -O3 
#
#
EX27_OBJ=ex27.o cg.o
EX28_OBJ=ex28.o cg.o
EX29_OBJ=ex29.o cg.o
EX29b_OBJ=ex29b.o cg_omp.o
#
TARGET= ex27 ex28 ex29 ex29b
#
all:	$(TARGET)	

ex27	:$(EX27_OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $(EX27_OBJ)

ex28	:$(EX28_OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $(EX28_OBJ)

ex29	:$(EX29_OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $(EX29_OBJ)

ex29b	:$(EX29b_OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $(EX29b_OBJ) -fopenmp

###########################################
cg.o:	cg.cpp cg.h
	$(CXX) $(CXXFLAGS) -c cg.cpp

cg_omp.o:	cg_omp.cpp cg_omp.h
	$(CXX) $(CXXFLAGS) -c cg_omp.cpp -fopenmp

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	rm -f *.o $(TARGET) 

