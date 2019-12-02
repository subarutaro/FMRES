CXX?=g++
CXXFLAGS+= -O2 -std=c++11
CXXFLAGS+= -fopenmp
INC= -I../../../src/

HED= ../../../src/array.hxx ../../../src/wham.hxx condition.hxx

all:	wham.out
wham.out:	main.cxx $(HED)
	$(CXX) $(CXXFLAGS) $(INC) main.cxx -o wham.out
clean:
	rm wham.out
