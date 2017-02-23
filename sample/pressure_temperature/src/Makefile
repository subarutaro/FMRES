CXX=g++
CXXFLAGS+= -std=c++11 -g
#CXXFLAGS+= -Wall
INC+= -I../../../src/

CUDAPATH:=/usr/local/cuda
NVCC=$(CUDAPATH)/bin/nvcc
NVCCFLAGS=
LIB+= -L$(CUDAPATH)/lib
INC+= -I$(CUDAPATH)/include

SRC= main.cu
HED= molecular_dynamics.cuh integrator.cuh constant.cuh
OUT= remd.out

all:	$(OUT)
$(OUT):	$(SRC) $(HED)
	$(NVCC) $(NVCCFLAGS) $(CXXFLAGS) $(INC) $(LIB) $(SRC) -o $(OUT)
clean:
	rm $(OUT)

