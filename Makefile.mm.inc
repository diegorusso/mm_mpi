CC      = gcc-5
CC_MPI  = mpicc
INCS    = -I$(INCLUDE)
CFLAGS  = -O3 -mcmodel=medium
LD      = $(CC)
LD_MPI  = $(CC_MPI)
LDFLAGS = $(CFLAGS)
LIBS    =  -lpthread
