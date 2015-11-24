CC      = gcc-5
CC_MPI  = mpicc
INCS    = -I /System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers/
CFLAGS  = -O3 -mcmodel=medium -fopenmp ${INCS}
LD      = $(CC)
LD_MPI  = $(CC_MPI)
LDFLAGS = $(CFLAGS) -llapack -lclapack -lcblas
LIBS    = -L /System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/ -lpthread