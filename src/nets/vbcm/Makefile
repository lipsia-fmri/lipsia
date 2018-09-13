
LDLIBS = -lgsl -lcblas -latlas -lviaio3.0 -lvia3.0 -lm  -fopenmp -lz
#CFLAGS  = -g
CFLAGS  += -fopenmp


PROG = vbcm

SRC = vbcm.c VoxelMap.c
OBJ=$(SRC:.c=.o)

${PROG}: $(OBJ)

clean:
	-rm -f *.o *~
