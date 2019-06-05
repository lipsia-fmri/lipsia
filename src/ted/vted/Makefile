
LDLIBS = -lgsl -lgslcblas -lviaio3.0 -lm -fopenmp -lz

#CFLAGS  = -g
CFLAGS  += -fopenmp

PROG = vted
SRC = vted.c VoxelMap.c EdgeDensity.c ZMatrix.c Histogram.c quantile.c Median.c

OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
