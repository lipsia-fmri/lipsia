
LDLIBS = -lfftw3 -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm

PROG = vpreprocess
SRC = vpreprocess.c Convolve.c FreqFilter.c GaussConv.c MinVals.c SpatialFilter.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
