
LDLIBS = -lfftw3 -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm -lz

PROG = vpreprocess
SRC = vpreprocess.c FreqFilter.c Gauss4d.c GaussFilter.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
