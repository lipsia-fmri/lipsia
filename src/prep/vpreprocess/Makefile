
LDLIBS = -lvia3.0 -lviaio3.0 -lgsl -lgslcblas -lm -lz

PROG = vpreprocess
SRC = vpreprocess.c Gauss4d.c GaussFilter.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
