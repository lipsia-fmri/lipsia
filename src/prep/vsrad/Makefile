
LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm -lz

PROG = vsrad
SRC = vsrad.c SRAD.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
