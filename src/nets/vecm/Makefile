LDLIBS = -lgsl -lgslcblas -lviaio3.0 -lm -fopenmp -lz
#CFLAGS  = -g
CFLAGS  += -fopenmp


PROG = vecm
SRC = vecm.c ECM_IO.c MatrixECM.c MatrixProjection.c Irreducible.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~ ${PROG}
