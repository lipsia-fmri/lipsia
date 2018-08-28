
LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm -fopenmp -lz
#CFLAGS  = -g
CFLAGS  += -fopenmp


PROG = vlisa_generic
SRC = vlisa_generic.c ../utils/BilateralFilter.c ../utils/Hotspot.c ../utils/FDR.c ../utils/quantile.c 
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~ ../utils/*.o ${PROG}
