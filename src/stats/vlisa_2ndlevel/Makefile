
LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm -fopenmp -lz
#CFLAGS  = -g
CFLAGS  += -fopenmp


PROG = vlisa_2ndlevel
SRC = vlisa_2ndlevel.c ../utils/BilateralFilter.c ../utils/Hotspot.c ../utils/FDR.c ../utils/ttest.c ../utils/quantile.c 2ndLevel.c pseudoinv.c Read2ndLevelDesign.c genperm.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~ ../utils/*.o ${PROG}
