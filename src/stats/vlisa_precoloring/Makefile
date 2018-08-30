LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm -fopenmp -lz
CFLAGS  += -fopenmp
#CFLAGS  += -g

PROG = vlisa_precoloring
SRC = vlisa_precoloring.c gauss.c pseudoinv.c GLM.c Design.c ReadData.c HemoModel.c Covariates.c \
../utils/Hotspot.c ../utils/quantile.c ../utils/FDR.c ../utils/BilateralFilter.c 


OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
