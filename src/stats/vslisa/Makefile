LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm -fopenmp -lz
CFLAGS  += -fopenmp
#CFLAGS  += -g

PROG = vslisa
SRC = vslisa.c gauss.c pseudoinv.c GLM.c ReadData.c HemoModel.c Covariates.c \
../utils/Hotspot.c ../utils/quantile.c ../utils/FDR.c ../utils/BilateralFilter.c 


OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
