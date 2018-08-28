LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm -fopenmp -lz
CFLAGS  += -fopenmp
#CFLAGS  += -g

PROG = vlisa_prewhitening
SRC = vlisa_prewhitening.c Design.c ReadData.c HemoModel.c Covariates.c \
WhiteGLM.c whitecov.c whitecov2.c dfs.c prewhite.c gsl_futils.c \
../utils/Hotspot.c ../utils/quantile.c ../utils/FDR.c ../utils/BilateralFilter.c 


OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
