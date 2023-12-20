LDLIBS = -lgsl -lgslcblas -lviaio3.0 -lvia3.0 -lm -lz -fopenmp
#CFLAGS  = -g
CFLAGS  += -fopenmp

PROG = vsml_statistics
SRC = vsml_statistics.c ReadFiles.c Folds.c PLS1.c Confounds.c Loop.c Fit.c SelProb.c Norms.c Performance.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~  ${PROG}
