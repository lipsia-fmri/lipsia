LDLIBS = -lgsl -lgslcblas -lviaio3.0 -lvia3.0 -lm -lz -fopenmp
#CFLAGS  = -g
CFLAGS  += -fopenmp

PROG = vsml
SRC = vsml.c ReadFiles.c PLS1.c Loop.c SelProb.c Norms.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~  ${PROG}
