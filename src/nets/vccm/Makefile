
LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm

PROG = vccm
SRC = vccm.c CCC.c Kendall.c OCCC.c
OBJ=$(SRC:.c=.o)
HDRS = 

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
