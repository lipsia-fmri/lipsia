
LDLIBS = -lvia3.0 -lviaio3.0 -lm -lgsl -lgslcblas 

PROG = vcuttrials
SRC = vcuttrials.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~ ${PROG}
