
LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm -lz


PROG = vreadconnectome
SRC = vreadconnectome.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~  ${PROG}
