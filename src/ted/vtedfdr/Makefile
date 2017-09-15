
LDLIBS = -lgsl -lgslcblas -lviaio3.0 -lm -lz

PROG = vtedfdr
SRC = vtedfdr.c

OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o ${PROG} *~
