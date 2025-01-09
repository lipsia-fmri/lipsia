LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm -lz
#CFLAGS  = -g


PROG = vrim
SRC = vrim.c Border3d.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~ ${PROG}
