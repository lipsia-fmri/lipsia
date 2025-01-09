LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm -lz
#CFLAGS  = -g


PROG = vmetric
SRC = vmetric.c EDist.c CDist3d.c CDT3d.c ../cylutils/Conv.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~ ${PROG}
