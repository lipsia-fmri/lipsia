LDLIBS = -lgsl -lgslcblas -lviaio3.0 -lvia3.0 -lm -lz
#CFLAGS += -g


PROG = vfdr
SRC = vfdr.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~
