LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm  -fopenmp -lz
#CFLAGS  = -g
#CFLAGS  += -fopenmp



PROG = vcylarim_probe
SRC = vcylarim_probe.c ../cylutils/Conv.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~ ${PROG}
