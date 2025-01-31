LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm  -fopenmp -lz
#CFLAGS  = -g
CFLAGS  += -fopenmp



PROG = vcylarim
SRC = vcylarim.c LayerGLM.c ../cylutils/EHist.c ../cylutils/Cylinders.c ../cylutils/Conv.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~ ${PROG}
