LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm  -fopenmp -lz
#CFLAGS  = -g
CFLAGS  += -fopenmp


PROG = vcylarim_roi
SRC = vcylarim_roi.c ROIstats.c ROIprint.c ROIbeta.c ../cylutils/Conv.c ../cylutils/Cylinders.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~ ${PROG}
