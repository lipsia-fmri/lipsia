
LDLIBS =  -lviaio3.0 -lm -lgsl -lgslcblas -lvia3.0 -lz
#CFLAGS = -g

PROG = vdetrend
SRC = vdetrend.c SubtractMovingAverage.c SubtractPolynomial.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~
