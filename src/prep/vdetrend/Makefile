
LDLIBS =  -lviaio3.0 -lm -lgsl -lgslcblas -lvia3.0
#CFLAGS = -g

PROG = vdetrend
SRC = vdetrend.c Detrend.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~
