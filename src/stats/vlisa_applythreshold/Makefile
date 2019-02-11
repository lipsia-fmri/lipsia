LDLIBS = -lgsl -lgslcblas -lvia3.0 -lviaio3.0 -lm -lz
#CFLAGS  += -g

PROG = vlisa_applythreshold
SRC = vlisa_applythreshold.c


OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
