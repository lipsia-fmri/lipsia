LDLIBS = -lgsl -lgslcblas -lviaio3.0 -lm -lz
#CFLAGS += -g

PROG = vapplymask
SRC = vapplymask.c
OBJ=$(SRC:.c=.o)
HDRS = 

${PROG}: ${OBJ}

clean:
	-rm -f *.o *~
