
LDLIBS = -lviaio3.0 -lm -lgsl -lgslcblas -lz -lm 
#CFLAGS += -g

PROG = vnifti
SRC = vnifti.c Vista2Nii.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
