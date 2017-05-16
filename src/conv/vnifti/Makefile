
LDLIBS = -lgsl -lgslcblas -lniftiio -lviaio3.0 -lvia3.0 -lm -lz
#CFLAGS = -g

PROG = vnifti
SRC = vnifti.c Vista2Nii.c Nii2Vista.c ByteSwap.c newfile_io.c
OBJ=$(SRC:.c=.o)

${PROG}: ${OBJ}

clean:
	-rm -f ${PROG} *.o *~
