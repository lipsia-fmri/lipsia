
SRC = nifti1_io.c znzlib.c
OBJ = $(SRC:.c=.o)

TARGET = ${LIPSIA_LIB}/libniftiio.a

default: $(TARGET)

$(TARGET): $(OBJ)
	ar rcs $(TARGET) $(OBJ)

clean:
	rm -f *.o *~ $(TARGET)
