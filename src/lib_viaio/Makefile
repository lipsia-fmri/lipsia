CFLAGS+=-fPIC

SRC = Alloc.c Attr.c ConvertC.c ConvertL.c ConvertR.c ConvertO.c Crop.c \
      Dictionary.c Error.c FileIO.c GetVersion.c Image.c ImageType.c\
      List.c Option.c PackData.c Param.c ReadPlain.c Type.c Stats.c \
      Fill.c stringtoken.c Edges.c EdgesType.c Volumes.c VolumesType.c \
      Flip.c Op.c Invert.c History.c BlockFileIO.c \
      GeoInfo.c slicetime.c Transpose.c \
      newfile_io.c ByteSwap.c Nii2Vista.c nifti1_io.c znzlib.c

OBJ = $(SRC:.c=.o)

STATIC = ${LIPSIA_LIB}/libviaio3.0.a
SHARED = ${LIPSIA_LIB}/libviaio3.0.so

default: $(STATIC) $(SHARED)

$(SHARED): $(OBJ)
	gcc -shared -o $@ $^

$(STATIC): $(OBJ)
	ar rcs $@ $^
clean:
	rm -f ${OBJ} $(STATIC) $(SHARED) *~
