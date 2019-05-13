
CFLAGS+=-fPIC
#LDLIBS=-lgsl -lviaio3.0
LDLIBS=

SRC = Aniso2d.c Aniso3d.c Bicubic.c Binarize.c Binmorph3d.c Border3d.c BorderPoint.c \
      CDT3d.c ChamferDist3d.c Contrast.c Convolve.c DeleteSmall.c Dist2d.c \
      EuclideanDist3d.c GaussFilter.c GenusLee.c GreyMorph3d.c Label2d.c Label3d.c Magnitude.c\
      MatrixInverse.c Median.c NNSample3d.c NNScale3d.c Pixel.c QuickMorph.c\
      Rotate2d.c RotationMatrix.c Sample2d.c Sample3d.c Scale2d.c Scale3d.c \
      SelectBig.c ShapeMoments.c Shear.c SimplePoint.c Skel2d.c Skel3d.c Smooth3d.c \
      Spline.c Thin3d.c Topoclass.c VCheckPlane.c VolumesOps.c VPoint_hpsort.c \
      Contrast.c RegistrationUtils.c Resample.c StatsConversions.c FreqFilter.c

OBJ = $(SRC:.c=.o)

STATIC = ${LIPSIA_LIB}/libvia3.0.a
SHARED = ${LIPSIA_LIB}/libvia3.0.so

default: $(STATIC) $(SHARED)

$(SHARED): $(OBJ)
	$(CC) -shared -o $@ $^ $(LDFLAGS) $(LDLIBS)

$(STATIC): $(OBJ)
	ar rcs $@ $^

clean:
	rm -f ${OBJ} $(STATIC) *~
