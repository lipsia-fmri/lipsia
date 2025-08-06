
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../cylutils/cyl.h"

double LaminarZMax(gsl_vector *y)
{
  size_t i;
  double zpos=0,zneg=0;
  for (i=0; i<y->size; i++) {
    if (y->data[i] > zpos) zpos = y->data[i];
    if (y->data[i] < zneg) zneg = y->data[i];
  }
  double z=zpos;
  if (fabs(zneg) > z) z = zneg;
  return z;
}
