/*! \file
  Select most frequent voxel label

The input image must be of type "ubyte" or "short".
The output image will be of type "bit" and contain
only those pixels which belong to the most frequent
class. The maximum number of classes allowed is currently
set to 8192.

\par Author:
Gabriele Lohmann, MPI-CBS
*/

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <stdio.h>
#include <stdlib.h>


#define FillTable(type) \
{ \
  type *pp = (type *) VImageData(src);	 \
  for (j=0; j<VImageNPixels(src); j++) { \
    i = (long)(*pp++); \
    table[i]++; \
  } \
}


VImage VDeleteSmall (VImage src,VImage dest,int msize)
{
  long i,j;
  int b,r,c,nbands,nrows,ncols;
  VRepnKind repn;

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);
  repn   = VPixelRepn(src);


  /* get largest label */
  double u=0,umax=0;
  for (b=0; b<VImageNBands(src); b++) {
    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {
	u = VGetPixel(src,b,r,c);
	if (u > umax) umax = u;
	if (u < 0) VSetPixel(src,b,r,c,0);
      }
    }
  }
  int maxlabel = (int)(umax + 2.0);
  size_t *table = (size_t *) VCalloc(maxlabel,sizeof(size_t));
  if (table == NULL) VError("VSelectBig: err allocating table");
  for (i=0; i<maxlabel; i++) table[i] = 0;


  /* get histogram of component sizes */
  switch (repn) {
  case VUByteRepn:
    FillTable(VUByte);
    break;
  case VShortRepn:
    FillTable(VShort);
    break; 
  case VUShortRepn:
    FillTable(VUShort);
    break;
  case VIntegerRepn:
    FillTable(VInteger);
    break;
  case VUIntegerRepn:
    FillTable(VShort);
    break;
  case VLongRepn:
    FillTable(VLong);
    break;
  case VULongRepn:
    FillTable(VULong);
    break;

  default:
    VError(" VDeleteSmall: illegal pixel repn");
  }


  /* create output image */
  dest = VSelectDestImage("VDeleteSmall",dest,nbands,nrows,ncols,VBitRepn);
  if (! dest) VError("Error creating destination image");
  VCopyImageAttrs (src, dest);
  VFillImage(dest,VAllBands,0);

  for (b=0; b<VImageNBands(src); b++) {
    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {
	i = (long)VGetPixel(src,b,r,c);
	if (table[i] >= msize && i > 0) VPixel(dest,b,r,c,VBit) = 1; 
      }
    }
  }
  return dest;
}
