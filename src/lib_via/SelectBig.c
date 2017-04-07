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



VImage VSelectBig (VImage src,VImage dest)
{
  long i,j,i0;
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
    VError(" VSelectBig: illegal pixel repn");
  }


  /* get largest component */
  long maxsize = 0;
  i0 = -1;
  for (i=1; i<maxlabel; i++) {
    if (table[i] > maxsize) {
      i0 = i;
      maxsize = table[i];
    }
  }
  if (i0 < 0)
    VError(" input image is zero.");
  VFree(table);


  /* create output image */
  dest = VSelectDestImage("VSelectBig",dest,nbands,nrows,ncols,VBitRepn);
  if (! dest) VError("Error creating destination image");
  VCopyImageAttrs (src, dest);
  VFillImage(dest,VAllBands,0);

  for (b=0; b<VImageNBands(src); b++) {
    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {
	i = (long)VGetPixel(src,b,r,c);
	if (i==i0) VPixel(dest,b,r,c,VBit) = 1; 
      }
    }
  }
  return dest;
}
