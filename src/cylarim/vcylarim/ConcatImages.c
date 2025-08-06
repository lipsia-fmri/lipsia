#include <viaio/Vlib.h>
#include <viaio/VImage.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>


VAttrList VConcatImages(VImage *src,VAttrList geolist,int n)
{
  int i,b,r,c;
  double u;
  
  int nslices = VImageNBands(src[0]);
  int nrows   = VImageNRows(src[0]);
  int ncols   = VImageNColumns(src[0]);
  VImage *dest = (VImage *)VCalloc(nslices,sizeof(VImage));
  for (b=0; b<nslices; b++) {
    dest[b] = VCreateImage(n,nrows,ncols,VFloatRepn);
    VCopyImageAttrs (src[0], dest[b]);
    VFillImage(dest[b],VAllBands,0);
  }

  for (i=0; i<n; i++) {
    for (b=0; b<nslices; b++) {
      for (r=0; r<nrows; r++) {
	for (c=0; c<ncols; c++) {
	  u = VGetPixel(src[i],b,r,c);
	  VSetPixel(dest[b],i,r,c,u);
	}
      }
    }
  }

  /* output */
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  for (b=0; b<nslices; b++) {
    VAppendAttr (out_list,"image",NULL,VImageRepn,dest[b]);
  }
  return out_list;
}


/* Function to strip directories and suffixes from a filename */
void strip_filename(char *path)
{
  if (path == NULL || *path == '\0') return;
  char *base = basename(path); 
  strcpy(path, base);
  char *dot = strrchr(path, '.');
  if (dot != NULL && dot != path) {
    *dot = '\0';
  }
}
