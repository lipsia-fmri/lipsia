/*
** filter: most freqent value in a local nieghbourhood
**
** G.Lohmann, Jan 2016
*/


/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>
#include <via/via.h>

/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))
extern VImage VBorderImage3d (VImage src,VImage dest);


#define GM_CSF 1
#define GM_WM 2
#define GM_INSIDE 3


VDictEntry ADJDict[] = {
  { "6", 0 },
  { "18", 1 },
  { "26", 2 },
  { NULL }
};


VImage XBorder(VImage src,VImage dest)
{
  int b,r,c,bb,rr,cc;
  int nslices = VImageNBands(src);
  int nrows = VImageNRows(src);
  int ncols = VImageNColumns(src);

  
  if (dest == NULL) dest = VCreateImageLike(src);
  VFillImage(dest,VAllBands,0);
  
  for (b=1; b<nslices-1; b++) {
    for (r=1; r<nrows-1; r++) {
      for (c=1; c<ncols-1; c++) {		       
	if (VPixel(src,b,r,c,VBit) == 0) continue;

	for (bb=b-1; bb<=b+1; bb++) {
	  for (rr=r-1; rr<=r+1; rr++) {
	    for (cc=c-1; cc<=c+1; cc++) {
	      if (VPixel(src,bb,rr,cc,VBit) == 0) VPixel(dest,b,r,c,VBit) = 1;
	    }
	  }
	}
      }
    }
  }
  return dest;
}

VImage XBinarize(VImage src,int id)
{
  int b,r,c,i;
  int nslices = VImageNBands(src);
  int nrows  = VImageNRows(src);
  int ncols  = VImageNColumns(src);

  VImage dest = VCreateImage(nslices,nrows,ncols,VBitRepn);
  VFillImage(dest,VAllBands,0);
  VCopyImageAttrs (src, dest);
  
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	i = (int)VGetPixel(src,b,r,c);
	if (i != id) continue;
	VPixel(dest,b,r,c,VBit) = 1;
      }
    }
  }
  return dest;
}

VImage VRim(VImage src,VImage border,int gm,int wm,int csf)
{ 
  int i,j,m,k,l;
  int b,r,c;
  int adjdef=1;
  int nslices = VImageNBands(src);
  int nrows  = VImageNRows(src);
  int ncols  = VImageNColumns(src);
  VRepnKind repn = VPixelRepn(src);
  if (repn != VUByteRepn && repn != VShortRepn && repn != VIntegerRepn)
    VError("input must be ubyte/short/int");


  fprintf(stderr," gm: %d,  wm: %d,  csf: %d\n",gm,wm,csf);
  
  VImage dest = VCreateImage(nslices,nrows,ncols,VUByteRepn);
  VFillImage(dest,VAllBands,0);
  VCopyImageAttrs (src, dest);

  for (b=1; b<nslices-1; b++) {
    for (r=1; r<nrows-1; r++) {
      for (c=1; c<ncols-1; c++) {

	i = VGetPixel(src,b,r,c);
	if (i==gm) VPixel(dest,b,r,c,VUByte) = GM_INSIDE;

	if (VPixel(border,b,r,c,VBit) == 0) continue;
		
	for (m=-1; m<=1; m++) {
	  for (k=-1; k<=1; k++) {
	    for (l=-1; l<=1; l++) {

	      if (adjdef == 0) {     /* 6 adjacency */
		if (ABS(m)+ABS(k)+ABS(l) > 1) continue;
	      }
	      j = VGetPixel(src,b+m,r+k,c+l);
	      if (j == csf) VPixel(dest,b,r,c,VUByte) = GM_CSF;
	      if (j == wm) VPixel(dest,b,r,c,VUByte) = GM_WM;
	    }
	  }
	}
      }
    }
  }
  return dest;
}


int main (int argc, char **argv)
{
  static VShort gm=2;
  static VShort wm=1;
  static VShort csf=0;
  static VOptionDescRec options[] = {
    {"gm", VShortRepn,1,(VPointer) &gm,VOptionalOpt,NULL,"GM"},
    {"wm", VShortRepn,1,(VPointer) &wm,VOptionalOpt,NULL,"WM"},
    {"csf", VShortRepn,1,(VPointer) &csf,VOptionalOpt,NULL,"CSF"},
  };
  VString in_filename=NULL;
  FILE *in_file, *out_file;
  char *prg=GetLipsiaName("vrim");
  fprintf (stderr, "%s\n", prg);

  
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);

  
  /* read image */
  VAttrList list = VReadAttrList(in_filename,0L,TRUE,FALSE);
  VImage src = VReadImage(list);
  if (src == NULL) VError(" no src found");
  VAttrList geolist = VGetGeoInfo(list);


  VImage tmp = XBinarize(src,(int)gm);
  /* VImage border = VBorderImage3d (tmp,NULL); */
  VImage border = XBorder(tmp,NULL);
  VImage dest = VRim(src,border,(int)gm,(int)wm,(int)csf);


  /* output */
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,dest);
  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
