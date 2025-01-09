/*
** Input: a segmentation image containing at least three classes:
**      grey matter, white matter, CSF
**
** Output:
**      metric image from two distamce transforms
**
** G.Lohmann, MPI-KYB,  May 2024
*/


/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern VImage VEDist3d(VImage src,VImage dest);
extern VImage VCDist3d(VImage src,VImage dest);
extern void VCDT3d(VImage src,VImage dest,int inside,int outside,int obstacle);

void XSetPixels(VImage src,VImage dest,int id,VBit val)
{
  int b,r,c,i;
  int nslices = VImageNBands(src);
  int nrows  = VImageNRows(src);
  int ncols  = VImageNColumns(src);
  
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	i = (int)VGetPixel(src,b,r,c);
	if (i != id) continue;
	VPixel(dest,b,r,c,VBit) = val;
      }
    }
  }
}


/* normalize to [0,1] using both distance transforms */
void XMetric(VImage dist1,VImage dist2)
{
  size_t i;
  double y;
  VFloat *p1 = VImageData(dist1);
  VFloat *p2 = VImageData(dist2);
  
  for (i=0; i<VImageNPixels(dist1); i++) {
   
    if (p1[i] > 0 && p2[i] > 0) {
      y = p1[i]/p2[i];
      p1[i] = (float)(y/(y+1.0));
    }
    if (p1[i] < 0) p1[i] = 0;
    if (p1[i] > 1) p1[i] = 1;
  }
}

VDictEntry TYPDict[] = {
  { "cdt", 0, 0,0,0,0 },
  { "chamfer", 1, 0,0,0,0  },
  { "euclidean", 2, 0,0,0,0  },
  { NULL, 0,0,0,0,0 }
};

int main (int argc, char **argv)
{
  static VShort gm=2;
  static VShort wm=1;
  static VShort csf=0;
  static VShort type = 0;
  static VOptionDescRec options[] = {
    {"gm", VShortRepn,1,(VPointer) &gm,VOptionalOpt,NULL,"GM"},
    {"wm", VShortRepn,1,(VPointer) &wm,VOptionalOpt,NULL,"WM"},
    {"csf", VShortRepn,1,(VPointer) &csf,VOptionalOpt,NULL,"CSF"},
    {"type", VShortRepn, 1, & type, VOptionalOpt,TYPDict,"Distance metric" },
  };
  VString in_filename=NULL;
  FILE *in_file, *out_file;
  VImage tmp = NULL;
  size_t i;
  char *prg=GetLipsiaName("vmetric");
  fprintf (stderr, "%s\n", prg);

  
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);

  
  /* read image */
  VAttrList list = VReadAttrList(in_filename,0L,TRUE,FALSE);
  VImage src = VReadImage(list);
  if (src == NULL) VError(" no src found");
  VAttrList geolist = VGetGeoInfo(list);

  int nslices = VImageNBands(src);
  int nrows  = VImageNRows(src);
  int ncols  = VImageNColumns(src);

  VImage dist1 = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(dist1,VAllBands,0);
  VCopyImageAttrs (src,dist1);

  VImage dist2 = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(dist2,VAllBands,0);
  VCopyImageAttrs (src,dist2);


  
  /* constrained distance transform */
  if (type == 0) {
    VCDT3d(src,dist1,gm,wm,0);
    VCDT3d(src,dist2,gm,csf,0);
    XMetric(dist1,dist2);
  }

  /* chamfer or euclidean DT */
  else {
    tmp = VCreateImage(nslices,nrows,ncols,VBitRepn);
    VFillImage(tmp,VAllBands,0);
    VCopyImageAttrs (src,tmp);

    /* distance transform 1 */
    fprintf(stderr," dist1...\n");
    XSetPixels(src,tmp,wm,1);
    XSetPixels(src,tmp,csf,0);
    XSetPixels(src,tmp,gm,0);
    VCDist3d(tmp,dist1);
    
    if (type == 1) VEDist3d(tmp,dist1);
    else if (type == 2) VCDist3d(tmp,dist1);
    

    /* distance transform 2 */
    fprintf(stderr," dist2...\n");
    VFillImage(tmp,VAllBands,1);
    XSetPixels(src,tmp,gm,0);
    XSetPixels(src,tmp,wm,0);
    if (type == 1) VEDist3d(tmp,dist2);
    else if (type == 2) VCDist3d(tmp,dist2);
  
    /* extract grey matter */
    VFillImage(tmp,VAllBands,0);
    XSetPixels(src,tmp,gm,1);
    
    VBit *pp = VImageData(tmp);
    VFloat *p1 = VImageData(dist1);
    VFloat *p2 = VImageData(dist2);
    for (i=0; i<VImageNPixels(tmp); i++) {
      if (pp[i] == 0) p1[i] = p2[i] = 0;
    }
    XMetric(dist1,dist2);
  }
  
  
  /* output */
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,dist1);
  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
