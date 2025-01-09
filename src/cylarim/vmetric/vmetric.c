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


extern VImage Convert2Repn(VImage src,VImage dest,VRepnKind repn);
extern void VCDT3d(VImage src,VImage dest,int inside,int outside,int obstacle);


/* normalize to [0,1] using both distance transforms */
void XMetric(VImage dist1,VImage dist2)
{
  size_t i;
  VFloat *p1 = VImageData(dist1);
  VFloat *p2 = VImageData(dist2);
  
  for (i=0; i<VImageNPixels(dist1); i++) {
   
    if (p1[i] > 0 || p2[i] > 0) {
      p1[i] = p1[i]/(p1[i] + p2[i]);
    }
    if (p1[i] < 0) p1[i] = 0;
    if (p1[i] > 1) p1[i] = 1;
  }
}


int main (int argc, char **argv)
{
  static VShort gm=3;
  static VShort wm=2;
  static VShort csf=1;
  static VOptionDescRec options[] = {
    {"gm", VShortRepn,1,(VPointer) &gm,VOptionalOpt,NULL,"GM"},
    {"wm", VShortRepn,1,(VPointer) &wm,VOptionalOpt,NULL,"WM"},
    {"csf", VShortRepn,1,(VPointer) &csf,VOptionalOpt,NULL,"CSF"},
  };
  VString in_filename=NULL;
  FILE *in_file, *out_file;
  char *prg=GetLipsiaName("vmetric");
  fprintf (stderr, "%s\n", prg);

  
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);

  
  /* read image */
  VImage rim=NULL;
  VAttrList list = VReadAttrList(in_filename,0L,TRUE,FALSE);
  VImage tmp = VReadImage(list);
  if (tmp == NULL) VError(" no src found");
  VAttrList geolist = VGetGeoInfo(list);
  rim = Convert2Repn(tmp,rim,VUByteRepn);
  VDestroyImage(tmp);

  
  int nslices = VImageNBands(rim);
  int nrows  = VImageNRows(rim);
  int ncols  = VImageNColumns(rim);

  VImage dist1 = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(dist1,VAllBands,0);
  VCopyImageAttrs (rim,dist1);

  VImage dist2 = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(dist2,VAllBands,0);
  VCopyImageAttrs (rim,dist2);

  
  /* constrained distance transform */
  VCDT3d(rim,dist1,gm,wm,0);
  VCDT3d(rim,dist2,gm,csf,0);
  XMetric(dist1,dist2);


  /* add rim pts to metric image */
  VFloat *px = VImageData(dist1);
  VUByte *pu = VImageData(rim);
  size_t i;
  for (i=0; i<VImageNPixels(rim); i++) {
    if (pu[i] == 1) px[i] = 1;
    if (pu[i] == 2) px[i] = 0.0001;
  }

  
  /* output */
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,dist1);
  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
