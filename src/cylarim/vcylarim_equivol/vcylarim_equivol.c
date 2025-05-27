/*
** Cylarim:  Cylinder-based Algorithm for Laminar fMRI Data Analysis
**
** G.Lohmann, MPI-KYB, Oct 2024
*/


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>

#include "../cylutils/cyl.h"

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

extern VImage Convert2Repn(VImage src,VImage dest,VRepnKind repn);
extern Cylinders *VCylinder(VImage rim,VImage metric,double radius,VBoolean);
extern void HistEqualize(Cylinders *,VImage,VImage,VImage,double);


void XWriteOutput(VImage image,VAttrList geolist,char *filename)
{
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,image);
  FILE *fp = fopen(filename,"w");
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
}


VImage Cylarim_equivol(VImage metric,VImage rim,double radius,double alpha)
{

  /* create cylinder data struct */
  Cylinders *cyl = VCylinder(rim,metric,radius,TRUE);

  
  /* ini output image */
  VImage wimage = VCreateImageLike(metric);
  VFillImage(wimage,VAllBands,0);

  
  /* equivolume correction */
  HistEqualize(cyl,wimage,metric,rim,alpha);

  return wimage;
}


int main (int argc, char **argv)
{
  static VString    rim_filename="";
  static VString    mask_filename="";
  static VBoolean   equivol = TRUE;
  static VFloat     radius = 2.0;
  static VFloat     alpha = -1;
  static VShort     nproc = 0;
  static VOptionDescRec options[] = {
    {"rim", VStringRepn,1,(VPointer) &rim_filename,VRequiredOpt,NULL,"rim image"},
    {"mask", VStringRepn,1,(VPointer) &mask_filename,VOptionalOpt,NULL,"mask image"},
    {"radius", VFloatRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Cylinder radius in mm"},
    {"alpha", VFloatRepn,1,(VPointer) &alpha,VOptionalOpt,NULL,"Alpha"},
    {"equivol", VBooleanRepn,1,(VPointer) &equivol,VOptionalOpt,NULL,"Whether to output equivol metric"},
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"},
   };
  VString in_filename=NULL;
  FILE *in_file, *out_file;
  size_t i;
  char *prg=GetLipsiaName("vcylarim_equivol");
  fprintf (stderr, " %s\n", prg);

  
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);

  
  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  if (nproc > 0 && nproc < num_procs) num_procs = nproc;
  if (num_procs > 64) num_procs=64;
  fprintf(stderr," using %d cores\n",(int)num_procs);
  omp_set_num_threads(num_procs);
#endif /* _OPENMP */

 
  
  /* metric image */
  VImage metric=NULL;  
  VAttrList xlist = VReadAttrList(in_filename,0L,TRUE,FALSE);
  if (xlist == NULL) VError(" error reading file %s",in_filename);
  VImage tmp = VReadImage(xlist);
  if (tmp == NULL) VError(" err reading %s",in_filename);
  if (VPixelRepn(tmp) != VFloatRepn && VPixelRepn(tmp) != VDoubleRepn)
    VError(" metric image should be 4-byte or 8-byte float");
  metric = Convert2Repn(tmp,metric,VFloatRepn);
  VDestroyImage(tmp);
  VAttrList geolist = VGetGeoInfo(xlist);


  
  /* rim image */
  VImage rim=NULL;
  VAttrList rlist = VReadAttrList(rim_filename,0L,TRUE,FALSE);
  if (rlist == NULL) VError(" error reading file %s",rim_filename);
  VImage rtmp = VReadImage(rlist);
  if (rtmp == NULL) VError(" err reading %s",rim_filename);
  if (VImageNPixels(rtmp) != VImageNPixels(metric))
    VError(" inconsistent image dimensions: rim vs metric");
  rim = Convert2Repn(rtmp,rim,VUByteRepn);
  VDestroyImage(rtmp);

  
  
  /* mask image */
  VFloat *px = VImageData(metric);
  VUByte *pu = VImageData(rim);
  VImage mask=NULL;
  if (strlen(mask_filename) > 1) {
    VAttrList mlist = VReadAttrList(mask_filename,0L,TRUE,FALSE);
    if (mlist == NULL) VError(" error reading file %s",mask_filename);
    VImage mtmp = VReadImage(mlist);
    if (mtmp == NULL) VError(" err reading %s",mask_filename);
    if (VImageNPixels(mtmp) != VImageNPixels(metric))
      VError(" inconsistent image dimensions: mask vs metric");
    mask = Convert2Repn(mtmp,mask,VBitRepn);
    VDestroyImage(mtmp);
  
    VBit *pm = VImageData(mask);  /* apply mask, if present */
    for (i=0; i<VImageNPixels(rim); i++) {
      if (pm[i] == 0 || pu[i] == 0) { px[i] = 0; pu[i] = 0; }
    }
  }
  

  /* add rim pts to metric image in case they are not included */
  for (i=0; i<VImageNPixels(rim); i++) {
    if (pu[i] == 1) px[i] = 1;
    if (pu[i] == 2) px[i] = 0.0001;
  }

  
  /* Main */
  VImage wimage = Cylarim_equivol(metric,rim,(double)radius,(double)alpha);


  
  /* write output */
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  char cbuf[64];
  sprintf(cbuf, "%f",radius);
  VAppendAttr(out_list,"radius",NULL,VStringRepn,(VString)cbuf);
  if (equivol) VAppendAttr(out_list,"equivol",NULL,VImageRepn,metric);
  else VAppendAttr(out_list,"wimage",NULL,VImageRepn,wimage);

  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
