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

#include "../cylutils/cyl.h"

extern void ROIprint(VImage zmap,VImage metric,VImage roi,int,VString filename);
extern VImage Convert2Repn(VImage src,VImage dest,VRepnKind repn);
  
int main (int argc, char **argv)
{
  static VString  roi_filename="";
  static VString  metric_filename="";
  static VString  txt_filename="";
  static VShort   nbins = 5;
  static VOptionDescRec options[] = {
    {"metric", VStringRepn,1,(VPointer) &metric_filename,VRequiredOpt,NULL,"metric image"},
    {"out", VStringRepn,1,(VPointer) &txt_filename,VRequiredOpt,NULL,"Output txt file"},
    {"roi", VStringRepn,1,(VPointer) &roi_filename,VRequiredOpt,NULL,"ROI image"},
    {"nbins", VShortRepn,1,(VPointer) &nbins,VOptionalOpt,NULL,"Number of bins"},
  };
  VString in_filename=NULL;
  char *prg=GetLipsiaName("vcylarim_roi");
  fprintf (stderr, " %s\n", prg);

  
  /* parse command line */
  VParseFilterCmdX(VNumber(options),options,argc,argv,&in_filename,NULL);

 
  /* read zmap  */
  VImage zmap=NULL;
  VAttrList zlist = VReadAttrList(in_filename,0L,TRUE,FALSE);
  VImage ztmp = VReadImage(zlist);
  if (ztmp == NULL) VError(" no zmap found");
  zmap = Convert2Repn(ztmp,zmap,VFloatRepn);
  VDestroyImage(ztmp);

  
  /* metric image */
  VImage metric=NULL;
  VAttrList xlist = VReadAttrList(metric_filename,0L,TRUE,FALSE);
  if (xlist == NULL) VError(" error reading file %s",metric_filename);
  VImage tmp = VReadImage(xlist);
  if (tmp == NULL) VError(" err reading %s",metric_filename);
  if (VImageNPixels(tmp) != VImageNPixels(zmap))
    VError(" inconsistent image dimensions: metric vs zmap");
  if (VPixelRepn(tmp) != VFloatRepn && VPixelRepn(tmp) != VDoubleRepn)
    VError(" metric image should be 4-byte or 8-byte float");
  metric = Convert2Repn(tmp,metric,VFloatRepn);
  VDestroyImage(tmp);

 
  /* ROI image */
  VAttrList mlist = VReadAttrList(roi_filename,0L,TRUE,FALSE);
  if (mlist == NULL) VError(" error reading file %s",roi_filename);
  VImage mtmp = VReadImage(mlist);
  if (mtmp == NULL) VError(" err reading %s",roi_filename);
  if (VImageNPixels(mtmp) != VImageNPixels(zmap))
    VError(" inconsistent image dimensions: ROI vs zmap");
  VImage roi = Convert2Repn(mtmp,NULL,VBitRepn);
  VDestroyImage(mtmp);

  
  /* main */
  ROIprint(zmap,metric,roi,(int)nbins,txt_filename);


  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
