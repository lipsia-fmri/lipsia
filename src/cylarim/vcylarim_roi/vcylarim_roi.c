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


extern VImage Convert2Repn(VImage src,VImage dest,VRepnKind repn);
extern void ROIprint(VImage zmap,VImage metric,VImage *,int,VImage roi);
extern void ROIstats(VImage *betaimage,VImage *covimage,VImage,int,int,VImage roi);
extern int BetaStats(VImage metric,VImage *betaimage,VImage roi);
extern void ROIbeta(VImage *betaimage,VImage metric,VImage roi,VString filename);
  
int main (int argc, char **argv)
{
  static VString  zmap_filename="";
  static VString  metric_filename="";
  static VString  rim_filename="";
  static VString  roi_filename="";
  static VString  txt_filename="";
  static VBoolean xprint=FALSE;
  static VOptionDescRec options[] = {
    {"zmap", VStringRepn,1,(VPointer) &zmap_filename,VRequiredOpt,NULL,"zmap"},
    {"metric", VStringRepn,1,(VPointer) &metric_filename,VRequiredOpt,NULL,"metric image"},
    {"rim", VStringRepn,1,(VPointer) &rim_filename,VRequiredOpt,NULL,"rim image"},
    {"roi", VStringRepn,1,(VPointer) &roi_filename,VRequiredOpt,NULL,"ROI image"},
    {"txt", VStringRepn,1,(VPointer) &txt_filename,VRequiredOpt,NULL,"Output txt file"},
    {"print",VBooleanRepn,1,(VPointer) &xprint,VOptionalOpt,NULL,"Print full"},
  };
  VString in_file=NULL;
  size_t n=0;
  char *prg=GetLipsiaName("vcylarim_roi");
  fprintf (stderr, " %s\n", prg);

  
  /* parse command line */
  VParseFilterCmdX(VNumber(options),options,argc,argv,&in_file,NULL);


  /* beta images */
  VLong tr=0L;
  VAttrList list = VReadAttrList(in_file,tr,TRUE,FALSE);
  if (list == NULL) VError(" error reading input file %s",in_file);
  VImage src;
  VAttrListPosn posn;
  int nbeta=0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (strcmp(VGetAttrName (& posn),"nbeta") == 0) {
      VGetAttrValue (& posn, NULL,VShortRepn, & nbeta);
    }
  }
  fprintf(stderr," nbeta: %d\n",nbeta);
  VImage *betaimage = (VImage *)VCalloc(nbeta,sizeof(VImage));
  
  n=0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    if (strcmp(VGetAttrName (& posn),"beta") == 0) {
      VGetAttrValue (& posn, NULL,VImageRepn, & src);
      if (n >= nbeta) VError(" too many beta images ");
      betaimage[n] = VCopyImage(src,NULL,VAllBands);
      n++;
    }
  }
  if (n < 1) VError(" no beta images found");


  /* covariance */
  int ncov=0;
  VImage *covimage = NULL;
  if (nbeta == 4) {
    ncov=6;
    covimage = (VImage *)VCalloc(ncov,sizeof(VImage));
    n=0;
    for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      if (strcmp(VGetAttrName (& posn),"covariance") == 0) {
	VGetAttrValue (& posn, NULL,VImageRepn, & src);
	if (n >= ncov) VError(" too many cov images ");
	covimage[n] = VCopyImage(src,NULL,VAllBands);
	n++;
      }
    }
  }


  
  /* read zmap  */
  VImage zmap=NULL;
  VAttrList zlist = VReadAttrList(zmap_filename,0L,TRUE,FALSE);
  VImage ztmp = VReadImage(zlist);
  if (ztmp == NULL) VError(" no zmap found");
  if (VPixelRepn(ztmp) != VFloatRepn && VPixelRepn(ztmp) != VDoubleRepn)
    VError(" zmap image should be 4-byte or 8-byte float");
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

  
  /* rim image */
  VImage rim=NULL;
  VAttrList rlist = VReadAttrList(rim_filename,0L,TRUE,FALSE);
  if (rlist == NULL) VError(" error reading file %s",rim_filename);
  VImage rtmp = VReadImage(rlist);
  if (rtmp == NULL) VError(" err reading %s",rim_filename);
  if (VImageNPixels(rtmp) != VImageNPixels(zmap))
    VError(" inconsistent image dimensions: rim vs zmap");
  rim = Convert2Repn(rtmp,rim,VUByteRepn);
  VDestroyImage(rtmp);

  
  /* add rim pts to metric image in case they are not included */
  VUByte *pu = VImageData(rim);
  VFloat *px = VImageData(metric);
  size_t i;
  for (i=0; i<VImageNPixels(rim); i++) {
    if (pu[i] == 1) px[i] = 1;
    if (pu[i] == 2) px[i] = 0.0001;
  }

  
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
  ROIbeta(betaimage,metric,roi,txt_filename);

  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
