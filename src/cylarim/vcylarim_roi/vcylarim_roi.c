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
extern void ROIplot(VImage *betaimage,VImage roi,VString filename);
  
int main (int argc, char **argv)
{
  static VString  roi_filename="";
  static VString  txt_filename="";
  static VOptionDescRec options[] = {
    {"roi", VStringRepn,1,(VPointer) &roi_filename,VRequiredOpt,NULL,"ROI image"},
    {"txt", VStringRepn,1,(VPointer) &txt_filename,VRequiredOpt,NULL,"Output txt file"},
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

 
  /* ROI image */
  VAttrList mlist = VReadAttrList(roi_filename,0L,TRUE,FALSE);
  if (mlist == NULL) VError(" error reading file %s",roi_filename);
  VImage mtmp = VReadImage(mlist);
  if (mtmp == NULL) VError(" err reading %s",roi_filename);
  if (VImageNPixels(mtmp) != VImageNPixels(betaimage[0]))
    VError(" inconsistent image dimensions: ROI vs betaimage");
  VImage roi = Convert2Repn(mtmp,NULL,VBitRepn);
  VDestroyImage(mtmp);

  
  /* main */
  ROIplot(betaimage,roi,txt_filename);

  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
