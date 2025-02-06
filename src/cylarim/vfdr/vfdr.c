/*
** B-H-FDR applied to zmap
**
** G.Lohmann, Jan 2018
*/
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>


VImage VFDR(VImage src,double alpha,VBoolean dependence) 
{
  size_t i=0,j=0,n=0;
  double z=0;


  /* count number of positive voxels */
  n = 0;
  VFloat *pp = VImageData(src);
  for (i = 0; i < VImageNPixels(src); i++) {
    if (pp[i] > 0) n++;
  }

  /* convert z to p */
  gsl_vector *pvals = gsl_vector_calloc(n);
  j = 0;
  for (i = 0; i < VImageNPixels(src); i++) {
    z = (double)pp[i];
    if (z > 0) {
      pvals->data[j] = gsl_cdf_ugaussian_Q(z);
      j++;
    }
  }
  
  /* sort p-vals */
  gsl_sort_vector(pvals);

  
  /* dependence correction factor, Benjamini–Hochberg–Yekutieli procedure */
  double c=1.0;
  if (dependence) {
    c=0;
    for (i=1; i<=n; i++) c += 1.0/(double)i;
  }
  double nc = c*(double)n;


  /* get FDR threshold */
  double q=0,pthr=2.0;
  for (i=n; i>0; i--) {
    q = alpha * ((double)i)/nc;
    if (pvals->data[i-1] <= q) {
      pthr = pvals->data[i-1];
      break;
    }
  }
  if (pthr > 1.0) VError(" no voxel above FDR threshold");

  
  /* convert p to z */
  double eps = DBL_EPSILON;
  double p1 = 1.0 - eps;
  if (pthr < eps) pthr = eps;
  if (pthr > p1) pthr = p1;
  double zthr = gsl_cdf_ugaussian_Qinv(pthr);
  fprintf(stderr," FDR  z-threshold= %.4f\n",zthr);

  
  /* fill dest image */
  VImage dest = VCreateImageLike(src);
  VFillImage(dest,VAllBands,0);
  
  VFloat *pq = VImageData(dest);
  j=0;
  for (i = 0; i < VImageNPixels(src); i++) {
    z = (double)(pp[i]);
    if (z > zthr) { pq[i] = z; j++; }
  }
  fprintf(stderr," voxels above FDR threshold: %lu of %lu,  ratio: %.4f\n",j,n,(double)j/(double)n);
  return dest;
}




int main (int argc, char *argv[])
{
  static VFloat  alpha = 0.05; 
  static VBoolean dependence=FALSE;
  static VOptionDescRec options[] = {
    {"alpha",VFloatRepn,1,(VPointer) &alpha,VOptionalOpt,NULL,"FDR significance level"},
    {"dependence",VBooleanRepn,1,(VPointer) &dependence,VOptionalOpt,NULL,"Dependence correction"},
  };
  FILE *in_file,*out_file;
  char *prg=GetLipsiaName("vfdr");
  fprintf (stderr, "%s\n", prg);
  gsl_set_error_handler_off();
  
  
  /* Parse command line arguments: */
  VParseFilterCmd (VNumber (options), options,argc,argv,&in_file,&out_file);

  
  /* Read the input file: */
  VAttrList list = VReadFile (in_file, NULL);
  if (! list) exit (1);
  VAttrList geolist = VGetGeoInfo(list);
  VImage src = VReadImage(list);
  if (VPixelRepn(src) != VFloatRepn) VError(" input zmap must be 4-byte float");


  /* main */
  VImage dest = VFDR(src,(double)alpha,dependence);

  
  /* output */
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,dest);
  VWriteFile (out_file, out_list);
  fprintf (stderr, "%s: done.\n",argv[0]);
  return 0;
}
