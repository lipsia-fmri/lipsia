/*
** Implementation of LISA algorithm
** for statistical inference of fMRI images
**
** generic plug-in.
** A file containing a list of all 3D permutation images must be supplied.
**
** G.Lohmann, April 2017
*/

#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>
#include <viaio/option.h>
#include <viaio/os.h>
#include <viaio/VImage.h>
#include <via/via.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

#define ABS(x) ((x) > 0 ? (x) : -(x))

extern void VIsolatedVoxels(VImage src,float threshold);
extern void VHistogram(gsl_histogram *histogram,VString filename);
extern void VCheckImage(VImage src);
extern void FDR(VImage src,VImage dest,gsl_histogram *nullhist,gsl_histogram *realhist,double alpha);
extern void ImageStats(VImage src,double *,double *,double *hmin,double *hmax);
extern void VBilateralFilter(VImage src,VImage dest,int radius,double var1,double var2,int);
extern double VImageVar(VImage src);
extern void  VImageCount(VImage src);
extern void  VGetHistRange(VImage src,double *hmin,double *hmax);
extern float VGetMode(VImage src);
extern void  VZScale(VImage src,float mode,float stddev);
extern void HistoUpdate(VImage,gsl_histogram *);




/* make sure all input images are in float and have the same number of pixels */
void CheckImageTypes(VImage zmap,VImage *permimages,int numperm)
{
  size_t npixels = VImageNPixels(zmap);
  int i;
  for (i=0; i<numperm; i++) {
    if (npixels != VImageNPixels(permimages[i]))
      VError(" inconsistent number of pixels in permutation image %d",i);
    if (VPixelRepn(permimages[i]) != VFloatRepn) VError(" perm image %d is not in float repn");
  }
}


int main (int argc, char *argv[])
{
  static VString  filename = "";
  static VFloat   alpha = 0.05;
  static VShort   radius = 2;
  static VFloat   rvar = 2.0;
  static VFloat   svar = 2.0;
  static VShort   numiter = 2;
  static VBoolean centering = FALSE;
  static VBoolean cleanup = TRUE;
  static VShort   nproc = 0;
  static VOptionDescRec  options[] = {
    {"permutations",VStringRepn,1,(VPointer) &filename,VRequiredOpt,NULL,"List of all permutation images"},
    {"alpha",VFloatRepn,1,(VPointer) &alpha,VOptionalOpt,NULL,"FDR significance level"},
    {"radius",VShortRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Bilateral parameter (radius in voxels)"},
    {"rvar",VFloatRepn,1,(VPointer) &rvar,VOptionalOpt,NULL,"Bilateral parameter (radiometric)"},
    {"svar",VFloatRepn,1,(VPointer) &svar,VOptionalOpt,NULL,"Bilateral parameter (spatial)"},
    {"filteriterations",VShortRepn,1,(VPointer) &numiter,VOptionalOpt,NULL,"Bilateral parameter (number of iterations)"},
    {"cleanup",VBooleanRepn,1,(VPointer) &cleanup,VOptionalOpt,NULL,"Whether to remove isolated voxels"},
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"},
  };
  FILE *out_file=NULL,*in_file=NULL;
  VString in_filename=NULL;
  char *prg_name=GetLipsiaName("vlisa_generic");
  fprintf (stderr, "%s\n", prg_name);
  gsl_set_error_handler_off ();


  /* Parse command line arguments and identify files: */
  VParseFilterCmdZ (VNumber (options), options, argc, argv,&in_file,&out_file,&in_filename);


  /* Read the input zmap file: */
  VAttrList list = VReadAttrListZ(in_file,in_filename,0L,TRUE,FALSE);
  VImage zmap1 = VReadImage(list);
  if (zmap1 == NULL) VError(" no input zmap image found");
  if (VPixelRepn(zmap1) != VFloatRepn) VError(" input pixel repn must be float");
  VAttrList geolist = VGetGeoInfo(list);


  /* Read permutation file containing a list of 3D images */
  VAttrList listperm = VReadAttrList(filename,0L,FALSE,FALSE);
  int numperm  = VAttrListNumImages(listperm);
  VImage *zmap = VAttrListGetImages(listperm,numperm);
  CheckImageTypes(zmap1,zmap,numperm);
  fprintf(stderr," Number of permutation images: %d\n",(int)numperm);


  /* estimate null variance to adjust radiometric parameter, use first 30 permutations */
  double zvar=0,hmin=0,hmax=0;
  float stddev=1.0;
  int nperm=0;
  if (numperm > 0) {
    int tstperm = 30;
    if (tstperm > numperm) tstperm = numperm;
    double varsum=0,nx=0;
    for (nperm = 0; nperm < tstperm; nperm++) {
      zvar = VImageVar(zmap[nperm]);
      varsum += zvar;
      nx++;
    }
    double meanvar = varsum/nx;
    stddev = sqrt(meanvar);
  }


  /* get non-permuted hotspot map */
  float mode=0;
  if (centering) mode = VGetMode(zmap1);
  if (numperm > 0) VZScale(zmap1,mode,stddev);

  VImage dst1 = VCreateImageLike(zmap1);
  VBilateralFilter(zmap1,dst1,(int)radius,(double)rvar,(double)svar,(int)numiter);



  /* ini histograms */
  VGetHistRange(dst1,&hmin,&hmax);
  size_t nbins = 20000;
  gsl_histogram *hist0 = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (hist0,hmin,hmax);
  gsl_histogram *histz = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (histz,hmin,hmax);
  HistoUpdate(dst1,histz);


  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  if (nproc > 0 && nproc < num_procs) num_procs = nproc;
  fprintf(stderr," using %d cores\n",(int)num_procs);
  omp_set_num_threads(num_procs);
#endif /* _OPENMP */

  /* do random permutations */
#pragma omp parallel for shared(zmap) schedule(dynamic)
  for (nperm = 0; nperm < numperm; nperm++) {
    if (nperm%20 == 0) fprintf(stderr," perm  %4d  of  %d\r",nperm,(int)numperm);

    float mode=0;
    if (centering) mode = VGetMode(zmap[nperm]);
    VZScale(zmap[nperm],mode,stddev);

    VImage dst = VCreateImageLike (zmap1);
    VBilateralFilter(zmap[nperm],dst,(int)radius,(double)rvar,(double)svar,(int)numiter);

    #pragma omp critical
    {
      HistoUpdate(dst,hist0);
    }
    VDestroyImage(dst);
  }


  /* apply fdr */
  VImage fdrimage = VCopyImage (dst1,NULL,VAllBands);
  if (numperm > 0) {
    FDR(dst1,fdrimage,hist0,histz,(double)alpha);

    if (cleanup && alpha < 1.0) {
      VIsolatedVoxels(fdrimage,(float)(1.0-alpha));
    }
  }


  /* output */
  VAttrList out_list = VCreateAttrList ();
  VHistory(VNumber(options),options,prg_name,&list,&out_list);
  VSetGeoInfo(geolist,out_list);
  VAppendAttr (out_list,"image",NULL,VImageRepn,fdrimage);
  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "\n%s: done.\n", argv[0]);
  exit(0);
}
