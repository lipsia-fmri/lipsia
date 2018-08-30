/*
** Implementation of LISA algorithm
** for statistical inference of fMRI images
**
** 2nd level inference (one-sample test)
**
** G.Lohmann, 2017
*/
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>
#include <viaio/option.h>
#include <viaio/os.h>
#include <viaio/VImage.h>

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

extern void   VIsolatedVoxels(VImage src,float threshold);
extern void   VHistogram(gsl_histogram *histogram,VString filename);
extern void   VCheckImage(VImage src);
extern void   FDR(VImage src,VImage dest,gsl_histogram *nullhist,gsl_histogram *realhist,double);
extern double ttest1(double *data1,int n);
extern void   ImageStats(VImage src,double *,double *,double *hmin,double *hmax);
extern void   VBilateralFilter(VImage src,VImage dest,int radius,double var1,double var2,int);
extern double VImageVar(VImage src);
extern void   VImageCount(VImage src);
extern void   VGetHistRange(VImage src,double *hmin,double *hmax);
extern void   VZScale(VImage src,float mode,float stddev);
extern float  VGetMode(VImage src);
extern double t2z(double,double);
extern void   HistoUpdate(VImage,gsl_histogram *);

/* generate permutation table */
int **genperm(gsl_rng *rx,int n,int numperm)
{
  int i,j;
  int **table = (int **) VCalloc(numperm,sizeof(int *));
  for (i = 0; i < numperm; i++) {
    table[i] = (int *) VCalloc(n,sizeof(int));
    for (j=0; j<n; j++) {
      table[i][j] = 0;
      if (gsl_ran_bernoulli (rx,(double)0.5) == 1) table[i][j] = 1;
    }
  }
  return table;
}



/* onesample t-test */
void OnesampleTest(VImage *src1,VImage dest,int *permtable,int n)
{
  int i,k,b,r,c,nslices,nrows,ncols;
  double nx,ave,var,u,t,z,df;
  double tiny=1.0e-8;
  extern void avevar(double *data,int n,double *a,double *v);

  nslices = VImageNBands(src1[0]);
  nrows   = VImageNRows(src1[0]);
  ncols   = VImageNColumns(src1[0]);
  VFillImage(dest,VAllBands,0);
  double *data = (double *) VCalloc(n,sizeof(double));

  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	k = 0;
	for (i=0; i<n; i++) {
	  u = (double)VPixel(src1[i],b,r,c,VFloat);
	  if (ABS(u) > tiny) {
	    if (permtable[i] > 0) u = -u;
	    data[k] = u;
	    k++;
	  }
	}
	if (k < n-2) continue;
	avevar(data,k,&ave,&var);
	if (var < tiny) continue;
	nx   = (double)k;
	t    = sqrt(nx) * ave/sqrt(var);
	df   = nx - 1.0;
	z    = t2z(t,df);
	if (t < 0) z = -z;
	VPixel(dest,b,r,c,VFloat) = z;
      }
    }
  }
  VFree(data);
}



int main (int argc, char *argv[])
{
  static VArgVector in_files1;
  static VString  out_filename="";
  static VFloat   alpha = 0.05;
  static VShort   radius = 2;
  static VFloat   rvar = 2.0;
  static VFloat   svar = 2.0;
  static VShort   numiter = 2;
  static VShort   numperm = 5000;
  static VLong    seed = 99402622;
  static VBoolean centering = FALSE;
  static VBoolean cleanup = TRUE;
  static VShort   nproc = 0;
  static VOptionDescRec options[] = {
    {"in", VStringRepn, 0, & in_files1, VRequiredOpt, NULL,"Input files" },
    {"out", VStringRepn, 1, & out_filename, VRequiredOpt, NULL,"Output file" },
    {"alpha",VFloatRepn,1,(VPointer) &alpha,VOptionalOpt,NULL,"FDR significance level"},
    {"perm",VShortRepn,1,(VPointer) &numperm,VOptionalOpt,NULL,"Number of permutations"},
    {"seed",VLongRepn,1,(VPointer) &seed,VOptionalOpt,NULL,"Seed for random number generation"},
    {"radius",VShortRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Bilateral parameter (radius in voxels)"},
    {"rvar",VFloatRepn,1,(VPointer) &rvar,VOptionalOpt,NULL,"Bilateral parameter (radiometric)"},
    {"svar",VFloatRepn,1,(VPointer) &svar,VOptionalOpt,NULL,"Bilateral parameter (spatial)"},
    {"filteriterations",VShortRepn,1,(VPointer) &numiter,VOptionalOpt,NULL,"Bilateral parameter (number of iterations)"},
    {"cleanup",VBooleanRepn,1,(VPointer) &cleanup,VOptionalOpt,NULL,"Whether to remove isloated voxels"},
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"},
  };

  FILE *fp=NULL;
  VString in_filename;
  VAttrList list1=NULL,out_list=NULL,geolist=NULL;
  VImage *src1;
  int i,nimages,npix=0;
  char *prg_name=GetLipsiaName("vlisa_onesample");
  fprintf (stderr, "%s\n", prg_name);


  /* parse command line */
  if (! VParseCommand (VNumber (options), options, & argc, argv)) {
    VReportUsage (argv[0], VNumber (options), options, NULL);
    exit (EXIT_FAILURE);
  }
  if (argc > 1) {
    VReportBadArgs (argc, argv);
    exit (EXIT_FAILURE);
  }


  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  if (nproc > 0 && nproc < num_procs) num_procs = nproc;
  fprintf(stderr," using %d cores\n",(int)num_procs);
  omp_set_num_threads(num_procs);
#endif /* _OPENMP */


  /* images  */
  nimages = in_files1.number;
  src1 = (VImage *) VCalloc(nimages,sizeof(VImage));
  for (i = 0; i < nimages; i++) {
    in_filename = ((VString *) in_files1.vector)[i];
    list1   = VReadAttrList(in_filename,0L,TRUE,FALSE);
    src1[i] = VReadImage(list1);
    if (src1[i] == NULL) VError(" no input image found");
    if (VPixelRepn(src1[i]) != VFloatRepn) VError(" input pixel repn must be float");
    if (i == 0) npix = VImageNPixels(src1[i]);
    else if (npix != VImageNPixels(src1[i])) VError(" inconsistent image dimensions");

    /* use geometry info from 1st file */
    if (geolist == NULL) geolist = VGetGeoInfo(list1);
  }



  /* ini random permutations */
  size_t n = nimages;
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);
  int nperm=0;
  int **permtable = genperm(rx,(int)n,(int)numperm);


  /* estimate null variance based on first 30 permutations */
  double hmin=0,hmax=0;
  float stddev=1.0;
  if (numperm > 0) {
    int tstperm = 30;
    if (tstperm > numperm) tstperm = numperm;
    double varsum=0,nx=0;

#pragma omp parallel for shared(src1,permtable) schedule(dynamic)
    for (nperm = 0; nperm < tstperm; nperm++) {
      VImage zmap = VCreateImageLike(src1[0]);
      OnesampleTest(src1,zmap,permtable[nperm],nimages);
#pragma omp critical
      {
	varsum += VImageVar(zmap);
	nx++;
      }
      VDestroyImage(zmap);
    }
    double meanvar = varsum/nx;
    stddev = sqrt(meanvar);
  }



  /* no permutation */
  int *nopermtable = (int *) VCalloc(n,sizeof(int));
  VImage dst1  = VCreateImageLike (src1[0]);
  VImage zmap1 = VCreateImageLike(src1[0]);
  OnesampleTest(src1,zmap1,nopermtable,nimages);

  if (numperm == 0) {
    double z = VImageVar(zmap1);
    stddev = (float)(sqrt(z)); /* update stddev */
  }
  float mode=0;
  if (centering) mode = VGetMode(zmap1);
  if (numperm > 0) VZScale(zmap1,mode,stddev);
  VBilateralFilter(zmap1,dst1,(int)radius,(double)(rvar),(double)svar,(int)numiter);


  /* ini histograms */
  VGetHistRange(dst1,&hmin,&hmax);
  size_t nbins = 20000;
  gsl_histogram *hist0 = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (hist0,hmin,hmax);
  gsl_histogram *histz = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (histz,hmin,hmax);
  HistoUpdate(dst1,histz);



  /* random permutations */
#pragma omp parallel for shared(src1,permtable) schedule(dynamic)
  for (nperm = 0; nperm < numperm; nperm++) {
    if (nperm%20 == 0) fprintf(stderr," perm  %4d  of  %d\r",nperm,(int)numperm);

    VImage zmap = VCreateImageLike(src1[0]);
    VImage dst  = VCreateImageLike (zmap);
    OnesampleTest(src1,zmap,permtable[nperm],nimages);
    float mode=0;
    if (centering) mode = VGetMode(zmap);
    VZScale(zmap,mode,stddev);
    VBilateralFilter(zmap,dst,(int)radius,(double)(rvar),(double)svar,(int)numiter);

#pragma omp critical
    {
      HistoUpdate(dst,hist0);
    }
    VDestroyImage(dst);
    VDestroyImage(zmap);
  }


  /* apply fdr */
  VImage fdrimage = VCopyImage (dst1,NULL,VAllBands);
  if (numperm > 0) {
    FDR(dst1,fdrimage,hist0,histz,(double)alpha);

    if (cleanup && alpha < 1.0) {
      VIsolatedVoxels(fdrimage,(float)(1.0-alpha));
    }
    if (alpha < 1.0) VImageCount(fdrimage);
  }



  /* write output to disk */
  out_list = VCreateAttrList ();
  VHistory(VNumber(options),options,prg_name,&list1,&out_list);
  VSetGeoInfo(geolist,out_list);
  VAppendAttr (out_list,"image",NULL,VImageRepn,fdrimage);
  fp = VOpenOutputFile (out_filename, TRUE);
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
  fprintf (stderr, "\n%s: done.\n", argv[0]);
  exit(0);
}
