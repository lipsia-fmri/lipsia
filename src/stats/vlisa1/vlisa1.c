/*
** Implementation of LISA algorithm
** for statistical inference of fMRI images
**
** 2nd level inference applied to one sample
**
** G.Lohmann, 2017
*/

#include "viaio/Vlib.h"
#include "viaio/file.h"
#include "viaio/mu.h"
#include "viaio/option.h"
#include "viaio/os.h"
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

#define Median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))
#define ABS(x) ((x) > 0 ? (x) : -(x))

extern double t2z(double,double);
extern void VIsolatedVoxels(VImage src,VImage,float threshold);
extern void VHistogram(gsl_histogram *histogram,VString filename);
extern void VCheckImage(VImage src);
extern void FDR(VImage src,VImage dest,double alpha,gsl_histogram *nullhist,gsl_histogram *realhist,VString filename);
extern double ttest1(double *data1,int n);
extern void ImageStats(VImage src,double *,double *,double *hmin,double *hmax);
extern void Hotspot(VImage src,VImage dst,VImage tmp,int type,int radius,double var1,double var2,int numiter);


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



/* update histogram */
void HistoUpdate(VImage src1,gsl_histogram *hist)
{
  float u,tiny = 1.0e-6;
  size_t i;
  float xmin = gsl_histogram_min (hist);
  float xmax = gsl_histogram_max (hist);

  float *pp1 = VImageData(src1);
  for (i=0; i<VImageNPixels(src1); i++) {
    u = *pp1++;
    if (ABS(u) < tiny) continue;
    if (u > xmax) u = xmax-tiny;
    if (u < xmin) u = xmin+tiny;
    gsl_histogram_increment (hist,u);
  }
}

void OnesampleTest(VImage *src1,VImage dest,int *permtable,int n)
{
  int i,k,b,r,c,nslices,nrows,ncols;
  double s1,s2,nx,mean,var,u,t;
  double tiny=1.0e-8;

  nslices = VImageNBands(src1[0]);
  nrows   = VImageNRows(src1[0]);
  ncols   = VImageNColumns(src1[0]);
  VFillImage(dest,VAllBands,0);

  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	k = 0;
	s1 = s2 = 0;
	for (i=0; i<n; i++) {
	  u = (double)VPixel(src1[i],b,r,c,VFloat);
	  if (ABS(u) > tiny) {
	    if (permtable[i] > 0) u = -u;
	    s1 += u;
	    s2 += u*u;
	    k++;
	  }
	}
	if (k < n-2) continue;
	nx = (double)k;
	mean = s1/nx;
	var = (s2 - nx * mean * mean) / (nx - 1.0);
	t  = sqrt(nx) * mean/sqrt(var);
	VPixel(dest,b,r,c,VFloat) = t;
      }
    }
  }
}


int main (int argc, char *argv[])
{
  static VArgVector in_files1;
  static VString  out_filename="";
  static VFloat   alpha = 0.05; 
  static VShort   radius = 2;
  static VString  fdrfilename= "";
  static VFloat   rvar = 1.9;
  static VFloat   svar = 2.0;
  static VShort   numiter = 2;
  static VShort   numperm = 2000;
  static VLong    seed = 99402622;
  static VBoolean cleanup = TRUE;
  static VShort   nproc = 10;
  static VOptionDescRec options[] = {
    {"in", VStringRepn, 0, & in_files1, VRequiredOpt, NULL,"Input files" },
    {"out", VStringRepn, 1, & out_filename, VOptionalOpt, NULL,"Output file" },
    {"alpha",VFloatRepn,1,(VPointer) &alpha,VOptionalOpt,NULL,"FDR significance level"},
    {"perm",VShortRepn,1,(VPointer) &numperm,VOptionalOpt,NULL,"Number of permutations"},    
    {"seed",VLongRepn,1,(VPointer) &seed,VOptionalOpt,NULL,"Seed for random number generation"},    
    {"radius",VShortRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Neighbourhood radius in voxels"},
    {"rvar",VFloatRepn,1,(VPointer) &rvar,VOptionalOpt,NULL,"Bilateral parameter (radiometric)"},
    {"svar",VFloatRepn,1,(VPointer) &svar,VOptionalOpt,NULL,"Bilateral parameter (spatial)"},
    {"numiter",VShortRepn,1,(VPointer) &numiter,VOptionalOpt,NULL,"Number of iterations in bilateral filter"},  
    {"cleanup",VBooleanRepn,1,(VPointer) &cleanup,VOptionalOpt,NULL,"Whether to apply cleanup"},      
    {"filename",VStringRepn,1,(VPointer) &fdrfilename,VOptionalOpt,NULL,"Name of output fdr txt-file"},    
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"},
  };

  FILE *fp=NULL;
  VStringConst in_filename;
  VAttrList list1=NULL,out_list=NULL,geolist=NULL;
  VAttrListPosn posn;
  VImage src=NULL,*src1;
  int i,nimages,npix=0;
  char *prg_name=GetLipsiaName("vlisa1");
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
    src1[i] = NULL;
    in_filename = ((VStringConst *) in_files1.vector)[i];
    fp = VOpenInputFile (in_filename, TRUE);
    list1 = VReadFile (fp, NULL);
    if (! list1)  VError("Error reading image");
    fclose(fp);

    /* use geometry info from 1st file */
    if (geolist == NULL) geolist =  VGetGeoInfo(list1);


    for (VFirstAttr (list1, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      VGetAttrValue (& posn, NULL, VImageRepn, & src);
      if (VPixelRepn(src) != VFloatRepn) continue;

      if (i == 0) npix = VImageNPixels(src);
      else if (npix != VImageNPixels(src)) VError(" inconsistent image dimensions");

      src1[i] = src;
      break;
    }
    if (src1[i] == NULL) VError(" no image found in %s",in_filename);
  }


  /* random permutations */
  size_t n = nimages;
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);
  int nperm=0;
  int **permtable = genperm(rx,(int)n,(int)numperm);


  /* no permutation */
  int *nopermtable = (int *) VCalloc(n,sizeof(int));
  VImage tmpimage = VCreateImageLike(src1[0]);
  VImage dst1 = VCreateImageLike (src1[0]);
  VImage zmap1 = VCreateImageLike(src1[0]);
  VFillImage(dst1,VAllBands,0);
  VFillImage(zmap1,VAllBands,0);
  OnesampleTest(src1,zmap1,nopermtable,nimages);


  /* adjust radiometric parameter for bilateral filter if needed */
  double vmin=0,vmax=0,hmin=0,hmax=0,gave=0,gvar=0,globalvar=0;
  double xrvar = rvar;
  ImageStats(zmap1,&gave,&globalvar,&vmin,&vmax); 
  xrvar *= globalvar;


  /* get non-permuted hotspot map */
  int type = 0;
  if (cleanup) type = 1;
  Hotspot(zmap1,dst1,tmpimage,(int)type,(int)radius,(double)xrvar,(double)svar,(int)numiter);
  VDestroyImage(tmpimage);


  /* ini histograms */
  size_t nbins = 10000;
  ImageStats(dst1,&gave,&gvar,&vmin,&vmax);
  hmin = vmin - 4.0;
  hmax = vmax + 4.0;
  if (hmin < -12.0) hmin = -12.0;
  if (hmax > 12.0) hmax = 12.0;

  gsl_histogram *hist0 = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (hist0,hmin,hmax);
  gsl_histogram *histz = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (histz,hmin,hmax);
  HistoUpdate(dst1,histz);
  type = 0;  /* no cleanup for permutations */


  /* random permutations */
#pragma omp parallel for shared(src1,permtable) schedule(dynamic)
  for (nperm = 0; nperm < numperm; nperm++) {
    if (nperm%20 == 0) fprintf(stderr," perm  %4d  of  %d\r",nperm,(int)numperm);

    VImage tmpimage = NULL;
    VImage zmap = VCopyImage(src1[0],NULL,VAllBands);
    VFillImage(zmap,VAllBands,0);
    OnesampleTest(src1,zmap,permtable[nperm],nimages);
    VImage dst = VCreateImageLike (zmap);
    tmpimage = VCreateImageLike (zmap);
    Hotspot(zmap,dst,tmpimage,(int)type,(int)radius,(double)rvar,(double)svar,(int)numiter);

#pragma omp critical 
    {
      HistoUpdate(dst,hist0);
    }

    VDestroyImage(dst);
    VDestroyImage(tmpimage);
    VDestroyImage(zmap);
  }
  

  /* apply fdr */
  VImage fdrimage = VCopyImage (dst1,NULL,VAllBands);
  if (numperm > 0) FDR(dst1,fdrimage,(double)alpha,hist0,histz,fdrfilename);


  /* remove isolated voxels */
  VImage fdrimage_clean = VCopyImage(fdrimage,NULL,VAllBands);
  if (cleanup && alpha < 1.0)
    VIsolatedVoxels(fdrimage,fdrimage_clean,(float)(1.0-alpha));


  /*
  ** output
  */
  out_list = VCreateAttrList ();
  VHistory(VNumber(options),options,prg_name,&list1,&out_list);
  VSetGeoInfo(geolist,out_list);
  VAppendAttr (out_list,"image",NULL,VImageRepn,fdrimage_clean);
  fp = VOpenOutputFile (out_filename, TRUE);
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
  fprintf (stderr, "\n%s: done.\n", argv[0]);
  exit(0);
}
