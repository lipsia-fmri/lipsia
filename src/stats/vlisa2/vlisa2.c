/*
** Implementation of LISA algorithm
** for statistical inference of fMRI images
**
** 2nd level inference, two-sample ttest, paired t-test
**
** G.Lohmann, MPI-KYB, 2018
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
#include <gsl/gsl_histogram2d.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

#define ABS(x) ((x) > 0 ? (x) : -(x))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))


extern double ttest2(double *data1,double *data2,int n1,int n2);
extern double xtest2(double *data1,double *data2,int n);
extern double paired_ttest(double *data1,double *data2,int n);
extern double welchtest(double *data1,double *data2,int n1,int n2);
extern void   VIsolatedVoxels(VImage src,float threshold);
extern void   VHistogram(gsl_histogram *histogram,VString filename);
extern void   VCheckImage(VImage src);
extern void   FDR(VImage src,VImage dest,gsl_histogram *nullhist,gsl_histogram *realhist,double);
extern double ttest1(double *data1,int n);
extern void   ImageStats(VImage src,double *,double *,double *hmin,double *hmax);
extern void   VBilateralFilter(VImage src,VImage dest,int radius,double var1,double var2,int);
extern double VImageVar(VImage src);
extern void   VGetHistRange(VImage src,double *hmin,double *hmax);
extern void   VZScale(VImage src,float mode,float stddev);
extern float  VGetMode(VImage src);
extern void   HistoUpdate(VImage,gsl_histogram *);


/* generate permutation table */
int **genperm(long seed,int n1,int n2,int numperm,int testtype)
{
  int i,j,n=n1+n2;
  int **table = (int **) VCalloc(numperm,sizeof(int *));
  int *base = (int *) VCalloc(n,sizeof(int));
  for (j=0; j<n; j++) base[j] = j;

  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);


  /* for each permutation:  */
  for (i = 0; i < numperm; i++) {

    /* unpaired test */
    if (testtype == 0 || testtype == 2) {
      table[i] = (int *) VCalloc(n,sizeof(int));
      gsl_ran_shuffle (rx,base,n,sizeof(int));
      for (j=0; j<n; j++) {
	table[i][j] = base[j];
      }
    }

    /* paired test */
    if (testtype == 1) {
      table[i] = (int *) VCalloc(n1,sizeof(int));
      for (j=0; j<n1; j++) {
	table[i][j] = 0;
	if (gsl_ran_bernoulli (rx,(double)0.5) == 1) table[i][j] = 1;
      }
    }
  }
  VFree(base);
  return table;
}



void TTest(VImage *src1,VImage *src2,int *permtable,VImage dest,int n1,int n2,int testtype)
{
  int i,b,r,c,nslices,nrows,ncols;
  double u=0,v=0,z=0;
  int k1=0,k2=0,k=0;
  int n = n1+n2;
  int min1 = n1-2;
  int min2 = n2-2;

  gsl_set_error_handler_off();

  nslices = VImageNBands(src1[0]);
  nrows   = VImageNRows(src1[0]);
  ncols   = VImageNColumns(src1[0]);
  VFillImage(dest,VAllBands,0);


  /* alloc */
  double *data1 = (double *) VCalloc(n1,sizeof(double));
  double *data2 = (double *) VCalloc(n2,sizeof(double));
  double *data = (double *) VCalloc(n,sizeof(double));


  /* for every voxel */
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	/* unpaired test */
	if (testtype == 0 || testtype == 2) {

	  k = 0;
	  for (i=0; i<n1; i++) {
	    data[k++] = VPixel(src1[i],b,r,c,VFloat);
	  }
	  for (i=0; i<n2; i++) {
	    data[k++] = VPixel(src2[i],b,r,c,VFloat);
	  }

	  k1 = 0;
	  for (i=0; i<n1; i++) {
	    u = data[permtable[i]];
	    if (fabs(u) > 0) {
	      data1[k1++] = u;
	    }
	  }
	  if (k1 < min1) continue;

	  k2 = 0;
	  for (i=0; i<n2; i++) {
	    u = data[permtable[i+n1]];
	    if (fabs(u) > 0) {
	      data2[k2++] = u;
	    }
	  }
	  if (k2 < min2) continue;
	}


	/* paired test */
	if (testtype == 1) {
	  k=0;
	  for (i=0; i<n1; i++) {
	    if (permtable[i] == 0) {  /* no flip */
	      u = VPixel(src1[i],b,r,c,VFloat);
	      v = VPixel(src2[i],b,r,c,VFloat);
	    }
	    else {   /* flip */
	      u = VPixel(src2[i],b,r,c,VFloat);
	      v = VPixel(src1[i],b,r,c,VFloat);
	    }
	    if (fabs(u) > 0 && fabs(v) > 0) {
	      data1[k] = u;
	      data2[k] = v;
	      k++;
	    }
	  }
	  if (k < min1) continue;
	}



	/* t-test */
	z = 0;
	switch(testtype) {
	case 0 :
	  z = ttest2(data1,data2,k1,k2);     /* pooled variance */
	  break;
	case 1 :
	  z = paired_ttest(data1,data2,k);   /* paired */
	  break;
	case 2 :
	  z = welchtest(data1,data2,k1,k2);  /* unequal variance */
	  break;
	default:
	  VError(" unknown testtype");
	}
	VPixel(dest,b,r,c,VFloat) = z;
      }
    }
  }
  VFree(data);
  VFree(data1);
  VFree(data2);
}


VDictEntry TSTDict[] = {
  { "ttest", 0, 0,0,0,0 },
  { "paired", 1, 0,0,0,0 },
  { "welch", 2, 0,0,0,0  },
  { NULL, 0,0,0,0,0 }
};


int main (int argc, char *argv[])
{
  static VArgVector in_files1;
  static VArgVector in_files2;
  static VString  out_filename="";
  static VFloat   alpha = 0.05; 
  static VShort   testtype = 0; 
  static VShort   radius = 2;
  static VFloat   rvar = 2.0;
  static VFloat   svar = 2.0;
  static VShort   numiter = 2;
  static VShort   numperm = 2000;
  static VLong    seed = 99402622;
  static VBoolean centering = FALSE;
  static VBoolean cleanup = TRUE;
  static VShort   nproc = 0;
  static VOptionDescRec options[] = {
    {"in1", VStringRepn, 0, & in_files1, VRequiredOpt, NULL,"Input files 1" },
    {"in2", VStringRepn, 0, & in_files2, VRequiredOpt, NULL,"Input files 2" },
    {"out", VStringRepn, 1, & out_filename, VRequiredOpt, NULL,"Output file" },
    {"alpha",VFloatRepn,1,(VPointer) &alpha,VOptionalOpt,NULL,"FDR significance level"},
    {"perm",VShortRepn,1,(VPointer) &numperm,VOptionalOpt,NULL,"Number of permutations"},
    {"test",VShortRepn,1,(VPointer) &testtype,VOptionalOpt,TSTDict,"type of test"},      
    {"seed",VLongRepn,1,(VPointer) &seed,VOptionalOpt,NULL,"Seed for random number generation"},
    {"radius",VShortRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Bilateral parameter (radius in voxels)"},
    {"rvar",VFloatRepn,1,(VPointer) &rvar,VOptionalOpt,NULL,"Bilateral parameter (radiometric)"},
    {"svar",VFloatRepn,1,(VPointer) &svar,VOptionalOpt,NULL,"Bilateral parameter (spatial)"},
    {"filteriterations",VShortRepn,1,(VPointer) &numiter,VOptionalOpt,NULL,"Bilateral parameter (number of iterations)"},
    {"cleanup",VBooleanRepn,1,(VPointer) &cleanup,VOptionalOpt,NULL,"Whether to apply cleanup"},      
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"},
  };

  FILE *fp=NULL;
  VString in_filename,str1,str2;
  VAttrList list1=NULL,list2=NULL,out_list=NULL,geolist=NULL;
  int i,nimages1=0,nimages2=0,npix=0;
  char *prg_name=GetLipsiaName("vlisa2");
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


  /* number of inputs */
  nimages1 = in_files1.number;
  nimages2 = in_files2.number;
  fprintf(stderr," nimages= %d  %d\n",nimages1,nimages2);


  /* images 1  */
  VImage *src1 = (VImage *) VCalloc(nimages1,sizeof(VImage));
  for (i = 0; i < nimages1; i++) {
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


  /* images 2  */
  VImage *src2 = (VImage *) VCalloc(nimages2,sizeof(VImage));
  for (i = 0; i < nimages2; i++) {
    in_filename = ((VString *) in_files2.vector)[i];
    list2   = VReadAttrList(in_filename,0L,TRUE,FALSE);
    src2[i] = VReadImage(list2);
    if (src2[i] == NULL) VError(" no input image found");
    if (VPixelRepn(src2[i]) != VFloatRepn) VError(" input pixel repn must be float");
    if (i == 0) npix = VImageNPixels(src2[i]);
    else if (npix != VImageNPixels(src2[i])) VError(" inconsistent image dimensions");
  }

  if (testtype == 1) {  /* print filenames to terminal */
    for (i = 0; i < nimages1; i++) {
      str1 = ((VString *) in_files1.vector)[i];
      str2 = ((VString *) in_files2.vector)[i];
      fprintf(stderr," %3d:  %s   %s\n",i,str1,str2);
    }
  }
  else {
    fprintf(stderr," Group 1:\n");
    for (i = 0; i < nimages1; i++) {
      str1 = ((VString *) in_files1.vector)[i];
      fprintf(stderr," %3d:  %s\n",i,str1);
    }
    fprintf(stderr,"\n Group 2:\n");
    for (i = 0; i < nimages2; i++) {
      str2 = ((VString *) in_files2.vector)[i];
      fprintf(stderr," %3d:  %s\n",i,str2);
    }
  }


  /* random permutations */
  size_t n = nimages1+nimages2;
  int nperm=0;
  int **permtable = genperm((long)seed,(int)nimages1,(int)nimages2,(int)numperm,(int)testtype);
  int *nopermtable = (int *) VCalloc(n,sizeof(int));
  for (i=0; i<n; i++) nopermtable[i]=i;
  if (testtype == 1) for (i=0; i<n; i++) nopermtable[i]=0;


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
      TTest(src1,src2,permtable[nperm],zmap,nimages1,nimages2,(int)testtype); 
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
  VImage dst1  = VCreateImageLike (src1[0]);
  VImage zmap1 = VCreateImageLike(src1[0]);
  VFillImage(zmap1,VAllBands,0);
  TTest(src1,src2,nopermtable,zmap1,nimages1,nimages2,(int)testtype);

  if (numperm == 0) {
    double z = VImageVar(zmap1);
    stddev = (float)(sqrt(z)); /* update stddev */
  }
  float mode=0;
  if (centering) mode = VGetMode(zmap1); 
  if (numperm > 0) VZScale(zmap1,mode,stddev);
  VBilateralFilter(zmap1,dst1,(int)radius,(double)rvar,(double)svar,(int)numiter);


  /* ini histograms */
  VGetHistRange(dst1,&hmin,&hmax);
  /* fprintf(stderr," Histogram range:  [%.3f, %.3f]\n",hmin,hmax); */
  size_t nbins = 10000;
  gsl_histogram *hist0 = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (hist0,hmin,hmax);
  gsl_histogram *histz = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (histz,hmin,hmax);
  HistoUpdate(dst1,histz);


#pragma omp parallel for shared(src1,src2,permtable) schedule(dynamic)
  for (nperm = 0; nperm < numperm; nperm++) {
    if (nperm%20 == 0) fprintf(stderr," perm  %4d  of  %d\r",nperm,(int)numperm);
 
    VImage zmap = VCreateImageLike(src1[0]);
    VImage dst  = VCreateImageLike (zmap);
    TTest(src1,src2,permtable[nperm],zmap,nimages1,nimages2,(int)testtype);    

    float mode=0;
    if (centering) mode = VGetMode(zmap);
    VZScale(zmap,mode,stddev);
    VBilateralFilter(zmap,dst,(int)radius,(double)rvar,(double)svar,(int)numiter);


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
  }


  /* write output to disk */
  out_list = VCreateAttrList ();
  VHistory(VNumber(options),options,prg_name,&list1,&out_list);
  VSetGeoInfo(geolist,out_list);
  VAppendAttr (out_list,"image",NULL,VImageRepn,fdrimage);
  fp = VOpenOutputFile (out_filename, TRUE);
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
  fprintf (stderr, "\n");
  fprintf (stderr, "%s: done.\n", argv[0]);
  exit(0);
}
