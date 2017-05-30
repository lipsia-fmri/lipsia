/*
** Task-related edge density, main algorithm
**
** Ref:  Lohmann et al (2016) PLoS One, 11(6):e0158185
** 
** G.Lohmann, MPI-KYB, Jan 2015
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_float.h>

#include "viaio/Vlib.h"
#include "viaio/VImage.h"
#include "viaio/mu.h"
#include "viaio/option.h"


#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))

extern VImage *VImagePointer(VAttrList list,int *nt);
extern void   VReadImagePointer(VAttrList list,VImage *src);
extern VImage VoxelMap(VImage mask,size_t *nvoxels);
extern void   VDataMatrix(VImage *src,int first,int len,VImage map,gsl_matrix_float *X);
extern void   VCheckMatrix(gsl_matrix_float *X);
extern long   VMaskCoverage(VAttrList list,VImage mask);
extern float  ZMatrix(gsl_histogram *histogram,gsl_matrix_float *SNR1,gsl_matrix_float *SNR2,int n1,int n2,
		      VImage roi,VImage map,VImage,int,float elength,float quantile,int,int);
extern void   GetSNR(gsl_matrix_float **X1,gsl_matrix_float **X2,int *,int n,gsl_matrix_float *SNR,int);
extern void   GetMedian(gsl_matrix_float **X1,gsl_matrix_float **X2,int *,int n,gsl_matrix_float *SNR,int);
extern void   VPrintHistogram(gsl_histogram *histogram,int numperm,VString filename);
extern void   HistoUpdate(float *A,size_t nvox,gsl_histogram *hist);

extern size_t EdgeDensity(float *C,int *I,int *J,size_t nedges,
			  gsl_histogram *histogram,gsl_matrix_float *SNR1,gsl_matrix_float *SNR2,int n1,int n2,
			  VImage roi,VImage map,VImage,int,float,float,int,int,int,float,int);

extern size_t EstimateEdges(float *E,int *I,int *J,size_t fulldim,
			    gsl_histogram *TedHist,gsl_matrix_float *SNR1,gsl_matrix_float *SNR2,int n1,int n2,
			    VImage roi,VImage map,VImage mapimage,int adjdef,float elength,float zthreshold,float,int);


/* generate permutation table */
int genperm(gsl_rng *rx,int *table,int n,int numperm)
{
  int s=0;
  int ks = 0;
  int kc=3;
  if (numperm < 10) {
    if (n%2==0) kc = n/2;
    else kc = n/2-1;
  }
  int iter=0;
  while ((ks < kc || ks > n-kc)) {
    if (iter > 100) VError(" not possible to generate random permutations, perhaps not enough input data ?");
    ks=0;
    for (s=0; s<n; s++) {
      table[s] = 0;
      if (gsl_ran_bernoulli (rx,(double)0.5) == 1) {
	table[s] = 1;
	ks++;
      }
    }
    iter++;
  }
  return ks;
}


void VReleaseStorage(VAttrList list) 
{
  VImage src;
  VAttrListPosn posn;

  for(VFirstAttr(list, & posn); VAttrExists(& posn); VNextAttr(& posn)) {
    if(VGetAttrRepn(& posn) != VImageRepn) continue;
    VGetAttrValue(& posn, NULL, VImageRepn, & src);
    VDestroyImage(src);
    src = NULL;
  }
}


/*  Check mask coverage */
VAttrList VCheckMaskCoverage(VStringConst *in_filenames,size_t n,VImage mask)
{
  size_t i;  
  VAttrList list=NULL;
  VAttrList geolist = NULL;
  FILE *fp=NULL;

  for (i=0; i<n; i++) {
    fp = VOpenInputFile (in_filenames[i], TRUE);
    list = VReadFile (fp, NULL);
    if (! list)  VError("Error reading image");
    fclose(fp);

    /* read geometry info, if unknown */
    if (geolist == NULL) geolist = VGetGeoInfo(list);

    long count = VMaskCoverage(list,mask);
    if (count > 100) fprintf(stderr," incomplete mask_coverage  %s:  %ld\n",in_filenames[i],count);
    VReleaseStorage(list);
    list = NULL;
  }
  return geolist;
}

gsl_matrix_float **VReadImageData(VStringConst *in_filenames,VImage *src,VImage map,size_t n,size_t nvox,size_t first,size_t len)
{
  size_t i;
  VAttrList list=NULL;
  FILE *fp=NULL;

  /* allocate data struct */
  gsl_matrix_float **X = VCalloc(n,sizeof(gsl_matrix_float));
  for (i=0; i<n; i++) {
    X[i] = gsl_matrix_float_calloc(nvox,len);
    if (!X[i]) VError(" err allocating data matrix");
  }
  /* read data */
  for (i=0; i<n; i++) {
    fp = VOpenInputFile (in_filenames[i], TRUE);
    list = VReadFile (fp, NULL);
    if (! list)  VError("Error reading image");
    fclose(fp);
    VReadImagePointer(list,src);
    VDataMatrix(src,(int)first,(int)len,map,X[i]);
    VCheckMatrix(X[i]);
    VReleaseStorage(list);
    list = NULL;
  }
  return X;
}




VDictEntry ADJDict[] = {
  { "6", 0 },
  { "18", 1 },
  { "26", 2 },
  { NULL }
};

VDictEntry TypeDict[] = {
  { "SNR", 0 },
  { "median", 1 },
  { NULL }
};

VDictEntry MetricDict[] = {
  { "pearson", 0 },
  { "spearman", 1 },
  { NULL }
};


int main (int argc,char *argv[])
{
  static VArgVector in_files1;
  static VArgVector in_files2;
  static VString  out_filename = "";
  static VString  mask_filename = "";
  static VString  roi_filename = "";
  static VString  hist_filename= "";
  static VShort   first = 0;
  static VShort   len = 0;
  static VLong    seed = 99402622;
  static VFloat   qthreshold = 0.99;
  static VFloat   elength = 5;
  static VFloat   noise_cutoff = 0.01;
  static VShort   adjdef = 2;
  static VShort   numperm = 0;
  static VShort   type = 0;
  static VShort   metric = 0;
  static VShort   nproc = 10;
  static VShort   step = 2;
  static VOptionDescRec  options[] = {
    {"in1", VStringRepn, 0, & in_files1, VRequiredOpt, NULL,"Input files 1" },
    {"in2", VStringRepn, 0, & in_files2, VRequiredOpt, NULL,"Input files 2" },
    {"out", VStringRepn, 1, & out_filename, VOptionalOpt, NULL,"Output file" },
    {"mask", VStringRepn, 1, & mask_filename, VRequiredOpt, NULL,"Mask file" },
    {"roi", VStringRepn, 1, & roi_filename, VOptionalOpt, NULL,"ROI file" },
    {"histogram",VStringRepn,1,(VPointer) &hist_filename,VRequiredOpt,NULL,"Output histogram filename"},
    {"permutations",VShortRepn,1,(VPointer) &numperm,VOptionalOpt,NULL,"Number of permutations"},
    {"qthreshold",VFloatRepn,1,(VPointer) &qthreshold,VOptionalOpt,NULL,"Initial quantile threshold"}, 
    {"adj",VShortRepn,1,(VPointer) &adjdef,VOptionalOpt,ADJDict,"Definition of adjacency"},
    {"seed",VLongRepn,1,(VPointer) &seed,VOptionalOpt,NULL,"Seed for random number generator"},
    {"first",VShortRepn,1,(VPointer) &first,VOptionalOpt,NULL,"First timepoint to use within trial"},
    {"len",VShortRepn,1,(VPointer) &len,VOptionalOpt,NULL,"Number of timepoints to use, '0' to use all"},
    {"type",VShortRepn,1,(VPointer) &type,VOptionalOpt,TypeDict,"Type of trial average"},
    {"metric",VShortRepn,1,(VPointer) &metric,VOptionalOpt,MetricDict,"Correlation metric"},
    {"edgelength",VFloatRepn,1,(VPointer) &elength,VOptionalOpt,NULL,"Minimal edge length in voxels"},   
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"},
    /*     {"noisecutoff",VFloatRepn,1,(VPointer) &noise_cutoff,VOptionalOpt,NULL,"estimate of noisy edges"}, */
    /*     {"step",VShortRepn,1,(VPointer) &step,VOptionalOpt,NULL,"first iteration ZMatrix step size"} */
  };
  FILE *fp=NULL;
  VAttrList list=NULL,list1=NULL,geolist=NULL;
  VAttrListPosn posn;
  VImage mask=NULL;
  size_t nvox=0,i=0;


  /* parse command line */
  if (! VParseCommand (VNumber (options), options, & argc, argv)) {
    VReportUsage (argv[0], VNumber (options), options, NULL);
    exit (EXIT_FAILURE);
  }
  if (argc > 1) {
    VReportBadArgs (argc, argv);
    exit (EXIT_FAILURE);
  }


  /* whether to produce output file, or only compute histogram */
  VBoolean histonly = FALSE;
  if (strlen(out_filename) < 2) histonly = TRUE;



  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  int jprocs = num_procs;
  if (nproc > 0 && nproc < num_procs) jprocs = nproc;
  fprintf(stderr," using %d cores of %d\n",(int)jprocs,(int)num_procs);
  omp_set_num_threads(jprocs);
#endif /* _OPENMP */



  /* input filenames */
  size_t n1 = (size_t)in_files1.number;
  size_t n2 = (size_t)in_files2.number;
  if (n1 != n2) VError(" n1 != n2, %d %d",n1,n2);

  VStringConst *in_filenames1 = (VStringConst *) VCalloc(n1,sizeof(VStringConst));
  for (i=0; i<n1; i++) {
    in_filenames1[i] = ((VStringConst *) in_files1.vector)[i];
  }
  VStringConst *in_filenames2 = (VStringConst *) VCalloc(n2,sizeof(VStringConst));
  for (i=0; i<n2; i++) {
    in_filenames2[i] = ((VStringConst *) in_files2.vector)[i];
  }

  

  /* read mask */
  fp = VOpenInputFile (mask_filename, TRUE);
  list1 = VReadFile (fp, NULL);
  if (! list1) VError("Error reading mask file");
  fclose(fp);

  for (VFirstAttr (list1, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & mask);
    if (VPixelRepn(mask) == VFloatRepn || VPixelRepn(mask) == VDoubleRepn) {
      mask = NULL;
      continue;
    }
  }
  if (mask == NULL) VError(" no mask found");
  int nslices = VImageNBands(mask);
  int nrows = VImageNRows(mask);
  int ncols = VImageNColumns(mask);


  /* read ROI */
  VAttrList list2 = NULL;
  VImage roi = NULL;
  if (strlen(roi_filename) > 1) {
    fp = VOpenInputFile (roi_filename, TRUE);
    list2 = VReadFile (fp, NULL);
    if (! list2) VError("Error reading roi file");
    fclose(fp);

    for (VFirstAttr (list2, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      VGetAttrValue (& posn, NULL,VImageRepn, & roi);
      if (VPixelRepn(roi) == VFloatRepn || VPixelRepn(roi) == VDoubleRepn) {
	roi = NULL;
	continue;
      }
    }
    if (roi == NULL) VError(" no roi found");
  }

  /* Check mask coverage */
  VAttrList geo=NULL;
  geo = VCheckMaskCoverage(in_filenames1,n1,mask);
  geo = VCheckMaskCoverage(in_filenames2,n2,mask);
  if (geolist == NULL) geolist = geo;



  /* voxel map */
  VImage map = VoxelMap(mask,&nvox);


  /* get image dimensions */
  fp = VOpenInputFile (in_filenames1[0], TRUE);
  list = VReadFile (fp, NULL);
  if (! list)  VError("Error reading image");
  fclose(fp);

  int nt=0;
  VImage *src = VImagePointer(list,&nt);

  if (first >= nt || first < 0) VError(" illegal value, first= %d, nt= %d",first,nt);
  if (len <= 0) len = nt-first;
  if (first + len >= nt) len = nt-first;
  if (len < 2) VError(" len= %d",len);
  fprintf(stderr," image: %d x %d x %d, nt: %d, nvox: %ld\n",(int)nslices,(int)nrows,(int)ncols,nt,nvox);

  /* read image data */
  gsl_matrix_float **X1 = VReadImageData(in_filenames1,src,map,n1,nvox,first,len);
  gsl_matrix_float **X2 = VReadImageData(in_filenames2,src,map,n2,nvox,first,len);

 
  /* voxel addresses */
  VImage mapimage = VCreateImage(nslices,nrows,ncols,VIntegerRepn);
  VFillImage(mapimage,VAllBands,(VInteger)(-1));
  for (i=0; i<nvox; i++) {
    int b = VPixel(map,0,0,i,VShort);
    int r = VPixel(map,0,1,i,VShort);
    int c = VPixel(map,0,2,i,VShort);
    VPixel(mapimage,b,r,c,VInteger) = (VInteger)i;
  }


  /* tmp storage for SNR time courses */
  int dim=0;
  if (type < 3) dim = len;
  else if (type == 3) dim = len*n1;
  else VError(" illegal type %d",type);
  gsl_matrix_float *SNR1 = gsl_matrix_float_calloc(nvox,dim);
  gsl_matrix_float *SNR2 = gsl_matrix_float_calloc(nvox,dim);


  /* ini zhist */
  size_t hbins = 10000;
  double hmin = 0.0,hmax = 1.001;
  gsl_histogram *TedHist = gsl_histogram_alloc (hbins);
  gsl_histogram_set_ranges_uniform (TedHist,hmin,hmax);
  gsl_histogram_reset(TedHist);


  /* fuer z-threshold */
  size_t nbins = 25000;
  hmin = -1.001,hmax = 1.001;
  gsl_histogram *ZvalHist = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (ZvalHist,hmin,hmax);
  gsl_histogram_reset(ZvalHist);


  /* ini random */
  size_t n = n1;
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);
  int *table = (int *) VCalloc(n,sizeof(int));
  for (i=0; i<n; i++) table[i] = 0;


  /* allocate output structure */
  size_t nedges=0;
  size_t fulldim = (nvox*(nvox+1)/2);
  double qx = 1.0-qthreshold;
  qx += 0.001;
  size_t nedges_estimated = (size_t) (qx*(double)fulldim);


  /* ini matrices */
  float *E=NULL;   /* matrix of edge densities */
  int *I=NULL;     /* row indices in matrix E */
  int *J=NULL;     /* col indices in matrix E */

  if (histonly == FALSE && noise_cutoff < 0.00001) {
    E = VCalloc(nedges_estimated,sizeof(float));  /* matrix of edge densities */
    I = VCalloc(nedges_estimated,sizeof(int));    /* row indices in matrix E */
    J = VCalloc(nedges_estimated,sizeof(int));    /* col indices in matrix E */
  }


  /* main loop */
  size_t s;
  int ks=0,nperm=0;
  int startperm=0;
  if (numperm > 0) startperm = 1;


  for (nperm = startperm; nperm <= numperm; nperm++) {

    fprintf(stderr,"\n perm %3d: ",nperm);
    if (nperm > 0) {
      ks = genperm(rx,table,n,(int)numperm);
    }
    for (s=0; s<n; s++) {
      fprintf(stderr,"%d",table[s]);
    }
    fprintf(stderr,"  %3d\n",ks);


    /* get SNR data */
    if (type == 0) {
      GetSNR(X1,X2,table,n,SNR1,metric);
      GetSNR(X2,X1,table,n,SNR2,metric);
    }
    if (type == 1) {
      GetMedian(X1,X2,table,n,SNR1,metric);
      GetMedian(X2,X1,table,n,SNR2,metric);
    }


    /* corr matrices */
    float zthr = ZMatrix(ZvalHist,SNR1,SNR2,n1,n2,roi,map,mapimage,(int)adjdef,
			 (float)elength,(float)qthreshold,(int)step,(int)metric);
    
    /* estimate number of truly needed edges, estimate noise */
    if (noise_cutoff > 0.0 && nperm == numperm && histonly == FALSE) {
      size_t old_estimate = nedges_estimated;
      nedges_estimated = EstimateEdges(E,I,J,old_estimate,TedHist,SNR1,SNR2,n1,n2,roi,map,mapimage,
				       (int)adjdef,elength,zthr,noise_cutoff,(int)metric);
      
      E = VCalloc(nedges_estimated,sizeof(float));  /* matrix of edge densities */
      I = VCalloc(nedges_estimated,sizeof(int));    /* row indices in matrix E */
      J = VCalloc(nedges_estimated,sizeof(int));    /* col indices in matrix E */
    }


    /* edge densities */
    fprintf(stderr," Computing edge densities...\n");
    nedges = EdgeDensity(E,I,J,nedges_estimated,TedHist,SNR1,SNR2,n1,n2,roi,map,mapimage,(int)adjdef,
			 elength,zthr,nperm,numperm,(int)1,noise_cutoff,(int)metric);
  }

  
  /* free storage */
  for (i=0; i<n; i++) {
    gsl_matrix_float_free(X1[i]);
    gsl_matrix_float_free(X2[i]);
  }
  gsl_matrix_float_free(SNR1);
  gsl_matrix_float_free(SNR2);
  VFree(table);


  /* print histogram */
  VPrintHistogram(TedHist,(int)numperm,hist_filename);

  if (histonly == TRUE) {
    fprintf (stderr," %s: done.\n", argv[0]);
    exit(0);
  }


  /* write edge density matrix in triplet sparse format */
  size_t edim = nedges*3;
  float filesize1 = (float)(edim*sizeof(float))/(1000.0*1000.0);
  float filesize2 = (float)(edim*sizeof(float))/(1000.0*1000.0*1000.0);
  fprintf(stderr," nedges= %ld,  filesize= %.3f MByte  (%.3f GByte)\n",nedges,filesize1,filesize2);


  /* output */
  VAttrList elist = VCreateAttrList();
  VAttrList ilist = VCreateAttrList();
  VAttrList jlist = VCreateAttrList();
  VBundle edgedens = VCreateBundle ("data",elist,nedges*sizeof(float),(VPointer)E);
  VBundle rowindex = VCreateBundle ("data",ilist,nedges*sizeof(int),(VPointer)I);
  VBundle colindex = VCreateBundle ("data",jlist,nedges*sizeof(int),(VPointer)J);

  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"EdgeDensity",NULL,VBundleRepn,(VBundle)edgedens);
  VAppendAttr(out_list,"RowIndex",NULL,VBundleRepn,(VBundle)rowindex);
  VAppendAttr(out_list,"ColIndex",NULL,VBundleRepn,(VBundle)colindex);
  VAppendAttr(out_list,"nedges",NULL,VLongRepn,(VLong)nedges);
  VAppendAttr(out_list,"map",NULL,VImageRepn,map);
  VAppendAttr(out_list,"qthreshold",NULL,VFloatRepn,(VFloat)qthreshold);

  FILE *fp_out = VOpenOutputFile (out_filename, TRUE);
  if (! VWriteFile (fp_out, out_list)) exit (1);
  fclose(fp_out);
  fprintf (stderr," %s: done.\n\n", argv[0]);
  exit(0);
}
