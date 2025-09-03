/*
** Cylarim:  Cylinder-based Algorithm for Laminar fMRI Data Analysis
**
** G.Lohmann, MPI-KYB, Oct 2024
*/


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>

#include "../cylutils/cyl.h"

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

extern void GetResolution(VImage src,gsl_vector *reso);
extern VImage Convert2Repn(VImage src,VImage dest,VRepnKind repn);
extern Cylinders *VCylinder(VImage rim,VImage metric,double radius);
extern void HistEqualize(Cylinders *,VImage,VImage,VImage);

extern double LaminarMedian(gsl_vector *,gsl_vector *,size_t,long,gsl_vector *,gsl_vector *);
extern double LaminarNBins(gsl_vector *,gsl_vector *,gsl_vector *,int);
extern double LaminarConvex(gsl_vector *,gsl_vector *,gsl_vector *,int);
extern double LaminarPeak(gsl_vector *,gsl_vector *,gsl_vector *,int);
extern double LaminarFit(gsl_vector *,gsl_vector *,gsl_vector *);
extern double LaminarLinear(gsl_vector *,gsl_vector *,gsl_vector *);
extern double LaminarZMax(gsl_vector *y);
extern VAttrList VConcatImages(VImage *src,VAttrList geolist,int n);
extern void strip_filename(char *path);


void XWriteOutput(VImage image,VAttrList geolist,char *filename)
{
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,image);
  FILE *fp = fopen(filename,"w");
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
}

/* count number of data points in current cylinder */
size_t NumPoints(VImage zmap,VImage metric,Cylinders *cyl,size_t cid)
{
  size_t m = cyl->addr[cid]->size;
  size_t i,k;
  int b,r,c;
  double w=0,z=0;
  size_t n=0;
  for (i=0; i<m; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    w = VPixel(metric,b,r,c,VFloat);
    if (w < 0.05 || w > 0.95) continue;
    z = VPixel(zmap,b,r,c,VFloat);
    if (fabs(z) < TINY) continue;
    n++;
  }
  return n;
}


/* fill vectors as prep for laminar statistics */
int PrepStats(VImage zmap,VImage metric,Cylinders *cyl,size_t cid,gsl_vector *y,gsl_vector *mvec)
{
  size_t m = cyl->addr[cid]->size;
  size_t i,j,k;
  int b,r,c;
  double w=0,z=0;

  int verbose=0;
  j=0;
  for (i=0; i<m; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    w = VPixel(metric,b,r,c,VFloat);
    if (w < 0.05 || w > 0.95) continue;
    z = VPixel(zmap,b,r,c,VFloat);
    if (fabs(z) < TINY) continue;
    mvec->data[j] = w;
    y->data[j] = z;
    j++;
  }
  return verbose;
}


void WriteToFile(VImage *betaimage,VAttrList geolist,int nbeta,
		 float radius,VBoolean equivol,VLong numperm,VLong seed,
		 VString name,VString out_filename)
{
  char cbuf[5000];
  VAttrList out_list = NULL;
  if (nbeta > 1) out_list = VConcatImages(betaimage,geolist,nbeta);
  else {
    out_list = VCreateAttrList();
    VSetGeoInfo(geolist,out_list);
    VAppendAttr(out_list,"image",NULL,VImageRepn,betaimage[0]);
  }
  if (numperm > 0) {
    VPrependAttr(out_list,"numperm",NULL,VLongRepn,numperm);
    VPrependAttr(out_list,"seed",NULL,VLongRepn,seed);
  }
  memset(cbuf,'\0',sizeof(cbuf));
  sprintf(cbuf,"%.4f",radius);
  VPrependAttr(out_list,"radius",NULL,VStringRepn,cbuf);
  if (equivol==TRUE) VPrependAttr(out_list,"equivolume",NULL,VStringRepn,(VString)"true");
  VPrependAttr(out_list,"name",NULL,VStringRepn,name);  

  memset(cbuf,'\0',sizeof(cbuf));
  sprintf(cbuf,"%s_%s.v",out_filename,name);
  fprintf(stderr," writing %s\n",cbuf);
  FILE *fp = VOpenOutputFile (cbuf, TRUE);
  if (! VWriteFile (fp,out_list)) exit (1);
  fclose(fp);
}

void WriteStats(Cylinders *cyl,size_t k,VImage rim,gsl_vector *beta,VImage *coeff)
{
  size_t i,j;
  for (i=0; i<cyl->addr[k]->size; i++) {
    size_t l = cyl->addr[k]->data[i];
    int b = gsl_matrix_int_get(cyl->xmap,l,0);
    int r = gsl_matrix_int_get(cyl->xmap,l,1);
    int c = gsl_matrix_int_get(cyl->xmap,l,2);

#pragma omp critical
    {
      for (j=0; j<beta->size; j++) {
	VPixel(coeff[j],b,r,c,VFloat) += beta->data[j];
      }
    }
  }
}


/* normalize images */
void WNormalizeImage(VImage *coeff,int ncoeff,VImage wimage,VImage rim)
{
  size_t i,j;
  VFloat *pw = VImageData(wimage);
  VUByte *pu = VImageData(rim);
  VFloat *pb = NULL;
  float wmin = 0.01;
  for (j=0; j<ncoeff; j++) {
    pb = VImageData(coeff[j]);
    for (i=0; i<VImageNPixels(wimage); i++) {
      if (pw[i] > wmin) pb[i] /= pw[i];
      else pb[i] = 0;
      if (gsl_finite((double)pb[i])==0) pb[i] = 0;
      if (pu[i] != 3) pb[i] = 0;  /* exclude rim border voxels */
    }
  }
}



VImage Cylarim(VImage zmap,VImage metric,VImage rim,double radius,
	       VBoolean equivol,size_t numperm,long seed,
	       VBoolean x3bins,int nbins,VBoolean xpeak,VBoolean xconvex,VBoolean xlinear,VBoolean xR2,
	       VImage *bin3_coeff,VImage *bin3_perm,VImage *nbin_coeff,VImage *peak_coeff,
	       VImage convex_coeff,VImage linear_coeff,
	       VImage *R2_coeff)
{
  size_t k;
  
  /* create cylinder data struct */
  Cylinders *cyl = VCylinder(rim,metric,radius);

  
  /* equivolume correction */
  VImage wimage = VCreateImageLike(zmap);
  VFillImage(wimage,VAllBands,0);
  if (equivol) {
    fprintf(stderr," Equivolume...\n");
    HistEqualize(cyl,wimage,metric,rim);
  }
  
  /* LayerStats */  
  size_t progress=0;
  size_t np = (size_t)((float)cyl->numcylinders/100.0);
  if (np < 1) np = 1;
  VFillImage(wimage,VAllBands,0);
  VImage zabsimage = VCreateImageLike(wimage);
  VFillImage(zabsimage,VAllBands,0);

 
#pragma omp parallel for shared(progress)
  for (k=0; k<cyl->numcylinders; k++) {
    size_t i;
    if (k%np==0) fprintf(stderr," LaminarAnalysis:  %7.3f\r",(float)(progress)/(float)cyl->numcylinders); 
    progress++;

    size_t n = NumPoints(zmap,metric,cyl,k);
    if (n < NPTS) continue;
    gsl_vector *y = gsl_vector_calloc(n);
    gsl_vector *mvec = gsl_vector_calloc(n);
    int verbose = PrepStats(zmap,metric,cyl,k,y,mvec);
    double tss = gsl_stats_tss(y->data,1,y->size);
    if (tss < TINY) {
      gsl_vector_free(y);
      gsl_vector_free(mvec);
      continue;
    }
    
    if (x3bins) {
      gsl_vector *beta = gsl_vector_calloc(n3bins);
      gsl_vector *zval = gsl_vector_calloc(p3bins);
      if (LaminarMedian(y,mvec,numperm,seed,beta,zval) > 0) {
	WriteStats(cyl,k,rim,beta,bin3_coeff);
	WriteStats(cyl,k,rim,zval,bin3_perm);
      }
      gsl_vector_free(beta);
      gsl_vector_free(zval);
    }
    if (nbins > 0) {
      gsl_vector *beta = gsl_vector_calloc(nbins);
      if (LaminarNBins(y,mvec,beta,verbose) > 0) {
	WriteStats(cyl,k,rim,beta,nbin_coeff);
      }
      gsl_vector_free(beta);
    }
    if (xpeak) {
      gsl_vector *beta = gsl_vector_calloc(npeak);
      if (LaminarPeak(y,mvec,beta,verbose) > 0) {
	WriteStats(cyl,k,rim,beta,peak_coeff);
      }
      gsl_vector_free(beta);
    }

    if (xconvex) {
      gsl_vector *beta = gsl_vector_calloc(1);
      if (LaminarConvex(y,mvec,beta,verbose) > 0) {
	WriteStats(cyl,k,rim,beta,&convex_coeff);
      }
      gsl_vector_free(beta);
    }

    if (xlinear) {
      gsl_vector *beta = gsl_vector_calloc(1);
      if (LaminarLinear(y,mvec,beta) > 0) {
	WriteStats(cyl,k,rim,beta,&linear_coeff);
      }
      gsl_vector_free(beta);
    }

    if (xR2) {
      gsl_vector *beta = gsl_vector_calloc(nR2);
      if (LaminarFit(y,mvec,beta) > 0) {
	WriteStats(cyl,k,rim,beta,R2_coeff);
      }
      gsl_vector_free(beta);
    }
    
    /* max absolute z-value */
    double zmax = LaminarZMax(y);
    
    gsl_vector_free(y);
    gsl_vector_free(mvec);


#pragma omp critical
    {    
      for (i=0; i<cyl->addr[k]->size; i++) {
	size_t l = cyl->addr[k]->data[i];
	int b = gsl_matrix_int_get(cyl->xmap,l,0);
	int r = gsl_matrix_int_get(cyl->xmap,l,1);
	int c = gsl_matrix_int_get(cyl->xmap,l,2);
	VPixel(zabsimage,b,r,c,VFloat) += zmax;
	VPixel(wimage,b,r,c,VFloat) += 1.0;
      }
    }
  }

  
  /* normalize images */
  if (x3bins) {
    WNormalizeImage(bin3_coeff,(int)n3bins,wimage,rim);
    WNormalizeImage(bin3_perm,(int)p3bins,wimage,rim);
  }
  if (nbins > 0) {
    WNormalizeImage(nbin_coeff,(int)nbins,wimage,rim);
  }
  if (xpeak) {
    WNormalizeImage(peak_coeff,(int)npeak,wimage,rim);
  }
  if (xconvex) {
    WNormalizeImage(&convex_coeff,(int)1,wimage,rim);
  }
  if (xlinear) {
    WNormalizeImage(&linear_coeff,(int)1,wimage,rim);
  }
  if (xR2) {
    WNormalizeImage(R2_coeff,(int)nR2,wimage,rim);
  }
  WNormalizeImage(&zabsimage,(int)1,wimage,rim);
  return zabsimage;
}

int main (int argc, char **argv)
{
  static VString    in_filename="";
  static VString    out_filename="";
  static VString    metric_filename="";
  static VString    rim_filename="";
  static VString    mask_filename="";
  static VFloat     radius = 2.0;
  static VBoolean   equivol = FALSE;
  static VBoolean   x3bins = FALSE;
  static VShort     xnbins = 0;
  static VBoolean   xconvex = FALSE;
  static VBoolean   xlinear = FALSE;
  static VBoolean   xpeak = FALSE;
  static VBoolean   xR2 = FALSE;
  static VBoolean   xzabs = FALSE;
  static VLong      seed = 5555;
  static VShort     numperm = 1000;
  static VShort     nproc = 0;
  static VOptionDescRec options[] = {
    {"in", VStringRepn, 1, & in_filename, VRequiredOpt, NULL,"Input file" },
    {"out", VStringRepn, 1, & out_filename, VRequiredOpt, NULL,"Output file" },
    {"metric", VStringRepn,1,(VPointer) &metric_filename,VRequiredOpt,NULL,"metric image"},
    {"rim", VStringRepn,1,(VPointer) &rim_filename,VRequiredOpt,NULL,"rim image"},
    {"mask", VStringRepn,1,(VPointer) &mask_filename,VOptionalOpt,NULL,"mask image"},
    {"radius", VFloatRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Cylinder radius in mm"},
    {"equivol", VBooleanRepn,1,(VPointer) &equivol,VOptionalOpt,NULL,"Equivolume correction"},
    {"3bins", VBooleanRepn,1,(VPointer) &x3bins,VOptionalOpt,NULL,"Compute 3bins"},
    {"nbins", VShortRepn,1,(VPointer) &xnbins,VOptionalOpt,NULL,"Binning using N bins"},
    {"peak", VBooleanRepn,1,(VPointer) &xpeak,VOptionalOpt,NULL,"Compute peaks and troughs"},
    {"concave", VBooleanRepn,1,(VPointer) &xconvex,VOptionalOpt,NULL,"Compute concave up"},
    {"linear", VBooleanRepn,1,(VPointer) &xlinear,VOptionalOpt,NULL,"Compute linear fit"},
    {"R2", VBooleanRepn,1,(VPointer) &xR2,VOptionalOpt,NULL,"Compute model fits (R^2)"},
    {"maxabs", VBooleanRepn,1,(VPointer) &xzabs,VOptionalOpt,NULL,"Compute max abs"},
    {"nperm", VShortRepn, 1, & numperm, VOptionalOpt,NULL,"Number of permutations" },
    {"seed", VLongRepn, 1, & seed, VOptionalOpt, NULL,"Seed for random number generator" },
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"},
   };
  size_t i;
  char *prg=GetLipsiaName("vcylarim");
  fprintf (stderr, " %s\n", prg);

  
  /* parse command line */
  if (! VParseCommand (VNumber (options), options, & argc, argv)) {
    VReportUsage (argv[0], VNumber (options), options, NULL);
    exit (EXIT_FAILURE);
  }
  if (argc > 1) {
    VReportBadArgs (argc, argv);
    exit (EXIT_FAILURE);
  }
  if (radius < 0.0001) VError(" radius must be positive");
  if (!x3bins && xnbins==0 && !xpeak && !xlinear && !xconvex && ! xR2 && !xzabs)
    VError(" no output will be generated, all compute options are set to 'false' ");
  if (xnbins > 100) VError(" nbins > 100 does not make sense");
			    

  
  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  if (nproc > 0 && nproc < num_procs) num_procs = nproc;
  fprintf(stderr," using %d cores\n",(int)num_procs);
  omp_set_num_threads(num_procs);
#endif
  
  /* read zmap  */
  VImage zmap=NULL;
  VAttrList zlist = VReadAttrList(in_filename,0L,TRUE,FALSE);
  VImage ztmp = VReadImage(zlist);
  if (ztmp == NULL) VError(" no zmap found");
  VAttrList geolist = VGetGeoInfo(zlist);
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

  
  
  /* mask image */
  VFloat *px = VImageData(metric);
  VFloat *pz = VImageData(zmap);
  VUByte *pu = VImageData(rim);
  VImage mask=NULL;
  if (strlen(mask_filename) > 1) {
    VAttrList mlist = VReadAttrList(mask_filename,0L,TRUE,FALSE);
    if (mlist == NULL) VError(" error reading file %s",mask_filename);
    VImage mtmp = VReadImage(mlist);
    if (mtmp == NULL) VError(" err reading %s",mask_filename);
    if (VImageNPixels(mtmp) != VImageNPixels(zmap))
      VError(" inconsistent image dimensions: mask vs zmap");
    mask = Convert2Repn(mtmp,mask,VBitRepn);
    VDestroyImage(mtmp);
  
    VBit *pm = VImageData(mask);  /* apply mask, if present */
    for (i=0; i<VImageNPixels(rim); i++) {
      if (pm[i] == 0 || pu[i] == 0) { px[i] = 0; pu[i] = 0; pz[i] = 0; }
    }
  }

  /* exclude voxels not covered by zmap */
  for (i=0; i<VImageNPixels(zmap); i++) {
    if (fabs(pz[i]) < TINY) { pu[i] = 0; px[i] = 0; }
  }
  

  /* add rim pts to metric image in case they are not included */
  for (i=0; i<VImageNPixels(rim); i++) {
    if (pu[i] == 1) px[i] = 1;
    if (pu[i] == 2) px[i] = 0.0001;
  }



  /* alloc 3bins */
  VImage *bin3_coeff = NULL;
  VImage *bin3_perm = NULL;
  if (x3bins) {
    bin3_coeff = (VImage *)VCalloc(n3bins,sizeof(VImage));
    for (i=0; i<n3bins; i++) {
      bin3_coeff[i] = VCreateImageLike(zmap);
      VFillImage(bin3_coeff[i],VAllBands,0);
    }
    bin3_perm = (VImage *)VCalloc(p3bins,sizeof(VImage));
    for (i=0; i<p3bins; i++) {
      bin3_perm[i] = VCreateImageLike(zmap);
      VFillImage(bin3_perm[i],VAllBands,0);
    }
  }

  /* alloc nbins */
  VImage *nbin_coeff = NULL;
  if (xnbins > 0) {
    nbin_coeff = (VImage *)VCalloc((int)xnbins,sizeof(VImage));
    for (i=0; i<xnbins; i++) {
      nbin_coeff[i] = VCreateImageLike(zmap);
      VFillImage(nbin_coeff[i],VAllBands,0);
    }
  }
  
  /* alloc peak */
  VImage *peak_coeff = NULL;
  if (xpeak) {
    peak_coeff = (VImage *)VCalloc(npeak,sizeof(VImage));
    for (i=0; i<npeak; i++) {
      peak_coeff[i] = VCreateImageLike(zmap);
      VFillImage(peak_coeff[i],VAllBands,0);
    }
  }


  /* alloc convex */
  VImage convex_coeff = NULL;
  if (xconvex) {
    convex_coeff = VCreateImageLike(zmap);
    VFillImage(convex_coeff,VAllBands,0);
  }

  /* alloc linear ramp */
  VImage linear_coeff = NULL;
  if (xlinear) {
    linear_coeff = VCreateImageLike(zmap);
    VFillImage(linear_coeff,VAllBands,0);
  }

  /* alloc R2 */
  VImage *R2_coeff = NULL;
  if (xR2) {
    R2_coeff = (VImage *)VCalloc(nR2,sizeof(VImage));
    for (i=0; i<nR2; i++) {
      R2_coeff[i] = VCreateImageLike(zmap);
      VFillImage(R2_coeff[i],VAllBands,0);
    }
  }

  
  /* Main */
  VImage zcover = Cylarim(zmap,metric,rim,(double)radius,equivol,(size_t)numperm,(long)seed,
			  x3bins,(int)xnbins,xpeak,xconvex,xlinear,xR2,
			  bin3_coeff,bin3_perm,nbin_coeff,peak_coeff,convex_coeff,linear_coeff,R2_coeff);

  /* write output images */
  fprintf(stderr,"\n");
  strip_filename(out_filename);

  if (x3bins) {
    WriteToFile(bin3_coeff,geolist,n3bins,
		(float)radius,(VBoolean)equivol,(VLong)0,(VLong)0,
		"3bins_coeff",out_filename);
    WriteToFile(bin3_perm,geolist,p3bins,
		(float)radius,(VBoolean)equivol,(VLong)numperm,(VLong)seed,
		"3bins_zvals",out_filename);
  }
  if (xnbins > 0) {
    WriteToFile(nbin_coeff,geolist,(int)xnbins,
		(float)radius,(VBoolean)equivol,(VLong)0,(VLong)0,
		"nbins_coeff",out_filename);
  }
  if (xpeak) {
    WriteToFile(peak_coeff,geolist,npeak,
		(float)radius,(VBoolean)equivol,(VLong)0,(VLong)0,
		"peak",out_filename);
  }
  if (xconvex) {
    WriteToFile(&convex_coeff,geolist,1,
		(float)radius,(VBoolean)equivol,(VLong)0,(VLong)0,
		"convex",out_filename);
  }

  if (xlinear) {
    WriteToFile(&linear_coeff,geolist,1,
		(float)radius,(VBoolean)equivol,(VLong)0,(VLong)0,
		"linear",out_filename);
  }
  
  if (xR2) {
    WriteToFile(R2_coeff,geolist,nR2,
		(float)radius,(VBoolean)equivol,(VLong)0,(VLong)0,
		"R2",out_filename);
  }
  if (xzabs) {
    WriteToFile(&zcover,geolist,1,(float)radius,(VBoolean)equivol,(VLong)0,(VLong)0,"zabs",out_filename);
  }
    
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
