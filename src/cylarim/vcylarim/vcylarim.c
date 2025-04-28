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
extern int LaminarGLM(gsl_vector *,gsl_vector *,size_t,long,gsl_vector *,gsl_vector *);
extern int LaminarMean(gsl_vector *,gsl_vector *,size_t,long,gsl_vector *,gsl_vector *);
extern int LaminarTrimmedMean(gsl_vector *,gsl_vector *,size_t,long,gsl_vector *,gsl_vector *);
extern int LaminarMedian(gsl_vector *,gsl_vector *,size_t,long,gsl_vector *,gsl_vector *);
extern VImage VCylCover(VImage,VImage,VImage,Cylinders *cyl);

    
void XWriteOutput(VImage image,VAttrList geolist,char *filename)
{
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,image);
  FILE *fp = fopen(filename,"w");
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
}


/* fill vectors as prep for laminat statistics */
int PrepStats(VImage zmap,VImage metric,Cylinders *cyl,size_t cid,gsl_vector *y,gsl_vector *mvec,size_t n)
{
  size_t i,k,nzero=0;
  int b,r,c;
  for (i=0; i<n; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    mvec->data[i] = VPixel(metric,b,r,c,VFloat);
    y->data[i] = VPixel(zmap,b,r,c,VFloat);
    if (fabs(y->data[i]) < TINY) nzero++;
  }
  if (nzero > n-1) return -1;  /* z-values are all zero */
  return 1;
}



VImage Cylarim(VImage zmap,VImage metric,VImage rim,double radius,
	       VBoolean equivol,int type,size_t numperm,long seed,VImage *betaimage,VImage *zvalimage)
{
  size_t i,j,k;

  
  /* create cylinder data struct */
  Cylinders *cyl = VCylinder(rim,metric,radius);  

  
  /* equivolume correction */
  VImage wimage = VCreateImageLike(zmap);
  VFillImage(wimage,VAllBands,0);
  if (equivol) {
    fprintf(stderr," Equivolume...\n");
    HistEqualize(cyl,wimage,metric,rim);
  }
  
  /* LayerGLM */  
  size_t progress=0;
  size_t np = (size_t)((float)cyl->numcylinders/100.0);
  if (np < 1) np = 1;
  int nbeta=3,nzval=3;
  VFillImage(wimage,VAllBands,0);
  
#pragma omp parallel for shared(progress)
  for (k=0; k<cyl->numcylinders; k++) {
    size_t i,j;
    if (k%np==0) fprintf(stderr," LaminarStats:  %7.3f\r",(float)(progress)/(float)cyl->numcylinders);
    progress++;

    size_t n = cyl->addr[k]->size;
    if (n < 20) continue; /* cylinder must be big enough for sufficient stats */
    int rtcode = 0;
    gsl_vector *beta = gsl_vector_calloc(nbeta);
    gsl_vector *zval = gsl_vector_calloc(nzval);
    gsl_vector *y = gsl_vector_calloc(n);
    gsl_vector *mvec = gsl_vector_calloc(n);
    if (PrepStats(zmap,metric,cyl,k,y,mvec,n) < 0) goto skip;
    
    switch (type) {
    case 0:
      rtcode = LaminarMedian(y,mvec,numperm,seed,beta,zval);
      break;
    case 1:
      rtcode = LaminarMean(y,mvec,numperm,seed,beta,zval);
      break;
    case 2:
      rtcode = LaminarGLM(y,mvec,numperm,seed,beta,zval);
      break;
    default:
      VError(" unknown type %d",type);
    }
    if (rtcode < 0) goto skip;
       
    for (i=0; i<cyl->addr[k]->size; i++) {
      size_t l = cyl->addr[k]->data[i];
      int b = gsl_matrix_int_get(cyl->xmap,l,0);
      int r = gsl_matrix_int_get(cyl->xmap,l,1);
      int c = gsl_matrix_int_get(cyl->xmap,l,2);
      if (VPixel(rim,b,r,c,VUByte) != 3) continue;

#pragma omp critical
      {
	VPixel(wimage,b,r,c,VFloat) += 1.0;
	
	for (j=0; j<nbeta; j++) {
	  VPixel(betaimage[j],b,r,c,VFloat) += beta->data[j];
	}
	
	for (j=0; j<nzval; j++) {
	  VPixel(zvalimage[j],b,r,c,VFloat) += zval->data[j];
	}
      }
    }
  skip: ;
    gsl_vector_free(y);
    gsl_vector_free(mvec);
    gsl_vector_free(beta);
    gsl_vector_free(zval);
  }

  
  /* normalize betaimages */
  VFloat *pw = VImageData(wimage);
  VFloat *pb = NULL;
  float wmin = 0.01;
  for (j=0; j<nbeta; j++) {
    pb = VImageData(betaimage[j]);
    for (i=0; i<VImageNPixels(wimage); i++) {
      if (pw[i] > wmin) pb[i] /= pw[i];
      else pb[i] = 0;
    }
  }

  /* normalize zval images */
  for (j=0; j<nzval; j++) {
    pb = VImageData(zvalimage[j]);
    for (i=0; i<VImageNPixels(wimage); i++) {
      if (pw[i] > wmin) pb[i] /= pw[i];
      else pb[i] = 0;
    }
  }

  VImage zcover = VCylCover(zmap,rim,wimage,cyl);
  return zcover;
}


VDictEntry TypDict[] = {
  { "median", 0, 0,0,0,0  },
  { "mean", 1, 0,0,0,0  },
  { "glm", 2, 0,0,0,0  },
  { NULL, 0,0,0,0,0 }
};

int main (int argc, char **argv)
{
  static VString    metric_filename="";
  static VString    rim_filename="";
  static VString    mask_filename="";
  static VFloat     radius = 2.0;
  static VBoolean   equivol = FALSE;
  static VShort     type = 0;
  static VLong      seed = 5555;
  static VShort     numperm = 1000;
  static VShort     nproc = 0;
  static VOptionDescRec options[] = {
    {"metric", VStringRepn,1,(VPointer) &metric_filename,VRequiredOpt,NULL,"metric image"},
    {"rim", VStringRepn,1,(VPointer) &rim_filename,VRequiredOpt,NULL,"rim image"},
    {"mask", VStringRepn,1,(VPointer) &mask_filename,VOptionalOpt,NULL,"mask image"},
    {"radius", VFloatRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Cylinder radius in mm"},
    {"equivol", VBooleanRepn,1,(VPointer) &equivol,VOptionalOpt,NULL,"Equivolume correction"},
    {"type", VShortRepn, 1, & type, VOptionalOpt,TypDict,"Type of model" },
    {"nperm", VShortRepn, 1, & numperm, VOptionalOpt,NULL,"Number of permutations" },
    {"seed", VLongRepn, 1, & seed, VOptionalOpt, NULL,"Seed for random number generator" },
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"},
   };
  VString in_filename=NULL;
  FILE *in_file, *out_file;
  size_t i;
  char *prg=GetLipsiaName("vcylarim");
  fprintf (stderr, " %s\n", prg);

  
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);
  if (radius < 0.0001) VError(" radius must be positive");
  if (type < 0 || type > 2) VError(" unknown type %d",type);
  if (seed < 0) VError(" seed must be non-negative");
  
  
  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  if (nproc > 0 && nproc < num_procs) num_procs = nproc;
  fprintf(stderr," using %d cores\n",(int)num_procs);
  omp_set_num_threads(num_procs);
#endif
  
  clock_t xstart = clock();

  fprintf(stderr," type:  %s\n",TypDict[type].keyword);

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
  

  /* add rim pts to metric image in case they are not included */
  for (i=0; i<VImageNPixels(rim); i++) {
    if (pu[i] == 1) px[i] = 1;
    if (pu[i] == 2) px[i] = 0.0001;
  }


  /* alloc output images */
  int nbeta=3,nzval=3;
  VImage *betaimage = (VImage *)VCalloc(nbeta,sizeof(VImage));
  for (i=0; i<nbeta; i++) {
    betaimage[i] = VCreateImageLike(zmap);
    VFillImage(betaimage[i],VAllBands,0);
  }
  VImage *zvalimage = (VImage *)VCalloc(nzval,sizeof(VImage));
  for (i=0; i<nzval; i++) {
    zvalimage[i] = VCreateImageLike(zmap);
    VFillImage(zvalimage[i],VAllBands,0);
  }
 
  
  /* Main */
  VImage zcover = Cylarim(zmap,metric,rim,(double)radius,equivol,(int)type,(size_t)numperm,(long)seed,betaimage,zvalimage);
  
  clock_t xend = clock();
  double time_spent = (double)(xend - xstart) / CLOCKS_PER_SEC;
  printf(" CPU time used: %.2f secs,  %.2f min,  %.2f std\n", time_spent,time_spent/60.0,time_spent/3600.0);

  
  /* write output */
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  char cbuf[64];
  sprintf(cbuf, "%f",radius);
  VAppendAttr(out_list,"type",NULL,VStringRepn,(VString)TypDict[type].keyword);
  VAppendAttr(out_list,"radius",NULL,VStringRepn,(VString)cbuf);
  if (equivol==TRUE) VAppendAttr(out_list,"equivolume",NULL,VStringRepn,(VString)"true");
  else VAppendAttr(out_list,"equivolume",NULL,VStringRepn,(VString)"false");
  VAppendAttr(out_list,"nbeta",NULL,VShortRepn,(VShort)nbeta);
  VAppendAttr(out_list,"seed",NULL,VLongRepn,(VLong)seed);
  for (i=0; i<nbeta; i++) VAppendAttr(out_list,"beta",NULL,VImageRepn,betaimage[i]);
  for (i=0; i<nzval; i++) VAppendAttr(out_list,"zvalimage",NULL,VImageRepn,zvalimage[i]);
  if (zcover != NULL) VAppendAttr(out_list,"cover",NULL,VImageRepn,zcover);

  if (! VWriteFile (out_file, out_list)) exit (1);

  
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
