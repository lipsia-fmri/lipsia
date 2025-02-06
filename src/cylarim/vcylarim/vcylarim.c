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
#include <gsl/gsl_errno.h>

#include "../cylutils/cyl.h"

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

extern void GetResolution(VImage src,gsl_vector *reso);
extern VImage Convert2Repn(VImage src,VImage dest,VRepnKind repn);
extern Cylinders *VCylinder(VImage rim,VImage metric,double radius,VBoolean,double);
extern void HistEqualize(Cylinders *,VImage,VImage,VImage);
extern double PermGLM(VImage,VImage,Cylinders *,size_t,size_t,double,gsl_vector *,gsl_vector *);
extern double LayerMean(VImage,VImage,Cylinders *,size_t,double,gsl_vector *,gsl_vector *);


void XWriteOutput(VImage image,VAttrList geolist,char *filename)
{
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,image);
  FILE *fp = fopen(filename,"w");
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
}


void Cylarim(VImage zmap,VImage metric,VImage rim,double radius,size_t numperm,
	     VBoolean equivol,VImage *betaimage,VImage *zvalimage,VBoolean addrim,double voxel_scale)
{
  size_t i,j,k;


  /* create cylinder data struct */
  Cylinders *cyl = VCylinder(rim,metric,radius,addrim,voxel_scale);  

  
  /* equivolume correction */
  VImage wimage = VCreateImageLike(zmap);
  VFillImage(wimage,VAllBands,0);
  if (equivol) {
    fprintf(stderr," Equivolume...\n");
    HistEqualize(cyl,wimage,metric,rim);
  }


  /* LayerGLM */  
  VFillImage(wimage,VAllBands,0);

  
  size_t progress=0;
  size_t np = (size_t)((float)cyl->numcylinders/100.0);
  if (np < 1) np = 1;
  int nbeta=4,nzval=3;

  
#pragma omp parallel for shared(progress)
  for (k=0; k<cyl->numcylinders; k++) {
    size_t i,j;
    if (k%np==0) fprintf(stderr," beta:  %7.3f\r",(float)(progress)/(float)cyl->numcylinders);
    progress++;

    gsl_vector *beta = gsl_vector_calloc(nbeta);
    gsl_vector *zval = gsl_vector_calloc(nzval);
    double w = PermGLM(zmap,metric,cyl,k,numperm,voxel_scale,beta,zval);
    
    for (i=0; i<cyl->addr[k]->size; i++) {
      size_t l = cyl->addr[k]->data[i];
      int b = gsl_matrix_int_get(cyl->xmap,l,0);
      int r = gsl_matrix_int_get(cyl->xmap,l,1);
      int c = gsl_matrix_int_get(cyl->xmap,l,2);

#pragma omp critical
      {
	VPixel(wimage,b,r,c,VFloat) += w;
	for (j=0; j<nbeta; j++) {
	  VPixel(betaimage[j],b,r,c,VFloat) += w*beta->data[j];
	}
	for (j=0; j<nzval; j++) {
	  VPixel(zvalimage[j],b,r,c,VFloat) += w*zval->data[j];
	}
      }
    }
    gsl_vector_free(beta);
    gsl_vector_free(zval);
  }
  fprintf(stderr,"                       \n");

  
  /* normalize betaimages */
  VFloat *pb = NULL;
  VFloat *pw = VImageData(wimage);
  for (j=0; j<nbeta; j++) {
    pb = VImageData(betaimage[j]);
    for (i=0; i<VImageNPixels(wimage); i++) {
      if (pw[i] > 0.001) pb[i] /= pw[i];
      else pb[i] = 0;
    }
  }


  /* normalize zval images */
  for (j=0; j<nzval; j++) {
    pb = VImageData(zvalimage[j]);
    for (i=0; i<VImageNPixels(wimage); i++) {
      if (pw[i] > 0.001) pb[i] /= pw[i];
      else pb[i] = 0;
    }
  }
}

int main (int argc, char **argv)
{
  static VString    metric_filename="";
  static VString    rim_filename="";
  static VString    mask_filename="";
  static VFloat     radius = 2.0;
  static VShort     numperm = 1;
  static VBoolean   equivol = FALSE;
  static VBoolean   addrim = TRUE;
  static VShort     nproc = 0;
  static VOptionDescRec options[] = {
    {"metric", VStringRepn,1,(VPointer) &metric_filename,VRequiredOpt,NULL,"metric image"},
    {"rim", VStringRepn,1,(VPointer) &rim_filename,VRequiredOpt,NULL,"rim image"},
    {"mask", VStringRepn,1,(VPointer) &mask_filename,VOptionalOpt,NULL,"mask image"},
    {"radius", VFloatRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Cylinder radius in mm"},
    {"perm", VShortRepn,1,(VPointer) &numperm,VOptionalOpt,NULL,"Number of permutations"},
    {"equivol", VBooleanRepn,1,(VPointer) &equivol,VOptionalOpt,NULL,"Equivolume correction"},
    {"addrim", VBooleanRepn,1,(VPointer) &addrim,VOptionalOpt,NULL,"Include rim"},
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"},
   };
  VString in_filename=NULL;
  FILE *in_file, *out_file;
  size_t i;
  char *prg=GetLipsiaName("vcylarim");
  fprintf (stderr, " %s\n", prg);

  
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);
  if (numperm < 1) VError(" number of permutations must be positive");
  if (radius < 0.00001) VError(" radius must be positive");


  
  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  if (nproc > 0 && nproc < num_procs) num_procs = nproc;
  if (num_procs > 64) num_procs=64;
  fprintf(stderr," using %d cores\n",(int)num_procs);
  omp_set_num_threads(num_procs);
#endif
  

  /* read zmap  */
  VImage zmap=NULL;
  VAttrList zlist = VReadAttrList(in_filename,0L,TRUE,FALSE);
  VImage ztmp = VReadImage(zlist);
  if (ztmp == NULL) VError(" no zmap found");
  if (VPixelRepn(ztmp) != VFloatRepn && VPixelRepn(ztmp) != VDoubleRepn)
    VError(" zmap image should be 4-byte or 8-byte float");
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
  int nbeta=4,nzval=3;
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
  double voxel_scale=1.0;
 
  
  /* Main */
  Cylarim(zmap,metric,rim,(double)radius,(size_t)numperm,equivol,betaimage,zvalimage,addrim,voxel_scale);

  
  /* write output */
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  char cbuf[64];
  sprintf(cbuf, "%f",radius);
  VAppendAttr(out_list,"radius",NULL,VStringRepn,(VString)cbuf);
  if (equivol==TRUE) VAppendAttr(out_list,"equivolume",NULL,VStringRepn,(VString)"true");
  else VAppendAttr(out_list,"equivolume",NULL,VStringRepn,(VString)"false");
  VAppendAttr(out_list,"nbeta",NULL,VShortRepn,(VShort)nbeta);
  for (i=0; i<nbeta; i++) VAppendAttr(out_list,"beta",NULL,VImageRepn,betaimage[i]);
  for (i=0; i<nzval; i++) VAppendAttr(out_list,"zvalimage",NULL,VImageRepn,zvalimage[i]);

  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
