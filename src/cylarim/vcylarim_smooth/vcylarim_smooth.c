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
extern int LaminarSmooth(VImage zmap,VImage metric,Cylinders *cyl,size_t cid,VImage dest,VImage ndest);
      
void XWriteOutput(VImage image,VAttrList geolist,char *filename)
{
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,image);
  FILE *fp = fopen(filename,"w");
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
}


VImage Cylarim_Smooth(VImage zmap,VImage metric,VImage rim,double radius,VBoolean equivol)
{
  size_t i,k;

  
  /* create cylinder data struct */
  Cylinders *cyl = VCylinder(rim,metric,radius);

  
  /* equivolume correction */
  VImage wimage = VCreateImageLike(zmap);
  VFillImage(wimage,VAllBands,0);
  if (equivol) {
    fprintf(stderr," Equivolume...\n");
    HistEqualize(cyl,wimage,metric,rim);
  }

  VImage dest = VCreateImageLike(zmap);
  VFillImage(dest,VAllBands,0);
  VFillImage(wimage,VAllBands,0);
  
  for (k=0; k<cyl->numcylinders; k++) {
    LaminarSmooth(zmap,metric,cyl,k,dest,wimage);
  }

  VFloat *pz = VImageData(zmap);
  VFloat *pd = VImageData(dest);
  VFloat *pw = VImageData(wimage);
  for (i=0; i<VImageNPixels(dest); i++) {
    if (pw[i] > 0.001) pd[i] /= pw[i];
    else pd[i] = pz[i];
  }
  return dest;
}


int main (int argc, char **argv)
{
  static VString    metric_filename="";
  static VString    rim_filename="";
  static VString    mask_filename="";
  static VFloat     radius = 1.0;
  static VBoolean   equivol = FALSE;
  static VShort     nproc = 0;
  static VOptionDescRec options[] = {
    {"metric", VStringRepn,1,(VPointer) &metric_filename,VRequiredOpt,NULL,"metric image"},
    {"rim", VStringRepn,1,(VPointer) &rim_filename,VRequiredOpt,NULL,"rim image"},
    {"mask", VStringRepn,1,(VPointer) &mask_filename,VOptionalOpt,NULL,"mask image"},
    {"radius", VFloatRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Cylinder radius in mm"},
    {"equivol", VBooleanRepn,1,(VPointer) &equivol,VOptionalOpt,NULL,"Equivolume correction"},
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"},
   };
  VString in_filename=NULL;
  FILE *in_file, *out_file;
  size_t i;
  char *prg=GetLipsiaName("vcylarim_smooth");
  fprintf (stderr, " %s\n", prg);

  
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);
  if (radius < 0.0001) VError(" radius must be positive");
  
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


  /* Main */
  VImage dest = Cylarim_Smooth(zmap,metric,rim,(double)radius,equivol);
  
  /* write output */
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  char cbuf[64];
  sprintf(cbuf, "%f",radius);
  VAppendAttr(out_list,"radius",NULL,VStringRepn,(VString)cbuf);
  if (equivol==TRUE) VAppendAttr(out_list,"equivolume",NULL,VStringRepn,(VString)"true");
  else VAppendAttr(out_list,"equivolume",NULL,VStringRepn,(VString)"false");
  VAppendAttr(out_list,"zmap_smoothed",NULL,VImageRepn,dest);
  if (! VWriteFile (out_file, out_list)) exit (1);

  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
