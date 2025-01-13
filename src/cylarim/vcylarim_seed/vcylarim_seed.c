/*
** Cylam:  Cylinder-based Algorithm for Laminar fMRI Data Analysis
**
** G.Lohmann, MPI-KYB, Oct 2024
*/


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>

#include "../cylutils/cyl.h"

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

extern VImage Convert2Repn(VImage src,VImage dest,VRepnKind repn);
extern Cylinders *VCylinder(VImage rim,VImage metric,double radius);
extern void HistEqualize(Cylinders *,VImage,VImage);
extern void GetResolution(VImage src,gsl_vector *reso);

void ApplySeedMask(VImage metric,VImage rim,int b0,int r0,int c0,int wn)
{
  int b,r,c,d,rad2=wn*wn;

  for (b=0; b<VImageNBands(rim); b++) {
    for (r=0; r<VImageNRows(rim); r++) {
      for (c=0; c<VImageNColumns(rim); c++) {
	d = (b0-b)*(b0-b) + (r0-r)*(r0-r) + (c-c0)*(c-c0);
	if (d > rad2) {
	  VPixel(rim,b,r,c,VUByte) = 0;
	  VPixel(metric,b,r,c,VFloat) = 0;
	}
      }
    }
  }
}


VImage VFindCyl(VImage zmap,VImage metric,VImage rim,double radius,VBoolean equivol,VBoolean single,
		int type,VString txt_filename,int b1,int r1,int c1,int b2,int r2,int c2)
{
  size_t i,j,k,l,n;
  double s=0,nx=0,d1=0,d2=0,x[3],xx[3],y[3];
  FILE *fp = NULL;

  int nslices = VImageNBands(rim);
  int nrows = VImageNRows(rim);
  int ncols = VImageNColumns(rim);
  if (c1 < 0 || c1 >= ncols) VError(" illegal column address (%d), must be < %d",c1,ncols);
  if (r1 < 0 || r1 >= nrows) VError(" illegal row address (%d), must be < %d",r1,nrows);
  if (b1 < 0 || b1 >= nslices) VError(" illegal slice address (%d), must be < %d",b1,nslices);
  if (b2 >= nslices) VError(" illegal slice address (%d), must be < %d",b2,nslices);
  if (r2 >= nslices) VError(" illegal row address (%d), must be < %d",r2,nrows);
  if (c2 >= nslices) VError(" illegal column address (%d), must be < %d",c2,ncols);

  /* mask enclosing seed point */
  gsl_vector *reso = gsl_vector_calloc(3);
  GetResolution(rim,reso);
  double xreso = gsl_vector_min(reso); 
  int wn = (int)((2.0*radius+5.0)/xreso+0.5);
  ApplySeedMask(metric,rim,b1,r1,c1,wn);
  if (b2 >= 0) ApplySeedMask(metric,rim,b2,r2,c2,wn);

  /* for (i=0; i<3; i++) reso->data[i] = 1; */
  x[0] = reso->data[0]*(double)c1;
  x[1] = reso->data[1]*(double)r1;
  x[2] = reso->data[2]*(double)b1;

  xx[0] = xx[1] = xx[2] = -1;
  if (b1 >= 0) {
    xx[0] = reso->data[0]*(double)c2;
    xx[1] = reso->data[1]*(double)r2;
    xx[2] = reso->data[2]*(double)b2;
  }
  
  /* get cylinder struct */
  Cylinders *cyl = VCylinder(rim,metric,radius);

  
  /* equivolume correction */
  VImage wimage = VCreateImageLike(zmap);
  VFillImage(wimage,VAllBands,0);
  if (equivol) {
    fprintf(stderr," equivolume...\n");
    HistEqualize(cyl,wimage,metric);
  }

  
  /* find cylinder ID */
  double *table = (double *)VCalloc(cyl->numcylinders,sizeof(double));
  int *match = (int *)VCalloc(cyl->numcylinders,sizeof(int));
  i=n=0;
  
  for (k=0; k<cyl->numcylinders; k++) {
    s=nx=0;
    for (l=0; l<cyl->addr[k]->size; l++) {
      j = cyl->addr[k]->data[l];
      int b = gsl_matrix_int_get(cyl->xmap,j,0);
      int r = gsl_matrix_int_get(cyl->xmap,j,1);
      int c = gsl_matrix_int_get(cyl->xmap,j,2);
      if (b==b1 && r==r1 && c==c1) match[i++] = k;

      y[0] = reso->data[0]*(double)c;
      y[1] = reso->data[1]*(double)r;
      y[2] = reso->data[2]*(double)b;
      
      d1 = (x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]);
      s += exp(-d1);
      
      if (b2 >= 0) {
	d2 = (xx[0]-y[0])*(xx[0]-y[0]) + (xx[1]-y[1])*(xx[1]-y[1]) + (xx[2]-y[2])*(xx[2]-y[2]);
	s += exp(-d2);
      }
      nx++;
    }
    if (nx > 0.1) table[k] = s/nx;
  }
  n = i;

  
  /* find cylinder that best fits the seed voxel */
  size_t k0 = 0;
  double zmax = -1;
  for (k=0; k<cyl->numcylinders; k++) {
    if (table[k] > zmax) {
      zmax = table[k];
      k0 = k;
    }
  }
  if (zmax < 0) VError(" no cylinder found");


  /* map to dest image */
  VImage dest = VCreateImageLike(zmap);
  VFillImage(dest,VAllBands,0);
  VFloat u=0,z=0,w=0;

  /* map a single cylinder that best matches the seed voxel */
  if (single) {
    if (strlen(txt_filename) > 0) fp = fopen(txt_filename,"w");
    for (l=0; l<cyl->addr[k0]->size; l++) {
      j = cyl->addr[k0]->data[l];
      int b = gsl_matrix_int_get(cyl->xmap,j,0);
      int r = gsl_matrix_int_get(cyl->xmap,j,1);
      int c = gsl_matrix_int_get(cyl->xmap,j,2);
      z = VPixel(zmap,b,r,c,VFloat);
      u = VPixel(metric,b,r,c,VFloat);

      if (fp != NULL) fprintf(fp," %f %f\n",u,z);  
      w=0;
      if (type == 0) w = u;
      if (type == 1) w = z;
      if (type == 2) w = 1;
      VPixel(dest,b,r,c,VFloat) = w;
    }
    if (fp != NULL) fclose(fp);      
  }

  /* map all cylinders that match the seed voxel */
  else {
    fprintf(stderr," number of cylinders found: %lu\n",n);
    VFillImage(wimage,VAllBands,0);
    for (i=0; i<n; i++) {
      k = match[i];
      for (l=0; l<cyl->addr[k]->size; l++) {
	j = cyl->addr[k]->data[l];
	int b = gsl_matrix_int_get(cyl->xmap,j,0);
	int r = gsl_matrix_int_get(cyl->xmap,j,1);
	int c = gsl_matrix_int_get(cyl->xmap,j,2);
	z = VPixel(zmap,b,r,c,VFloat);
	u = VPixel(metric,b,r,c,VFloat);
	w=0;
	if (type == 0) w = u;
	if (type == 1) w = z;
	if (type == 2) w = 1;
	VPixel(dest,b,r,c,VFloat) += w;
	VPixel(wimage,b,r,c,VFloat) += 1;
      }
    }
    if (type < 2) {
      VFloat *pw = VImageData(wimage);
      VFloat *px = VImageData(dest);
      for (i=0; i<VImageNPixels(dest); i++) {
	if (pw[i] > 0.001) { px[i] /= pw[i]; }
      }
    }
  }

  return dest;
}


typedef struct SpointStruct{
  VShort x;
  VShort y;
  VShort z;
} SPoint;

VDictEntry TypDict[] = {
  { "metric", 0, 0,0,0,0  },
  { "zmap", 1, 0,0,0,0  },
  { "cover", 2, 0,0,0,0  },
  { NULL, 0,0,0,0,0 }
};

int main (int argc, char **argv)
{
  static SPoint     addr1;
  static SPoint     addr2;
  static VString    metric_filename="";
  static VString    rim_filename="";
  static VString    txt_filename="";
  static VFloat     radius = 2.0;
  static VShort     type = 0;
  static VBoolean   equivol = FALSE;
  static VBoolean   single = TRUE;
  static VOptionDescRec options[] = {
    {"seed1",VShortRepn,3,(VPointer) &addr1,VRequiredOpt,NULL,"Seed 1 (x,y,z)"},
    {"seed2",VShortRepn,3,(VPointer) &addr2,VOptionalOpt,NULL,"Seed 2 (x,y,z)"},
    {"metric", VStringRepn,1,(VPointer) &metric_filename,VRequiredOpt,NULL,"metric image"},
    {"rim", VStringRepn,1,(VPointer) &rim_filename,VRequiredOpt,NULL,"rim image"},
    {"report", VStringRepn,1,(VPointer) &txt_filename,VOptionalOpt,NULL,"output txt-file"},
    {"radius", VFloatRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Cylinder radius in mm"},
    {"equivol", VBooleanRepn,1,(VPointer) &equivol,VOptionalOpt,NULL,"Equivolume correction"},
    {"one", VBooleanRepn,1,(VPointer) &single,VOptionalOpt,NULL,"Whether to show only one cylinder"},
    {"type", VShortRepn,1,(VPointer) &type,VOptionalOpt,TypDict,"Output type"},
   };
  VString in_filename=NULL;
  FILE *in_file, *out_file;
  size_t i;
  char *prg=GetLipsiaName("vcylarim_seed");
  fprintf (stderr, " %s\n", prg);

  
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);

  
  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  if (num_procs > 64) num_procs=64;
  fprintf(stderr,"using %d cores\n",(int)num_procs);
  omp_set_num_threads(num_procs);
#endif /* _OPENMP */

  

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


  /* add rim pts to metric image in case they are not included */
  VUByte *pu = VImageData(rim);
  VFloat *px = VImageData(metric);
  for (i=0; i<VImageNPixels(rim); i++) {
    if (pu[i] == 1) px[i] = 1;
    if (pu[i] == 2) px[i] = 0.0001;
  }

  
  /* voxel address  */
  int b1=-1,r1=-1,c1=-1;
  b1 = addr1.z;
  r1 = addr1.y;
  c1 = addr1.x;
  fprintf(stderr," seed voxel 1:  %d  %d  %d \n",c1,r1,b1);

  int b2=-1,r2=-1,c2=-1;
  if (b2 >= 0) {
    b2 = addr2.z;
    r2 = addr2.y;
    c2 = addr2.x;
    fprintf(stderr," seed voxel 2:  %d  %d  %d \n",c2,r2,b2);
  }

  /* main */
  VImage dest = VFindCyl(zmap,metric,rim,(double)radius,equivol,single,(int)type,txt_filename,b1,r1,c1,b2,r2,c2);
  

  /* output */
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,dest); 
  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
