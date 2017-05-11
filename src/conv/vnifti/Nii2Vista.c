/*
** read nifti-1 to vista format
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include <gsl/gsl_math.h>

#include <nifti/nifti1.h>
#include <nifti/nifti1_io.h>

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

#define TINY 1.0e-10
#define ABS(x) ((x) > 0 ? (x) : -(x))


extern int isnan(double);
extern int isinf(double);
extern void VByteSwapNiftiHeader(nifti_1_header *hdr);
extern void VByteSwapData(char *data,size_t ndata,size_t nsize);



#define xgetval(data,index,Datatype) \
{ \
  Datatype *k = (Datatype *)(&data[index]); \
  u = (double)(*k); \
}


double VGetValue(char *data,size_t index,int datatype)
{
  double u=0;
  switch(datatype) {
  case DT_BINARY:
    xgetval(data,index,VBit);
    break;
  case DT_UNSIGNED_CHAR:
    xgetval(data,index,VUByte);
    break;
  case DT_SIGNED_SHORT:
    xgetval(data,index,short);
    break;
  case DT_SIGNED_INT:
    xgetval(data,index,int);
    break;
  case DT_FLOAT:
    xgetval(data,index,float);
    break;
  case DT_DOUBLE:
    xgetval(data,index,double);
    break;
  case DT_INT8:
    xgetval(data,index,VSByte);
    break;
  case DT_UINT16:
    xgetval(data,index,unsigned short);
    break;
  case DT_UINT32:
    xgetval(data,index,unsigned int);
    break;
  case DT_INT64:
    xgetval(data,index,long);
    break;
  case DT_UINT64:
    xgetval(data,index,unsigned long);
    break;
  default:
    VError(" unknown datatype %d",datatype);
  }
  return u;
}


/* get image statistics for re-scaling parameters */
void VDataStats(char *data,size_t ndata,size_t nsize,int datatype,double *xmin,double *xmax)
{
  size_t i,n;
  double zmin = VRepnMaxValue(VDoubleRepn);
  double zmax = VRepnMinValue(VDoubleRepn);

  n=0;
  for (i=0; i<ndata; i+= nsize) {
    double u = VGetValue(data,i,datatype);
    if (fabs(u) < TINY) continue;
    if (u < zmin) zmin = u;
    if (u > zmax) zmax = u;
    n++;
  }
  if (n < 1) VError(" no non-zero data points found");

  *xmin = (zmin+TINY);
  *xmax = (zmax-TINY);
}


/* 3D data */
void Nii2Vista3D(char *data,size_t nsize,size_t nslices,size_t nrows,size_t ncols,
		 int datatype,VRepnKind dst_repn,VBoolean do_scaling,double xmin,double xmax,
		 VString voxelstr,VShort tr,VAttrList out_list)
{
  size_t slice,row,col;
  VImage image = VCreateImage(nslices,nrows,ncols,dst_repn);
  if (!image) VError(" err allocating output image");
  VFillImage(image,VAllBands,0);
  VSetAttr(VImageAttrList(image),"voxel",NULL,VStringRepn,voxelstr);
  if (tr > 0) VSetAttr(VImageAttrList(image),"repetition_time",NULL,VShortRepn,tr);
  double umin  = VPixelMinValue(image);
  double umax  = VPixelMaxValue(image);
  size_t nrnc = nrows*ncols;
  double scalefactor = umax;

  /* check if scaling needed */
  VBoolean noscale = FALSE;
  if (datatype == DT_FLOAT && dst_repn == VFloatRepn) noscale = TRUE;
  if (datatype == DT_SIGNED_SHORT && dst_repn == VShortRepn) noscale = TRUE;
  if (datatype == DT_UINT16 && dst_repn == VUShortRepn) noscale = TRUE;
  if (datatype == DT_SIGNED_INT && dst_repn == VIntegerRepn) noscale = TRUE;
  if (datatype == DT_UINT32 && dst_repn == VUIntegerRepn) noscale = TRUE;
  if (datatype == DT_UNSIGNED_CHAR && dst_repn == VUByteRepn) noscale = TRUE;

  /* no scaling... */
  /*  klappt nicht bei isinf, isnan
  if (noscale && do_scaling == FALSE) {
    fprintf(stderr," no scaling\n");
    VImageData(image) = data;
  }
  */
  if (noscale == TRUE && do_scaling == FALSE) {
    /* fprintf(stderr," no scaling\n"); */
    for (slice=0; slice<nslices; slice++) {
      for (row=0; row<nrows; row++) {
	for (col=0; col<ncols; col++) {
	  const size_t index = (col + row*ncols + slice*nrnc)*nsize;
	  double u = VGetValue(data,index,datatype);
	  if (gsl_isnan(u) || gsl_isinf(u)) u=0;
	  VSetPixel(image,slice,row,col,u);
	}
      }
    }
  }


  /* scaling... */
  else {
    /* fprintf(stderr," scaling,  scalefactor: %f\n",scalefactor); */
    for (slice=0; slice<nslices; slice++) {
      for (row=0; row<nrows; row++) {
	for (col=0; col<ncols; col++) {
	  const size_t index = (col + row*ncols + slice*nrnc)*nsize;

	  double u = VGetValue(data,index,datatype);
	  if (gsl_isnan(u) || gsl_isinf(u)) u = 0;
	  if (do_scaling) {
	    if (ABS(u) > TINY) u = scalefactor * (u-xmin)/(xmax-xmin);
	  }
	  if (u < umin) u = umin;
	  if (u > umax) u = umax;
	  VSetPixel(image,slice,row,col,u);
	}
      }
    }
  }

  /* check of isnan, isinf */
  VAppendAttr(out_list,"image",NULL,VImageRepn,image);
}


/* 4D time series data */
void Nii2Vista4D(char *data,size_t nsize,size_t nslices,size_t nrows,size_t ncols,size_t nt,
		 int datatype,VBoolean do_scaling,double xmin,double xmax,
		 VString voxelstr,double *slicetime,VShort tr,VAttrList out_list)
{
  size_t slice,row,col,ti;
  size_t nrnc = nrows*ncols;
  size_t npix = nrnc*nslices;

  double umin = 0;
  double umax = VRepnMaxValue(VShortRepn);

  VImage *dst = (VImage *) VCalloc(nslices,sizeof(VImage));

  for (slice=0; slice<nslices; slice++) {
    dst[slice] = VCreateImage(nt,nrows,ncols,VShortRepn);
    if (!dst[slice]) VError(" err allocating image");
    VFillImage(dst[slice],VAllBands,0);

    VSetAttr(VImageAttrList(dst[slice]),"voxel",NULL,VStringRepn,voxelstr);
    if (tr > 0) VSetAttr(VImageAttrList(dst[slice]),"repetition_time",NULL,VShortRepn,tr);
    if (slicetime != NULL)
      VSetAttr(VImageAttrList(dst[slice]),"slice_time",NULL,VShortRepn,(VShort)slicetime[slice]);

    for (ti=0; ti<nt; ti++) {
      for (row=0; row<nrows; row++) {
	for (col=0; col<ncols; col++) {
	  const size_t index = (col + row*ncols + slice*nrnc + ti*npix)*nsize;
	  double u = VGetValue(data,index,datatype);
	  if (do_scaling) {
	    if (ABS(u) > TINY) {
	      u = umax * (u-xmin)/(xmax-xmin);
	    }
	  }
	  if (u < umin) u = umin;
	  if (u > umax) u = umax;
	  VPixel(dst[slice],ti,row,col,VShort) = (VShort)u;
	}
      }
    }
    VAppendAttr(out_list,"image",NULL,VImageRepn,dst[slice]);
  }
}




/* copy nifti header infos to geolist in vista header */
double *VGetNiftiHeader(VAttrList geolist,nifti_1_header hdr,VLong tr)
{

  /* units in lipsia: mm and sec */
  char xyzt = hdr.xyzt_units;
  int spaceunits = XYZT_TO_SPACE(xyzt);
  int timeunits  = XYZT_TO_TIME(xyzt);
  double xscale  = 1.0;
  double tscale  = 1.0;
  if (spaceunits == NIFTI_UNITS_MICRON) xscale=1000.0;
  if (timeunits == NIFTI_UNITS_SEC) tscale=1000.0;

 
  /* dim info */
  VSetAttr(geolist,"dim_info",NULL,VShortRepn,(VShort)hdr.dim_info);


  /* dim */
  /* hdr.dim[0]==3 means dim=3,  hdr.dim[0]==4 means dim=4 (timeseries) */
  float *E = VCalloc(8,sizeof(float));
  VAttrList elist = VCreateAttrList();
  VBundle ebundle = VCreateBundle ("bundle",elist,8*sizeof(float),(VPointer)E);
  VSetAttr(geolist,"dim",NULL,VBundleRepn,ebundle);
  int i;
  for (i=0; i<8; i++) E[i] = hdr.dim[i];
  for (i=5; i<8; i++) E[i] = 0;



  /* pixdim */  
  float *D = VCalloc(8,sizeof(float));
  for (i=0; i<8; i++) D[i] = hdr.pixdim[i];
  for (i=1; i<=3; i++) D[i] *= xscale; 
  D[4] *= tscale;
  if (tr > 0) D[4] = (double)tr;
  for (i=5; i<8; i++) D[i] = 0; 
  VAttrList dlist = VCreateAttrList();
  VBundle dbundle = VCreateBundle ("bundle",dlist,8*sizeof(float),(VPointer)D);
  VSetAttr(geolist,"pixdim",NULL,VBundleRepn,dbundle);



  /* qform */
  size_t dim=6;
  float *Q = VCalloc(dim,sizeof(float));
  VAttrList qlist = VCreateAttrList();
  VBundle qbundle = VCreateBundle ("bundle",qlist,dim*sizeof(float),(VPointer)Q);
  Q[0] = hdr.quatern_b;
  Q[1] = hdr.quatern_c;
  Q[2] = hdr.quatern_d;
  Q[3] = hdr.qoffset_x;
  Q[4] = hdr.qoffset_y;
  Q[5] = hdr.qoffset_z;
  VSetAttr(geolist,"qform_code",NULL,VShortRepn,(VShort)hdr.qform_code);
  VSetAttr(geolist,"qform",NULL,VBundleRepn,qbundle);



  /* sform */
  VImage sform = VCreateImage(1,4,4,VFloatRepn);
  VFillImage (sform,VAllBands,0);
  int j;
  for (j=0; j<4; j++) {
    VPixel(sform,0,0,j,VFloat) = hdr.srow_x[j];
    VPixel(sform,0,1,j,VFloat) = hdr.srow_y[j];
    VPixel(sform,0,2,j,VFloat) = hdr.srow_z[j];
  }
  VSetAttr(geolist,"sform_code",NULL,VShortRepn,(VShort)hdr.sform_code);
  VSetAttr(geolist,"sform",NULL,VImageRepn,sform);


  /* MRI encoding directions */
  int freq_dim=0,phase_dim=0,slice_dim=0;
  if (hdr.dim_info != 0) {
    freq_dim  = DIM_INFO_TO_FREQ_DIM (hdr.dim_info );
    phase_dim = DIM_INFO_TO_PHASE_DIM(hdr.dim_info );
    slice_dim = DIM_INFO_TO_SLICE_DIM(hdr.dim_info );
    VSetAttr(geolist,"freq_dim",NULL,VShortRepn,(VShort)freq_dim);
    VSetAttr(geolist,"phase_dim",NULL,VShortRepn,(VShort)phase_dim);
    VSetAttr(geolist,"slice_dim",NULL,VShortRepn,(VShort)slice_dim);
  }
  if (slice_dim == 0) return NULL;

  /* slicetiming information */
  int slice_start = hdr.slice_start;
  int slice_end = hdr.slice_end;
  int slice_code = hdr.slice_code;
  double slice_duration = hdr.slice_duration;
  if (slice_code == 0 || slice_duration < 1.0e-10) return NULL;
  fprintf(stderr," slice duration: %.2f ms\n", slice_duration);
  VSetAttr(geolist,"slice_start",NULL,VShortRepn,(VShort)slice_start);
  VSetAttr(geolist,"slice_end",NULL,VShortRepn,(VShort)slice_end);
  VSetAttr(geolist,"slice_code",NULL,VShortRepn,(VShort)slice_code);
  VSetAttr(geolist,"slice_duration",NULL,VFloatRepn,(VFloat)slice_duration);

  int nslices = (int)E[3];
  double *slicetimes = (double *) VCalloc(nslices,sizeof(double));
  VGetSlicetimes(slice_start,slice_end,slice_code,slice_duration,slicetimes);
  return slicetimes;
}



VAttrList Nifti1_to_Vista(char *hdr_file, char *data_file,VLong tr,VBoolean do_scaling)
{
  nifti_1_header hdr;

  /* open and read header */
  FILE *fp = fopen(hdr_file,"r");
  if (fp == NULL) {
    fprintf(stderr, "\nError opening header file %s\n",hdr_file);
    exit(1);
  }
 
  size_t ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  if (ret != 1) {
    fprintf(stderr, "\nError reading header file %s\n",hdr_file);
    exit(1);
  }
  fclose(fp);
  if ((strncmp(hdr.magic,"ni1\0",4) != 0) && (strncmp(hdr.magic,"n+1\0",4) != 0))
    VError(" not a nifti-1 file, magic number is %s",hdr.magic);


  int swap = NIFTI_NEEDS_SWAP(hdr);
  if (swap == 1) {
    VByteSwapNiftiHeader(&hdr);
  }

  /* get data type */
  VRepnKind dst_repn=VShortRepn;
  int datatype = (int)hdr.datatype;

  switch(datatype) {
  case DT_UNKNOWN:
    VError(" unknown data type");
    break;
  case DT_BINARY:
    dst_repn = VBitRepn;
    break;
  case DT_UNSIGNED_CHAR:
    dst_repn = VUByteRepn;
    break;
  case DT_SIGNED_SHORT:
    dst_repn = VShortRepn;
    break;
  case DT_SIGNED_INT:
    dst_repn = VIntegerRepn;
    break;
  case DT_FLOAT:
    dst_repn = VFloatRepn;
    break;
  case DT_DOUBLE:
    dst_repn = VDoubleRepn;
    break;
  case DT_INT8:
    dst_repn = VSByteRepn;
    break;
  case DT_UINT16:
    dst_repn = VUShortRepn;
    break;
  case DT_UINT32:
    dst_repn = VUIntegerRepn;
    break;
  case DT_INT64:
    dst_repn = VLongRepn;
    break;
  case DT_UINT64:
    dst_repn = VULongRepn;
    break;
  default:
    VError(" unknown repn %d",datatype);
  }



  /* open the datafile, jump to data offset */
  fp = fopen(data_file,"r");
  if (fp == NULL) {
    VError("Error opening data file %s",data_file);
  }
  ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);
  if (ret != 0) {
    VError("Error doing fseek in data file %s",(long)(hdr.vox_offset), data_file);
  }

  /* number of values stored at each time point */
  if (hdr.dim[5] > 1) VError("data type not supported, dim[5]= %d\n",hdr.dim[5]);


  /* image size */
  size_t dimtype = (size_t)hdr.dim[0];
  size_t ncols   = (size_t)hdr.dim[1];
  size_t nrows   = (size_t)hdr.dim[2];
  size_t nslices = (size_t)hdr.dim[3];
  size_t nt      = (size_t)hdr.dim[4];

  if (nt == 1) {
    hdr.dim[0] = 3;
    dimtype = 3;
  }
  if (nt > 1) dimtype = 4;
  fprintf(stderr," dim: %ldD,  image matrix= %ld %ld %ld,  nt= %ld\n",dimtype,ncols,nrows,nslices,nt);

  if (dimtype == 4) {
    dst_repn = VShortRepn; 
    if (hdr.datatype != DT_SIGNED_SHORT) do_scaling = TRUE;
    VWarning(" functional data will be rescaled to 16 bit representation");
  }

  /* fill data container */
  size_t bytesize = 8;
  if (dst_repn == VBitRepn) bytesize = 1;
  size_t nsize   = hdr.bitpix/bytesize;
  size_t npixels = nslices * nrows * ncols;
  size_t ndata   = nt * npixels * nsize;
  char *data = (char *) VMalloc(ndata);
  if (!data) VError(" error allocating data container");
  ret = fread(data,1,ndata,fp);
  if (ret != ndata) {
    fprintf(stderr," fread error, ret= %ld:\n",ret);
    VError("Error reading from %s",data_file);
  }
  fclose(fp);


  /* byte swap image data, if needed */
  if (swap == 1) {
    VByteSwapData(data,ndata,nsize);
  }


  /* functional data must be VShortRepn, rescale if needed */
  double xmin=0,xmax=0;
  if (do_scaling) {
    /* fprintf(stderr," scaling...\n"); */
    size_t npixels = nslices * nrows * ncols;
    size_t ndata   = nt * npixels * nsize;
    VDataStats(data,ndata,nsize,datatype,&xmin,&xmax);
    fprintf(stderr," data range: [%f, %f]\n",xmin,xmax);
  }



  /* repetition time (may be wrong in some cases) */
  char xyzt = hdr.xyzt_units;
  int tcode = XYZT_TO_TIME(xyzt);
  float factor = 1.0;
  if (tcode == NIFTI_UNITS_MSEC) factor = 1.0;
  if (tcode == NIFTI_UNITS_SEC) factor = 1000.0;
  if (nt > 1)  {
    if (tr == 0) tr = (short)(factor*hdr.pixdim[4]);
    fprintf(stderr," nt=%ld,  TR= %ld\n",nt,tr);
    if (tr < 1) VWarning(" implausible TR (%d ms), use parameter '-tr' to set correct TR",tr);
  }

  
  /* voxel reso */
  int buflen=512;
  char *voxelstr = (char *) VCalloc(buflen,sizeof(char));
  memset(voxelstr,0,buflen);
  sprintf(voxelstr,"%f %f %f",hdr.pixdim[1],hdr.pixdim[2],hdr.pixdim[3]);
  fprintf(stderr," voxel: %f %f %f\n",hdr.pixdim[1],hdr.pixdim[2],hdr.pixdim[3]);


  /* geometry information */
  VAttrList geolist = VCreateAttrList();
  double *slicetime = VGetNiftiHeader(geolist,hdr,tr);
  

  /* read nii image into vista attrlist */
  VAttrList out_list = VCreateAttrList();
  VAppendAttr(out_list,"geoinfo",NULL,VAttrListRepn,geolist);

  if (nt > 1 || dimtype == 4) {
    Nii2Vista4D(data,nsize,nslices,nrows,ncols,nt,datatype,do_scaling,xmin,xmax,voxelstr,slicetime,(VShort)tr,out_list);
  }
  else {
    Nii2Vista3D(data,nsize,nslices,nrows,ncols,datatype,dst_repn,do_scaling,xmin,xmax,voxelstr,(VShort)tr,out_list);
  }

  return out_list;
}
