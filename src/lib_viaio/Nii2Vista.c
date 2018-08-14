/*
** read nifti-1 to vista format
**
** G.Lohmann, MPI-KYB, 2016
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
#include <ctype.h>

#include <gsl/gsl_math.h>

#include <nifti/nifti1.h>
#include <nifti/nifti1_io.h>

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

#define TINY 1.0e-10
#define ABS(x) ((x) > 0 ? (x) : -(x))

extern int isinff(float x);
extern int isnanf(float x);

extern char *VReadGzippedData(char *filename,size_t *len);
extern char *VReadUnzippedData(char *filename,VBoolean nofail,size_t *size);
extern char *VReadDataContainer(char *filename,VBoolean nofail,size_t *size);
extern FILE *VReadInputFile (char *filename,VBoolean nofail);
extern int  CheckGzip(char *filename);

extern void VByteSwapNiftiHeader(nifti_1_header *hdr);
extern void VByteSwapData(char *data,size_t ndata,size_t nsize);
extern float strtof(const char *str, char **endptr);


/* read repetition time from description field instead of pixdim  */
float TRfromString(const char *str)
{
  float f_tr   = -1;    /* return value on failure */
  float factor = 1000;  /* expect seconds as default */
  char* c_tr;

  if( ( c_tr = strstr( str, "TR=" ) ) || ( c_tr = strstr( str, "TR =" ) ) ) {
      while( *c_tr++ != '=' )
    ;

      while( *c_tr && isspace( *c_tr ) )
	c_tr++;

      f_tr = strtof( c_tr, &c_tr );
      
      if( f_tr > 0.0 ) {   
	while( *c_tr && isspace( *c_tr ) )
	  c_tr++;
	
	if( strncasecmp( c_tr, "ms", 2 ) == 0 )
	  factor = 1;
      }
  }
  return f_tr * factor;
}



/* set dimension in geo info */
void VUpdateGeoinfo(VAttrList geolist,int dimtype,VLong tr)
{
  double *D = VGetGeoDim(geolist,NULL);
  D[0] = (double)dimtype;
  VSetGeoDim(geolist,D);

  double *E = VGetGeoPixdim(geolist,NULL);
  if (tr > 0) {
    E[4] = (double)tr;
    VSetGeoPixdim(geolist,E);
  }
} 


/*
 *  VIniImage
 *
 *  Allocate Image structure and fill voxel values with databuffer.
 *  Returns a pointer to the image if successful, zero otherwise.
 *  Useful for converting from 4D nifti to list of 3D vista images
 */
VImage VIniImage (int nbands, int nrows, int ncolumns, VRepnKind pixel_repn, char *databuffer)
{
  size_t row_size = ncolumns * VRepnSize (pixel_repn);
  size_t row_index_size = nbands * nrows * sizeof (char *);
  size_t band_index_size = nbands * sizeof (char **);
  size_t pixel_size=0;
  VImage image;
  int band=0, row=0;
  char *p = NULL;

  /* Check parameters: */
  if (nbands < 1) VError("VIniImage: Invalid number of bands: %d", (int) nbands);
  if (nrows < 1)  VError ("VIniImage: Invalid number of rows: %d", (int) nrows);
  if (ncolumns < 1) VError ("VIniImage: Invalid number of columns: %d",(int) ncolumns);


#define AlignUp(v, b) ((((v) + (b) - 1) / (b)) * (b))

  /* Initialize VImage data struct */ 
  pixel_size = VRepnSize (pixel_repn);
  p = VMalloc(AlignUp (sizeof (VImageRec) + row_index_size + band_index_size, pixel_size));

  image = (VImage) p;
  image->nbands = nbands;
  image->nrows  = nrows;
  image->ncolumns = ncolumns;
  image->flags    = VImageSingleAlloc;
  image->pixel_repn = pixel_repn;
  image->attributes = VCreateAttrList ();
  image->band_index = (VPointer **) (p += sizeof (VImageRec));
  image->row_index  = (VPointer *) (p += band_index_size);
  image->data = (VPointer) AlignUp((long) (databuffer+4), pixel_size);

  image->nframes = nbands;
  image->nviewpoints = image->ncolors = image->ncomponents = 1;

  /* Initialize the indices: */
  for (band = 0; band < nbands; band++)
    image->band_index[band] = image->row_index + band * nrows;
  for (row = 0, p = image->data; row < nbands * nrows; row++, p += row_size)
    image->row_index[row] = p;

  return image;

#undef AlignUp
}


#define xgetval(data,index,Datatype) \
{ \
  Datatype *k = (Datatype *)(&data[index]); \
  u = (float)(*k); \
}


float VGetValue(char *data,size_t index,int datatype)
{
  float u=0;
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
  if (isinff(u) || isnanf(u)) u = 0;
  return u;
}


void VCleanData(VImage src)
{
  int b,r,c;
  double u=0;
  for (b=0; b<VImageNBands(src); b++) {
    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {
	u = VGetPixel(src,b,r,c);
	if (gsl_isinf(u) || gsl_isnan(u)) {
	  u = 0;
	  VSetPixel(src,b,r,c,u);
	}
      }
    }
  }
}


/* get image statistics for re-scaling parameters */
void VDataStats(char *data,size_t ndata,size_t nsize,int datatype,float *xmin,float *xmax)
{
  size_t i,n;
  float u=0;
  float zmin = VRepnMaxValue(VFloatRepn);
  float zmax = VRepnMinValue(VFloatRepn);

  n=0;
  for (i=0; i<ndata; i+= nsize) {
    u = VGetValue(data,i,datatype);
    if (fabs(u) < TINY) continue;
    if (u < zmin) zmin = u;
    if (u > zmax) zmax = u;
    n++;
  }
  if (n < 1) VError(" no non-zero data points found");

  *xmin = (zmin+TINY);
  *xmax = (zmax-TINY);
}



/* list of 3D images */
void Nii2Vista3DList(char *data,size_t nsize,size_t nslices,size_t nrows,size_t ncols,size_t nt,
		     VRepnKind dst_repn,int datatype,float scl_slope,float scl_inter,
		     VString voxelstr,VLong tr,VAttrList out_list)
{
  int slice,row,col;
  float u=0;
  size_t i;
  size_t add=0;      
  size_t nrnc = nrows*ncols;
  size_t npix = nrnc*nslices;
  size_t ndata = nt * npix * nsize;  


  if (nsize == 1) add=4;
  if (nsize == 2) add=2;
  if (nsize == 4) add=1;

  VImage *dst = (VImage *) VCalloc(nt,sizeof(VImage));
  if (nt > 1) VAppendAttr(out_list,"nimages",NULL,VLongRepn,(VLong) nt);

  for (i=0; i<nt; i++) {

    /* rescale if needed */
    if (fabs(scl_slope) > 0 || fabs(scl_inter) > 0) {

      VImage tmp = VCreateImage(nslices,nrows,ncols,dst_repn);
      VSetAttr(VImageAttrList(tmp),"voxel",NULL,VStringRepn,voxelstr);
      if (tr > 0) VSetAttr(VImageAttrList(tmp),"repetition_time",NULL,VLongRepn,(VLong)tr);
      for (slice=0; slice<nslices; slice++) {
	for (row=0; row<nrows; row++) {
	  for (col=0; col<ncols; col++) {
	    const size_t src_index = (col + row*ncols + slice*nrnc + i*npix)*nsize;
	    if (src_index >= ndata || src_index < 0) continue;
	    u = VGetValue(data,src_index,datatype);
	    u = scl_slope*u + scl_inter;
	    VSetPixel(tmp,slice,row,col,(double)u);	    
	  }
	}
      }
      VAppendAttr(out_list,"image",NULL,VImageRepn,tmp);
    }

    /* otherwise just copy */
    else {
      const size_t index = (i*npix-add)*nsize;
      dst[i] = VIniImage(nslices,nrows,ncols,dst_repn,&data[index]);
      VSetAttr(VImageAttrList(dst[i]),"voxel",NULL,VStringRepn,voxelstr);
      if (tr > 0) VSetAttr(VImageAttrList(dst[i]),"repetition_time",NULL,VLongRepn,(VLong)tr);
      VAppendAttr(out_list,"image",NULL,VImageRepn,dst[i]);
      VCleanData(dst[i]);
    }
  }
}




/* 4D time series data */
void Nii2Vista4D(char *data,size_t nsize,size_t nslices,size_t nrows,size_t ncols,size_t nt,
		 VRepnKind dst_repn,int datatype,VBoolean do_scaling,float scl_slope,float scl_inter,
		 VString voxelstr,double *slicetime,VLong tr,VAttrList out_list)
{
  size_t slice,row,col,ti;
  size_t nrnc = nrows*ncols;
  size_t npix = nrnc*nslices;
  size_t ndata = nt * npix * nsize;
  size_t add=0;
  float u=0;

  /* rescale to 16bit integer if needed */
  float xmin=0,xmax=0,umin=0,umax=0;
  if (do_scaling && dst_repn != VShortRepn) {
    dst_repn = VShortRepn;
    VDataStats(data,ndata,nsize,datatype,&xmin,&xmax);
    fprintf(stderr," data range: [%f, %f]\n",xmin,xmax);
    umin = 0;
    umax = VRepnMaxValue(VShortRepn);
  }


  VImage *dst = (VImage *) VCalloc(nslices,sizeof(VImage));

  for (slice=0; slice<nslices; slice++) {
    dst[slice] = VCreateImage(nt,nrows,ncols,dst_repn);
    if (!dst[slice]) VError(" err allocating image");
    VFillImage(dst[slice],VAllBands,0);

    VSetAttr(VImageAttrList(dst[slice]),"voxel",NULL,VStringRepn,voxelstr);
    if (tr > 0) VSetAttr(VImageAttrList(dst[slice]),"repetition_time",NULL,VLongRepn,(VLong)tr);
    if (slicetime != NULL)
      VSetAttr(VImageAttrList(dst[slice]),"slice_time",NULL,VShortRepn,(VShort)slicetime[slice]);

    for (ti=0; ti<nt; ti++) {
      for (row=0; row<nrows; row++) {
	for (col=0; col<ncols; col++) {
	  const size_t src_index = (col + row*ncols + slice*nrnc + ti*npix + add)*nsize;
	  if (src_index >= ndata || src_index < 0) continue;

	  if (!do_scaling) {
	    u = VGetValue(data,src_index,datatype);
	    if (fabs(scl_slope) > 0 || fabs(scl_inter) > 0) {
	      u = scl_slope*u + scl_inter;
	    }
	    if (VPixelRepn(dst[slice]) == VShortRepn) {
	      VPixel(dst[slice],ti,row,col,VShort) = u;
	    }
	    else if (VPixelRepn(dst[slice]) == VFloatRepn) {
	      VPixel(dst[slice],ti,row,col,VFloat) = u;
	    }
	    else {
	      VSetPixel(dst[slice],ti,row,col,(double)u);
	    }
	  }
	  else {
	    u = VGetValue(data,src_index,datatype);
	    u = umax * (u-xmin)/(xmax-xmin);	    
	    if (u < umin) u = umin;
	    if (u > umax) u = umax;
	    VPixel(dst[slice],ti,row,col,VShort) = (VShort)u;
	  }
	}
      }
    }
    VCleanData(dst[slice]);
    VAppendAttr(out_list,"image",NULL,VImageRepn,dst[slice]);
  }
}



/* copy nifti header infos to geolist in vista header */
double *VGetNiftiHeader(VAttrList geolist,nifti_1_header hdr,VLong tr)
{
  /* if (hdr.sizeof_hdr != 348) VError(" incorrect header size %lu",hdr.sizeof_hdr); */

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


  /* override if TR info is in description field */
  float xtr = TRfromString(hdr.descrip);
  if (xtr > 0) {
    D[4] = (double)(xtr*tscale);
    fprintf(stderr," reading TR from description field, TR= %.4f ms\n",D[4]);
  }
  if (fabs(D[0]) < 0.0001) D[0] = 1.0;  /* if not specified, assume pixdim[0] = 1 */


  for (i=5; i<8; i++) D[i] = 0; 
  VAttrList dlist = VCreateAttrList();
  VBundle dbundle = VCreateBundle ("bundle",dlist,8*sizeof(float),(VPointer)D);
  VSetAttr(geolist,"pixdim",NULL,VBundleRepn,dbundle);


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

  VShort qform_code = hdr.qform_code;
  if (hdr.sform_code==0 && qform_code==0) qform_code = 1;    /* if both codes are unspecified, assume scanner coord */
  VSetAttr(geolist,"qform_code",NULL,VShortRepn,(VShort)qform_code);
  VSetAttr(geolist,"qform",NULL,VBundleRepn,qbundle);

  
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



VAttrList Nifti1_to_Vista(char *databuffer,VLong tr,VBoolean attrtype,VBoolean do_scaling,VBoolean *ok)
{

  /* read header */
  nifti_1_header hdr;
  memcpy(&hdr,databuffer,MIN_HEADER_SIZE);
  if ((strncmp(hdr.magic,"ni1\0",4) != 0) && (strncmp(hdr.magic,"n+1\0",4) != 0))
    VError(" not a nifti-1 file, magic number is %s",hdr.magic);


  int swap = NIFTI_NEEDS_SWAP(hdr);
  if (swap == 1) {
    VByteSwapNiftiHeader(&hdr);
  }

  /* get data type */
  VRepnKind dst_repn = DT_UNKNOWN;
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
    VError(" unknown data type %d",datatype);
  }

  /* apply scaling if needed */
  float scl_slope = (float)hdr.scl_slope;
  float scl_inter = (float)hdr.scl_inter;
  if (fabs(scl_slope-1.0) < 1.0e-6 && fabs(scl_inter) < 1.0e-6) goto skip;
  if (fabs(scl_slope) > 0 || fabs(scl_inter) > 0) {
    dst_repn = VFloatRepn;
    do_scaling = FALSE;  /* no automatic scaling */  
    fprintf(stderr," apply scaling, slope: %f, intercept= %f\n",hdr.scl_slope,hdr.scl_inter);
  }
 skip: ;
  /*
  dst_repn = VShortRepn;
  scl_slope = 1.0;
  scl_inter = 0;
  */



  /* number of values stored at each time point */
  if (hdr.dim[5] > 1) VError("data type not supported, dim[5]= %d\n",hdr.dim[5]);


  /* image size */
  /* size_t dimtype = (size_t)hdr.dim[0]; */
  size_t ncols   = (size_t)hdr.dim[1];
  size_t nrows   = (size_t)hdr.dim[2];
  size_t nslices = (size_t)hdr.dim[3];
  size_t nt      = (size_t)hdr.dim[4];


  /* fill data container */
  size_t bytesize = 8;
  if (dst_repn == VBitRepn) bytesize = 1;
  size_t nsize   = hdr.bitpix/bytesize;
  size_t npixels = nslices * nrows * ncols;
  size_t ndata   = nt * npixels * nsize;
  size_t vox_offset = (size_t)hdr.vox_offset;
  size_t startdata = MIN_HEADER_SIZE;
  if (vox_offset > 0) startdata = vox_offset;
  char *data = &databuffer[startdata];


  /* byte swap image data, if needed */
  if (swap == 1) {
    VByteSwapData(data,ndata,nsize);
  }


  /* repetition time (may be wrong in some cases) */
  char xyzt = hdr.xyzt_units;
  int tcode = XYZT_TO_TIME(xyzt);
  float factor = 1.0;
  float xtr=0;
  if (tcode == NIFTI_UNITS_MSEC) factor = 1.0;
  if (tcode == NIFTI_UNITS_SEC) factor = 1000.0;

  if (nt > 1)  {
    if (tr == 0) tr = (short)(factor*hdr.pixdim[4]);

    /* override if TR info is in description field */
    xtr = TRfromString(hdr.descrip);
    if (xtr > 0) {
      tr = (short) (xtr*factor);
      fprintf(stderr," reading TR from description field, TR= %ld ms\n",tr);
    }
    /* fprintf(stderr," nt=%ld,  TR= %ld milliseconds\n",nt,tr); */
    if (tr < 50) VError(" implausible TR (%d ms), use 'vnifti' for data conversion to specify TR on the command line",tr);
  }

  
  /* voxel reso */
  int blen=512;
  char *voxelstr = (char *) VCalloc(blen,sizeof(char));
  memset(voxelstr,0,blen);
  sprintf(voxelstr,"%f %f %f",hdr.pixdim[1],hdr.pixdim[2],hdr.pixdim[3]);


  /* geometry information */
  VAttrList geolist = VCreateAttrList();
  double *slicetime = VGetNiftiHeader(geolist,hdr,tr);


  /* read nii image into vista attrlist */
  VAttrList out_list = VCreateAttrList();
  VAppendAttr(out_list,"geoinfo",NULL,VAttrListRepn,geolist);


  (*ok) = FALSE;
  if (nt <= 1) {      /* output one 3D image */
    Nii2Vista3DList(data,nsize,nslices,nrows,ncols,nt,dst_repn,datatype,scl_slope,scl_inter,voxelstr,tr,out_list);
    VUpdateGeoinfo(geolist,(int)3,tr);
  }
  else if (attrtype == FALSE) {        /* output list of 3D images */
    Nii2Vista3DList(data,nsize,nslices,nrows,ncols,nt,dst_repn,datatype,scl_slope,scl_inter,voxelstr,tr,out_list);
    VUpdateGeoinfo(geolist,(int)3,tr);
  }
  else if (attrtype == TRUE && nt > 1) {   /* output one 4D image */
    Nii2Vista4D(data,nsize,nslices,nrows,ncols,nt,dst_repn,datatype,do_scaling,scl_slope,scl_inter,voxelstr,slicetime,tr,out_list);
    VUpdateGeoinfo(geolist,(int)4,tr);
    (*ok) = TRUE;
  }

  return out_list;
}
