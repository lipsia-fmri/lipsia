/*
** read/write nifti-1 to vista format
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

#include <nifti/nifti1.h>
#include <nifti/nifti1_io.h>

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352


/* write nifti1 header  */
void VWriteNiftiHeader(FILE *fp,nifti_1_header *hdr)
{
  size_t ret = fwrite(hdr, MIN_HEADER_SIZE, 1, fp);
  if (ret != 1) VError("Error writing header");

  /*  write extender pad */
  nifti1_extender pad;
  pad.extension[0] = 0;
  pad.extension[1] = 0;
  pad.extension[2] = 0;
  pad.extension[3] = 0;
  ret = fwrite(&pad, 4, 1, fp);
  if (ret != 1) VError("Error writing header file extension pad");
}


#define xsetval(datatype) \
{ \
  datatype *A = (datatype *) VCalloc(ndata,nsize); \
  size_t ret=0,slice,row,col,ti;		   \
  size_t index = 0; \
  for (ti=0; ti<nt; ti++) {		  \
    for (slice=0; slice<nslices; slice++) { \
      for (row=0; row<nrows; row++) { \
	for (col=0; col<ncols; col++) {	  \
	  double val = VPixel(src[slice],ti,row,col,datatype); \
	  A[index] = (datatype)val; \
	  index++; \
	} \
      } \
    } \
  } \
  ret = fwrite(A,nsize,ndata,fp); \
  if (ret != ndata) VError("Error writing data"); \
  VFree(A); \
}

void VWriteData(VImage *src,int nslices,int nrows,int ncols,int nt,
		int datatype,size_t ndata,size_t nsize,FILE *fp)
{
  switch(datatype) {
  case DT_BINARY:
    xsetval(VBit);
    break;
  case DT_UNSIGNED_CHAR:
    xsetval(VUByte);
    break;
  case DT_SIGNED_SHORT:
    xsetval(short);
    break;
  case DT_SIGNED_INT:
    xsetval(int);
    break;
  case DT_FLOAT:
    xsetval(float);
    break;
  case DT_DOUBLE:
    xsetval(double);
    break;
  case DT_INT8:
    xsetval(VSByte);
    break;
  case DT_UINT16:
    xsetval(unsigned short);
    break;
  case DT_UINT32:
    xsetval(unsigned int);
    break;
  case DT_INT64:
    xsetval(long);
    break;
  case DT_UINT64:
    xsetval(unsigned long);
    break;
  default:
    VError(" unknown datatype %d",datatype);
  }
}

  /* get data type */
int GetDataType(VRepnKind repn)
{
  int datatype = 0;
  switch (repn) {
  case VBitRepn:
    datatype = DT_BINARY;
    break;
  case VUByteRepn:
    datatype = DT_UNSIGNED_CHAR;
    break;
  case VShortRepn:
    datatype = DT_INT16;
    break;
  case VUShortRepn:
    datatype = DT_UINT16;
    break;
  case VIntegerRepn:
    datatype = DT_INT32;
    break;
  case VUIntegerRepn:
    datatype = DT_UINT32;
    break;
  case VLongRepn:
    datatype = DT_INT64;
    break;
  case VULongRepn:
    datatype = DT_UINT64;
    break;
  case VFloatRepn:
    datatype = DT_FLOAT;
    break;
  case VDoubleRepn:
    datatype = DT_DOUBLE;
    break;
  default:
    VError(" unknown datatype");
  }
  return datatype;
}

void Vista2Nii3D(VAttrList list,nifti_1_header *hdr,FILE *fp)
{
  VImage src = VReadImage(list);
  if (src == NULL) VError(" no input image found");
  size_t ncols = VImageNColumns(src);
  size_t nrows = VImageNRows(src);
  size_t nslices = VImageNBands(src);


  /* data type */
  VRepnKind repn = VPixelRepn(src);
  int datatype = GetDataType(repn);
  

  /* get voxel size */
  VAttrList geolist = VGetGeoInfo(list);
  double *pixdim = (double *)VCalloc(8,sizeof(double));
  if (geolist != NULL) pixdim = VGetGeoPixdim(geolist,pixdim);
  VString str=NULL;
  VFloat resx=1.0,resy=1.0,resz=1.0;
  if (VGetAttr (VImageAttrList (src), "voxel", NULL,
		VStringRepn, (VPointer) & str) == VAttrFound) {
    sscanf(str,"%f %f %f",&resx,&resy,&resz);
    pixdim[1] = resx;
    pixdim[2] = resy;
    pixdim[3] = resz;
  }


  /* write header infos */
  hdr->pixdim[1] = pixdim[1];
  hdr->pixdim[2] = pixdim[2];
  hdr->pixdim[3] = pixdim[3];
  hdr->pixdim[4] = pixdim[4];
  /* hdr->xyzt_units = SPACE_TIME_TO_XYZT(NIFTI_UNITS_MM,NIFTI_UNITS_MSEC); */

  hdr->dim[0] = 3;
  hdr->dim[1] = VImageNColumns(src);
  hdr->dim[2] = VImageNRows(src);
  hdr->dim[3] = VImageNBands(src);
  hdr->dim[4] = 1;
  hdr->dim[5] = 0;
  hdr->dim[6] = 0;
  hdr->dim[7] = 0;
  hdr->datatype = datatype;
  hdr->bitpix  = VPixelPrecision(src);


  /* write header */
  VWriteNiftiHeader(fp,hdr);


  /* write to disk */
  size_t ndata = ncols*nrows*nslices;
  size_t bytesize = 8;
  if (repn == VBitRepn) bytesize = 1;
  size_t nsize = VPixelPrecision(src)/bytesize;
  size_t ret = fwrite(VImageData(src),nsize,ndata,fp);
  if (ret != ndata) VError("Error writing data");
}


void Vista2Nii4D(VAttrList list,nifti_1_header *hdr,FILE *fp)
{  
  int nrows=0,ncols=0,nt=0;
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  VImageDimensions(src,nslices,&nt,&nrows,&ncols);


  /* read repetition time, voxel size */
  VAttrList geolist = VGetGeoInfo(list);
  double *pixdim = (double *)VCalloc(8,sizeof(double));
  if (geolist != NULL) pixdim = VGetGeoPixdim(geolist,pixdim);
  VString str=NULL;
  VFloat resx=1.0,resy=1.0,resz=1.0,tr=1.0;
  if (VGetAttr (VImageAttrList (src[0]), "repetition_time", NULL,
		VFloatRepn, (VPointer) & tr) == VAttrFound) {
    pixdim[4] = tr;
  }
  if (VGetAttr (VImageAttrList (src[0]), "voxel", NULL,
		VStringRepn, (VPointer) & str) == VAttrFound) {
    sscanf(str,"%f %f %f",&resx,&resy,&resz);
    pixdim[1] = resx;
    pixdim[2] = resy;
    pixdim[3] = resz;
  }
  fprintf(stderr," image dims: %d x %d x %d,  nt: %d\n",ncols,nrows,nslices,nt);
  fprintf(stderr," voxel resolution: %.4f %.4f %.4f,  TR= %.4f secs\n",resx,resy,resz,tr/1000.0);


  /* data type */
  VRepnKind repn = VPixelRepn(src[0]);
  int datatype = GetDataType(repn);


  /* image dims */
  size_t bytesize = 8;
  hdr->dim[0] = 4;
  hdr->dim[1] = ncols;
  hdr->dim[2] = nrows;
  hdr->dim[3] = nslices;
  hdr->dim[4] = nt;
  hdr->dim[5] = 0;
  hdr->dim[6] = 0;
  hdr->dim[7] = 0;
  hdr->datatype = datatype;
  hdr->bitpix   = VPixelPrecision(src[0]);

  hdr->pixdim[1] = pixdim[1];
  hdr->pixdim[2] = pixdim[2];
  hdr->pixdim[3] = pixdim[3];
  hdr->pixdim[4] = pixdim[4];
  /*   hdr->xyzt_units = SPACE_TIME_TO_XYZT(NIFTI_UNITS_MM,NIFTI_UNITS_MSEC); */


  /* write header */
  VWriteNiftiHeader(fp,hdr);

  size_t nsize  = hdr->bitpix/bytesize;
  size_t nrnc = (size_t)nrows*(size_t)ncols;
  size_t npix = nrnc*(size_t)nslices;
  size_t ndata = (size_t)nt * npix;

  /* write to disk */
  VWriteData(src,nslices,nrows,ncols,nt,datatype,ndata,nsize,fp);
}



void Vista_to_Nifti1(VAttrList list,VString filename)
{
  VBundle bundle=NULL;
  VImage sform=NULL;
  size_t j;
  nifti_1_header hdr;


  /* get geoinfo */
  VAttrList geolist = VGetGeoInfo(list);
  if (geolist == NULL) VError(" no geoinfo found");


  /* quaternion */
  VShort qform_code = 0;
  VShort sform_code = 0;
  VGetAttr (geolist, "qform_code", NULL,VShortRepn, (VPointer) & qform_code);
  VGetAttr (geolist, "sform_code", NULL,VShortRepn, (VPointer) & sform_code);

  if (sform_code != 0) {
    if (VGetAttr (geolist, "sform", NULL,VImageRepn, (VPointer) & sform) != VAttrFound)
      VError(" sform not found");
  }

  float *Q = NULL;
  if (qform_code != 0) {
    if (VGetAttr (geolist, "qform", NULL,VBundleRepn, (VPointer) & bundle) != VAttrFound)
      VError(" qform not found");
    Q = bundle->data;
  }

  if (VGetAttr (geolist, "pixdim", NULL,VBundleRepn, (VPointer) & bundle) != VAttrFound)
    VError(" pixdim not found");
  float *P = bundle->data;


  /* set some header fields */
  bzero((void *)&hdr, sizeof(hdr));
  hdr.sizeof_hdr = MIN_HEADER_SIZE;
  hdr.pixdim[0] = P[0];
  hdr.pixdim[1] = P[1];
  hdr.pixdim[2] = P[2];
  hdr.pixdim[3] = P[3];
  hdr.pixdim[4] = P[4];
  hdr.pixdim[5] = P[5];
  hdr.pixdim[6] = 0;
  hdr.pixdim[7] = 0;

  hdr.qform_code = qform_code;
  if (Q != NULL) {
    hdr.quatern_b = Q[0];
    hdr.quatern_c = Q[1];
    hdr.quatern_d = Q[2];
    hdr.qoffset_x = Q[3];
    hdr.qoffset_y = Q[4];
    hdr.qoffset_z = Q[5];
  }

  hdr.sform_code = sform_code;
  if (sform_code != 0) {
    for (j=0; j<4; j++) {
      hdr.srow_x[j] = VPixel(sform,0,0,j,VFloat);
      hdr.srow_y[j] = VPixel(sform,0,1,j,VFloat);
      hdr.srow_z[j] = VPixel(sform,0,2,j,VFloat);
    }
  }

  hdr.vox_offset = (float) NII_HEADER_SIZE;
  hdr.scl_slope = 0.0;
  hdr.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_MSEC;


  /* slicetimes */
  VShort sbuf=0;
  VFloat fbuf=0;
  if (VGetAttr (geolist, "slice_start", NULL,VShortRepn, (VPointer) & sbuf) == VAttrFound){
    hdr.slice_start = sbuf;
  }
  if (VGetAttr (geolist, "slice_end", NULL,VShortRepn, (VPointer) & sbuf) == VAttrFound){
    hdr.slice_end = sbuf;
  }
  if (VGetAttr (geolist, "slice_code", NULL,VShortRepn, (VPointer) & sbuf) == VAttrFound){
    hdr.slice_code = sbuf;
  }
  if (VGetAttr (geolist, "slice_duration", NULL,VFloatRepn, (VPointer) & fbuf) == VAttrFound){
    hdr.slice_duration = fbuf;
  }

  /* magic number */
  strncpy(hdr.magic, "n+1\0", 4);


  /* open output file */
  FILE *fp = fopen(filename,"w");
  if (fp == NULL) VError("Error opening file %s",filename);


  /* write image data */  
  int nimages = VAttrListNumImages(list);
  if (nimages < 2) Vista2Nii3D(list,&hdr,fp);
  else Vista2Nii4D(list,&hdr,fp);

  /* close */
  fclose(fp);
  return;
}
