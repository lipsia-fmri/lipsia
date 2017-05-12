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



void Vista2Nii3D(VAttrList list,float *D,nifti_1_header *hdr,FILE *fp)
{  
  VAttrListPosn posn;
  VImage src=NULL;
  size_t ncols   = (size_t) D[1];
  size_t nrows   = (size_t) D[2];
  size_t nslices = (size_t) D[3];
  fprintf(stderr," image dims: %ld %ld %ld\n",ncols,nrows,nslices);

  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src);
    break;
  }
  if (src == NULL) VError(" no image found");


  /* get data type */
  VRepnKind repn = VPixelRepn(src);
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

  /* header infos */
  hdr->dim[0] = 3;
  hdr->dim[1] = D[1];
  hdr->dim[2] = D[2];
  hdr->dim[3] = D[3];
  hdr->dim[4] = 1;
  hdr->dim[5] = 0;
  hdr->dim[6] = 0;
  hdr->dim[7] = 0;
  hdr->datatype  = datatype;
  hdr->bitpix    = VPixelPrecision(src);


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


void Vista2Nii4D(VAttrList list,float *D,nifti_1_header *hdr,FILE *fp)
{
  size_t ncols   = (size_t) D[1];
  size_t nrows   = (size_t) D[2];
  size_t nslices = (size_t) D[3];
  size_t nt      = (size_t) D[4];
  fprintf(stderr," image dims: %ld %ld %ld, nt= %ld\n",ncols,nrows,nslices,nt);

  VAttrListPosn posn;
  VImage *src = (VImage *) VCalloc(nslices,sizeof(VImage));
  size_t i = 0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src[i]);
    if (VPixelRepn(src[i]) != VShortRepn) continue;
    i++;
  }

  /* image dims */
  hdr->dim[0] = 4;
  hdr->dim[1] = D[1];
  hdr->dim[2] = D[2];
  hdr->dim[3] = D[3];
  hdr->dim[4] = D[4];
  hdr->dim[5] = 0;
  hdr->dim[6] = 0;
  hdr->dim[7] = 0;
  hdr->datatype = DT_SIGNED_SHORT;
  hdr->bitpix   = VPixelPrecision(src[0]);

  size_t nrnc  = nrows*ncols;
  size_t npix  = nrnc*nslices;
  size_t ndata = (npix*nt);
  VShort *data = (VShort *) VCalloc(ndata,sizeof(VShort));
  if (!data) VError(" err allocating data ");

  size_t slice,row,col,ti;
  for (slice=0; slice<nslices; slice++) {
    for (ti=0; ti<nt; ti++) {
      for (row=0; row<nrows; row++) {
	for (col=0; col<ncols; col++) {
	  VShort val = (VShort) VPixel(src[slice],ti,row,col,VShort);
	  const size_t index = col + row*ncols + slice*nrnc + ti*npix;
	  data[index] = val;
	}
      }
    }
  }

  /* write header */
  VWriteNiftiHeader(fp,hdr);

  /* write to disk */
  size_t bytesize = 8;
  size_t nsize = VPixelPrecision(src[0])/bytesize;
  size_t ret = fwrite(data,nsize,ndata,fp);
  if (ret != ndata) VError("Error writing data");
  VFree(data);
}



void Vista_to_Nifti1(VAttrList list,VString filename)
{
  VImage sform=NULL;
  size_t j;
  nifti_1_header hdr;


  VBundle bundle;
  VAttrList geolist = VGetGeoInfo(list);
  if (geolist == NULL) VError(" no geoinfo found");


  /* dim */
  if (VGetAttr (geolist, "dim", NULL,VBundleRepn, (VPointer) & bundle) != VAttrFound)
    VError("dim not found");
  float *D = bundle->data;
  long nt = D[4];

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
  if (nt < 2) Vista2Nii3D(list,D,&hdr,fp);
  else Vista2Nii4D(list,D,&hdr,fp);

  /* close */
  fclose(fp);
  return;
}
