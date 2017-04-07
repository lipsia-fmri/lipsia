/*
** read and write geometry information
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <viaio/Vlib.h>
#include <viaio/mu.h>


void VCopyGeoInfo(VAttrList src_list,VAttrList dst_list)
{
  VAttrListPosn src_posn,dst_posn;
  VAttrList geolist=NULL;
  VString str=NULL;

  for (VFirstAttr (src_list, & src_posn); VAttrExists (& src_posn); VNextAttr (& src_posn)) {
    if (strncmp(VGetAttrName(&src_posn), "geoinfo", strlen("geoinfo")) == 0 ) {
      VGetAttrValue (& src_posn, NULL, VAttrListRepn, & geolist);
      break;
    }
  }
  if (geolist != NULL) {
    VExtractAttr(dst_list, "geoinfo", NULL, VAttrListRepn, &str, FALSE);
    if (VLookupAttr (dst_list,"history", &dst_posn) == TRUE) {
      VInsertAttr (&dst_posn,TRUE,"geoinfo",NULL,VAttrListRepn,geolist);
    }
    else
      VPrependAttr(dst_list,"geoinfo",NULL,VAttrListRepn,geolist);
  }
}



void VSetGeoInfo(VAttrList geolist,VAttrList dst_list)
{
  if (geolist == NULL) return;
  VString str=NULL;
  VAttrListPosn dst_posn;
  VExtractAttr(dst_list, "geoinfo", NULL, VAttrListRepn, &str, FALSE);
  if (VLookupAttr (dst_list,"history", &dst_posn) == TRUE) {
    VInsertAttr (&dst_posn,TRUE,"geoinfo",NULL,VAttrListRepn,geolist);
  }
  else
    VPrependAttr(dst_list,"geoinfo",NULL,VAttrListRepn,geolist);
}


VAttrList VGetGeoInfo(VAttrList list)
{
  VAttrList geolist=NULL;
  VAttrListPosn posn;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VAttrListRepn) continue;
    if (strcmp(VGetAttrName(&posn), "geoinfo") == 0 ||
	strcmp(VGetAttrName(&posn), "imageinfo") == 0 ) {
      VGetAttrValue (& posn, NULL, VAttrListRepn, & geolist);
      break;
    }
  }
  return geolist;
}


double *VGetGeoPixdim(VAttrList geolist,double *pixdim)
{
  int i,n=8;
  VBundle bundle=NULL;
  if (pixdim == NULL) pixdim = (double *) VCalloc(n,sizeof(double));
  for (i=0; i<n; i++) pixdim[i] = 1.0;
  if (VGetAttr (geolist, "pixdim", NULL,VBundleRepn, (VPointer) & bundle) == VAttrFound) {
    float *D = bundle->data;
    for (i=0; i<n; i++) pixdim[i] = (double)D[i];
  }
  else
    VError(" pixdim information not found in geolist");
  return pixdim;
}

double *VGetGeoDim(VAttrList geolist,double *dim)
{
  int i,n=8;
  VBundle bundle=NULL;
  if (dim == NULL) dim = (double *) VCalloc(n,sizeof(double));
  if (VGetAttr (geolist, "dim", NULL,VBundleRepn, (VPointer) & bundle) == VAttrFound) {
    float *D = bundle->data;
    for (i=0; i<n; i++) dim[i] = (double)D[i];
  }
  else
    VError(" dim information not found in geolist");
  return dim;
}


void VSetGeoDim(VAttrList geolist,double *dim)
{
  if (dim == NULL) return;
  int i,n=8;
  VBundle bundle=NULL;
  if (VGetAttr (geolist, "dim", NULL,VBundleRepn, (VPointer) & bundle) == VAttrFound) {
    float *D = bundle->data;
    for (i=0; i<n; i++) D[i] = dim[i];
  }
  else {
    VAttrList dlist = VCreateAttrList();
    float *D = (float *) VCalloc(8,sizeof(float));
    for (i=0; i<8; i++) D[i] = (float)dim[i];
    bundle = VCreateBundle ("bundle",dlist,8*sizeof(float),(VPointer)D);
    VSetAttr(geolist,"dim",NULL,VBundleRepn,bundle);
  }
}


void VSetGeoPixdim(VAttrList geolist,double *pixdim)
{
  if (pixdim == NULL) return;
  int i,n=8;
  VBundle bundle=NULL;
  if (VGetAttr (geolist, "pixdim", NULL,VBundleRepn, (VPointer) & bundle) == VAttrFound) {
    float *D = bundle->data;
    for (i=0; i<n; i++) D[i] = pixdim[i];
  }
  else {
    VAttrList dlist = VCreateAttrList();
    float *D = (float *) VCalloc(8,sizeof(float));
    for (i=0; i<8; i++) D[i] = (float)pixdim[i];
    bundle = VCreateBundle ("bundle",dlist,8*sizeof(float),(VPointer)D);
    VSetAttr(geolist,"pixdim",NULL,VBundleRepn,bundle);
  }
}


double *VGetGeoSlicetime(VAttrList geolist,double *slicetime,int nslices)
{
  VShort sbuf=0;
  VFloat fbuf=0;
  int slice_start=-1,slice_end=-1,slice_code=-1;
  double slice_duration=-1.0;

  if (VGetAttr (geolist, "slice_start", NULL,VShortRepn, (VPointer) & sbuf) == VAttrFound){
    slice_start = sbuf;
  }
  if (VGetAttr (geolist, "slice_end", NULL,VShortRepn, (VPointer) & sbuf) == VAttrFound){
    slice_end = sbuf;
  }
  if (VGetAttr (geolist, "slice_code", NULL,VShortRepn, (VPointer) & sbuf) == VAttrFound){
    slice_code = sbuf;
  }
  if (VGetAttr (geolist, "slice_duration", NULL,VFloatRepn, (VPointer) & fbuf) == VAttrFound){
    slice_duration = fbuf;
  }
  if (slice_code <= 0) return NULL;
  if (slice_duration <= 0) return NULL;
  if (slice_start < 0) slice_start=0;
  if (slice_end < 0) slice_end=nslices-1;


  if (slice_end >= nslices) {
    VWarning(" inconsistent slicetime info, nslices=%d, slice_end=%d",nslices,slice_end);
    slice_end=nslices-1;
  }
  if (slicetime == NULL) slicetime = (double *) VCalloc(nslices,sizeof(double));

  VGetSlicetimes(slice_start,slice_end,slice_code,slice_duration,slicetime);
  return slicetime;
}
