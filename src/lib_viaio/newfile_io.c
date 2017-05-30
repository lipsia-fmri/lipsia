/*
** file IO update for lipsia programs
** lipsia programs can read files with extensions .v  .v.gz  .nii  .nii.gz
** lipsia outputs files with extensions .v
**
** G.Lohmann, MPI-KYB, May 2017
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>
#include <errno.h>

#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>

#include <zlib.h>

#define N 1048576  /* one megabyte */

extern VAttrList Nifti1_to_Vista(char *databuffer,VLong tr,VBoolean,VBoolean,VBoolean *);


/* check if gzipped */
int CheckGzip(char *filename)
{
  int n = strlen(filename);
  if (n < 3) return 0;
  if (filename[n-1] == 'z' && filename[n-2] == 'g' && filename[n-3] == '.') return 1;
  else return 0;
}


/* read gzipped file into data container, unpack one MByte at a time */
char *VReadGzippedData(char *filename,size_t *len)
{
  int i=0,err=0,ulen=0;
  char *buf = (char *) VCalloc(N,sizeof(char));
  gzFile *fp = (gzFile *)gzopen(filename,"rb");
  if (fp == NULL) VError("Error: Failed to gzopen %s\n", filename);

  char *data = (char *)malloc(sizeof(char));
  /* int gflag = 0; */
  size_t n=0,m=0;
  while(!gzeof(fp)) {
    ulen = gzread(fp,buf,(unsigned int)N);
    if (ulen == Z_NULL) VError(" gzread error, %s",gzerror(fp, &err));

    /* check gzip magic number */
    /*
    if (gflag == 0) {
      if ((data[0] == 0x1f) && (data[1] == 0x8b)) gflag = 1;
      else VError(" not in gzip format, magic number incorrect, %x %x ",data[0],data[1]);
    }
    */

    m = n;
    n += (size_t)ulen;
    data = (char *)realloc(data,(size_t)n);
    if (data == NULL) {
      VError(" gzread, error in realloc");
    }
    for (i=0; i<ulen; i++) {
      data[i+m] = buf[i];
    }
    m = n;
  }
  VFree(buf);

  (*len) = n;
  return data;
}


char *VReadUnzippedData(char *filename,VBoolean nofail,size_t *size)
{
  struct stat st;
  stat(filename, &st);
  size_t buflen = st.st_size;
  char *data = (char *) malloc(buflen);
  FILE *fp = VOpenInputFile (filename,nofail);
  if (fp == NULL) VError("Error opening file %s\n", filename);
  int ret = fread(data,1,buflen,fp);
  if (ret != buflen) 
    VError("Error reading from %s, fread error= %ld\n",filename,ret);
  fclose(fp);
  (*size) = buflen;
  return data;
}


char *VReadDataContainer(char *filename,VBoolean nofail,size_t *size)
{
  char *databuffer=NULL;
  size_t buflen=0;

  if (CheckGzip(filename)) {
    databuffer = VReadGzippedData(filename,&buflen);
  }
  else {
    databuffer = VReadUnzippedData(filename,nofail,&buflen);
  }
  (*size) = buflen;
  return databuffer;
}


FILE *VOpenStream(char *databuffer,size_t size)
{
  extern FILE *fmemopen(void *buf, size_t size, const char *mode);
  FILE *fs = fmemopen (databuffer,size,"r");
  if (!fs) VError(" fmemopen failed, %s",strerror(errno));
  return fs;
}



/* get type of format (.nii or .v) */
int getformat(char *filename)
{
  char *last=NULL;
  int n = strlen(filename);

  int type = -1;
  if (n > 2) {
    last = &filename[n-2];    
    if (strncmp(last,".v",2) == 0) type = 0;
  }
  if (n > 5) {
    last = &filename[n-5];
    if (strncmp(last,".v.gz",5) == 0) type = 0;
  }
  if (n > 4) {
    last = &filename[n-4];
    if (strncmp(last,".nii",4) == 0) type = 1;
  }
  if (n > 7) {
    last = &filename[n-7];
    if (strncmp(last,".nii.gz",7) == 0) type = 1;
  }
  return type;
}


/* read data into vista attrlist */
VAttrList VReadAttrList(VString filename,VLong tr,VBoolean attrtype,VBoolean do_scaling)
{
  VAttrList list = NULL;

  /* get format type */
  int type = getformat(filename);
  if (type != 0 && type != 1) VError( "illegal format type %d",type);


  /* read data buffer */
  size_t bufsize=0;
  char *databuffer = VReadDataContainer(filename,TRUE,&bufsize);


  /* vista format */
  if (type == 0) { 
    FILE *fp = VOpenStream(databuffer,bufsize);
    list = VReadFile (fp,NULL);
    if (! list) VError(" can't read input file");
    fclose(fp);
  }

  /* nifti format */
  else if (type == 1) {
    VBoolean ok = FALSE;
    list = Nifti1_to_Vista(databuffer,tr,attrtype,do_scaling,&ok);
    if (ok) VFree(databuffer);
  }

  else {
    VError(" illegal input format");
  }

  return list;
}


/* read first image from file */
VImage VReadImage(VAttrList list)
{
  VAttrListPosn posn;
  VImage image=NULL;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & image);
    break;
  }
  return image;
}


int VAttrListNumImages(VAttrList list)
{
  VAttrListPosn posn;
  VImage tmp=NULL;
  int n=0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & tmp);
    n++;
  }
  return n;
}

VImage *VAttrListGetImages(VAttrList list,int n)
{
  VAttrListPosn posn;
  VImage *images = (VImage *) VCalloc(n,sizeof(VImage));
  int i=0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    if (i >= n) break;
    VGetAttrValue (& posn, NULL,VImageRepn, & images[i]);
    i++;
  }
  return images;
}


void VImageDimensions(VImage *src,int nimages,int *bands,int *rows,int *cols)
{
  int i=0;
  int nbands = VImageNBands(src[0]);
  int nrows = VImageNRows(src[0]);
  int ncols = VImageNColumns(src[0]);

  for (i=1; i<nimages; i++) {
    if (VImageNBands(src[i]) != nbands) VError(" inconsistent num bands");
    if (VImageNRows(src[i]) != nrows) VError(" inconsistent num rows");
    if (VImageNColumns(src[i]) != ncols) VError(" inconsistent num columns");
  }
  *bands = nbands;
  *rows = nrows;
  *cols = ncols;
}
