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

#define N 5242880  /* five megabyte */


extern VAttrList Nifti1_to_Vista(char *databuffer,VLong tr,VBoolean,VBoolean,VBoolean *);
extern int CheckGzip(char *filename);
extern FILE *VOpenStream(char *databuffer,size_t size);


/* read gzipped file into data container, unpack one MByte at a time */
char *XReadGzippedData(char *filename,size_t *len)
{
  int i=0,err=0,ulen=0;
  char *buf = (char *) VCalloc(N,sizeof(char));
  gzFile fp = (gzFile) gzopen(filename,"rb");
  if (fp == NULL) VError("Error: Failed to gzopen %s\n", filename);

  char *data = (char *)malloc(sizeof(char));

  size_t n=0,m=0;
  while(!gzeof(fp)) {
    ulen = gzread(fp,buf,(unsigned int)N);
    if (ulen == Z_NULL) VError(" gzread error, %s",gzerror(fp, &err));

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


char *XReadUnzippedData(char *filename,VBoolean nofail,size_t *size)
{
  struct stat st;
  stat(filename, &st);
  size_t buflen = st.st_size;
  char *data = (char *) malloc(buflen);
  FILE *fp = VOpenInputFile (filename,nofail);
  if (fp == NULL) VError("Error opening file %s\n", filename);
  size_t ret = fread(data,1,buflen,fp);
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
    databuffer = XReadGzippedData(filename,&buflen);
  }
  else {
    databuffer = XReadUnzippedData(filename,nofail,&buflen);
  }
  (*size) = buflen;
  return databuffer;
}



/* get type of format (.nii or .v) */
int xgetformat(char *filename)
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
  int type = xgetformat(filename);
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
