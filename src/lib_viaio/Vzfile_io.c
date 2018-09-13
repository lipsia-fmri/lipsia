/*
** file IO update for lipsia programs
** lipsia programs can read files with extensions .v  .v.gz  .nii  .nii.gz
** lipsia outputs files with extensions .v
**
** G.Lohmann, MPI-KYB, May 2017, updated Feb 2018
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


/* check if gzipped */
int CheckGzip(char *filename)
{
  int n = strlen(filename);
  if (n < 3) return 0;
  if (filename[n-1] == 'z' && filename[n-2] == 'g' && filename[n-3] == '.') return 1;
  else return 0;
}

/* check if input comes from a pipe */
int CheckPipe(char *filename)
{
  if (filename[0]=='-') return 1;
  else return 0;
}


/* read gzipped file into data container, unpack one MByte at a time */
char *VReadGzippedDataZ(gzFile *fp,size_t *len)
{
  int i=0,err=0,ulen=0;
  char *buf = (char *) VCalloc(N,sizeof(char));
  char *data = (char *)malloc(sizeof(char));
  
  if (fp==0) VError(" gzip open error, fp=0");

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
  gzclose(fp);
  VFree(buf);

  (*len) = n;
  return data;
}


char *VReadUnzippedDataZ(FILE *fp,VBoolean nofail,size_t buflen)
{
  char *data = (char *) malloc(buflen);
  size_t ret = fread(data,1,buflen,fp);
  if (ret != buflen) 
    VError(" fread error= %ld\n",ret);
  fclose(fp);
  return data;
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
  if (n==1) {
    if (filename[0]=='-') type = 0;  /* pipe */
  }
  if (n > 2) {
    last = &filename[n-2];    
    if (strncmp(last,".v",2) == 0) type = 0;
  }
  if (n > 5) {
    last = &filename[n-5];
    if (strncmp(last,".v.gz",5) == 0) type = 1;
  }
  if (n > 4) {
    last = &filename[n-4];
    if (strncmp(last,".nii",4) == 0) type = 2;
  }
  if (n > 7) {
    last = &filename[n-7];
    if (strncmp(last,".nii.gz",7) == 0) type = 3;
  }
  return type;
}


/*
** read data into vista attrlist 
** differs from 'VReadAttrList' in that it supports pipes.
** therefore, the first parameter (*in_file) has to be different.
*/
VAttrList VReadAttrListZ(FILE *in_file,VString in_filename,
			 VLong tr,VBoolean attrtype,VBoolean do_scaling)
{
  VAttrList list = NULL;
  FILE *fp=NULL;
  size_t bufsize=0;
  char *databuffer = NULL;
  VBoolean ok = FALSE;
  gzFile gz_file;

  struct stat st;
  if (stat(in_filename, &st) != 0) {
    VError(" Error opening %s,  %s",in_filename,strerror(errno));
  }
  size_t buflen = st.st_size;

  int format_type = getformat(in_filename);


  switch (format_type) {

  case 0:    /* vista format, unzipped or pipe */
    list = VReadFile (in_file,NULL);
    fclose(in_file);
    break;

  case 1:   /* vista format, gzipped */
    gz_file = gzopen(in_filename,"rb");
    if (gz_file == NULL) VError("Error: Failed to gzopen %s\n", in_filename);
    databuffer = VReadGzippedDataZ(gz_file,&bufsize);
    fp = VOpenStream(databuffer,bufsize);
    list = VReadFile (fp,NULL);
    fclose(fp);
    break;

  case 2:  /* nifti format, unzipped */
    databuffer = VReadUnzippedDataZ(in_file,TRUE,buflen);
    list = Nifti1_to_Vista(databuffer,tr,attrtype,do_scaling,&ok);
    if (! list) VError(" can't read input file");
    if (ok) VFree(databuffer);
    break;

  case 3:  /* nifti format, gzipped */
    gz_file = gzopen(in_filename,"rb");
    if (gz_file == NULL) VError("Error: Failed to gzopen %s\n", in_filename);
    databuffer = VReadGzippedDataZ(gz_file,&bufsize);
    list = Nifti1_to_Vista(databuffer,tr,attrtype,do_scaling,&ok);
    if (ok) VFree(databuffer);
    break;

  default:
    VError(" ZZZ,  illegal input format %d",format_type);
  }
  if (! list) VError(" can't read input file");

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
    if (VPixelRepn(images[0]) != VPixelRepn(images[i])) VError(" inconsistent pixel repn");
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
