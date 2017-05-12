#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
#include "zlib.h"


#define N 1024*32

void zerr(int ret)
{
  switch (ret) {
  case Z_ERRNO:
    fputs("gzip error, Z_ERRNO\n", stderr);
    break;
  case Z_STREAM_ERROR:
    fputs("gzip error, invalid compression level\n", stderr);
    break;
  case Z_DATA_ERROR:
    fputs("gzip error, invalid or incomplete deflate data\n", stderr);
    break;
  case Z_MEM_ERROR:
    fputs("gzip error, out of memory\n", stderr);
    break;
  case Z_VERSION_ERROR:
    fputs("gzip error, zlib version mismatch\n", stderr);
    break;
  default:
     fputs("gzip error, some error\n", stderr);
     break;
  }
}


/* read gzipped file into data container */
char *ReadGzippedData(char *filename,size_t *len)
{
  size_t i=0,n=0,m=N;
  char buf[N];
  gzFile *fp = (gzFile *)gzopen(filename,"rb");
  if (fp == NULL) ("Error: Failed to gzopen %s\n", filename);

  /* get file size */
  size_t buflen=0;
  while(!gzeof(fp)) {
    int rtcode = gzread(fp,buf,m);
    if (rtcode == Z_NULL) zerr(rtcode);
    buflen += (size_t)rtcode;
  }
  gzrewind(fp);

  /* read into buffer */
  char *data = (char *) malloc(buflen);
  n=0;
  while(!gzeof(fp)) {
    int rtcode = gzread(fp,buf,m);
    if (rtcode == Z_NULL) zerr(rtcode);
    for (i=0; i<m; i++) {
      if (i+n < buflen) data[i+n] = buf[i];
    }
    n += m;
  }
  
  (*len) = buflen;
  return data;
}




/* read unzipped file into data container */
char *ReadUnzippedData(char *filename,size_t *len)
{
  /* get file size */
  struct stat st;
  stat(filename, &st);
  size_t buflen = st.st_size;

  /* read from file */
  FILE *fp = fopen(filename,"rb");
  if (fp == NULL) ("Error: Failed to open %s\n", filename);
  char *data = (char *) malloc(buflen);
  int ret = fread(data,1,buflen,fp);
  if (ret != buflen) 
    VError("Error reading from %s, fread error= %ld\n",filename,ret);
  fclose(fp);
  
  (*len) = buflen;
  return data;
}
