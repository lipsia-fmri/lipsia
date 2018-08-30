/*
** read nuisance covariates file into a matrix
**
** G.Lohmann, MPI-KYB, 2018
*/

#include <viaio/Vlib.h>
#include <viaio/mu.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


extern int VStringToken (char *,char *,int,int);
extern int line_empty(char *buf,int len);
extern int CheckBuffer(char *buf,int len);
extern int VistaFormat(char *buf,int len);
extern int test_ascii(int val);


gsl_matrix *VReadCovariates(VString cfile,VBoolean normalize)
{
  float x=0;
  int i=0,j=0,nrows=0,ncols=0,len=10000,tlen=80;
  char *buf = (char *)VCalloc(len,sizeof(char));
  char *token = (char *)VCalloc(tlen,sizeof(char));

  FILE *fp = fopen(cfile,"r");
  if (fp == NULL) VError(" err opening covariates file %s",cfile);
  

  /* get matrix dimensions */
  while (!feof(fp)) {

    for (j=0; j<len; j++) buf[j] = '\0';
    if (fgets(buf,len,fp) == NULL) break;
    if (VistaFormat(buf,len) > 1) VError(" File %s must be a text file",cfile);
    if (CheckBuffer(buf,len) < 1) continue;
    if (! test_ascii((int)buf[0])) VError(" File %s: line %d begins with an illegal character (%c)",cfile,nrows,buf[0]);

    j = 0;
    while (VStringToken(buf,token,j,tlen)) {
      sscanf(token,"%f",&x);
      j++;
    }
    if (j > ncols) ncols = j;
    nrows++;
  }
  rewind(fp);
  fprintf(stderr," Nuisance covariates:  %d x %d\n",nrows,ncols);
  if (ncols < 1) VError(" Nuisance covariates: no columns in file %s",cfile);


  /* read into matrix */
  gsl_matrix *dest = gsl_matrix_calloc(nrows,ncols);
  i=0;
  while (!feof(fp)) {
    for (j=0; j<len; j++) buf[j] = '\0';
    if (fgets(buf,len,fp) == NULL) break;
    if (CheckBuffer(buf,len) < 1) continue;

    j = 0;
    while (VStringToken(buf,token,j,tlen)) {
      if (j >= ncols) break;
      sscanf(token,"%f",&x);
      gsl_matrix_set(dest,i,j,(double)x);
      j++;
    }
    i++;
    if (i >= nrows) break;
  }
  if (!normalize) return dest;


  /* normalize */
  double u=0,sum1=0,sum2=0,nx=0,mean=0,sigma=0,tiny=1.0e-6;
  for (i=0; i<dest->size2; i++) {
    sum1 = sum2 = nx = 0;
    for (j=0; j<dest->size1; j++) {
      u = gsl_matrix_get(dest,j,i);
      sum1 += u;
      sum2 += u*u;
      nx++;
    }
    mean = sum1/nx;
    sigma = sqrt((sum2 - nx * mean * mean) / (nx - 1.0));
    if (sigma < tiny) VError(" Nuisance covariates:  column %d is constant",i);

    for (j=0; j<dest->size1; j++) {
      u = gsl_matrix_get(dest,j,i);
      u = (u-mean)/sigma;
      gsl_matrix_set(dest,j,i,u);
    }
  }
  return dest;
}

