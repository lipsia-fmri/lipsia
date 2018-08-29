/*
** construct a design file from an ASCII file
** no hemodynamic modelling
**
** G.Lohmann, Dec 2006
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

#define LEN  10000   /* max number of characters per line in file */
#define NCOV   200   /* max number of additional covariates     */
#define M       50

gsl_matrix *VReadCovariates(VString cfile,VBoolean normalize)
{
  FILE *fp=NULL;
  char buf[LEN],token[80];
  fp = fopen(cfile,"r");
  if (fp == NULL) VError(" err opening covariates file %s",cfile);

  /* get matrix dimensions */
  float x=0;
  int i=0,j=0,nrows=0,ncols=0;
  while (!feof(fp)) {
    for (j=0; j<LEN; j++) buf[j] = '\0';
    fgets(buf,LEN,fp);
    if (buf[0] == '#' || buf[0] == '%') continue;
    if (strlen(buf) < 1) continue;

    j = 0;
    while (VStringToken(buf,token,j,M)) {
      sscanf(token,"%f",&x);
      if (j >= NCOV) VError(" too many additional covariates (%d), max is %d",j,NCOV-1);
      j++;
    }
    if (j > ncols) ncols = j;
    nrows++;
  }
  rewind(fp);
  /* nrows--; */
  fprintf(stderr," covariates,  nrows: %d, ncols: %d\n",nrows,ncols);


  /* read into matrix */
  gsl_matrix *dest = gsl_matrix_calloc(nrows,ncols);
  i=0;
  while (!feof(fp)) {
    for (j=0; j<LEN; j++) buf[j] = '\0';
    fgets(buf,LEN,fp);
    if (buf[0] == '#' || buf[0] == '%') continue;
    if (strlen(buf) < 1) continue;

    j = 0;
    while (VStringToken(buf,token,j,M)) {
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
  double u=0,sum1=0,sum2=0,nx=0,mean=0,sigma=0,tiny=1.0e-8;
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
    if (sigma < tiny) VError(" i=%d, sigma= %f",i,sigma);

    for (j=0; j<dest->size1; j++) {
      u = gsl_matrix_get(dest,j,i);
      u = (u-mean)/sigma;
      gsl_matrix_set(dest,j,i,u);
    }
  }
  return dest;
}

