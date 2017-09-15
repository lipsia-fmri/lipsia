/*
** GLM (general linear modeling)
**
** G.Lohmann
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>


typedef struct TrialStruct {
  int   id;
  float onset;
  float duration;
  float height;
} Trial;


#define ABS(x) ((x) > 0 ? (x) : -(x))
#define SQR(x) ((x)*(x))

extern float kth_smallest(float *a, size_t n, size_t k);
#define Median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

extern void VectorConvolve(const double *src,gsl_vector *dst,gsl_vector *kernel);
extern gsl_matrix *PseudoInv(gsl_matrix *A,gsl_matrix *B);
extern void printmat(gsl_matrix *R,char *str);
extern void printvec(gsl_vector *x,char *str);


/*
** general linear model (GLM)
*/
void VGLM(gsl_matrix *Data,gsl_matrix *X,gsl_matrix *XInv,gsl_vector *con,VImage map,VImage zmap)
{
  int i;
  int m = Data->size2;
  int n = con->size;
  gsl_set_error_handler_off();
  gsl_vector *y = gsl_vector_calloc (m);
  gsl_vector *beta = gsl_vector_calloc (n);


  /* compute pseudoinverse */
  XInv = PseudoInv(X,XInv);


  /* main loop */
  VFillImage(zmap,VAllBands,0);
  size_t nvox=0;
  for (nvox=0; nvox < Data->size1; nvox++) {

    y->data = gsl_matrix_ptr(Data,nvox,0);
    gsl_blas_dgemv(CblasNoTrans,1.0,XInv,y,0.0,beta);

    /* contrast image */
    double sum=0;
    for (i=0; i<beta->size; i++) {
      sum += (beta->data[i]*con->data[i]);
    }
    int b = VPixel(map,0,0,nvox,VShort);
    int r = VPixel(map,0,1,nvox,VShort);
    int c = VPixel(map,0,2,nvox,VShort);
    VPixel(zmap,b,r,c,VFloat) = sum;
  }

  gsl_vector_free(beta);
  gsl_vector_free(y);
}





int test_ascii(int val)
{
  if (val >= 'a' && val <= 'z') return 1;
  if (val >= 'A' && val <= 'Z') return 1;
  if (val >= '0' && val <= '9') return 1;
  if (val ==  ' ') return 1;
  if (val == '\0') return 1;
  if (val == '\n') return 1;
  if (val == '\r') return 1;
  if (val == '\t') return 1;
  if (val == '\v') return 1;
  return 0;
}



/* parse design file */
Trial *ReadDesign(VString designfile,int *numtrials,int *nevents)
{
  FILE *fp=NULL;
  int  i,j,id,len=2014;
  float onset=0,duration=0,height=0;
  char *buf = (char *)VCalloc(len,sizeof(char));

  fp = fopen(designfile,"r");
  if (!fp) VError(" error opening design file %s",designfile);
  int ntrials = 0;
  while (!feof(fp)) {
    for (j=0; j<len; j++) buf[j] = '\0';
    if (fgets(buf,len,fp) == NULL) break;
    if (strlen(buf) < 2) continue;
    if (buf[0] == '%' || buf[0] == '#') continue;
    if (! test_ascii((int)buf[0])) VError(" input file must be a text file");
    ntrials++;
  }
  *numtrials = ntrials;

  Trial *trial = (Trial *) VCalloc(ntrials+1,sizeof(Trial));

  i = (*nevents) = 0;
  rewind(fp);
  while (!feof(fp)) {
    for (j=0; j<len; j++) buf[j] = '\0';
    if (fgets(buf,len,fp) == NULL) break;
    if (strlen(buf) < 2) continue;
    if (buf[0] == '%' || buf[0] == '#') continue;
    if (! test_ascii((int)buf[0])) VError(" input file must be a text file");

    /* remove non-alphanumeric characters */
    for (j=0; j<strlen(buf); j++) {
      if (buf[j] == '\v') buf[j] = ' '; /* remove tabs */
      if (buf[j] == '\t') buf[j] = ' ';
    }
    if (sscanf(buf,"%d %f %f %f",&id,&onset,&duration,&height) != 4)
      VError(" line %d: illegal input format",i+1);
    if (duration < 0.5 && duration >= -0.0001) duration = 0.5;
    trial[i].id       = id;
    trial[i].onset    = onset;
    trial[i].duration = duration;
    trial[i].height   = height;
    i++;
    if (id > (*nevents)) (*nevents) = id;
  }
  (*nevents) = (*nevents) + 1;
  fclose(fp);
  return trial;
}



Trial *CopyTrials(Trial *trial,int numtrials)
{
  int i;
  Trial *newtrial = (Trial *) VCalloc(numtrials,sizeof(Trial));
  for (i=0; i<numtrials; i++) {
    newtrial[i].id       = trial[i].id;
    newtrial[i].onset    = trial[i].onset;
    newtrial[i].duration = trial[i].duration;
    newtrial[i].height   = trial[i].height;
  }
  return newtrial;
}




/* concatenate trials from all runs */
Trial *ConcatenateTrials(Trial **trial,int *numtrials,float *run_duration,int dlists,int sumtrials)
{
  int i,j;
  Trial *alltrials = (Trial *) VCalloc(sumtrials,sizeof(Trial));
  
  int ii=0;
  float add=0;
  for (i=0; i<dlists; i++) {
    for (j=0; j<numtrials[i]; j++) {
      alltrials[ii].id       = trial[i][j].id;
      alltrials[ii].onset    = trial[i][j].onset + add;
      alltrials[ii].duration = trial[i][j].duration;
      alltrials[ii].height   = trial[i][j].height;
      ii++;
    }
    add += run_duration[i];
  }
  return alltrials;
}

