/*
** Read design file
**
** G.Lohmann, MPI-KYB 2018
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


extern void printmat(gsl_matrix *R,char *str);
extern void printvec(gsl_vector *x,char *str);

int test_ascii(int val)
{
  if (val >= 'a' && val <= 'z') return 1;
  if (val >= 'A' && val <= 'Z') return 1;
  if (val >= '0' && val <= '9') return 1;
  if (val ==  ' ') return 1;
  if (val ==  '-') return 1;
  if (val ==  '+') return 1;
  if (val ==  '.') return 1;
  if (val == '\0') return 1;
  if (val == '\n') return 1;
  if (val == '\r') return 1;
  if (val == '\t') return 1;
  if (val == '\v') return 1;
  return 0;
}


int line_empty(char *buf,int len)
{
  int i;
  for (i=0; i<len; i++) {
    if (buf[i] != ' ' && buf[i] != '\0' && buf[i] != '\n') return 0;
  }
  return 1;
}

int CheckBuffer(char *buf,int len)
{
  int j;
  if(strlen(buf) < 1) return 0;
  if (buf[0] == '%' || buf[0] == '#' || buf[0] == '/') return 0;  /* comment */
  for (j=0; j<len; j++) {
    if (buf[j] == '\v') buf[j] = ' '; /* remove tabs */
    if (buf[j] == '\t') buf[j] = ' ';
  }
  if (line_empty(buf,len) > 0) return 0;
  return 1;
}


int VistaFormat(char *buf,int len)
{
  if (strncmp(buf,"V-data",6) == 0) return 1;
  return 0;
}


/* parse design file */
Trial *ReadDesign(VString designfile,int *numtrials,int *nevents)
{
  FILE *fp=NULL;
  int  i,j,id,len=4096;
  float onset=0,duration=0,height=0;
  char *buf = (char *)VCalloc(len,sizeof(char));

  fp = fopen(designfile,"r");
  if (!fp) VError(" error opening design file %s",designfile);
  int ntrials=0;
  while (!feof(fp)) {
    for (j=0; j<len; j++) buf[j] = '\0';
    if (fgets(buf,len,fp) == NULL) break;
    if (VistaFormat(buf,len) > 1) VError(" Design file must be a text file");
    if (CheckBuffer(buf,len) < 1) continue;
    if (! test_ascii((int)buf[0])) VError(" File %s: line %d begins with an illegal character (%c)",designfile,ntrials,buf[0]);
    ntrials++;    
  }
  *numtrials = ntrials;

  Trial *trial = (Trial *) VCalloc(ntrials+1,sizeof(Trial));

  i = (*nevents) = 0;
  rewind(fp);
  while (!feof(fp)) {
    for (j=0; j<len; j++) buf[j] = '\0';
    if (fgets(buf,len,fp) == NULL) break;
    if (CheckBuffer(buf,len) < 1) continue;

    int nch = sscanf(buf,"%d %f %f %f",&id,&onset,&duration,&height);
    if (nch > 0 && nch != 4) VError(" line %d: illegal input format",i+1);
    if (duration < 0.5 && duration >= -0.0001) duration = 0.5;
    trial[i].id       = id;
    trial[i].onset    = onset;
    trial[i].duration = duration;
    trial[i].height   = height;
    i++;
    if (id > (*nevents)) (*nevents) = id;
  }
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



/* check if trial labels are complete and positive */
void CheckTrialLabels(Trial *trial,int numtrials)
{
  int i,j;
  int maxlabel = -1;
  for (i=0; i<numtrials; i++) {
    j = trial[i].id;
    if (j > maxlabel) maxlabel = j;
    if (j < 1) VError(" trial[%d] has illegal label %d",i,j);
  }

  int *tmp = (int *) VCalloc(maxlabel+1,sizeof(int));
  for (i=0; i<numtrials; i++) {
    j = trial[i].id;
    tmp[j]++;
  }
  for (i=1; i<=maxlabel; i++) {
    if (tmp[i] < 1) VError(" Label %d missing in design file",i);
    if (tmp[i] < 4) VWarning(" Label %d has too few trials (%d) for random permutations",i,tmp[i]);
  }
  VFree(tmp);
}

