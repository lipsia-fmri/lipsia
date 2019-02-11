/*
**  parse second level design file
**
** G.Lohmann
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

extern int VStringToken(char *, char *, int, int);

#define LEN  10000


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
  if (val == '-') return 1;
  if (val == '+') return 1;
  if (val == '.') return 1;
  return 0;
}



int line_empty(char *buf,int len)
{
  int i;
  for (i=0; i<len; i++) {
    if (buf[i] != ' ' && buf[i] != '\n' && buf[i] != 0) return 0;
  }
  return 1;
}


int CheckBuffer(char *buf,int len)
{
  int j;

  if(strlen(buf) < 1) return 0;
  if (buf[0] == '%' || buf[0] == '#' || buf[0] == '/' || buf[0] == '\n') return 0;

  /* remove tabs */
  for (j=0; j<len; j++) {
    if (buf[0] == '\t') buf[j] = ' ';
  }
  if (line_empty(buf,len) > 0) return 0;
  return 1;
}


int VistaFormat(char *buf,int len)
{
  if (strncmp(buf,"V-data",6) == 0) return 1;
  return 0;
}



/*
** read a 2nd level design file
*/
gsl_matrix *XRead2ndLevel(VString filename) 
{
  int i, j, nrows, ncols;
  double val;
  char buf[LEN], token[32];

  fprintf(stderr," Reading %s\n",filename);

  FILE *fp = VOpenInputFile (filename, TRUE);
  if (!fp) VError(" error opening %s",filename);

  nrows = ncols = 0;
  while(!feof(fp)) {
    memset(buf, 0, LEN);
    if(!fgets(buf, LEN, fp)) continue;
    if (VistaFormat(buf,LEN) > 0) VError(" Design file must be a text file");
    if (CheckBuffer(buf,LEN) < 1) continue;
    if (! test_ascii((int)buf[0])) VError(" Design file must be a text file");

    j = 0;
    while(VStringToken(buf, token, j, 30)) {
      if(!sscanf(token, "%lf", &val))
	VError("illegal text string in design file: %s", buf);
      j++;
    }

    if (ncols == 0) ncols = j;
    if (j < 1) continue;
    else if(ncols != j)
      VError(" inconsistent number of columns in row %d", nrows + 1);
    nrows++;
  }
  rewind(fp);


  /* fill design matrix */
  gsl_matrix *X = gsl_matrix_calloc(nrows,ncols);

  i = 0;
  while(!feof(fp)) {
    memset(buf, 0, LEN);
    if(!fgets(buf, LEN, fp))  continue;    
    if (CheckBuffer(buf,LEN) < 1) continue;

    j = 0;
    while(VStringToken(buf, token, j, 30)) {
      sscanf(token, "%lf", &val);
      gsl_matrix_set(X,i,j,val);
      j++;
    }
    if (j < 1) continue;
    i++;
  }
  fclose(fp);

  return X;
}


/* read exchangeability file */
void XReadExchange(VString filename,int *exchange,int n)
{
  int len=1024;
  char token[32];
  char *buf = (char *)VCalloc(len,sizeof(char));
  FILE *fp = VOpenInputFile (filename, TRUE);
  if (!fp) VError(" error opening %s",filename);
  fprintf(stderr," Reading %s\n",filename);

  int nrows = 0;
  while(!feof(fp)) {
    memset(buf, 0, len);
    if(!fgets(buf, len, fp)) break;
    if (VistaFormat(buf,LEN) > 0) VError(" Exchangeability file must be a text file");
    if (CheckBuffer(buf,LEN) < 1) continue;
    if (! test_ascii((int)buf[0])) VError(" Exchangeability file must be a text file");
    nrows++;    
  }
  rewind(fp);

  if (nrows != n) 
    VError(" Number of rows in exchangeability file (%d) does not match number of rows in design file (%d)",nrows,n);

  int i=0,j=0,k=0;
  while(!feof(fp)) {
    memset(buf, 0, len);
    if(!fgets(buf, len, fp))  break;
    if (CheckBuffer(buf,LEN) < 1) continue;

    k = 0;
    while(VStringToken(buf, token, k, 30)) {
      if (!sscanf(token, "%d", &j))
	VError("illegal text string in exchangeability file: %s", buf);
      if (j < 0) VError(" Values in exchangeability file must be non-negative");
      if (j > n) VError(" Illegal value (%d) in exchangeability file, values must be <= %d",j,n);
      exchange[i] = j;
      k++;
    }
    if (k > 1) VError(" Exchangeability file must have exactly one column");
    i++;
  }
  VFree(buf);
  fclose(fp);
}

/*
int main()
{
  int i,j;
  gsl_matrix *X = XRead2ndLevel("bla.mat");


  fprintf(stderr,"X:\n");
  for (i=0; i<X->size1; i++) {
    for (j=0; j<X->size2; j++) {
      fprintf(stderr," %.2f",gsl_matrix_get(X,i,j));
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
  exit(0);
}
*/
