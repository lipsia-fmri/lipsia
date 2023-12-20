/*
** read preprocessed data into a data format compatible with the SML-algorithm
** The input data must be in csv-format, the output is in vista-format.
** The parameter '-ncomponents' specifies the number of parcels/components of the atlas/decomposition used.
**
** Example:
**   vreadconnectome -in connectome_subj007.csv -ncomponents 100 -out connectome_subj007.v
**
**  G.Lohmann, MPI-KYB, 2023
*/

#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>
#include <viaio/VImage.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc,char *argv[])
{
  static VString in_filename = "";
  static VLong  xcomponents = 100;
  static VOptionDescRec  options[] = {
    {"in",VStringRepn,1,(VPointer) &in_filename,VRequiredOpt,NULL,"Input file"},
    {"ncomponents",VLongRepn,1,(VPointer) &xcomponents,VRequiredOpt,NULL,"components"},
  };
  FILE *out_file;
  size_t i;
  float u=0;

  
  /* parse command line */
  VParseFilterCmd (VNumber (options),options,argc,argv,NULL,&out_file);

  size_t ncomponents = (size_t) xcomponents;  
  size_t dim = ncomponents*(ncomponents-1)/2;
  size_t len  = 30*dim;
  char *buf  = (char *) VCalloc(len,sizeof(char));
  char *token = "";


  /*
  FILE *fq = fopen(in_filename,"w");
  for (i=0; i<dim-1; i++) {
    fprintf(fq,"%.3f,",(float)i);
  }
  fprintf(fq,"%.3f\n",(float)i);
  fclose(fq);
  exit(0);
  */

  
  FILE *fp = fopen(in_filename,"r");
  if (!fp) VError(" err opening %s",in_filename);
  memset(buf,0,len);
  if (fgets(buf,len,fp) == FALSE) VError(" err reading input file");

  float *A = (float *)VCalloc(dim,sizeof(float));
  token = strtok(buf,",");
  for (i=0; i<dim; i++) {
    if (token == NULL) VError(" file has only %lu entries. It should have %lu entries",i,dim);
    sscanf(token,"%f",&u);
    token = strtok(NULL,",");
    A[i] = u;
  }
  fclose(fp);


  VAttrList alist = VCreateAttrList();
  VBundle corr = VCreateBundle ("data",alist,dim*sizeof(float),(VPointer)A);
  VAttrList out_list = VCreateAttrList();
  VAppendAttr(out_list,"nedges",NULL,VLongRepn,(VLong)dim);
  VAppendAttr(out_list,"Correlation",NULL,VBundleRepn,(VBundle)corr);
  VAppendAttr(out_list,"ncomponents",NULL,VLongRepn,(VLong)ncomponents);   
  if (! VWriteFile (out_file, out_list)) exit (1);
  
  fprintf (stderr, "%s: done.\n", argv[0]);
  exit(0);
}
