/*
** Detrending
**
** G.Lohmann, 2018, MPI-KYB
*/
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

extern void VSubtractPolynomial(VAttrList list,VShort type,VShort i0);
extern void VSubtractMovingAverage(VAttrList list,float window,int);
extern void VDetrend(VAttrList list,VShort type,VShort del);

VDictEntry TypeDict[] = {
  { "demean", 0 },
  { "linear", 1 },
  { "cubic", 2 },
  { "order5", 3 },
  { NULL }
};


int main(int argc, char *argv[]) 
{
  static VShort del = 0;
  static VShort type = 0;  
  static VBoolean linear=TRUE;
  static VFloat window = 100;
  static VOptionDescRec  options[] = {
    {"type", VShortRepn, 1, (VPointer) &type, VOptionalOpt,TypeDict, "Type of detrending"},   
    {"linear", VBooleanRepn, 1, (VPointer) &linear, VOptionalOpt,NULL, "Whether to first subtract linear drift"},   
    {"window", VFloatRepn, 1, (VPointer) &window, VOptionalOpt,NULL, "Window size in seconds (only for demeaning)"},
    {"del", VShortRepn, 1, (VPointer) &del, VOptionalOpt, NULL, "Number of initial timepoints to ignore"},
  };
  VString in_file=NULL;
  FILE *out_file = NULL;
  char *prg=GetLipsiaName("vdetrend");
  fprintf(stderr, "%s\n", prg);

  VParseFilterCmdX(VNumber(options), options, argc, argv, &in_file, &out_file);
 

  /* read the file */ 
  VLong tr=0L;
  VAttrList list = VReadAttrList(in_file,tr,TRUE,FALSE);
  if (list == NULL) VError(" error reading input file %s",in_file);


  /* perform detrending */
  if (linear && type != 1) VSubtractPolynomial(list,(int)1,(int)del);

  if (type == 0) {
    VSubtractMovingAverage(list,(float)window,(int)del);
  }
  else {
    VSubtractPolynomial(list,type,(int)del);
  }

  /* output */
  VHistory(VNumber(options),options,prg,&list,&list);
  if(! VWriteFile(out_file, list)) exit(1);
  fprintf(stderr, "%s: done.\n", argv[0]);
  exit(0);
}
