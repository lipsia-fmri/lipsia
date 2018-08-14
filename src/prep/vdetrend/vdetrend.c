#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

extern void VDetrend(VAttrList list,VFloat minval,VShort type,VShort del);

VDictEntry TypeDict[] = {
  { "linear", 0 },
  { "cubic", 1 },
  { "order5", 2 },
  { NULL }
};


int main(int argc, char *argv[]) 
{
  static VShort del = 1;
  static VShort type = 0;
  static VFloat minval = 0;
  static VOptionDescRec  options[] = {
    {"del", VShortRepn, 1, (VPointer) &del, VOptionalOpt, NULL, "Number of initial timepoints to ignore"},
    {"type", VShortRepn, 1, (VPointer) &type, VOptionalOpt,TypeDict, "Type of detrending"},
    {"minval", VFloatRepn, 1, (VPointer) &minval, VOptionalOpt, NULL, "Signal threshold"}
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
  VDetrend(list,minval,type,del);


  /* output */
  VHistory(VNumber(options),options,prg,&list,&list);
  if(! VWriteFile(out_file, list)) exit(1);
  fprintf(stderr, "%s: done.\n", argv[0]);
  exit(0);
}
