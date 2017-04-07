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
  FILE *in_file = NULL, *out_file = NULL;
  VAttrList list = NULL;
  char *prg=GetLipsiaName("vdetrend");
  fprintf(stderr, "%s\n", prg);

  VParseFilterCmd(VNumber(options), options, argc, argv, &in_file, &out_file);
 

  /* read the file */
  if(!(list = VReadFile(in_file, NULL))) exit(1);
  fclose(in_file);

  /* perform detrending */
  VDetrend(list,minval,type,del);


  /* output */
  VHistory(VNumber(options),options,prg,&list,&list);
  if(! VWriteFile(out_file, list)) exit(1);
  fprintf(stderr, "%s: done.\n", argv[0]);
  exit(0);
}
