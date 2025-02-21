/*
** filter: most freqent value in a local nieghbourhood
**
** G.Lohmann, Jan 2016
*/

/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern VImage Convert2Repn(VImage src,VImage dest,VRepnKind repn);


int main (int argc, char **argv)
{
  static VString rim_filename="";
  static VOptionDescRec options[] = {
   {"rim", VStringRepn,1,(VPointer) &rim_filename,VRequiredOpt,NULL,"rim image"},
  };
  VString in_filename=NULL;
  FILE *in_file, *out_file;
  char *prg=GetLipsiaName("vdelrim");
  fprintf (stderr, "%s\n", prg);

  
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);

  /* read input qimage */
  VAttrList list = VReadAttrList(in_filename,0L,TRUE,FALSE);
  VImage src = VReadImage(list);
  if (src == NULL) VError(" no src found");

  /* rim image */
  VImage rim=NULL;
  VAttrList rlist = VReadAttrList(rim_filename,0L,TRUE,FALSE);
  if (rlist == NULL) VError(" error reading file %s",rim_filename);
  VImage rtmp = VReadImage(rlist);
  if (rtmp == NULL) VError(" err reading %s",rim_filename);
  if (VImageNPixels(rtmp) != VImageNPixels(src))
    VError(" inconsistent image dimensions: rim vs src");
  rim = Convert2Repn(rtmp,rim,VUByteRepn);
  VDestroyImage(rtmp);

  int b,r,c,k;
  for (b=0; b<VImageNBands(src); b++) {
    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {
	k = VPixel(rim,b,r,c,VUByte);
	if (k != 3) VSetPixel(src,b,r,c,0);
      }
    }
  }

  /* output */
  if (! VWriteFile (out_file, list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
