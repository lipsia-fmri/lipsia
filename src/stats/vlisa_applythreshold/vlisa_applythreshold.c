/*
** apply FDR threshold to a map produced by LISA statistical inference
** This program is useful if no FDR-threshold was applied in lisa inference (-alpha 1.0)
**
** G.Lohmann, MPI-KYB, Nov 2-18
*/

/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/option.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void VApplyThreshold(VImage src,float alpha)
{
  size_t i=0,npos=0;
  float beta=1.0-alpha;

  /* apply threshold */
  VFloat *pp = VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    float u = *pp;
    if (u < beta) (*pp) = 0;
    else npos++;
    pp++;
  }
  fprintf(stderr," number of voxels passing threshold %.3f:  %lu\n",alpha,npos);
}


int main(int argc, char *argv[]) 
{
  static VFloat alpha = 0.05;
  static VOptionDescRec  options[] = {
    {"threshold", VFloatRepn, 1, (VPointer) &alpha, VOptionalOpt, NULL, "FDR threshold"},
  };
  FILE *in_file, *out_file;
  VAttrList list = NULL;
  VAttrListPosn posn;
  VImage src=NULL;
  char *prg_name = GetLipsiaName("vlisa_applythreshold");
  fprintf (stderr, "%s\n", prg_name);


  VParseFilterCmd(VNumber(options), options, argc, argv, &in_file, &out_file);
  if (alpha < 0 || alpha > 1) VError(" parameter 'alpha' must be in range [0,1]");
  
  if(!(list = VReadFile(in_file, NULL))) exit(1);
  fclose(in_file);

  for(VFirstAttr(list, & posn); VAttrExists(& posn); VNextAttr(& posn)) {
    if(VGetAttrRepn(& posn) != VImageRepn) continue;
    VGetAttrValue(& posn, NULL, VImageRepn, & src);
    if(VPixelRepn(src) != VFloatRepn) continue;
    VApplyThreshold(src,(float)alpha);
  }
  if (src == NULL) VError(" no input image found");

  VHistory(VNumber(options), options, prg_name, &list, &list);
  if(! VWriteFile(out_file, list)) exit(1);
  fprintf(stderr, "%s: done.\n", argv[0]);
  return 0;
}
