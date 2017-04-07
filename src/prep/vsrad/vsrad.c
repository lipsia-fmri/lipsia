/*
** Speckle reducing anisotropic diffusion
**
** Lit: Y.Yu, S.T. Acton, 
**  IEEE Trans Image Proc., Vol. 11, No. 11, Nov 2002
**
**
**  G.Lohmann, MPI-CBS, 2008
*/


#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern VImage VSRAD(VImage,VImage,VShort,VShort,VFloat);


void
DeleteBorder(VImage src)
{
  int b,r,c;

  r = VImageNRows(src)-1;
  for (b=0; b<VImageNBands(src); b++) {
    for (c=0; c<VImageNColumns(src); c++) {
      VSetPixel(src,b,0,c,0);
      VSetPixel(src,b,r,c,0);
    }
  }

  c = VImageNColumns(src)-1;
  for (b=0; b<VImageNBands(src); b++) {
    for (r=0; r<VImageNRows(src); r++) {
      VSetPixel(src,b,r,0,0);
      VSetPixel(src,b,r,c,0);
    }
  }
}



int 
main (int argc,char *argv[])
{
  static VShort numiter = 7;
  static VShort type    = 0;
  static VFloat rho     = 0.1666;
  static VOptionDescRec  options[] = {
    {"iter",VShortRepn,1,(VPointer) &numiter,
       VOptionalOpt,NULL,"number of iterations"},
    {"type",VShortRepn,1,(VPointer) &type,
       VOptionalOpt,NULL,"type of diffusion function (0 or 1)"},
    {"rho",VFloatRepn,1,(VPointer) &rho,
       VOptionalOpt,NULL,"speckle parameter"}
  };
  FILE *in_file,*out_file;
  VAttrList list=NULL;
  VAttrListPosn posn;
  VImage src=NULL,dest=NULL;
  int n,nimages;
  char *prg=GetLipsiaName("vsrad");
  fprintf (stderr, "%s\n", prg);


  VParseFilterCmd (VNumber (options),options,argc,argv,&in_file,&out_file);

  if (rho < 0 || rho > 1) VError(" rho must be in [0,1]");

  if (! (list = VReadFile (in_file, NULL))) exit (1);
  fclose(in_file);

  nimages = 0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) == VImageRepn) nimages++;
  }

  n = 0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src);
    fprintf(stderr," image %4d of %d\r",n,nimages);
    dest = VSRAD(src,dest,numiter,type,rho);
    DeleteBorder(dest);
    src = VCopyImage(dest,src,VAllBands);
    VSetAttrValue (& posn, NULL,VImageRepn,src);
    n++;
  }

  VHistory(VNumber(options),options,prg,&list,&list);
  if (! VWriteFile (out_file, list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
