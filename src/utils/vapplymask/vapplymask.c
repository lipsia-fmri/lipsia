/*
** apply mask
**
** G.Lohmann, MPI-KYB, May 2018
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "viaio/Vlib.h"
#include "viaio/VImage.h"
#include "viaio/mu.h"
#include "viaio/option.h"



void ApplyMask3d(VImage src,VImage mask)
{
  int b,r,c;
  double u;
  int nslices = VImageNBands(src);
  int nrows = VImageNRows(src);
  int ncols = VImageNColumns(src);
  
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	u = VGetPixel(mask,b,r,c);
	if (u < TINY) VSetPixel(src,b,r,c,0.0);
      }
    }
  }
}


void ApplyMask4d(VAttrList list,VImage mask)
{
  int b,r,c,j;
  double u;

  int nrows=0,ncols=0,nt=0;
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  VImageDimensions(src,nslices,&nt,&nrows,&ncols);

  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	u = VGetPixel(mask,b,r,c);
	if (u < TINY) {
	  for (j=0; j<nt; j++) {
	    VSetPixel(src[b],j,r,c,0.0);
	  }
	}
      }
    }
  }
}



int main (int argc,char *argv[])
{ 
  static VString filename="";
  static VOptionDescRec  options[] = {   
    {"mask", VStringRepn,1,(VPointer) &filename,VRequiredOpt,NULL,"Mask file"},
  };
  FILE *out_file=NULL;
  VString in_file=NULL;
  char *prg = GetLipsiaName("vapplymask");


  /* Parse command line arguments and identify files: */
  VParseFilterCmdX (VNumber (options), options, argc,argv,&in_file,&out_file);


  /* Read mask image */
  VAttrList mlist = VReadAttrList(filename,0L,TRUE,FALSE);
  if (mlist == NULL) VError(" error reading mask file");
  VImage mask = VReadImage(mlist);


  /* read image data */
  VAttrList list = VReadAttrList(in_file,0L,TRUE,FALSE);
  if (list == NULL) VError(" error reading input file %s",in_file);
  int nimages = VAttrListNumImages(list);



  /* main process */
  if (nimages < 2) {  
    VImage src = VReadImage(list);  
    ApplyMask3d(src,mask);
  }
  else {
    ApplyMask4d(list,mask);
  }


  /* output */
  VHistory(VNumber(options),options,prg,&list,&list);
  VWriteFile (out_file,list);
  fclose(out_file);
  fprintf (stderr, "%s: done.\n", argv[0]);
  exit(0);
}
