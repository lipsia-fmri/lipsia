/*
** get cortical rim  using a segmentation into GM/WM/CSF
**
** G.Lohmann, MPI-KYB, May 2025
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


#define GM_CSF 1
#define GM_WM 2
#define GM_INSIDE 3


VImage VRim(VImage src,int gm,int wm,int csf)
{ 
  int j,k,n0,n1;
  int b,r,c,bb,rr,cc;
  int nslices = VImageNBands(src);
  int nrows  = VImageNRows(src);
  int ncols  = VImageNColumns(src);

  fprintf(stderr," gm: %d,  wm: %d,  csf: %d\n",gm,wm,csf);

  VImage dest = VCreateImage(nslices,nrows,ncols,VUByteRepn);
  VFillImage(dest,VAllBands,0);
  VCopyImageAttrs (src, dest);

  for (b=1; b<nslices-1; b++) {
    for (r=1; r<nrows-1; r++) {
      for (c=1; c<ncols-1; c++) {

	k = (int)VGetPixel(src,b,r,c);
	if (k != gm) continue;
	
	n0=n1=0;
	for (bb=b-1; bb<=b+1; bb++) {
	  for (rr=r-1; rr<=r+1; rr++) {
	    for (cc=c-1; cc<=c+1; cc++) {

	      /* 6-adjacency */
	      if (fabs(b-bb) + fabs(r-rr) + fabs(c-cc) > 1) continue;
	      
	      j = (int)VGetPixel(src,bb,rr,cc);
	      if (j == csf) n0++;
	      if (j == wm) n1++;
	    }
	  }
	}
	k = GM_INSIDE;
	if (n0 > 0) k = GM_CSF;
	if (n1 > n0) k = GM_WM;
	VPixel(dest,b,r,c,VUByte) = k;
      }
    }
  }
  return dest;
}


int main (int argc, char **argv)
{
  static VShort gm=2;
  static VShort wm=1;
  static VShort csf=0;
  static VOptionDescRec options[] = {
    {"gm", VShortRepn,1,(VPointer) &gm,VOptionalOpt,NULL,"GM"},
    {"wm", VShortRepn,1,(VPointer) &wm,VOptionalOpt,NULL,"WM"},
    {"csf", VShortRepn,1,(VPointer) &csf,VOptionalOpt,NULL,"CSF"},
  };
  VString in_filename=NULL;
  FILE *in_file, *out_file;
  char *prg=GetLipsiaName("vrim");
  fprintf (stderr, "%s\n", prg);

  
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);

  
  /* read image */
  VAttrList list = VReadAttrList(in_filename,0L,TRUE,FALSE);
  VImage src = VReadImage(list);
  if (src == NULL) VError(" no src found");
  VAttrList geolist = VGetGeoInfo(list);

  /* mark cortical rim */
  VImage dest = VRim(src,(int)gm,(int)wm,(int)csf);


  /* output */
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,dest);
  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
