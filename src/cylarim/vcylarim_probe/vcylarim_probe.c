/*
** cylarim: cylinder-based laminar fMRI analysis
**
** G.Lohmann, MPI-KYB, Nov 2024
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


float Conjunction(float u,float v)
{
  float z=0;
  if (u*v > 0) {
    if (u > 0 && v > 0) {
      if (v < u) z = v;
      else z = u;
    }
    if (u < 0 && v < 0) {
      if (v > u) z = v;
      else z = u;
    }
  }
  return z;
}

void VImageConjunction(VImage src1,VImage src2,VImage dest)
{
  size_t i;
  size_t n=VImageNPixels(src1);
  VFloat *p1 = VImageData(src1);
  VFloat *p2 = VImageData(src2);
  VFloat *pd = VImageData(dest);
  for (i=0; i<n; i++) {
    pd[i] = Conjunction(p1[i],p2[i]);
  }
}

void VImageDiff(VImage src1,VImage src2,VImage dest)
{
  size_t i;
  size_t n=VImageNPixels(src1);
  VFloat *p1 = VImageData(src1);
  VFloat *p2 = VImageData(src2);
  VFloat *pd = VImageData(dest);
  for (i=0; i<n; i++) {
    pd[i] = (p1[i]-p2[i]);
  }
}


void VImageMaxIndex(VImage *src,int nimages,VImage dest)
{
  int b,r,c,j,k;
  double u,xmax,s=0;
  int nslices = VImageNBands(dest);
  int nrows   = VImageNRows(dest);
  int ncols   = VImageNColumns(dest);

  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	k=-1;
	xmax = -99999;
	s=0;
	for (j=0; j<nimages; j++) {
	  u = VGetPixel(src[b],j,r,c);
	  s += fabs(u);
	  if (u > xmax) {
	    xmax = u;
	    k = j+1;
	  }
	}
	if (k > 0 && s > TINY) VSetPixel(dest,b,r,c,(float)(k));
      }
    }
  }
}

void VImageMinIndex(VImage *src,int nimages,VImage dest)
{
  int b,r,c,j,k;
  double u,xmin,s;
  int nslices = VImageNBands(dest);
  int nrows   = VImageNRows(dest);
  int ncols   = VImageNColumns(dest);

  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	k=-1;
	s=0;
	xmin = 99999;
	for (j=0; j<nimages; j++) {
	  u = VGetPixel(src[b],j,r,c);
	  s += fabs(u);
	  if (u < xmin) {
	    xmin = u;
	    k = j+1;
	  }
	}
	if (k > 0 && s > TINY) VSetPixel(dest,b,r,c,(float)(k));
      }
    }
  }
}



void VImageInvert(VImage src)
{
  size_t i;
  size_t n=VImageNPixels(src);
  VFloat *pp = VImageData(src);
  for (i=0; i<n; i++) {
    pp[i] = -pp[i];
  }
}

void VImageZMax(VImage *src,int nimages,VImage dest)
{
  int j,b,r,c;
  int nslices = VImageNBands(dest);
  int nrows   = VImageNRows(dest);
  int ncols   = VImageNColumns(dest);
  float u=0,z=0,zpos,zneg;
 
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	z=zpos=zneg=0;
	for (j=0; j<nimages; j++) {
	  u = VGetPixel(src[b],j,r,c);
	  if (u > zpos) zpos = u;
	  if (u < zneg) zneg = u;
	}
	z = zpos;
	if (fabs(zneg) > z) z = zneg;
	VSetPixel(dest,b,r,c,z);
      }
    }
  }
}

void SelectSlice(VImage *src,int slice,VImage dest)
{
  int b,r,c;
  double u;
  int nslices = VImageNBands(dest);
  int nrows   = VImageNRows(dest);
  int ncols   = VImageNColumns(dest);
  if (slice >= nslices) VError(" slice %d too big, max is: %d",slice,nslices);
  if (slice < 0) VError(" slice id must be non-negative");
  
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	u = VGetPixel(src[b],slice,r,c);
	VSetPixel(dest,b,r,c,u);
      }
    }
  }
}


VDictEntry TypDict[] = {
  { "volume", 0, 0,0,0,0 },
  { "top_d", 1, 0,0,0,0 },
  { "top_m", 2, 0,0,0,0 },
  { "top_s", 3, 0,0,0,0 },
  { "bottom_d", 4, 0,0,0,0 },
  { "bottom_m", 5, 0,0,0,0 },
  { "bottom_s", 6, 0,0,0,0 },
  { "max_id", 7, 0,0,0,0 },
  { "min_id", 8, 0,0,0,0 },
  { "maxabs", 9, 0,0,0,0 },
  { NULL, 0,0,0,0,0 }
};

int main(int argc, char *argv[])
{
  static VShort slice = 0;
  static VShort type = 0;
  static VShort xinvert = FALSE;
  static VOptionDescRec  options[] = {
    {"volume", VShortRepn,1,(VPointer) &slice,VOptionalOpt,NULL,"Volume to be selected"},
    {"xinvert", VBooleanRepn,1,(VPointer) &xinvert,VOptionalOpt,NULL,"Whether to invert the selected volume"},
    {"type", VShortRepn,1,(VPointer) &type,VOptionalOpt,TypDict,"Output type"},
  };
  FILE *out_file=NULL;
  VString in_file=NULL;
  char *prg=GetLipsiaName("vcylarim_probe");
  fprintf (stderr, " %s\n", prg);

  
  /* parse command line */
  VParseFilterCmdX(VNumber(options),options,argc,argv,&in_file,&out_file);

  
  /* read the file */
  VLong tr=0L;
  VAttrList list = VReadAttrList(in_file,tr,TRUE,FALSE);
  if (list == NULL) VError(" error reading input file %s",in_file);

  
  /* get name */
  VString cbuf;
  /*
    int i,j;
  char cbuf1[30],cbuf2[30];
  memset(cbuf1,0,sizeof(cbuf1));
  memset(cbuf2,0,sizeof(cbuf2));
  ntype=0,stype=0;
  */
  if (VGetAttr (list, "name", NULL,VStringRepn, (VPointer) & cbuf) == VAttrFound) {
    fprintf(stderr," name: %s\n",cbuf);
    if (strcmp(cbuf,"zabs")==0)  VError(" nothing to be done on a 'zabs' image");
    /*
    i=0;
    while (cbuf[i] != '_') { cbuf1[i]=cbuf[i]; i++; }
    i++;
    j=0;
    while (i<strlen(cbuf)) {
      cbuf2[j++]=cbuf[i++];
    }
    if (strcmp(cbuf1,"3bins")==0) ntype = 0;
    else if (strcmp(cbuf1,"peak")==0) ntype = 1;
    else if (strcmp(cbuf1,"shape")==0) ntype = 2;
    else if (strcmp(cbuf1,"R2")==0) ntype = 3;

    if (strcmp(cbuf2,"coeffs")==0) stype = 0;
    else if (strcmp(cbuf2,"stats")==0) stype = 1;
    */
  }
  if (type > 1 && type < 7 && strcmp(cbuf,"3bins_zvals") != 0 )
    VError(" type '%s' only available for '3bins_zvals', not for '%s' ",TypDict[type].keyword,cbuf);
  if (type == 9 && strcmp(cbuf,"3bins_coeff") != 0 )
    VError(" type '%s' only available for '3bins_coeff', not for '%s' ",TypDict[type].keyword,cbuf);

  
  /* ini functional data struct  */
  int nt,nrows,ncols;
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  VImageDimensions(src,nslices,&nt,&nrows,&ncols);  
  int nimages=nt;

  /* ini dest image */
  VImage dest = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VCopyImageAttrs (src[0],dest);
  VFillImage(dest,VAllBands,0);
  VImage atmp = VCreateImageLike(dest);
  VImage btmp = VCreateImageLike(dest);
  VImage ctmp = VCreateImageLike(dest);
  VFillImage(atmp,VAllBands,0);
  VFillImage(btmp,VAllBands,0);
  VFillImage(ctmp,VAllBands,0);

  /* cases */
  switch (type) {
  case 0:
    if (slice < 0 || slice >= nimages) VError(" slice %d not available in '%s' ",slice,cbuf);
    SelectSlice(src,(int)slice,dest);
    if (xinvert) {
      VImageInvert(dest);
      fprintf(stderr," invert selected slice %d\n",slice);
    }
    break;
    
  case 1:   /* top_d */
    SelectSlice(src,(int)0,atmp);
    SelectSlice(src,(int)1,btmp);
    VImageConjunction(atmp,btmp,dest);
    break;
    
  case 2:   /* top_m */
    SelectSlice(src,(int)0,atmp);
    SelectSlice(src,(int)2,btmp);
    VImageInvert(atmp);
    VImageConjunction(atmp,btmp,dest);
    break;
    
  case 3:    /* top_s */
    SelectSlice(src,(int)1,atmp);
    SelectSlice(src,(int)2,btmp);
    VImageInvert(atmp);
    VImageInvert(btmp);
    VImageConjunction(atmp,btmp,dest);
    break;

  case 4:   /* bottom_d, deep is smallest */
    SelectSlice(src,(int)0,atmp);
    VImageInvert(atmp);
    SelectSlice(src,(int)1,btmp);
    VImageInvert(btmp);
    VImageConjunction(atmp,btmp,dest);
    break;

  case 5:   /* bottom_m, middle is smallest */
    SelectSlice(src,(int)0,atmp);
    SelectSlice(src,(int)2,btmp);
    VImageInvert(btmp);
    VImageConjunction(atmp,btmp,dest);
    break;

  case 6:   /* bottom_s, superficial is smallest */
    SelectSlice(src,(int)1,atmp);
    SelectSlice(src,(int)2,btmp);
    VImageConjunction(atmp,btmp,dest);
    break;
  
  case 7:
    VImageMaxIndex(src,nimages,dest);
    break;
  case 8:
    VImageMinIndex(src,nimages,dest);
    break;
  case 9:
    VImageZMax(src,nimages,dest);
    break;
  default: ;
    VError(" unknown type %d",type);
  }
 
  
  /* output */
  VAttrList out_list = VCreateAttrList();
  VAttrList geolist = VGetGeoInfo(list);
  VSetGeoInfo(geolist,out_list);
  VHistory(VNumber(options),options,prg,&list,&out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,dest);
  VWriteFile (out_file, out_list);
  fprintf (stderr, "%s: done.\n", argv[0]);
  exit(0);
}
