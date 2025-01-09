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

/* 
** approximation: convert t to p values 
*/
double Xt2z(double t,double df)
{
  double z=0,u;

  u = df*log(1.0+t*t/df)*(1.0-0.5/df);
  if (u <= 0) return 0;
  z = sqrt(u);
  if (t < 0) z = -z;
  return z;
}


double gmin(double x,double y)
{
  if (x < y) return x;
  else return y;
}


/* contrast value */
double getcontrast(gsl_vector *beta,gsl_vector *bcov,double edf,double zthr,int type)
{
  double z=0;
  double bx0=beta->data[0];
  double bx1=beta->data[1];
  double bx2=beta->data[2];
  
  double c00=bcov->data[0];
  double c01=bcov->data[1];
  double c02=bcov->data[2];
  double c11=bcov->data[3];
  double c12=bcov->data[4];
  double c22=bcov->data[5];

  double t01 = (bx0 - bx1)/sqrt(c00 + c11 - 2.0*c01);
  double t02 = (bx0 - bx2)/sqrt(c00 + c22 - 2.0*c02);
  double t12 = (bx1 - bx2)/sqrt(c11 + c22 - 2.0*c12);

  double z01 = Xt2z(t01,edf);
  double z02 = Xt2z(t02,edf);
  double z12 = Xt2z(t12,edf);


  z=0;
  switch (type) {

  case 0: z = bx0; break;
  case 1: z = bx1; break;
  case 2: z = bx2; break;

  case 3: z = 2.0*bx0 - bx1 - bx2; break;
  case 4: z = 2.0*bx1 - bx0 - bx2; break;
  case 5: z = 2.0*bx2 - bx1 - bx0; break;

  case 6: z = gmin(z01,z02); break;    /* deep>middle && deep>superficial */
  case 7: z = gmin(-z01,z12); break;   /* middle>deep && middle>superficial */
  case 8: z = gmin(-z12,-z02); break;  /* superficial>deep && superficial>middle */
    
  case 9:  z = gmin(-z01,-z02);  break;  /* deep<middle && deep<superficial */
  case 10: z = gmin(z01,-z12);  break;   /* middle<deep && middle<superficial */
  case 11: z = gmin(z02,z12);   break;   /* superficial<deep && superficial<middle */

  case 12:  /* max ID */
    z=0;
    if (z01 > zthr && z02 > zthr) z = 1;
    if (-z01 > zthr && z12 > zthr) z = 2;
    if (-z02 > zthr && -z12 > zthr) z = 3;
    break;

  case 13:  /* min ID */
    z=0;
    if (z01 < zthr && z02 < zthr) z = 1;
    if (-z01 < zthr && z12 < zthr) z = 2;
    if (-z02 < zthr && -z12 < zthr) z = 3;
    break;
    
  default:
    VError(" unknown type %d",type);
  }
  return z;
}


VDictEntry TypDict[] = {
  { "deep", 0, 0,0,0,0  },
  { "middle", 1, 0,0,0,0  },
  { "superficial", 2, 0,0,0,0  },

  { "xdeep", 3, 0,0,0,0  },
  { "xmiddle", 4, 0,0,0,0  },
  { "xsuperficial", 5, 0,0,0,0  },
  
  { "cdeep", 6, 0,0,0,0  },
  { "cmiddle", 7, 0,0,0,0  },
  { "csuperficial", 8, 0,0,0,0  },
 
  { "notdeep", 9, 0,0,0,0  },
  { "notmiddle", 10, 0,0,0,0  },
  { "notsuperficial", 11, 0,0,0,0  },

  { "max", 12, 0,0,0,0 },
  { "min", 13, 0,0,0,0 },
  { NULL, 0,0,0,0,0 }
};


int main(int argc, char *argv[])
{
  static VShort type = 0;
  static VFloat zthr = 0;
  static VOptionDescRec  options[] = {
    {"type", VShortRepn,1,(VPointer) &type,VOptionalOpt,TypDict,"Output type"},
    {"zthr", VFloatRepn,1,(VPointer) &zthr,VOptionalOpt,NULL,"z-threshold"},
  };
  FILE *out_file=NULL;
  VString in_file=NULL;
  char *prg=GetLipsiaName("vcylarim_stats");
  fprintf (stderr, " %s\n", prg);

  
  /* parse command line */
  VParseFilterCmdX(VNumber(options),options,argc,argv,&in_file,&out_file);

  
  /* read beta images */
  VLong tr=0L;
  VAttrList list = VReadAttrList(in_file,tr,TRUE,FALSE);
  if (list == NULL) VError(" error reading input file %s",in_file);
  VImage src;
  VAttrListPosn posn;
  int nbeta=0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (strcmp(VGetAttrName (& posn),"nbeta") == 0) {
      VGetAttrValue (& posn, NULL,VShortRepn, & nbeta);
    }
  }
  fprintf(stderr," nbeta: %d\n",nbeta);
  VImage *betaimage = (VImage *)VCalloc(nbeta,sizeof(VImage));
  
  int n=0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    if (strcmp(VGetAttrName (& posn),"beta") == 0) {
      VGetAttrValue (& posn, NULL,VImageRepn, & src);
      if (n >= nbeta) VError(" too many beta images ");
      betaimage[n] = VCopyImage(src,NULL,VAllBands);
      n++;
    }
  }
  if (n < 1) VError(" no beta images found");

  int ncov=0;
  VImage *covimage = NULL;
  VImage edfimage = NULL;
  if (nbeta == 4) {
    ncov=6;
    covimage = (VImage *)VCalloc(ncov,sizeof(VImage));
    n=0;
    for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      if (strcmp(VGetAttrName (& posn),"covariance") == 0) {
	VGetAttrValue (& posn, NULL,VImageRepn, & src);
	if (n >= ncov) VError(" too many cov images ");
	covimage[n] = VCopyImage(src,NULL,VAllBands);
	n++;
      }
      if (strcmp(VGetAttrName (& posn),"edf") == 0) {
	VGetAttrValue (& posn, NULL,VImageRepn, & src);
	edfimage = VCopyImage(src,NULL,VAllBands);
      }
    }
  }


  /* alloc output image */
  VImage dest = VCreateImageLike(betaimage[0]);
  VFillImage(dest,VAllBands,0);


  /* beta contrasts */
  size_t i,j;
  size_t npixels = VImageNPixels(betaimage[0]);
  gsl_vector *beta = gsl_vector_calloc(nbeta);
  gsl_vector *bcov=NULL;
  if (ncov>0) bcov = gsl_vector_calloc(ncov);

  VFloat *pb,*pc;
  VFloat *pp = VImageData(dest);
  VFloat *pe = VImageData(edfimage);
  
  for (i=0; i<npixels; i++) {

    double s=0;
    for (j=0; j<nbeta; j++) {
      pb = VImageData(betaimage[j]);
      beta->data[j] = (double)(pb[i]);
      s += fabs(beta->data[j]);
    }
    if (s < TINY) continue;
    
    for (j=0; j<ncov; j++) {
      pc = VImageData(covimage[j]);
      bcov->data[j] = (double)(pc[i]);
    }
    pp[i] = (float)getcontrast(beta,bcov,(double)pe[i],(double)zthr,(int)type);
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
