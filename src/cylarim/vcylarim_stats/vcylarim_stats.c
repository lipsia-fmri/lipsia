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



double gmin(double x,double y)
{
  if (x < y) return x;
  else return y;
}

double gmax(double x,double y)
{
  if (x > y) return x;
  else return y;
}


/* contrast value */
double getcontrast(gsl_vector *beta,gsl_vector *zval,double zthr,int type)
{
  double z=0;
  double bx0=beta->data[0];
  double bx1=beta->data[1];
  double bx2=beta->data[2];

  double z01 = zval->data[0];
  double z02 = zval->data[1];
  double z12 = zval->data[2];

  double z10 = -z01;
  double z20 = -z02;
  double z21 = -z12;

  
  z=0;
  switch (type) {

  case 0: z = bx0; break;
  case 1: z = bx1; break;
  case 2: z = bx2; break;

  case 3: z = 2.0*bx0 - bx1 - bx2; break;
  case 4: z = 2.0*bx1 - bx0 - bx2; break;
  case 5: z = 2.0*bx2 - bx1 - bx0; break;

  case 6: z = gmin(z01,z02); break;    /* deep>middle && deep>superficial */
  case 7: z = gmin(z10,z12); break;    /* middle>deep && middle>superficial */
  case 8: z = gmin(z21,z20); break;    /* superficial>deep && superficial>middle */
    
  case 9:  z = gmax(z10,z20);  break;   /* deep<middle or deep<superficial */
  case 10: z = gmax(z01,z21);  break;   /* middle<deep or middle<superficial */
  case 11: z = gmax(z02,z12);  break;   /* superficial<deep or superficial<middle */

    
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
  
  { "zdeep", 6, 0,0,0,0  },
  { "zmiddle", 7, 0,0,0,0  },
  { "zsuperficial", 8, 0,0,0,0  },
 
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
  int nbeta = 3;
  VImage *betaimage = (VImage *)VCalloc(nbeta,sizeof(VImage));
  
  int n=0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    if (strcmp(VGetAttrName (& posn),"beta") == 0 && n < nbeta) {
      VGetAttrValue (& posn, NULL,VImageRepn, & src);
      betaimage[n] = VCopyImage(src,NULL,VAllBands);
      n++;
    }
  }
  if (n < 1) VError(" no beta images found");

  
  VImage *zvalimage = NULL;   /* zval images for permutation tests */
  size_t nzval=3;
  zvalimage = (VImage *)VCalloc(nzval,sizeof(VImage));
  n=0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    if (strcmp(VGetAttrName (& posn),"zvalimage") == 0) {
      VGetAttrValue (& posn, NULL,VImageRepn, & src);
      if (n >= nzval) VError(" too many zval images ");
      zvalimage[n] = VCopyImage(src,NULL,VAllBands);
      n++;
    }
  }


  /* alloc output image */
  VImage dest = VCreateImageLike(betaimage[0]);
  VFillImage(dest,VAllBands,0);


  /* beta contrasts */
  size_t i,j;
  size_t npixels = VImageNPixels(betaimage[0]);
  gsl_vector *beta = gsl_vector_calloc(nbeta);
  gsl_vector *zval = gsl_vector_calloc(nzval);

  VFloat *pb,*pz;
  VFloat *pp = VImageData(dest);
  
  for (i=0; i<npixels; i++) {

    double s=0;
    for (j=0; j<nbeta; j++) {
      pb = VImageData(betaimage[j]);
      beta->data[j] = (double)(pb[i]);
      s += fabs(beta->data[j]);
    }
    if (s < TINY) continue;
    
    for (j=0; j<nzval; j++) {
      pz = VImageData(zvalimage[j]);
      zval->data[j] = (double)(pz[i]);
    }
    pp[i] = (float)getcontrast(beta,zval,(double)zthr,(int)type);
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
