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

double xsigmoid(double x,double z)
{
  double scale = 1.0/(1.0 + exp(-(x-3.0)));
  return z * scale;
}


/* contrast value */
double getcontrast(gsl_vector *beta,gsl_vector *zval,int type)
{
  double z=0,s=0;

  /* means or beta estimates */
  double bx0=beta->data[0];
  double bx1=beta->data[1];
  double bx2=beta->data[2];

  /* twosample t-tests */
  double z01 = zval->data[0];
  double z02 = zval->data[1];
  double z12 = zval->data[2];

  double z10 = -z01;
  double z20 = -z02;
  double z21 = -z12;


  z=0;
  switch (type) {

  /* means or beta estimates */
  case 0: z = bx0; break;
  case 1: z = bx1; break;
  case 2: z = bx2; break;

  /* twosample t-tests */
  case 3: z = z01; break;
  case 4: z = z02; break;
  case 5: z = z12; break;
    
  /* twosample t-tests, signflipped */
  case 6: z = z10; break;
  case 7: z = z20; break;
  case 8: z = z21; break;

  /* conjunctions */
  case 9: z = gmin(z01,z02); break;    /* deep>middle && deep>superficial */
  case 10: z = gmin(z10,z12); break;    /* middle>deep && middle>superficial */
  case 11: z = gmin(z21,z20); break;    /* superficial>deep && superficial>middle */

  /* max ID */
  case 12:  
    z=0;
    if (bx0 > bx1 && bx0 > bx2) z = 1;
    if (bx1 > bx0 && bx1 > bx2) z = 2;
    if (bx2 > bx1 && bx2 > bx0) z = 3;
    break;

  /* min ID */
  case 13:
    z=0;
    if (bx0 < bx1 && bx0 < bx2) z = 1;
    if (bx1 < bx0 && bx1 < bx2) z = 2;
    if (bx2 < bx1 && bx2 < bx0) z = 3;
    break;

    /* max beta values */
  case 14:
    z = bx0;
    if (bx1 > z) z = bx1;
    if (bx2 > z) z = bx2;
    break;

    /* min beta values */
  case 15: 
    z = bx0;
    if (bx1 < z) z = bx1;
    if (bx2 < z) z = bx2;
    break;

    /* max absolute beta-values */
  case 16: 
    s = 1;
    z = fabs(bx0);
    if (bx0 < 0) s = -1;
    if (fabs(bx1) > z) { z = fabs(bx1); if (bx1 < 0) s = -1; }
    if (fabs(bx2) > z) { z = fabs(bx2); if (bx2 < 0) s = -1; }
    z *= s;
    break;

  default:
    VError(" unknown type %d",type);
  }

  return z;
}


VDictEntry TypDict[] = {
  { "d", 0, 0,0,0,0  },
  { "m", 1, 0,0,0,0  },
  { "s", 2, 0,0,0,0  },

  { "d-m", 3, 0,0,0,0  },
  { "d-s", 4, 0,0,0,0  },
  { "m-s", 5, 0,0,0,0  },
  
  { "m-d", 6, 0,0,0,0  },
  { "s-d", 7, 0,0,0,0  },
  { "s-m", 8, 0,0,0,0  },
  
  { "top_d", 9, 0,0,0,0  },
  { "top_m", 10, 0,0,0,0  },
  { "top_s", 11, 0,0,0,0  },

  { "max_id", 12, 0,0,0,0 },
  { "min_id", 13, 0,0,0,0 },

  { "max", 14, 0,0,0,0 },
  { "min", 15, 0,0,0,0 },
  { "maxabs", 16, 0,0,0,0 },
  
  { "zabs", 17, 0,0,0,0 },

  { NULL, 0,0,0,0,0 }
};


int main(int argc, char *argv[])
{
  static VShort type = 0;
  static VOptionDescRec  options[] = {
    {"type", VShortRepn,1,(VPointer) &type,VOptionalOpt,TypDict,"Output type"},
  };
  FILE *out_file=NULL;
  VString in_file=NULL;
  VImage dest = NULL;
  char *prg=GetLipsiaName("vcylarim_stats");
  fprintf (stderr, " %s\n", prg);

  
  /* parse command line */
  VParseFilterCmdX(VNumber(options),options,argc,argv,&in_file,&out_file);
  if (type > 17 || type < 0) VError(" unknown type %d",type);
 

  /* read input images */
  VLong tr=0L;
  VAttrList list = VReadAttrList(in_file,tr,TRUE,FALSE);
  if (list == NULL) VError(" error reading input file %s",in_file);
  VImage src;
  VAttrListPosn posn;


  /* maxabs z-values per cylinder, can be used to create a mask (not layer-specific) */
  if (type == 17) {
    for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      if (strcmp(VGetAttrName (& posn),"cover") == 0) {
	VGetAttrValue (& posn, NULL,VImageRepn, & src);
	dest = VCopyImage(src,dest,VAllBands);
	goto ende;
      }
    }
  }


  /* beta images */
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


  /* zval images for permutation tests */
  VImage *zvalimage = NULL;
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
  dest = VCreateImageLike(betaimage[0]);
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
    
    pp[i] = (float)getcontrast(beta,zval,(int)type);
  }

  
  /* output */
 ende: ;
  VAttrList out_list = VCreateAttrList();
  VAttrList geolist = VGetGeoInfo(list);
  VSetGeoInfo(geolist,out_list);
  VHistory(VNumber(options),options,prg,&list,&out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,dest);
  VWriteFile (out_file, out_list);
  fprintf (stderr, "%s: done.\n", argv[0]);
  exit(0);
}
