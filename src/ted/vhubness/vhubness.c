/*
** hubness map 
** this program expects as input a list of edges produced by "vted" 
** 
** Ref: Lohmann et al (2016) PLoS One, 11(6):e0158185
**
** G.Lohmann, May 2015
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#define SQR(x) ((x)*(x))
#define ABS(x) ((x) > 0 ? (x) : -(x))

VImage VEdge2Image(VBundle ebundle,VBundle ibundle,VBundle jbundle,long nedges,VImage pointmap,VImage roi,VFloat xmin,VFloat xmax)
{
  float *E = ebundle->data;
  int *I = ibundle->data;
  int *J = jbundle->data;
 

  /* alloc dest image */
  int nrows   = (int)VPixel(pointmap,0,3,1,VShort);  /* nrows */
  int ncols   = (int)VPixel(pointmap,0,3,2,VShort);  /* ncols */
  int nslices = (int)VPixel(pointmap,0,3,0,VShort);  /* nslices */
  int nvox    = VImageNColumns(pointmap);
  VImage dest = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);
  VCopyImageAttrs (pointmap,dest);
  fprintf(stderr," nvox: %d\n",nvox);


  /* map to image */
  long e=0;
  double nx=0,mx=0;
  for (e=0; e<nedges; e++) {
    int i = I[e];
    int j = J[e];
    float z = E[e];
    mx++;

    /* apply threshold */
    if (z < xmin || z > xmax) continue;

    /* map to image */
    int bi = VPixel(pointmap,0,0,i,VShort);
    int ri = VPixel(pointmap,0,1,i,VShort);
    int ci = VPixel(pointmap,0,2,i,VShort);

    int bj = VPixel(pointmap,0,0,j,VShort);
    int rj = VPixel(pointmap,0,1,j,VShort);
    int cj = VPixel(pointmap,0,2,j,VShort);

    if (roi) {
      double xi = VGetPixel(roi,bi,ri,ci);
      double xj = VGetPixel(roi,bj,rj,cj);
      if (xi < 0.5 && xj < 0.5) continue;
    }

    float u = VPixel(dest,bi,ri,ci,VFloat);
    float v = VPixel(dest,bj,rj,cj,VFloat);
    VPixel(dest,bi,ri,ci,VFloat) = u+1.0;
    VPixel(dest,bj,rj,cj,VFloat) = v+1.0;
    nx++;
  }
  if (nx < 1) VError(" no voxels found");
  fprintf(stderr," nx: %.0lf,  mx: %.0lf,  %lf\n",nx,mx,nx/mx);

  return dest;
}  


int main (int argc,char *argv[])
{
  static VFloat   xmax  = 1.0e+12;
  static VFloat   xmin  = -1.0e+12;
  static VString  filename = "";
  static VOptionDescRec  options[] = {  
    {"max",VFloatRepn,1,(VPointer) &xmax,VOptionalOpt,NULL,"max threshold"},
    {"min",VFloatRepn,1,(VPointer) &xmin,VOptionalOpt,NULL,"min threshold"},
    {"roi",VStringRepn,1,(VPointer) &filename,VOptionalOpt,NULL,"Region of interest mask"},
  };
  FILE *in_file, *out_file;
  VAttrList list=NULL,list1=NULL;
  VAttrListPosn posn;
  VString str;
  VBundle ebundle=NULL,ibundle=NULL,jbundle=NULL;
  VLong nedges=0;
  VImage dest=NULL,roi=NULL,map=NULL;
  char *prg=GetLipsiaName("vhubness");
  fprintf(stderr, "%s\n", prg);

  
  /* Parse command line arguments and identify files: */
  VParseFilterCmd (VNumber (options), options, argc, argv,& in_file, & out_file);


  /* Region of interest mask */
  if (strlen(filename) > 2) {
    FILE *fp = VOpenInputFile (filename, TRUE);
    list1 = VReadFile (fp, NULL);
    if (! list1) VError("Error reading ROI file");
    fclose(fp);
    for (VFirstAttr (list1, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      VGetAttrValue (& posn, NULL,VImageRepn, & roi);
      if (VPixelRepn(roi) != VBitRepn && VPixelRepn(roi) != VUByteRepn && VPixelRepn(roi) != VShortRepn) {
	roi = NULL;
	continue;
      }
    }
    if (roi == NULL) VError(" no ROI found");
  }
  
  
  /* Read the input file: */
  list = VReadFile (in_file, NULL);
  if (! list) exit (1);
  
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    str = VGetAttrName(&posn);

    if (strcmp(str,"EdgeDensity") == 0 && VGetAttrRepn (& posn) == VBundleRepn) {
      VGetAttrValue (& posn, NULL, VBundleRepn, & ebundle);
    }
    if (strcmp(str,"RowIndex") == 0 && VGetAttrRepn (& posn) == VBundleRepn) {
      VGetAttrValue (& posn, NULL, VBundleRepn, & ibundle);
    }
    if (strcmp(str,"ColIndex") == 0 && VGetAttrRepn (& posn) == VBundleRepn) {
      VGetAttrValue (& posn, NULL, VBundleRepn, & jbundle);
    }
    if (strcmp(str,"nedges") == 0) {
      VGetAttrValue (& posn, NULL, VLongRepn, & nedges);
    }
    if (strcmp(str,"map") == 0 && VGetAttrRepn (& posn) == VImageRepn) {
      VGetAttrValue (& posn, NULL, VImageRepn, & map);
    }
  }
  if (ibundle == NULL || jbundle == NULL || ebundle == NULL) VError(" data bundles missing");
  if (map == NULL) VError(" map not found");
  if (nedges == 0) VError(" no edges");
  

  /* map to image */
  dest = VEdge2Image(ebundle,ibundle,jbundle,(long)nedges,map,roi,xmin,xmax);
  

  /* update geoinfo */
  VAttrList out_list = VCreateAttrList();
  VAttrList geolist = VGetGeoInfo(list);
  if (geolist != NULL) {
    double *D = VGetGeoDim(geolist,NULL);
    D[0] = 3;
    D[4] = 1;
    VSetGeoDim(geolist,D);
    VSetGeoInfo(geolist,out_list);
  }


  /* Write out the results: */
  VAppendAttr(out_list,"image",NULL,VImageRepn,dest);
  VHistory(VNumber(options),options,prg,&list,&out_list);
  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
