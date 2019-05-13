/*
** CCM - connectivity consistency mapping
**
** G.Lohmann, July 2010
*/
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>
#include <viaio/option.h>
#include <viaio/os.h>
#include <viaio/VImage.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))

extern double VOCCC(gsl_matrix *,int,int);
extern double VICC(gsl_matrix *,int);
extern double VKendall_W(gsl_matrix *,int,int);
extern double SpearmanCorr(const float *data1,const float *data2, int n);

int metric = 1;
int verbose = 0;

VDictEntry TypeDict[] = {
  { "kendall", 0 },
  { "occc", 1 },
  { NULL }
};



double Correlation(const float *arr1,const float *arr2,int n)
{
  int i;
  double sx,sy,sxx,syy,sxy,rx;
  double tiny=1.0e-4;

  sxx = syy = sxy = sx = sy = 0;
  for (i=0; i<n; i++) {
    const double u = (double)arr1[i];
    const double v = (double)arr2[i];
    sx  += u;
    sy  += v;
    sxx += u*u;
    syy += v*v;
    sxy += u*v;
  }
  const double nx= n;
  const double u = nx*sxx - sx*sx;
  const double v = nx*syy - sy*sy;
  rx = 0;
  if (u*v > tiny)
    rx = (nx*sxy - sx*sy)/sqrt(u*v);
  return rx;
}


VImage VCCM(VImage **src, VImage mask, int nsubjects, int nslices,int first,int length,int type,int otype)
{
  VImage map=NULL,dest=NULL;
  int i,j,k,s,nvoxels,nt,len,last,b,r,c,bb,rr,cc,rad2,nrows,ncols;
  double u=0,v=0;
  gsl_matrix *data=NULL;
  gsl_matrix_float **mat=NULL;


  nrows = VImageNRows(mask);
  ncols = VImageNColumns(mask);


  /* count number of voxels, number of timesteps */
  nvoxels = nt = 0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VGetPixel(mask,b,r,c) < 0.001) continue;
	if (VImageNBands(src[0][b]) > nt)
	  nt = VImageNBands(src[0][b]);
	nvoxels++;
      }
    }
  }
  
  /* get shortest length of time series across subjects */
  nt = 99999;
  for (s=0; s<nsubjects; s++) {
    for (b=0; b<nslices; b++) {
      for (r=0; r<nrows; r++) {
	for (c=0; c<ncols; c++) {
	  if (VGetPixel(mask,b,r,c) < 0.001) continue;
	  if (VImageNBands(src[s][b]) < nt)
	    nt = VImageNBands(src[s][b]);
	}
      }
    }
  }


  /* get time steps to include */
  if (length < 1) length = nt-2;
  last = first + length -1;
  if (last >= nt) last = nt-1;
  if (first < 0) first = 1;

  nt = last - first + 1;
  if (nt < 2) VError(" not enough timesteps, nt= %d",nt);
  
  fprintf(stderr," first= %d, last= %d, ntimesteps= %d\n",
          (int)first,(int)last,(int)nt);
  fprintf(stderr," nsubjects= %d,  nvoxels= %d\n",nsubjects,nvoxels);
  fprintf(stderr," type= %s\n",TypeDict[type].keyword);
  


  /*
  ** store voxel addresses
  */
  map = VCreateImage(1,3,nvoxels,VIntegerRepn);
  if (map == NULL) VError(" error allocating addr map");
  VFillImage(map,VAllBands,0);

  i = 0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VGetPixel(mask,b,r,c) < 0.001) continue;
	if (i >= nvoxels) VError(" nvox= %d",i);
	VPixel(map,0,0,i,VInteger) = b;
	VPixel(map,0,1,i,VInteger) = r;
	VPixel(map,0,2,i,VInteger) = c;
	i++;
      }
    }
  }


  /*
  ** avoid casting to float, copy data to matrix
  */  
  mat = (gsl_matrix_float **) VCalloc(nsubjects,sizeof(gsl_matrix_float *));
  if (!mat) VError(" err allocating mat");

  for (s=0; s<nsubjects; s++) {
    mat[s] = gsl_matrix_float_calloc(nvoxels,nt);
    if (!mat[s]) VError(" err allocating mat[s]");

    for (i=0; i<nvoxels; i++) {      
      b = VPixel(map,0,0,i,VInteger);
      r = VPixel(map,0,1,i,VInteger);
      c = VPixel(map,0,2,i,VInteger);

      float *ptr = gsl_matrix_float_ptr(mat[s],i,0);
      int k;
      for (k=first; k<=last; k++) {
	*ptr++ = (float) VGetPixel(src[s][b],k,r,c);
      }
    }
  }

  /*
  ** for each seed voxel assess intersubject consistency 
  */
  dest = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);
  VCopyImageAttrs (src[0][0], dest);


  /*
  ** get correlation maps 
  */
  data = gsl_matrix_calloc(nsubjects,nvoxels);

  rad2 = 3*9;
  for (i=0; i<nvoxels; i++) {
    if (i%10 == 0) fprintf(stderr," i: %7d  of  %d\r",i,nvoxels);
    b = VPixel(map,0,0,i,VInteger);
    r = VPixel(map,0,1,i,VInteger);
    c = VPixel(map,0,2,i,VInteger);
        
    len = 0;
    for (s=0; s<nsubjects; s++) {
      const float *arr1 = gsl_matrix_float_ptr(mat[s],i,0);

      k = 0;
      for (j=0; j<nvoxels; j++) {
	if (i == j) continue;
	bb = VPixel(map,0,0,j,VInteger);
	rr = VPixel(map,0,1,j,VInteger);
	cc = VPixel(map,0,2,j,VInteger);
	/* exclude because of spatial smoothness: */
	if ((SQR(b-bb) + SQR(r-rr) + SQR(c-cc)) < rad2) continue; 

	const float *arr2 = gsl_matrix_float_ptr(mat[s],j,0);
	u = Correlation(arr1,arr2,nt);

	gsl_matrix_set(data,s,k,(double)u);
	k++;
      }
      len = k;
    }
    if (len >= data->size2) VError(" len= %d %d",len,data->size2);

    if (type == 0)
      v = VKendall_W(data,len,otype);
    else if (type == 1)
      v = VOCCC(data,len,verbose);
    else
      VError(" illegal type");

    b = VPixel(map,0,0,i,VInteger);
    r = VPixel(map,0,1,i,VInteger);
    c = VPixel(map,0,2,i,VInteger);
    VPixel(dest,b,r,c,VFloat) = v;
  }
  fprintf(stderr,"\n");
  return dest;
}


int main (int argc, char *argv[])
{
  static VArgVector in_files;
  static VString out_filename;
  static VString mask_filename;
  static VShort type   = 0; 
  static VShort otype  = 0; 
  static VShort first  = 2;
  static VShort length = 0;
  static VOptionDescRec options[] = {
    {"in", VStringRepn, 0, & in_files, VRequiredOpt, NULL,"Input files" },
    {"out", VStringRepn, 1, & out_filename, VRequiredOpt, NULL,"Output file" },
    {"mask",VStringRepn,1,(VPointer) &mask_filename,VRequiredOpt,NULL,"Mask file"},   
    {"first",VShortRepn,1,(VPointer) &first,VOptionalOpt,NULL,"First timestep to use"},
    {"length",VShortRepn,1,(VPointer) &length,VOptionalOpt,NULL,
     "Length of time series to use, '0' to use full length"},
    {"type",VShortRepn,1,(VPointer) &type,VOptionalOpt,TypeDict,"Type of metric"},
    {"gauss",VShortRepn,1,(VPointer) &otype,VOptionalOpt,NULL,"Use gaussian version of kendall"},
  };
  VString in_filename;
  VAttrList list=NULL,geolist=NULL;
  int b,r,c,i=0;
  float u;
  char *prg = GetLipsiaName("vccm");
  

  /*
  ** parse command line
  */
  if (! VParseCommand (VNumber (options), options, & argc, argv)) {
    VReportUsage (argv[0], VNumber (options), options, NULL);
    exit (EXIT_FAILURE);
  }
  if (argc > 1) {
    VReportBadArgs (argc, argv);
    exit (EXIT_FAILURE);
  }

  
  /*
  ** read mask
  */
  VAttrList mask_list = VReadAttrList(mask_filename,0L,TRUE,FALSE);
  if (mask_list == NULL) VError(" error reading %s",mask_filename);
  VImage mask = VReadImage(mask_list);
  if (mask == NULL) VError(" no mask found");



  /* 
  ** read input images  
  */
  int nslices=0,nt=0,nrows=0,ncols=0;
  int xslices=0,xnt=0,xrows=0,xcols=0;
  int nsubjects = (int)in_files.number;
  fprintf(stderr," nsubjects= %d\n",nsubjects);
  if (nsubjects < 2) VError(" not enough input files (%d), CCM should be used on >= 2 data sets",nsubjects);
  VImage **src = (VImage **) VCalloc(nsubjects,sizeof(VImage *));
  
  for (i=0; i<nsubjects; i++) {
    in_filename = ((VString *) in_files.vector)[i];
    list    = VReadAttrList(in_filename,0L,TRUE,FALSE);
    xslices = VAttrListNumImages(list);
    if (i==0) nslices = xslices;
    else if (nslices != xslices) VError(" inconsistent dimensions in input files");
    src[i] = VAttrListGetImages(list,nslices);
    if (i==0) VImageDimensions(src[i],nslices,&nt,&nrows,&ncols);
    else {
      VImageDimensions(src[i],xslices,&xnt,&xrows,&xcols);
      if (nrows != xrows) VError(" inconsistent dimensions in input files");
      if (ncols != xcols) VError(" inconsistent dimensions in input files");
      if (nt != xnt) VError(" inconsistent dimensions in input files");
    }

    /* use geometry info from 1st file */
    if (geolist == NULL) geolist = VGetGeoInfo(list);
  }


  /*
  **  CCM
  */
  VImage dest = VCCM(src,mask,nsubjects,nslices,(int)first,(int)length,(int)type,(int)otype);


  /* invert to get discordant map */
  VImage disc = VCreateImageLike(dest);
  VFillImage(disc,VAllBands,0);
  for (b=0; b<VImageNBands(dest); b++) {
    for (r=0; r<VImageNRows(dest); r++) {
      for (c=0; c<VImageNColumns(dest); c++) {
        if (VGetPixel(mask,b,r,c) < 1) continue;
        u = VPixel(dest,b,r,c,VFloat);
        VPixel(disc,b,r,c,VFloat) = 1-u;
      }
    }
  }

  
  /* 
  ** output
  */
  VAttrList out_list = VCreateAttrList();
  if (geolist != NULL) {
    double *D = VGetGeoDim(geolist,NULL);
    D[0] = 3;  /* 3D */
    D[4] = 1;  /* just one timestep */
    VSetGeoDim(geolist,D);
    VSetGeoInfo(geolist,out_list);
  }
  VAppendAttr (out_list,"concordant",NULL,VImageRepn,dest);
  VAppendAttr (out_list,"discordant",NULL,VImageRepn,disc);

  VHistory(VNumber(options),options,prg,&list,&out_list);
  FILE *fp = VOpenOutputFile (out_filename, TRUE);
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
  fprintf (stderr, "%s: done.\n", argv[0]);
  exit(0);
}
