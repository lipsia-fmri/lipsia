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



double
Correlation(const float *arr1,const float *arr2,int n)
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


VImage
VCCM(VImage **src, VImage mask, int nsubjects, int nslices,int first,int length,int type,int otype,VShort minval)
{
  VImage map=NULL,dest=NULL;
  int i,j,k,s,nvoxels,nt,len,last,b,r,c,bb,rr,cc,rad2,nrows,ncols;
  double u=0,v=0;
  gsl_matrix *data=NULL;
  gsl_matrix_float **mat=NULL;


  /* exclude voxels not covered by data */
  nrows = VImageNRows(mask);
  ncols = VImageNColumns(mask);

  for (s=0; s<nsubjects; s++) {
    for (b=0; b<nslices; b++) {
      if (VImageNRows(src[s][b]) < 3) {
	for (r=0; r<nrows; r++) {
	  for (c=0; c<ncols; c++) VSetPixel(mask,b,r,c,0);
	}
      }
      else {
	/*
	if (VImageNRows(src[s][b]) != nrows) VError(" inconsistent image dimensions, subj %d, nrows %d %d ",
						    s,nrows,VImageNRows(src[s][b]));
	if (VImageNColumns(src[s][b]) != ncols) 
	  VError(" inconsistent image dimensions, subj %d, ncols %d %d",s,ncols,VImageNColumns(src[s][b])); 
	*/
	for (r=0; r<nrows; r++) {
	  for (c=0; c<ncols; c++) {
	    if (VPixel(src[s][b],0,r,c,VShort) < minval) VSetPixel(mask,b,r,c,0);
	  }
	}
      }
    }
  }


  /* count number of voxels, number of timesteps */
  nvoxels = nt = 0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VGetPixel(mask,b,r,c) < 0.5) continue;
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
	  if (VGetPixel(mask,b,r,c) < 0.5) continue;
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
  
  fprintf(stderr,"# first= %d, last= %d, ntimesteps= %d\n",
          (int)first,(int)last,(int)nt);
  fprintf(stderr,"# nsubjects= %d,  nvoxels= %d\n",nsubjects,nvoxels);
  fprintf(stderr,"# type= %s\n",TypeDict[type].keyword);
  


  /*
  ** store voxel addresses
  */
  map = VCreateImage(1,3,nvoxels,VShortRepn);
  if (map == NULL) VError(" error allocating addr map");
  VFillImage(map,VAllBands,0);

  i = 0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VGetPixel(mask,b,r,c) < 0.5) continue;
	if (i >= nvoxels) VError(" nox= %d",i);
	VPixel(map,0,0,i,VShort) = b;
	VPixel(map,0,1,i,VShort) = r;
	VPixel(map,0,2,i,VShort) = c;
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
      b = VPixel(map,0,0,i,VShort);
      r = VPixel(map,0,1,i,VShort);
      c = VPixel(map,0,2,i,VShort);

      float *ptr = gsl_matrix_float_ptr(mat[s],i,0);
      int k;
      for (k=first; k<=last; k++) {
	*ptr++ = (float) VPixel(src[s][b],k,r,c,VShort);
      }
    }
  }

  /*
  ** for each seed voxel assess intersubject consistency 
  */
  dest = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);
  VCopyImageAttrs (src[0][0], dest);
  VSetAttr(VImageAttrList(dest),"modality",NULL,VStringRepn,"conimg");


  /*
  ** get correlation maps 
  */
  data = gsl_matrix_calloc(nsubjects,nvoxels);

  rad2 = 3*9;
  for (i=0; i<nvoxels; i++) {
    if (i%10 == 0) fprintf(stderr," i: %7d\r",i);
    b = VPixel(map,0,0,i,VShort);
    r = VPixel(map,0,1,i,VShort);
    c = VPixel(map,0,2,i,VShort);
        
    len = 0;
    for (s=0; s<nsubjects; s++) {
      const float *arr1 = gsl_matrix_float_ptr(mat[s],i,0);

      k = 0;
      for (j=0; j<nvoxels; j++) {
	if (i == j) continue;
	bb = VPixel(map,0,0,j,VShort);
	rr = VPixel(map,0,1,j,VShort);
	cc = VPixel(map,0,2,j,VShort);
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

    b = VPixel(map,0,0,i,VShort);
    r = VPixel(map,0,1,i,VShort);
    c = VPixel(map,0,2,i,VShort);
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
  static VShort minval = 0;
  static VOptionDescRec options[] = {
    {"in", VStringRepn, 0, & in_files, VRequiredOpt, NULL,"Input files" },
    {"mask",VStringRepn,1,(VPointer) &mask_filename,VRequiredOpt,NULL,"mask file"},   
    {"first",VShortRepn,1,(VPointer) &first,VOptionalOpt,NULL,"first timestep to use"},
    {"length",VShortRepn,1,(VPointer) &length,VOptionalOpt,NULL,
     "length of time series to use, '0' to use full length"},
    {"type",VShortRepn,1,(VPointer) &type,VOptionalOpt,TypeDict,"type"},
    {"gauss",VShortRepn,1,(VPointer) &otype,VOptionalOpt,NULL,"use gaussian version of kendall"},
    {"minval",VShortRepn,1,(VPointer) &minval,VOptionalOpt,NULL,"signal threshold"},
    {"out", VStringRepn, 1, & out_filename, VRequiredOpt, NULL,"Output file" }
  };
  FILE *fp=NULL,*fpo=NULL;
  VStringConst in_filename;
  VAttrList list=NULL,list1=NULL,out_list=NULL,geolist=NULL;
  VAttrListPosn posn;
  VImage xsrc,**src,mask=NULL,dest=NULL,disc=NULL;
  int i=0,j=0,nsubjects=0,nslices=0;
  int b,r,c;
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
  fp = VOpenInputFile (mask_filename, TRUE);
  list1 = VReadFile (fp, NULL);
  if (! list1) VError("Error reading mask file");
  fclose(fp);

  for (VFirstAttr (list1, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & mask);
    if (VPixelRepn(mask) != VBitRepn && VPixelRepn(mask) != VUByteRepn && VPixelRepn(mask) != VShortRepn) {
      mask = NULL;
      continue;
    }
  }
  if (mask == NULL) VError(" no mask found");



  /* 
  ** read images  
  */
  nsubjects = in_files.number;
  src = (VImage **) VCalloc(nsubjects,sizeof(VImage *));
  fprintf(stderr,"# nsubjects= %d\n",nsubjects);

  for (i=0; i<nsubjects; i++) {

    in_filename = ((VStringConst *) in_files.vector)[i];
    fprintf(stderr," %3d:  %s\n",i,in_filename);
    fp = VOpenInputFile (in_filename, TRUE);
    list = VReadFile (fp, NULL);
    if (! list)  VError("Error reading image");
    fclose(fp);

    if (geolist == NULL) geolist = VGetGeoInfo(list);

    j=0;
    for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      VGetAttrValue (& posn, NULL, VImageRepn, & xsrc);
      if (VPixelRepn(xsrc) != VShortRepn) continue;
      j++;
    }
    if (i == 0) nslices = j;
    else if (j != nslices) VError(" inconsistent number of slices %d %d",j,nslices);

    src[i] = (VImage *) VCalloc(nslices,sizeof(VImage));
    j=0;
    for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      VGetAttrValue (& posn, NULL, VImageRepn, & xsrc);
      if (VPixelRepn(xsrc) != VShortRepn) continue;
      if (i >= nsubjects || j >= nslices) VError(" i,j= %d %d",i,j);
      src[i][j] = xsrc;
      j++;
    }
  }


  /*
  ** do CCM
  */
  dest = VCCM(src,mask,nsubjects,nslices,(int)first,(int)length,(int)type,(int)otype,minval);


  /* invert to get discordant map */
  disc = VCreateImageLike(dest);
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
  out_list = VCreateAttrList();
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
  fpo = VOpenOutputFile (out_filename, TRUE);
  if (! VWriteFile (fpo, out_list)) exit (1);
  fclose(fpo);
  fprintf (stderr, "# %s: done.\n", argv[0]);
  exit(0);
}
