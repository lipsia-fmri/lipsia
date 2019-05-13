
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_wavelet2d.h>

extern void gsl_sort_vector(gsl_vector *);
#define ABS(x) ((x) > 0 ? (x) : -(x))



void WriteTest(gsl_matrix *data)
{
  FILE *fp=0;
  VAttrList alist;
  VImage tmp;
  double u;
  int i,j;

  tmp = VCreateImage(1,data->size1,data->size2,VFloatRepn);
  for (i = 0; i < data->size1; i++) {
    for (j = 0; j < data->size2; j++) {
      u = gsl_matrix_get(data,i,j);
      VPixel(tmp,0,i,j,VFloat) = u;
    }
  }
  alist = VCreateAttrList();
  VAppendAttr(alist,"image",NULL,VImageRepn,tmp);
  fp = fopen("qtest.v","w");
  VWriteFile (fp,alist);
  fclose(fp);
  exit(0);
}


double 
estimateSigma(gsl_matrix *data)
{
  int i,j,i2,j2,n,nn;
  gsl_vector *vec=NULL;
  double median=0,med0,u,xmax,tiny=1.0e-6;

  i2 = data->size1/2;
  j2 = data->size2/2;

  nn = 3 * (i2+1) * (j2+1);
  vec = gsl_vector_calloc(nn);


  /*
  ** first pass, median, no absolute values
  */
  xmax = VRepnMaxValue(VDoubleRepn);
  gsl_vector_set_all(vec,xmax);

  n = 0;
  for (i=0; i<i2; i++) {
    for (j=j2; j<data->size2; j++) {
      u = gsl_matrix_get(data,i,j);
      if (ABS(u) > tiny) gsl_vector_set(vec,n++,u);
    }
  }

  for (i=i2; i<data->size1; i++) {
    for (j=0; j<j2; j++) {
      u = gsl_matrix_get(data,i,j);
      if (ABS(u) > tiny) gsl_vector_set(vec,n++,u);
    }
  }

  for (i=i2; i<data->size1; i++) {
    for (j=j2; j<data->size2; j++) {
      u = gsl_matrix_get(data,i,j);
      if (ABS(u) > tiny) gsl_vector_set(vec,n++,u);
    }
  }
  
  gsl_sort_vector(vec);
  u  = (double)n/2.0;
  nn = (int)u;
  med0 = 0.5*(gsl_vector_get(vec,nn) + gsl_vector_get(vec,nn+1));

  /*
  ** median, 2nd pass
  */
  xmax = VRepnMaxValue(VDoubleRepn);
  gsl_vector_set_all(vec,xmax);

  n = 0;
  for (i=0; i<i2; i++) {
    for (j=j2; j<data->size2; j++) {
      u = gsl_matrix_get(data,i,j);
      if (ABS(u) <= tiny) continue;
      u = ABS(u - med0);
      gsl_vector_set(vec,n++,u);
    }
  }

  for (i=i2; i<data->size1; i++) {
    for (j=0; j<j2; j++) {
      u = gsl_matrix_get(data,i,j);
      if (ABS(u) <= tiny) continue;
      u = ABS(u - med0);
      gsl_vector_set(vec,n++,u);
    }
  }

  for (i=i2; i<data->size1; i++) {
    for (j=j2; j<data->size2; j++) {
      u = gsl_matrix_get(data,i,j);
      if (ABS(u) <= tiny) continue;
      u = ABS(u - med0);
      gsl_vector_set(vec,n++,u);
    }
  }
  
  gsl_sort_vector(vec);
  u  = (double)n/2.0;
  nn = (int)u;
  median = 0.5*(gsl_vector_get(vec,nn) + gsl_vector_get(vec,nn+1));

  return (median/0.455);
}




double sgn(double x)
{
  if (x < 0) return -1;
  else if (x > 0) return 1;
  else return 0;
}

void 
Filtering(gsl_matrix *data, int type, double intensity, double sigma)
{
  int i,j,n;
  double d,dd,alpha;

  n = data->size1;

  double tau = sigma * intensity * sqrt(2 * log10((double) (n*n)));

  for (i=0; i < data->size1; i++) {
    for (j=0; j < data->size2; j++) {
      d = gsl_matrix_get(data,i,j);


      /* Wiener-Filter */
      if (type == 0) {
	alpha = ((d * d) - intensity * (sigma * sigma)) / (d * d);
	if (alpha < 0) alpha = 0;
	d *= alpha;
      }

      /* Soft-Threshold */
      else if (type == 1) {
	if (ABS(d) < tau)
	  dd = 0;
	else
	  dd = sgn(d) * (ABS(d) - tau);
	d = dd;
      }

      /* Hard-Threshold */
      else if (type == 2) {
	if (ABS(d) < tau) d = 0;
      }

      gsl_matrix_set(data,i,j,d);
    }
  }
}

void
Shift(gsl_matrix *data,gsl_matrix *shiftdata,int ishift,int jshift)
{
  int i,j,ii,jj,n;
  double u;

  gsl_matrix_set_zero(shiftdata);

  n = data->size1;  
  for (i=0; i<n; i++) {
    ii = i+ishift;
    if (ii < 0 || ii >= n) continue;
    for (j=0; j<n; j++) {
      jj = j+jshift;
      if (jj < 0 || jj >= n) continue;
      u = gsl_matrix_get(data,i,j);
      gsl_matrix_set(shiftdata,ii,jj,u);
    }
  }
}


void
AddShift(gsl_matrix *data,gsl_matrix *shiftdata,int ishift,int jshift)
{
  int i,j,ii,jj,n;
  double u,v;
  n = data->size1;
  
  for (i=0; i<n; i++) {
    ii = i+ishift;
    if (ii < 0 || ii >= n) continue;
    for (j=0; j<n; j++) {
      jj = j+jshift;
      if (jj < 0 || jj >= n) continue;
      u = gsl_matrix_get(data,i,j);
      v = gsl_matrix_get(shiftdata,ii,jj);
      u = u+v;
      gsl_matrix_set(data,i,j,u);
    }
  }
}


VImage
VWavelets(VImage src,VImage dest,int filtertype,int wtype,int nshift,double level)
{
  static gsl_wavelet_workspace *work=NULL;
  static gsl_wavelet *w=NULL;
  static gsl_matrix *data=NULL,*tmp=NULL,*shiftdata=NULL;
  int nslices,nrows,ncols,i,j,slice,n,size,dim;
  double u,sum1,sum2,nx,mean,sig,sigma=0,avesig;
  extern VFloat VGetPixelValue (VImage,int);


  /*
  ** read image info
  */
  nslices = VImageNBands(src);
  nrows = VImageNRows(src);
  ncols = VImageNColumns(src);

  if (dest == NULL)
    dest = VCreateImage(nslices,nrows,ncols,VPixelRepn(src));
  VCopyImageAttrs (src, dest);
  VFillImage(dest,VAllBands,0);
  
  size = nrows;
  if (nrows < ncols) size = ncols;
  n = (int) ceil(log10((double) size) / log10(2.0));
  dim = (int) pow(2.0, (double) n);


  /*
  ** ini wavelets transform
  */
  if (work == NULL) {

    switch (wtype) {
    case 0:
      w = gsl_wavelet_alloc (gsl_wavelet_haar_centered,2);
      break;
    case 1:
      w = gsl_wavelet_alloc (gsl_wavelet_daubechies_centered,4);
      break;
    case 2:
      w = gsl_wavelet_alloc (gsl_wavelet_daubechies_centered,6);
      break;
    case 3:
      w = gsl_wavelet_alloc (gsl_wavelet_daubechies_centered,8);
      break;
    case 4:
      w = gsl_wavelet_alloc (gsl_wavelet_bspline_centered,103);
      break;
    case 5:
      w = gsl_wavelet_alloc (gsl_wavelet_bspline_centered,202);
      break;
    default: 
      VError("illegal wavelet type");
    }

    work = gsl_wavelet_workspace_alloc (dim);
    data = gsl_matrix_calloc(dim,dim);
    shiftdata = gsl_matrix_calloc(dim,dim);
    tmp = gsl_matrix_calloc(dim,dim);
  }
  

  avesig = 0;
  for (slice=0; slice<nslices; slice++) {

    /*
    ** read data
    */
    sum1 = sum2 = nx = 0;
    for (i=0; i<nrows; i++) {
      for (j=0; j<ncols; j++) {
	u = VGetPixel(src,slice,i,j);
	sum1 += u;
	sum2 += u*u;
	nx++;
      }
    }
    mean = sum1/nx;
    sig = sqrt((sum2 - nx * mean * mean) / (nx - 1.0));

    if (sig < 1.0e-3) {
      VFillImage(dest,slice,0);
      continue;
    }

    gsl_matrix_set_zero(data);
    for (i=0; i<nrows; i++) {
      for (j=0; j<ncols; j++) {
	u = VGetPixel(src,slice,i,j);
	u = (u-mean)/sig;
	gsl_matrix_set(data,i,j,u);
      }
    }
    gsl_matrix_memcpy(tmp,data);


    /*
    ** get sigma estimate (Donoho)
    */
    if (gsl_wavelet2d_nstransform_matrix_forward(w,data,work) != GSL_SUCCESS)
      VError(" gsl_wavelet2d_nstransform_matrix_forward failed");

    sigma = estimateSigma(data);
    /* fprintf(stderr," slice: %3d  %f\n",slice,sigma); */
    avesig += sigma;

    Filtering(data,filtertype,level,sigma);
    gsl_wavelet2d_nstransform_matrix_inverse(w,data,work);


    /*
    ** translation invariance, shift 4/8 algorithm
    */
    nx = 1;
    for (i=-nshift; i<=nshift; i++) {
      for (j=-nshift; j<=nshift; j++) {
	if (i == 0 && j == 0) continue;
	Shift(tmp,shiftdata,i,j);
	gsl_wavelet2d_nstransform_matrix_forward(w,shiftdata,work);
	Filtering(shiftdata,filtertype,level,sigma);
	gsl_wavelet2d_nstransform_matrix_inverse(w,shiftdata,work);
	AddShift(data,shiftdata,i,j);
	nx++;
      }
    }
    gsl_matrix_scale(data,1.0/nx);
    

    /*
    ** output
    */
    for (i=0; i<nrows; i++) {
      for (j=0; j<ncols; j++) {
	u = gsl_matrix_get(data,i,j);
	u = sig*u + mean;
	if (u < VPixelMinValue(dest)) u = VPixelMinValue(dest);
	if (u > VPixelMaxValue(dest)) u = VPixelMaxValue(dest);
	VSetPixel(dest,slice,i,j,u);
      }
    }
  }
  avesig /= (double)nslices;
  /* fprintf(stderr," kappa: %f\n",sqrt(avesig) * 1.75); */

  /*
  gsl_matrix_free(shiftdata);
  gsl_matrix_free(tmp);
  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);
  free (data);
  */
  return dest;
}


VDictEntry WDict[] = {
  { "haar", 0 },
  { "daub4", 1 },
  { "daub6", 2 },
  { "daub8", 3 },
  { "bspline103", 4 },
  { "bspline202", 5 },
  { NULL }
};


int main (int argc,char *argv[])
{ 
  static VShort ftype  = 0;
  static VShort wtype  = 0;
  static VShort nshift = 1;
  static VDouble level = 0.5;
  static VOptionDescRec  options[] = {
    {"wavelet",VShortRepn,1,(VPointer) &wtype,VOptionalOpt,WDict,"wavelet type"},
    {"filter",VShortRepn,1,(VPointer) &ftype,VOptionalOpt,NULL,"filter type"},
    {"shift",VShortRepn,1,(VPointer) &nshift,VOptionalOpt,NULL,"shift window size"},
    {"level",VDoubleRepn,1,(VPointer) &level,VOptionalOpt,NULL,"level"}
  };
  FILE *in_file,*out_file;
  VAttrList list=NULL;
  VAttrListPosn posn;
  VImage src=NULL,dest=NULL;
  int n=0,nimages=0;
  char *prg=GetLipsiaName("vdenoise");
  fprintf (stderr, "%s\n", prg);


  VParseFilterCmd (VNumber (options),options,argc,argv,&in_file,&out_file);


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
    dest = VWavelets(src,dest,(int)ftype,(int)wtype,(int)nshift,level);
    src = VCopyImage(dest,src,VAllBands);
    VSetAttrValue (& posn, NULL,VImageRepn,src);
    n++;
  }
  if (src == NULL) VError(" no input image found");


  VHistory(VNumber(options),options,prg,&list,&list);
  if (! VWriteFile (out_file, list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
