/*! \file
  Geometric transformations using cubic spline interpolation.


This file contains functions for cubic spline resampling
and 3D scaling. The transformation equation is:

   y = A(x-x0) + b

where x,x0,b,y are 1x3 vectors and A is a 3x3 matrix.
The vector x0 can be used to specify a position that
remains unchanged by the transformation.


\par Authors:
A.Hagert <hagert@cns.mpg.de>, 07.01.2003,
G.Lohmann, MPI-CBS
*/



/**************************************************/
/* Resample a 3D image cubic spline interpolation */
/*                                                */
/* A.Hagert <hagert@cns.mpg.de>, 07.01.2003       */
/*                                                */
/* bei 160x200x160 -> ca. 7-8 Minuten (gesamt)    */
/**************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>

#include <via/via.h>


#define ABS(x) ((x) > 0 ? (x) : -(x))
#define CEIL(x) ((int) ((x) + 0.9999999))

/* Hilfsfunktionen: */

/*****************************/
/* der kubische B-Spline ... */
/*****************************/

float beta3 (float x)    
{
  float r = 0.0;
  float y,y2;

  y  = ABS(x);
  y2 = 2 - y;

  if ((0 <= y) && (y < 1)) {
    r = 0.666666 - y*y + y*y*y*0.5;
  }
  else if ((1 <= y) && (y < 2)) {
    r = y2*y2*y2 / 6.0;
  }

  return r;
}

/**********************************/
/* komplette Spline-Funktion (3d) */
/**********************************/

float f (VImage field,float z,float y,float x)
{
  int nbands,nrows,ncols;
  int b, r, c, b_low, r_low, c_low, b_high, r_high, c_high;
  float temp = 0.0, t1 = 0.0, t2 = 0.0, t0=0;

  r_low = CEIL ((y-2)); /* ceil (y - (degree+1)/2); */
  c_low = CEIL ((x-2)); /* ceil (x - (degree+1)/2); */
  b_low = CEIL ((z-2)); /* ceil (z - (degree+1)/2); */

  nbands = VImageNBands(field);
  nrows  = VImageNRows(field);
  ncols  = VImageNColumns(field);

  r_high = r_low + 3; /* "3" ist der Grad der Interpolationspolynome */
  c_high = c_low + 3;
  b_high = b_low + 3;

  if (c_low < 0) c_low = 0;
  if (r_low < 0) r_low = 0;
  if (b_low < 0) b_low = 0;

  if (c_high >= ncols)  c_high = ncols-1;
  if (r_high >= nrows)  r_high = nrows-1;
  if (b_high >= nbands) b_high = nbands-1;
  
  for (b = b_low; b <= b_high; b++) {
    t1 = beta3 (z-b);
    for (r = r_low; r <= r_high; r++) {
      t2 = t1 * beta3 (y-r);
      for (c = c_low; c <= c_high; c++) {
	t0 = VPixel(field,b,r,c,VFloat);
	temp += t0 * t2 * beta3 (x-c);
      }
    }
  }
  
  if (temp < 0.0) temp = 0.0;
  if (temp > 255.0) temp = 255.0; 
  
  return (temp);
}

/*******************************************************/
/* Bestimmung der Koeffizienten c[k] - 1 dimensional ! */
/*******************************************************/

void compute_koeffizients (float *data, int anz)   
{
  int n;
  float z1;
  float c1[anz];
  float u;

  z1 = sqrt((double)3.0) - 2.0;

  /* Hinrichtung (wie brutal) */
  c1[0] = 0;
  u = 1;
  for (n = 0; n <= anz-1; n++) {  /* c1[0] - Startbedingung */
    c1[0] += data[n] * u;
    u *= z1;
  }
   
  for (n = 1; n <= anz-1; n++)  /* Filtern von links nach rechts */
    c1[n] = data[n] + z1 * c1[n-1];

  /* Rueckrichtung */

  data[anz-1] = z1 / (z1*z1 -1) * (c1[anz-1] + z1 * c1[anz-2]);

  for (n = anz-2; n >= 0; n--)  /* Filtern von rechts nach links */
    data[n] = z1 * (data[n+1] - c1[n]);
    
  for (n = 0; n <= anz-1; n++) data[n] *= 6.0;
}

/*******************************************************/
/* Hauptfunktion: Ermittlung des interpolierten Bildes */
/*******************************************************/
/*!
\fn VImage VCubicSplineSample3d (VImage src,VImage dest,VImage transform,
                                float b0,float r0,float c0,
			        int dst_nbands,int dst_nrows,int dst_ncolumns)
\brief Resample a 3D image using trilinear interpolation.

\param src   input image (any repn)
\param dest  output image (any repn)
\param transform  4x3 transformation image (float repn).
The first column of <transform> contains the translation vector.
The remaining three columns contains the 3x3 linear transformation matrix.
\param b0           slice address that remains fixed 
\param r0           row address that remains fixed 
\param c0           column address that remains fixed 
\param dst_nbands   number of output slices
\param dst_nrows    number of output rows
\param dst_ncolumns number of output columns
*/
VImage 
VCubicSplineSample3d (VImage src, VImage dest,VImage transform,
		      float b0,float r0,float c0,
		      int dest_nbands,int dest_nrows,int dest_ncolumns)
{
  VImage int_koeff=NULL;           /* Interpolationskoeffizienten (3d - Feld) */
  float *data;                     /* da kommen mal die Stuetzstellen rein */
  int b, r, c;                     /* Zaehler */
  VRepnKind src_repn;
  int src_nrows, src_ncolumns, src_nbands; 
  float a[3][3],ainv[3][3],detA;
  float bp, rp, cp, bx, rx, cx;
  float shift[3];
  float val;
  int i,j;

/* ================================================================================================ */

  /* V O R B E R E I T U N G E N */

  /* Daten einlesen, Speicher allokieren, ... */

  /* Bildgroesse */
  src_nrows    = VImageNRows (src);
  src_ncolumns = VImageNColumns (src);
  src_nbands   = VImageNBands (src);
  src_repn     = VPixelRepn (src);
  
  /* Allokierung des dynamischen 3D - Koeffizientenfeldes */
  int_koeff = VCreateImage(src_nbands,src_nrows,src_ncolumns,VFloatRepn);
  VFillImage(int_koeff,VAllBands,0);
  
  /* Ausgabebild erzeugen (noch leer) */
  dest = VSelectDestImage("VCubicSplineSample3d",dest,dest_nbands,dest_nrows,dest_ncolumns,src_repn);
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);
  
  /* fertig mit den V O R B E R E I T U N G E N */
  
/* ------------------------------------------------------------------------------------------------- */

  /* get inverse of transformation matrix : */

  /* get transformation matrix : */
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      a[i][j]  = VGetPixel(transform,0,i,j+1);
    }
  }

  /* get its inverse : */
  ainv[0][0] =  a[1][1]*a[2][2] - a[1][2]*a[2][1];
  ainv[1][0] = -a[1][0]*a[2][2] + a[1][2]*a[2][0];
  ainv[2][0] =  a[1][0]*a[2][1] - a[1][1]*a[2][0];

  ainv[0][1] = -a[0][1]*a[2][2] + a[0][2]*a[2][1];
  ainv[1][1] =  a[0][0]*a[2][2] - a[0][2]*a[2][0];
  ainv[2][1] = -a[0][0]*a[2][1] + a[0][1]*a[2][0];

  ainv[0][2] =  a[0][1]*a[1][2] - a[0][2]*a[1][1];
  ainv[1][2] = -a[0][0]*a[1][2] + a[0][2]*a[1][0];
  ainv[2][2] =  a[0][0]*a[1][1] - a[0][1]*a[1][0];


  /* determinant */
  detA = a[0][0]*ainv[0][0] + a[0][1]*ainv[1][0] + a[0][2]*ainv[2][0];
  if (detA == 0) VError(" VCubicSplineSample3d: transformation matrix is singular");

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      ainv[i][j] /= detA;
    }
  }

  for (i=0; i<3; i++) 
    shift[i] = VPixel(transform,0,i,0,VFloat);


  /* Ermittlung der Splinekoeffizienten ... */

  /* erstmal die horizontal waagerechten Abhaengigkeiten */
  data = (float *) VMalloc (src_ncolumns * sizeof (float));  /* Stuetzstellen */
  for (b = 0; b < src_nbands; b++) {
    for (r = 0; r < src_nrows; r++) {

      /* Einlesen der Stuetzstellen (1 Zeile) */
      for (c = 0; c < src_ncolumns; c++) {
	data[c] = VGetPixel(src, b, r, c);
      }

      /* Berechnen der Koeffizienten */
      compute_koeffizients (data, src_ncolumns);

      /* Eintragen in 3D-Koeffizientenfeld */
      for (c = 0; c < src_ncolumns; c++)
        VPixel(int_koeff,b,r,c,VFloat) = data[c];
    }
  }
  VFree (data);
  
  /* jetzt die horizontal senkrechten Abhaengigkeiten */
  data = (float *) VMalloc (src_nrows * sizeof (float));  /* Stuetzstellen */
  for (b = 0; b < src_nbands; b++) {
    for (c = 0; c < src_ncolumns; c++) {

      /* Einlesen der Stuetzstellen (1 Spalte) */
      for (r = 0; r < src_nrows; r++)
        data[r] = VPixel(int_koeff,b,r,c,VFloat);

      /* Berechnen der Koeffizienten */
      compute_koeffizients (data, src_nrows);
    
      /* Eintragen in 3D-Koeffizientenfeld */
      for (r = 0; r < src_nrows; r++)
        VPixel(int_koeff,b,r,c,VFloat) = data[r];
    }
  }
  VFree (data);
    
  /* jetzt die vertikal senkrechten Abhaengigkeiten */
  data = (float *) VMalloc (src_nbands * sizeof (float)); /* Stuetzstellen */
  for (r = 0; r < src_nrows; r++) {
    for (c = 0; c < src_ncolumns; c++) {

      /* Einlesen der Stuetzstellen (1 Spalte) */    
      for (b = 0; b < src_nbands; b++)
        data[b] = VPixel(int_koeff,b,r,c,VFloat);

      /* Berechnen der Koeffizienten */
      compute_koeffizients (data, src_nbands);
    
      /* Eintragen in 3D-Koeffizientenfeld */
      for (b = 0; b < src_nbands; b++)
        VPixel(int_koeff,b,r,c,VFloat) = data[b];
    }
  }
  VFree (data);
    
  /* ... Splinekoeffizienten liegen vor */


  /* Determines the value of each pixel in the destination image: */

  for (b = 0; b < dest_nbands; b++) {
    for (r = 0; r < dest_nrows; r++) {
      for (c = 0; c < dest_ncolumns; c++) {

	/* Ermittlung der transformierten Koordinaten */

	bx = (float) b - shift[0]; 
	rx = (float) r - shift[1]; 
	cx = (float) c - shift[2];

	bp = ainv[0][0] * bx + ainv[0][1] * rx + ainv[0][2] * cx;
	rp = ainv[1][0] * bx + ainv[1][1] * rx + ainv[1][2] * cx;
	cp = ainv[2][0] * bx + ainv[2][1] * rx + ainv[2][2] * cx;

	bp += b0;
	rp += r0;
	cp += c0;
        
	if (bp < -1.0 || bp > src_nbands) continue;
	if (rp < -1.0 || rp > src_nrows) continue;
	if (cp < -1.0 || cp > src_ncolumns) continue;
	
	val = f(int_koeff, bp, rp, cp);
	VSetPixel(dest, b, r, c, (VDouble) val);
      }
    }
  }
  
  VDestroyImage(int_koeff);
  return dest;
}



/*
** use cubic spline interpolation for scaling 
** G.Lohmann, April 2003
*/

/*!
\fn VImage VCubicSplineScale3d (VImage src, VImage dest,int dest_nbands,int dest_nrows,int dest_ncols,
                                float shift[3],float scale[3])	     
\brief 3D scaling using cubic spline interpolation
\param src   input image (any repn)
\param dest  output image (any repn)
\param dest_nbands number of output slices
\param dest_nrows  number of output rows
\param dest_ncols  number of output columns
\param shift[3]   translation vector (band,row,column)
\param scale[3]   scaling vector (band,row,column)
*/
VImage 
VCubicSplineScale3d (VImage src,VImage dest,int dest_nbands,int dest_nrows,int dest_ncols,
		     float shift[3],float scale[3])
{
  VImage transform=NULL;
  VFloat b0,r0,c0;

  dest = VSelectDestImage("VCubicSplineScale3d",dest,
			  dest_nbands,dest_nrows,dest_ncols,VPixelRepn(src));
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);

  transform = VCreateImage(1,3,4,VFloatRepn);
  VFillImage(transform,VAllBands,0);

  VPixel(transform,0,0,1,VFloat) = scale[0];
  VPixel(transform,0,1,2,VFloat) = scale[1];
  VPixel(transform,0,2,3,VFloat) = scale[2];
 

  VPixel(transform,0,0,0,VFloat) = shift[0];
  VPixel(transform,0,1,0,VFloat) = shift[1];
  VPixel(transform,0,2,0,VFloat) = shift[2];
  
  b0 = r0 = c0 = 0;

  dest = VCubicSplineSample3d (src,dest,transform,b0,r0,c0,
			       dest_nbands,dest_nrows,dest_ncols);

  VCopyImageAttrs (src, dest);
  return dest;
}
