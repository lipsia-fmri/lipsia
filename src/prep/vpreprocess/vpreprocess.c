/****************************************************************
 *
 * vpreprocess: vpreprocess.c
 *
 * Copyright (C) Max Planck Institute
 * for Human Cognitive and Brain Sciences, Leipzig
 *
 * Authors Gabriele Lohmann, Karsten Mueller <lipsia@cbs.mpg.de>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * $Id: vpreprocess.c 3181 2008-04-01 15:19:44Z karstenm $
 *
 *****************************************************************/

#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fftw3.h>

#define ABS(x) ((x) > 0 ? (x) : -(x))
extern void Gauss4d(VAttrList list,double *sigma);



int main(int argc, char *argv[]) 
{
  static VDouble  fwhm     =  0;
  static VFloat   high     =  0;
  static VFloat   low      =  0;
  static VBoolean stop     =  FALSE;
  static VFloat   sharp    =  0.8;
  static VFloat   minval   =  0;
  static VOptionDescRec  options[] = {
    {"fwhm", VDoubleRepn, 1, (VPointer) &fwhm, VOptionalOpt, NULL, "Spatial filter: FWHM in mm"},
    {"high", VFloatRepn, 1, (VPointer) &high, VOptionalOpt, NULL, "Temporal Filter: Cutoff for high pass/stop in seconds"},
    {"low",  VFloatRepn, 1, (VPointer) &low, VOptionalOpt, NULL, "Temporal Filter: Cutoff for low pass/stop in seconds"},
    {"stop",  VBooleanRepn, 1, (VPointer) &stop, VOptionalOpt, NULL, "Temporal Filter: Stop insted of pass filter"},
    {"minval", VFloatRepn, 1, (VPointer) &minval, VOptionalOpt, NULL, "Signal threshold"}
  };
  FILE *out_file = NULL,*in_file=NULL;
  VString in_filename=NULL;

  extern void VSpatialFilter(VAttrList, VDouble);
  extern void VFreqFilter(VAttrList, VFloat, VFloat, VBoolean, VFloat);
  char *prg=GetLipsiaName("vpreprocess");
  fprintf (stderr, "%s\n", prg);


  VParseFilterCmdZ(VNumber(options), options, argc, argv, &in_file, &out_file,&in_filename);

  if (fwhm < 0) VError("fwhm must be non-negative");
  if (low > high && high > 0) VError("low must be less than high");
  if (low < 0 || high < 0) VError("high and low must be non-negative");



  /* read the file */
  VLong tr=0L;
  VAttrList list = VReadAttrListZ(in_file,in_filename,tr,TRUE,FALSE);
  if (list == NULL) VError(" error reading input file %s",in_file);


  /* apply minval */
  if (minval > 0) VApplyMinval(list, (VFloat)minval);


  /* apply freq filter */
  if (low > 0 || high > 0) {
    VFreqFilter(list, high, low, stop, sharp);
  }


  /* apply spatial filter */
  if (fwhm > 0) {
    int i;
    double *sigma = (double *)VCalloc(3,sizeof(double));
    double *pixdim = (double *)VCalloc(8,sizeof(double));
    for (i=0; i<8; i++) pixdim[i] = 1.0;
    VAttrList geolist = VGetGeoInfo(list);
    if (geolist != NULL) pixdim = VGetGeoPixdim(geolist,pixdim);
    for (i=0; i<3; i++) sigma[i] = fwhm/2.355;
    for (i=0; i<3; i++) sigma[i] /= pixdim[i+1];
    fprintf(stderr," resolution: %.3f %.3f %.3f\n",pixdim[1],pixdim[2],pixdim[3]);
    fprintf(stderr," sigma: %.3f %.3f %.3f\n",sigma[0],sigma[1],sigma[2]);
    Gauss4d(list,sigma);
  }


  /* write output */
  VHistory(VNumber(options), options, prg, &list, &list);
  if(! VWriteFile(out_file, list)) exit(1);
  fprintf(stderr, "%s: done.\n", argv[0]);
  exit(0);
}
