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


int main(int argc, char *argv[]) {
  static VDouble  fwhm     =  0;
  static VFloat   high     =  0;
  static VFloat   low      =  0;
  static VBoolean stop     =  FALSE;
  static VFloat   sharp    =  0.8;
  static VShort   minval   =  0;
  static VOptionDescRec  options[] = {
    {"fwhm", VDoubleRepn, 1, (VPointer) &fwhm, VOptionalOpt, NULL, "Spatial filter: FWHM in mm"},
    {"high", VFloatRepn, 1, (VPointer) &high, VOptionalOpt, NULL, "Temporal Filter: Cutoff for high pass/stop in seconds"},
    {"low",  VFloatRepn, 1, (VPointer) &low, VOptionalOpt, NULL, "Temporal Filter: Cutoff for low pass/stop in seconds"},
    {"stop",  VBooleanRepn, 1, (VPointer) &stop, VOptionalOpt, NULL, "Temporal Filter: Stop insted of pass filter"},
    {"minval", VShortRepn, 1, (VPointer) &minval, VOptionalOpt, NULL, "Signal threshold"}
  };
  FILE *out_file = NULL;
  VString in_file=NULL;


  extern void VSpatialFilter(VAttrList, VDouble);
  extern void VApplyMinval(VAttrList, VShort);
  extern void VFreqFilter(VAttrList, VFloat, VFloat, VBoolean, VFloat);
  char *prg=GetLipsiaName("vpreprocess");
  fprintf (stderr, "%s\n", prg);


  VParseFilterCmdX(VNumber(options), options, argc, argv, &in_file, &out_file);

  if (fwhm < 0) VError("fwhm must be non-negative");
  if (low > high && high > 0) VError("low must be less than high");
  if (low < 0 || high < 0) VError("high and low must be non-negative");
  if (sharp < 0.01) VError(" sharp too small", sharp);


  /* read the file */
  VLong tr=0L;
  VAttrList list = VReadAttrList(in_file,tr,TRUE,FALSE);
  if (list == NULL) VError(" error reading input file %s",in_file);


  /* apply filtering */
  if (low > 0 || high > 0) VFreqFilter(list, high, low, stop, sharp);
  if (fwhm > 0) VSpatialFilter(list, fwhm);
  if (minval > 0) VApplyMinval(list, minval);


  /* write output */
  VHistory(VNumber(options), options, prg, &list, &list);
  if(! VWriteFile(out_file, list)) exit(1);
  fprintf(stderr, "%s: done.\n", argv[0]);
  exit(0);
}
