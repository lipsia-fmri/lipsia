/****************************************************************
 *
 * Copyright (C) Max Planck Institute
 * for Biological Cybernetics, Tuebingen, Germany
 *
 * Author: Gabriele Lohmann, 2015
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
 *****************************************************************/
/*
** read/write nifti-1 to/from vista format
*/
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>
#include <viaio/VImage.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

extern void Vista_to_Nifti1(VAttrList list,VString filename);
extern int getformat(char *filename);

int main(int argc,char *argv[])
{
  static VString out_filename = "";
  static VFloat xtr=0;
  static VBoolean attrtype = TRUE;
  static VBoolean do_scaling = FALSE;
  static VOptionDescRec  options[] = {
    {"out",VStringRepn,1,(VPointer) &out_filename,VRequiredOpt,NULL,"Output file"},
    {"tr",VFloatRepn,1,(VPointer) &xtr,VOptionalOpt,NULL,"Repetition time in seconds, use '0' to rely on header info"},
    {"attrtype",VBooleanRepn,1,(VPointer) &attrtype,VOptionalOpt,NULL,"Whether to output 4D data to lipsia 4D format"},
    {"scale",VBooleanRepn,1,(VPointer) &do_scaling,VOptionalOpt,NULL,"Whether to scale 4D data to 16bit integer"},
  };
  FILE *in_file=NULL;
  VString in_filename=NULL;
  char *prg_name=GetLipsiaName("vnifti");
  fprintf(stderr, "%s\n", prg_name);


  /* parse command line */
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,NULL,&in_filename);


  /* get TR */
  if (fabs(xtr) > 0 && xtr > 50) VError(" implausible TR (%.3f sec), TR must be specified in seconds",xtr);
  VLong tr=0;
  if (xtr > 0) tr = (VLong)(xtr*1000.0);


  /* get format types */
  int itype = getformat(in_filename);
  int otype = getformat(out_filename);


  /* read input data */
  VAttrList list = VReadAttrListZ(in_file,in_filename,(VLong)tr,attrtype,do_scaling);
  if (list == NULL) VError(" error reading input list");


  /* input from nifti-1, output to vista */
  if ((itype == 2 || itype == 3) && (otype == 0 || otype == 1)) {    
    FILE *fp_out = VOpenOutputFile (out_filename, TRUE);
    if (! VWriteFile (fp_out, list)) exit (1);
    fclose(fp_out);
  }


  /* vista to nifti-1 */
  else if ((itype < 2) && (otype == 2 || otype == 3)) {
    Vista_to_Nifti1(list,out_filename);
  }
  
  else if (itype == otype) {
    VError(" input equals output format");
  }
  
  else {
    VError(" illegal formats: %d %d\n",itype,otype);
  }

  fprintf (stderr, "%s: done.\n", argv[0]);
  exit(0);
}
