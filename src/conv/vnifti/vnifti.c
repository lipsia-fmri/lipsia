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
** read/write nifti-1 to vista format
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <nifti/nifti1.h>
#include <nifti/nifti1_io.h>

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

extern VAttrList Nifti1_to_Vista(char *hdr_file, char *data_file,VLong tr,VBoolean do_scaling);
extern void Vista_to_Nifti1(VAttrList list,VString filename);

/* get type of format (.nii or .v) */
int getformat(char *filename)
{
  int n = strlen(filename);
  int type = -1;
  if (filename[n-1] == 'i' && filename[n-2] == 'i' && filename[n-3] == 'n' && filename[n-4] == '.')
    type = 1;
  if (filename[n-1] == 'v' && filename[n-2] == '.') type = 0;
  return type;
}

int main(int argc,char *argv[])
{
  static VString in_filename = "";
  static VString out_filename = "";
  static VLong tr=0;
  static VBoolean scale = FALSE;
  static VOptionDescRec  options[] = {
    {"in",VStringRepn,1,(VPointer) &in_filename,VRequiredOpt,NULL,"input file"},
    {"out",VStringRepn,1,(VPointer) &out_filename,VRequiredOpt,NULL,"output file"},
    {"tr",VLongRepn,1,(VPointer) &tr,VOptionalOpt,NULL,"repetition time in milliseconds"},
    {"scale",VBooleanRepn,1,(VPointer) &scale,VOptionalOpt,NULL,"Whether to allow scaling"}
  };
  char *prg_name=GetLipsiaName("vnifti");
  fprintf(stderr, "%s\n", prg_name);

  /* parse command line */
  VParseFilterCmd (VNumber (options),options,argc,argv,NULL,NULL);


  /* get format types */
  int itype = getformat(in_filename);
  int otype = getformat(out_filename);


  /* nifti-1 to vista */
  if (itype == 1 && otype == 0) {
    VAttrList out_list = Nifti1_to_Vista(in_filename,in_filename,tr,scale);
    FILE *fp_out = VOpenOutputFile (out_filename, TRUE);
    if (! VWriteFile (fp_out, out_list)) exit (1);
    fclose(fp_out);
  }


  /* vista to nifti-1 */
  else if (itype == 0 && otype == 1) {
    FILE *fp = VOpenInputFile (in_filename, TRUE);
    VAttrList list = VReadFile (fp,NULL);
    if (! list) VError(" can't read input file");
    fclose(fp);
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
