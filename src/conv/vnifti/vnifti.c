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
#include <sys/stat.h>

#include <nifti/nifti1.h>
#include <nifti/nifti1_io.h>
#include <zlib.h>

extern char *VReadDataContainer(char *filename,VBoolean nofail,size_t *size);
extern FILE *VOpenStream(char *databuffer,size_t size);
extern void Vista_to_Nifti1(VAttrList list,VString filename);
extern VAttrList Nifti1_to_Vista(char *databuffer,VLong tr,VBoolean attrtype,VBoolean do_scaling,VBoolean *ok);


int main(int argc,char *argv[])
{
  static VString in_filename = "";
  static VString out_filename = "";
  static VLong tr=0;
  static VBoolean attrtype = TRUE;
  static VBoolean do_scaling = FALSE;
  static VOptionDescRec  options[] = {
    {"in",VStringRepn,1,(VPointer) &in_filename,VRequiredOpt,NULL,"Input file"},
    {"out",VStringRepn,1,(VPointer) &out_filename,VRequiredOpt,NULL,"Output file"},
    {"tr",VLongRepn,1,(VPointer) &tr,VOptionalOpt,NULL,"Repetition time in milliseconds, use '0' to rely on header info"},
    {"attrtype",VBooleanRepn,1,(VPointer) &attrtype,VOptionalOpt,NULL,"Whether to output 4D data to lipsia 4D format"},
    {"scale",VBooleanRepn,1,(VPointer) &do_scaling,VOptionalOpt,NULL,"Whether to scale 4D data to 16bit integer"},
  };
  char *prg_name=GetLipsiaName("vnifti");
  fprintf(stderr, "%s\n", prg_name);


  /* parse command line */
  VParseFilterCmd (VNumber (options),options,argc,argv,NULL,NULL);


  /* get format types */
  if (CheckGzip(out_filename)) VError(" output file cannot be gzipped");
  int itype = getformat(in_filename);
  if (itype < 0 || itype > 1) VError(" illegal input type");
  int otype = getformat(out_filename);
  if (otype < 0 || otype > 1) VError(" illegal output type");
  if (itype == otype) VError(" input equals output format");
  


  /* read data from input file */
  size_t bufsize=0;
  char *databuffer = VReadDataContainer(in_filename,TRUE,&bufsize);


  /* nifti-1 to vista */
  if (itype == 1 && otype == 0) {    
    VBoolean ok=FALSE;
    VAttrList out_list = Nifti1_to_Vista(databuffer,tr,attrtype,do_scaling,&ok);
    if (ok) VFree(databuffer);
    FILE *fp_out = VOpenOutputFile (out_filename, TRUE);
    if (! VWriteFile (fp_out, out_list)) exit (1);
    fclose(fp_out);
  }


  /* vista to nifti-1 */
  else if (itype == 0 && otype == 1) {
    FILE *fp = VOpenStream(databuffer,bufsize);
    VAttrList list = VReadFile (fp,NULL);
    if (! list) VError(" can't read input file");
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
