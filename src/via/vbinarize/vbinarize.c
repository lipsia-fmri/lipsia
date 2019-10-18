/****************************************************************
 *
 * Copyright (C) Max Planck Institute 
 * for Human Cognitive and Brain Sciences, Leipzig
 *
 * Author Gabriele Lohmann, 2004, <lipsia@cbs.mpg.de>
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
 * $Id: vbinarize.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *****************************************************************/

/*! \brief vbinarize -- threshold an image

\par Description
vbinarize thresholds an image, and produces a binary output image.

\par Usage

	<code>vbinarize</code>

	\param -in    input image
	\param -out   output image
	\param -min   lower threshold, values below <min> are set to zero.
	\param -max   upper threshold, values above <max> are set to zero.

\par Examples
  
\par Known bugs
none.

\file vbinarize.c
\author G.Lohmann, MPI-CBS, 2004
*/

#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <via/via.h>



int  main (int argc,char *argv[])
{  
  static VDouble xmin = -1.0e+16;
  static VDouble xmax = 1.0e+16;
  static VOptionDescRec  options[] = {
    {"min",VDoubleRepn,1,(VPointer) &xmin,VOptionalOpt,NULL,"lower threshold"},
    {"max",VDoubleRepn,1,(VPointer) &xmax,VOptionalOpt,NULL,"upper threshold"},
  };
  FILE *in_file,*out_file;
  VString in_filename=NULL;
  VAttrListPosn posn;
  VImage src=NULL,dest=NULL;
  char *prg=GetLipsiaName("vbinarize");
 
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);

  VAttrList list = VReadAttrListZ(in_file,in_filename,(VLong)0,TRUE,FALSE);
  if (list == NULL) VError(" error reading input list");

  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src);

    dest = VBinarizeImage (src,NULL,xmin,xmax);
    VSetAttrValue (& posn, NULL,VImageRepn,dest);
  }
  if (src == NULL) VError(" no input image found");


  VHistory(VNumber(options),options,prg,&list,&list);
  if (! VWriteFile (out_file, list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
