/****************************************************************
 *
 * vpreprocess: MinVals.c
 *
 * Copyright (C) Max Planck Institute 
 * for Human Cognitive and Brain Sciences, Leipzig
 *
 * <lipsia@cbs.mpg.de>
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
 * $Id: MinVals.c 3181 2008-04-01 15:19:44Z karstenm $
 *
 *****************************************************************/

#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <stdio.h>
#include <string.h>


void
VApplyMinval(VAttrList list,VShort minval)
{
  VImage src=NULL;
  VAttrListPosn posn;
  int r,c,j,n,m;
  double sum,nx,u;


  /* apply minval */
  m = 0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src);
    if (VPixelRepn(src) != VShortRepn) continue;
    if (VImageNRows(src) < 2) continue;

    n = VImageNBands(src);

    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {

	sum = nx = 0;
	for (j=0; j<n; j++) {
	  u = (double)VPixel(src,j,r,c,VShort);
	  sum += u;
	  nx++;
	}

	if (sum/nx < minval) {
	  for (j=0; j<n; j++) VPixel(src,j,r,c,VShort) = 0;
	  continue;
	}
	m++;
      }
    }
  }
  if (m < 1) VError("No voxels above threshold 'minval', use lower threshold");
}
