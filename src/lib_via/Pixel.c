/*
** Contrast enhancement
**
** Author:
**  G.Lohmann <lohmann@cns.mpg.de>, May 1998
*/

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>

#include <stdio.h>
#include <stdlib.h>
#include <via/via.h>


VFloat VReadPixel(VImage src,int b,int r,int c)
{
  switch(VPixelRepn(src)) {
  case VUByteRepn:
    return (VFloat) VPixel(src,b,r,c,VUByte);
  case VShortRepn:
    return (VFloat) VPixel(src,b,r,c,VShort);
  case VFloatRepn:
    return (VFloat) VPixel(src,b,r,c,VFloat);
  case VBitRepn:
    return (VFloat) VPixel(src,b,r,c,VBit);
  case VSByteRepn:
    return (VFloat) VPixel(src,b,r,c,VSByte);
  case VDoubleRepn:
    return (VFloat) VPixel(src,b,r,c,VDouble);
  default:
    ;
  }
  return 0.0f;
}


void VWritePixel(VImage src,int b,int r,int c,VFloat val)
{
  switch(VPixelRepn(src)) {
  case VUByteRepn:
    VPixel(src,b,r,c,VUByte) = val;
    return;
  case VShortRepn:
    VPixel(src,b,r,c,VShort) = val;
    return;
  case VFloatRepn:
    VPixel(src,b,r,c,VFloat) = val;
    return;
  case VBitRepn:
    VPixel(src,b,r,c,VBit) = val;
    return;
  case VSByteRepn:
    VPixel(src,b,r,c,VSByte) = val;
    return;
  case VDoubleRepn:
    VPixel(src,b,r,c,VDouble) = val;
    return;
  default:
    ;
  }
}





VFloat 
VGetPixelValue (VImage image, int i)
{
  VPointer p = (VPointer) VPixelIndex (image, i);

  switch (VPixelRepn (image)) {

  case VUByteRepn:
    return (VFloat) * (VUByte *) p;

  case VShortRepn:
    return (VFloat) * (VShort *) p;

  case VBitRepn:
    return (VFloat) * (VBit *) p;

  case VFloatRepn:
    return (VFloat) * (VFloat *) p;

  case VSByteRepn:
    return (VFloat) * (VFloat *) p;

  default:
    ;
  }
  return 0.0;
}


void 
VSetPixelValue (VImage image, int i, VFloat value)
{
  VPointer p = VPixelIndex (image, i);

  switch (VPixelRepn (image)) {

  case VUByteRepn:
    * (VUByte *) p = value;
  break;

  case VShortRepn:
    * (VShort *) p = value;
  break;

  case VBitRepn:
    * (VBit *) p = value;
  break;

  case VFloatRepn:
    * (VFloat *) p = value;
  break;

  case VSByteRepn:
    * (VFloat *) p = value;
  break;

  default:
    ;
  }
}

