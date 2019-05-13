/* From the standard C libaray: */
#include "viaio/Vlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *LipsiaVersion() 
{
  char *str="3.1.0";
  return str;
}


char *GetLipsiaName(char *prgname)
{
  int len1 = strlen(prgname);
  int len2 = strlen(LipsiaVersion());
  char *str = VCalloc(len1+len2+2,sizeof(char));
  sprintf(str,"%s V%s",prgname,LipsiaVersion());
  return str;
}
