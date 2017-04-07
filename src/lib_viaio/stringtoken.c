/*
** replaces 'strtok'
**
** G.Lohmann, Aug. 2005
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>


int VStringToken (char *buf,char *token,int index,int tlen)
{
  int i,j,len,count=0;

  len = 0;
  while (buf[len] != '\0' && buf[len] != '\n' && buf[len] != '\r') len++;

  i = count = 0;
  while (i < len) {

    j = 0;
    while (! isspace (buf[i]) && i < len && j < tlen) token[j++] = buf[i++];

    if (j > 0) {
      token[j] = '\0';
      if (count == index) return j;
      count++;
    }
    if (count > index) return 0;

    i++;
  }
  return 0;
}


/*
int
main(int argc, char *argv[])
{
  char *buf = "  , a 22 333xxx 4  555     zzz 1234567890  \r";
  char token[6];
  int i;
   

  i = 0;
  while (VToken(buf,token,i,6)) {
    fprintf(stderr," %3d:  %s\n",i,token);
    i++;
    if (i > 20) break;
  }
  exit(0);
}
*/
