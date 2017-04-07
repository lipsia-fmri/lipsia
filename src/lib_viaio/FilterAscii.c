/*
** remove spurious non-alphanumeric characters from an ASCII input file
**
** G.Lohmann, MPI-CBS, Jan 2005
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


void
VFilterAscii(FILE *fp_in,FILE *fp_out)
{
  int buf;

  rewind(fp_in);
  rewind(fp_out);
  while (!feof(fp_in)) {
    buf = fgetc(fp_in);
    if ((! isgraph(buf)) && (buf != '\n') && (buf != '\r')&& (buf != '\0')){
      buf = ' ';
    }
    if (buf == '\v') buf = ' '; /* remove tabs */
    if (buf == '\t') buf = ' ';
    fputc(buf,fp_out);
  }
  return;
}


/*
main()
{
  FILE *fp1,*fp2;
  int buf;

  fp1 = fopen("test1.txt","r");
  fp2 = fopen("test2.txt","w");
  VFilterAscii(fp1,fp2);
  fclose(fp1);
  fclose(fp2);
}
*/
