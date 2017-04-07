/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

/* From the standard C libaray: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define MAXPARAMLENGTH 500


char* VGetOptionValue (VOptionDescRec *option)
{
  char *ret;
  int n, i;
  char *vp;
  VDictEntry *dict;
  VLong ivalue;
  VDouble fvalue = 0.0;
  VStringConst svalue;
  VString ret1=NULL;

  ret = (char*)calloc(1000,sizeof(char));


  if (option->number == 0) {
    n = ((VArgVector *) option->value)->number;
    vp = (char *) ((VArgVector *) option->value)->vector;
  } else {
    n = option->number;
    vp = (char *) option->value;

  }
  for (i = 0; i < n; i++, vp += VRepnSize (option->repn)) {


    switch (option->repn) {
            
    case VBitRepn:
      ivalue = * (VBit *) vp;
      goto PrintLong;
            
    case VUByteRepn:
      ivalue = * (VUByte *) vp;
      goto PrintLong;
            
    case VSByteRepn:
      ivalue = * (VSByte *) vp;
      goto PrintLong;
            
    case VShortRepn:
      ivalue = * (VShort *) vp;
      goto PrintLong;

    case VUShortRepn:
      ivalue = * (VUShort *) vp;
      goto PrintLong;

    case VIntegerRepn:
      ivalue = * (VInteger *) vp;
      goto PrintLong;

    case VUIntegerRepn:
      ivalue = * (VUInteger *) vp;
      goto PrintLong;

    case VULongRepn:
      ivalue = * (VULong *) vp;
      goto PrintLong;
            
    case VLongRepn:
      ivalue = * (VLong *) vp;
    PrintLong: 
      sprintf (ret, "%ld", ivalue); 
      break;
            
    case VFloatRepn:
      fvalue = * (VFloat *) vp;
      goto PrintDbl;
            
    case VDoubleRepn:
      fvalue = * (VDouble *) vp;
    PrintDbl:   
      /* sprintf (ret, "%g", fvalue); */
      ret1 = VNewString(ret);
      sprintf (ret, "%s %g", ret1, fvalue);
      break;
            
    case VBooleanRepn:
      sprintf (ret, "%s", * (VBoolean *) vp ? "true" : "false");
      break;
    case VStringRepn:
      svalue = * (VString *) vp;
      if (! svalue)
        svalue = "(none)";
      else if (option->dict &&
               (dict = VLookupDictValue (option->dict, VDoubleRepn,
                                         svalue)))
        svalue = dict->keyword;
      sprintf (ret, "%s", svalue); 
      break;

    default:
      break;
    }
  }

  return ret;
}

char* 
VGetHistory(int noptions,VOptionDescRec *options,char *name) 
{
  int i,k, sw=0;
  char *history;
  char *item;
  char *ignore[] = { "in", "out", "ref", "contrast", NULL };

  history = (char*) malloc(sizeof(char*) * MAXPARAMLENGTH);
  history[0]='\0';
  strncat(history, name, strlen(name));
  strncat(history, " ", 1);

  for (i = 0; i < noptions; i++, options++) {
    sw=0;
    for (k=0; ignore[k] != NULL; k++) {
      if (strncmp(options->keyword, ignore[k], strlen(ignore[k])) == 0) sw=1;
    }
    if (sw==0) {
      item = VGetOptionValue (options);
      strncat(history, "-", 1);
      strncat(history, options->keyword, strlen(options->keyword));
      strncat(history, " ", 1);
      strncat(history, item, strlen(item));
      if (i+1<noptions) strncat(history, " ", 1);
    }
  }

  return history;
}


VAttrList VReadHistory(VAttrList *list)
{
  int i = 0;
  VAttrListPosn posn;
  VAttrList history_list=NULL;
  char history[]="history";

  for (VLastAttr((*list),&posn);VAttrExists(&posn);VPrevAttr(&posn)) {
    if (strncmp(VGetAttrName(&posn), history, strlen(history)) != 0 ) continue;
    if (VGetAttrRepn(&posn) == VAttrListRepn ) {
      if (i>0) VError("type mismatch while reading history");
      VGetAttrValue(&posn, NULL, VAttrListRepn, &history_list);
      break;
    }  
  }

  return history_list;
}

void
VPrependHistory(int noptions,VOptionDescRec *options,char *name,VAttrList *list)
{
  char *tok;
  char *newhistory;
  VString oname;

  /* Generate the new history entry */
  if ((newhistory = VGetHistory(noptions,options,name)) == NULL)
    VError("Error while building history string\n");

  tok = strtok(newhistory, " ");
  oname = malloc(strlen(tok)+1);
  strcpy(oname,tok);
  /* oname = (VStringConst)strdup(tok); */
  tok = strtok(NULL, "\0");

  /* Prepend history list */
  if ((*list) == NULL) (*list) = VCreateAttrList();
  VPrependAttr( (*list) ,oname, NULL, VStringRepn, tok); 
}

void  
VHistory(int noptions,VOptionDescRec *options,char *name,VAttrList *in_list,VAttrList *out_list) 
{
  VAttrList history_list=NULL;
  VAttrListPosn posn;
  char *history="history";

  /* Read history from list */
  history_list = VReadHistory(in_list);

  /* Prepend new history entry */
  VPrependHistory(noptions,options,name,&history_list);

  /* DELETE ANY history attributes in dest */
  for (VLastAttr((*out_list),&posn);VAttrExists(&posn);VPrevAttr(&posn)) {
    if (strncmp(VGetAttrName(&posn), history, strlen(history)) == 0 )
      VDeleteAttr(&posn);
  }

  /* Prepend history in dest */
  VPrependAttr( (*out_list),history,NULL,VAttrListRepn,history_list);
}
