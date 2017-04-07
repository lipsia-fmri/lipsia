/*
 *  $Id: Param.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains routines for parsing and reporting parameter values.
 */

/*
 *  Copyright 1994 University of British Columbia
 *
 *  Permission to use, copy, modify, distribute, and sell this software and its
 *  documentation for any purpose is hereby granted without fee, provided that
 *  the above copyright notice appears in all copies and that both that
 *  copyright notice and this permission notice appear in supporting
 *  documentation. UBC makes no representations about the suitability of this
 *  software for any purpose. It is provided "as is" without express or
 *  implied warranty.
 *
 *  Author: Arthur Pope, UBC Laboratory for Computational Intelligence
 */

/* From the Vista library: */
#include "viaio/Vlib.h"
#include "viaio/file.h"
#include "viaio/mu.h"
#include "viaio/option.h"
#include "viaio/os.h"

/* From the standard C library: */
#include <ctype.h>


/* Later in this file: */
static VBoolean LoadParamDefnFile (int, VOptionDescRec [], FILE *,
				   VPointer, VBoolean);
static VBoolean LoadVistaDataFile (int, VOptionDescRec [], FILE *,
				   VStringConst, VPointer, VBoolean);
static VOptionDescRec *LookupParamName (int, VOptionDescRec[], VStringConst);
static VReadFileFilterProc Filter;
static VOptionDescRec *MakeOptions (int, VOptionDescRec [], VPointer);

/* Buffer larger enough to contain any input line: */
typedef char Buffer[256];


/*
 *  VLoadParameters
 *
 *  Load values for a set of parameters from either a Vista data file
 *  or a parameter definition file.
 */
 
VBoolean VLoadParameters (int nparams, VOptionDescRec params[],
			  VStringConst filename, VStringConst object,
			  VPointer values, VBoolean ignore_unrecog)
{
  VOptionDescRec *p;
  FILE *f;
  VBoolean result = TRUE, unknown = FALSE;
  int i;

  /* Ensure each parameter is defined as taking a single value: */
  for (i = 0, p = params; i < nparams; i++, p++)
    if (p->number != 1) {
      VWarning ("VLoadParameters: Parameter %s takes multiple values",
		p->keyword);
      result = FALSE;
    }
  if (! result)
    return FALSE;

  /* Open the file: */
  if (! (f = VOpenInputFile (filename, FALSE)))
    return FALSE;

  /* If the filename ends in .vpd, its a parameter definition file: */
  i = strlen (filename);
  if (i > 4 && strcmp (& filename[i - 4], ".vpd") == 0)
    result = LoadParamDefnFile (nparams, params, f, values,
				ignore_unrecog);

  /* If the filename ends in .v, its a Vista data file: */
  else if (i > 2 && strcmp (& filename[i - 2], ".v") == 0)
    result = LoadVistaDataFile (nparams, params, f, object, values,
				ignore_unrecog);

  else {
    VWarning ("VLoadParameters: Suffix of %s not recognized", filename);
    unknown = TRUE;
    result = FALSE;
  }
	
  if (! unknown && ! result)
    VWarning ("VLoadParameters: Failed to process %s", filename);
    
  fclose (f);
  return result;
}


/*
 *  LoadParamDefnFile
 *
 *  Load parameter values from a parameter definition file.
 */

static VBoolean LoadParamDefnFile (int nparams, VOptionDescRec params[],
				   FILE *f, VPointer values,
				   VBoolean ignore_unrecog)
{
  VBoolean result = TRUE, name_encountered = FALSE;
  int i, line = 0;
  VOptionDescRec *p;
  VDictEntry *dict;
  VRepnKind repn;
  Buffer buf, name, value, comment;

  /* Parse each line: */
  while (fgets (buf, sizeof (buf), f)) {
    line++;
	    
    /* Skip comment lines and empty lines: */
    if (buf[0] == '#')
      continue;
    for (i = 0; buf[i] && isspace (buf[i]); i++) ;
    if (! buf[i])
      continue;
	    
    /* Skip the first non-comment line since it is supposed to contain the
       name of the parameter set: */
    if (! name_encountered) {
      name_encountered = TRUE;
      continue;
    }

    /* Extract the parameter name, type, value, and comment: */
    if (! VParseParamDefn (buf, name, & repn, value, comment)) {
    Error:	    VWarning ("VLoadParameters: Error at line %d", line);
    result = FALSE;
    continue;
    }
	    
    /* Look up the name: */
    if (! (p = LookupParamName (nparams, params, name))) {
      if (ignore_unrecog)
	continue;
      VWarning ("VLoadParameters: Parameter %s not recognized", name);
      goto Error;
    }
	    
    /* Ensure it's of the expected type: */
    if (repn != p->repn) {
      VWarning ("VLoadParameters: Expected %s to be a %s, not a %s",
		name, VRepnName (p->repn), VRepnName (repn));
      goto Error;
    }

    /* Choose a keyword dictionary, defaulting in the case of a Boolean
       parameter: */
    dict = p->dict ? p->dict :
      (p->repn == VBooleanRepn ? VBooleanDict : NULL);

    /* Decode and store the new parameter value: */
    if (! VDecodeAttrValue (value, dict, p->repn,
			    ((char *) values + (int) p->value)))
      goto Error;
	
    /* If the representation is a VString, VDecodeAttrValue may have just
       stored a pointer to value[], which is on the stack. In that case,
       clone value[]: */
    if (p->repn == VStringRepn && * (VStringConst *) p->value == value)
      * (VStringConst *) p->value = VNewString (value);
  }

  /* Check that we're at EOF: */
  if (result && ! (result &= feof (f)))
    VWarning ("LoadParameters: EOF expected at line %d", line);
  return result;
}


/*
 *  LoadVistaDataFile
 *
 *  Load parameter values from a specified object of a Vista data file.
 */

static VBoolean LoadVistaDataFile (int nparams, VOptionDescRec params[],
				   FILE *f, VStringConst object,
				   VPointer values, VBoolean ignore_unrecog)
{
  VAttrList list, sublist;
  VBoolean result = TRUE;
  VAttrListPosn posn;
  VOptionDescRec *p;
  VDictEntry *dict;

  /* Read the data files contents, except for binary data: */
  if (! (list = VReadFile (f, & Filter)))
    return FALSE;

  /* Fetch a sublist of attributes associated with the object having the
     specified name: */
  switch (VGetAttr (list, object, NULL, VAttrListRepn, & sublist)) {
	
  case VAttrFound:
    break;
	
  case VAttrMissing:
    goto Done;

  case VAttrBadValue:
    VWarning ("VLoadParameters: "
	      "The %s attribute's value isn't an attribute list", object);
    result = FALSE;
    goto Done;
  }

  /* For each attribute found in the sublist: */
  for (VFirstAttr (sublist, & posn);
       VAttrExists (& posn); VNextAttr (& posn)) {

    /* Look up the parameter definition coinciding with the attribute's
       name: */
    if (! (p = LookupParamName (nparams, params, VGetAttrName (& posn)))) {
      if (! ignore_unrecog) {
	VWarning ("VLoadParameters: Parameter %s not recognized",
		  VGetAttrName (& posn));
	result = FALSE;
      }
      continue;
    }

    /* Choose a keyword dictionary, defaulting in the case of a Boolean
       parameter: */
    dict = p->dict ? p->dict :
      (p->repn == VBooleanRepn ? VBooleanDict : NULL);

    /* Fetch and store the parameter's value: */
    if (! VGetAttrValue (& posn, dict, p->repn,
			 ((char *) values + (int) p->value))) {
      VWarning ("VLoadParameters: %s attribute has bad value",
		p->keyword);
      result = FALSE;
    }
  }

 Done:
  VDestroyAttrList (list);
  return result;
}


/*
 *  LookupParamName
 *
 *  Lookup a parameter name in a table of parameter definitions.
 */

static VOptionDescRec *LookupParamName (int nparams, VOptionDescRec params[],
					VStringConst name)
{
  for ( ; nparams > 0; nparams--, params++)
    if (strcmp (params->keyword, name) == 0)
      return params;
  return NULL;
}


/*
 *  Filter invoked by VReadFile to choose binary data for loading.
 */

static VBoolean Filter (VBundle b, VRepnKind repn)
{
  return FALSE;	/* don't want any */
}


/*
 *  VParseParamDefn
 *
 *  Parse a parameter definition of the form
 *	<name> <repn> <value> <comment>
 *  If the value is a string (enclosed in double quotes), it is returned
 *  with the double quotes stripped and any backslash escape characters
 *  removed.
 */

VBoolean VParseParamDefn (VStringConst buf, VString name, VRepnKind *repn,
			  VString value, VString comment)
{
  Buffer repn_buf;
  int i;
  char ch;

#define Repn(root) { "V" #root, V ## root ## Repn }
  static struct {
    VStringConst name;
    VRepnKind repn;
  } repn_names[] = {
    Repn (Bit), Repn (UByte), Repn (SByte), Repn (Short), Repn (Long),
    Repn (Float), Repn (Double), Repn (Boolean), Repn (String)
  };
#undef Repn

  /* Parse the parameter's name and representation type: */
  if (sscanf (buf, " %s %s %n", name, repn_buf, & i) != 2) {
  Error:	VWarning ("VParseParamDefn: Bad parameter definition");
  return FALSE;
  }
  buf += i;
  for (i = 0; i < VNumber (repn_names); i++)
    if (strcmp (repn_buf, repn_names[i].name) == 0)
      break;
  if (i == VNumber (repn_names)) {
    VWarning ("VParseParamDefn: Invalid representation: %s", repn_buf);
    return FALSE;
  }
  *repn = repn_names[i].repn;

  /* Parse the value: */
  if (*buf == '"') {
    buf++;
    while ((ch = *buf++) && ch != '"') {
      if (ch == '\\') {
	ch = *buf++;
	if (! ch) {
	  VWarning ("VParseParamDefn: Unterminated string");
	  return FALSE;
	}
      }
      *value++ = ch;
    }
    *value = 0;
  } else {
    if (sscanf (buf, "%s%n", value, & i) != 1)
      goto Error;
    buf += i;
  }

  /* The rest of the line, with leading and trailing space trimmed off,
     comprises the comment: */
  while (isspace (*buf))
    buf++;
  for (i = strlen (buf); i > 0 && isspace (buf[i - 1]); i--) ;
  strncpy (comment, buf, i);
  comment[i] = 0;

  return TRUE;
}


/*
 *  VParseParamOptions
 */
 
VBoolean VParseParamOptions (int nparams, VOptionDescRec params[],
			     int *argc, char **argv, VPointer values)
{
  VOptionDescRec *options = MakeOptions (nparams, params, values);
  VBoolean result;

  result = VParseCommand (nparams, options, argc, argv);
  VFree (options);
  return result;
}


/*
 *  VPrintParameters
 *
 *  Print the current values of a set of parameters.
 */
 
void VPrintParameters (FILE *f, int nparams, VOptionDescRec params[],
		       VPointer values)
{
  VOptionDescRec *options = MakeOptions (nparams, params, values);

  VPrintOptions (f, nparams, options);
  VFree (options);
}


/*
 *  VReportValidParamOptions
 */
 
void VReportValidParamOptions (int nparams, VOptionDescRec params[],
			       VPointer values)
{
  VOptionDescRec *options = MakeOptions (nparams, params, values);

  VReportValidOptions (nparams, options);
  VFree (options);
}


/*
 *  MakeOptions
 *
 *  Make a temporary copy of the option/parameter table with the value
 *  pointers set appropriately.
 */

static VOptionDescRec *MakeOptions (int nparams, VOptionDescRec params[],
				    VPointer values)
{
  VOptionDescRec *options = VMalloc (nparams * sizeof (VOptionDescRec));
  int i;

  VCopy (params, options, nparams);
  for (i = 0; i < nparams; i++)
    options[i].value = (char *) values + (int) (options[i].value);
  return options;
}
