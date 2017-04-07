/*
 *  $Id: FileIO.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains routines for reading attributes from files and writing
 *  attributes to files.
 */

/*
 *  Copyright 1993, 1994 University of British Columbia
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
#include "viaio/os.h"
#ifdef SupportUbcIff
#include "viaio/VImage.h"
#endif
#include "viaio/VList.h"

/* From the standard C library: */
#include <ctype.h>
#include <errno.h>

/* File identification string: */
VRcsId ("$Id: FileIO.c 3177 2008-04-01 14:47:24Z karstenm $");


/*
 *  Local definitions.
 */

/* Macro used in WriteFile, WriteAttrList, etc.: */
#define FailTest(put)	    if ((put) == EOF) goto Fail

/* Description of object with data block to be written later by VWriteFile: */
typedef struct {
    VAttrListPosn posn;		/* identify of object's attribute */
    VAttrList list;		/* attr list value referring to data */
    size_t length;		/* length of data block */
} DataBlock;

/* Local variables: */
static long offset;		/* current offset into file's binary data */
static VList data_list;		/* list of data blocks to write later */

/* Later in this file: */
static VBoolean ReadHeader (FILE *);
static VAttrList ReadAttrList (FILE *);
static char *ReadString (FILE *, char, VStringConst);
static VBoolean ReadDelimiter (FILE *);
static VBoolean ReadData (FILE *, VAttrList, VReadFileFilterProc *);
static VBoolean WriteAttrList (FILE *, VAttrList, int);
static VBoolean WriteAttr (FILE *, VAttrListPosn *, int);
static VBoolean WriteString (FILE *, const char *);
static VBoolean MySeek (FILE *, long);


/*
 *  VOpenInputFile, VOpenOutputFile
 *
 *  Open an input or output file, with "-" representing stdin or stdout.
 *  If nofail is TRUE, any failure is a fatal error.
 */

FILE *VOpenInputFile (VStringConst filename, VBoolean nofail)
{
    FILE *f;

    if (filename == NULL || strcmp (filename, "-") == 0)
	f = stdin;
    else if (! (f = fopen (filename, "r")))
	(nofail ? & VError : & VWarning)
	    ("Unable to open input file %s", filename);
    return f;
}

FILE *VOpenOutputFile (VStringConst filename, VBoolean nofail)
{
    FILE *f;

    if (filename == NULL || strcmp (filename, "-") == 0)
	f = stdout;
    else if (! (f = fopen (filename, "w")))
	(nofail ? & VError : & VWarning)
	    ("Unable to open output file %s", filename);
    return f;
}


/*
 *  VReadObjects
 *
 *  Read a Vista data file, extract object of a specified type, and return
 *  a vector of them plus a list of anything else found in the file.
 */

int VReadObjects (FILE *file, VRepnKind repn, VAttrList *attributes,
		  VPointer **objects)
{
    VAttrList list;
    VAttrListPosn posn;
    size_t i, nobjects = 0;
    VPointer *vector;

    /* Read the file's contents: */
    list = VReadFile (file, NULL);
    if (! list)
	return FALSE;

    /* Count the objects found: */
    for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn))
	nobjects += (VGetAttrRepn (& posn) == repn);
    if (nobjects == 0) {
	VWarning ("VReadObjects: No %s objects present in stream",
		  VRepnName (repn));
	VDestroyAttrList (list);
	return FALSE;
    }

    /* Allocate a vector of that many object pointers: */
    vector = VMalloc (nobjects * sizeof (VPointer));

    /* Extract the objects from the attribute list and place them in the
       vector: */
    for (VFirstAttr (list, & posn), i = 0; VAttrExists (& posn); )
	if (VGetAttrRepn (& posn) == repn) {
	    VGetAttrValue (& posn, NULL, repn, vector + i);
	    VDeleteAttr (& posn);
	    i++;
	} else VNextAttr (& posn);

    /* Return the objects and the remaining attributes: */
    *attributes = list;
    *objects = vector;
    return nobjects;
}


/*
 *  VReadFile
 *
 *  Read a Vista data file, returning an attribute list of its contents.
 */

VAttrList VReadFile (FILE *f, VReadFileFilterProc *filter)
{
    VAttrList list;
    size_t i;

    /* Ensure that the correct Vista data file header is there: */
    if (! ReadHeader (f))
	return NULL;

    /* Read all attributes in the file: */
    if (! (list = ReadAttrList (f)))
	return NULL;

    /* Swallow the delimiter and read the binary data following it: */
    offset = 0;
    if (! ReadDelimiter (f) || ! ReadData (f, list, filter)) {
	VDestroyAttrList (list);
	return NULL;
    }

    /* Now we should be at the end of the file: */
    i = fgetc (f);
    if (i != EOF) {
	ungetc (i, f);
	VWarning ("VReadFile: File continues beyond expected EOF");
    }
    return list;
}


/*
 *  ReadHeader
 *
 *  Read a Vista data file header.
 *  If the current version is 2 but the file is version 1, a warning message
 *  is produced but TRUE is still returned.
 */

static VBoolean ReadHeader (FILE *f)
{
    int version;

    if (fscanf (f, VFileHeader " %d", & version) != 1) {
	VWarning ("VReadFile: Vista data file header not found");
	return FALSE;
    }
	
    if (version == VFileVersion)
	return TRUE;

    if (version == 1 && VFileVersion == 2) {
	VWarning ("VReadFile: Obsolete data file -- pipe it thru v1to2");
	return TRUE;
    }

    VWarning ("VReadFile: Vista data file isn't version %d", VFileVersion);
    return FALSE;
}


/*
 *  ReadAttrList
 *
 *  Read a list of attributes from a stream.
 */

static VAttrList ReadAttrList (FILE *f)
{
    VAttrList sublist, list = VCreateAttrList ();
    VAttrRec *a;
    int ch = 0;
    size_t name_size;
    VBundle b;
    char buf[2], *str, name_buf[VMaxAttrNameLength + 1];

    /* Swallow a { marking the start of the attribute list: */
    if (fscanf (f, " %1s", buf) != 1 || buf[0] != '{') {
	VWarning ("VReadFile: Missing {");
	goto Fail;
    }
    
    /* For each attribute up to the next "}": */
    while (fscanf (f, " %255[^}: \t\n]", name_buf) == 1) {
	    name_size = strlen (name_buf);
	    
	    /* Read a : */
	    if (fscanf (f, " %1s", buf) != 1 || buf[0] != ':') {
		    VWarning ("VReadFile: Invalid %s attribute", name_buf);
		    goto Fail;
	    }
	    
	    do {
		    buf[0] = fgetc(f);
	    } while (buf[0] == ' ');
	    
	/* The first character of the value tells us whether its an attribute
	   list, quoted string, or unquoted string: */
	    if (buf[0] == '{') {
		    
		    /* The attribute value is another list of attributes: */
		    ungetc ('{', f);
		    if (! (sublist = ReadAttrList (f)))
			    goto Fail;
		    a = VMalloc (sizeof (VAttrRec) + name_size);
		    a->value = sublist;
		    a->repn = VAttrListRepn;
		    
	    } else if (buf[0] != '}' && buf[0] != '\n') {
		    
		    /* The value doesn't start with '{' -- parse a word or string: */
		    if (!(str = ReadString (f, buf[0], name_buf)))
			    goto Fail;
		    
		    while ((ch = fgetc (f)) && (ch == ' ' || ch == '\t')) ;
		    ungetc (ch, f);
		    
		    /* If the word is followed by an '{'... */
		    if (ch == '{') {
			    
			    /* ...then it's a typed value -- the word is it's type name
			       and the { is the start of it's attribute list value. */
			    b = VCreateBundle (str, NULL, 0, NULL);
			    if (! (sublist = ReadAttrList (f))) {
				    VFree (b);
				    goto Fail;
			    }
			    b->list = sublist;
			    a = VMalloc (sizeof (VAttrRec) + name_size);
			    a->repn = VBundleRepn;
			    a->value = b;
			    
		    } else {
			    
			    /* ...otherwise store it as a simple string value: */
			    a = VMalloc (sizeof (VAttrRec) + name_size + strlen (str) + 1);
			    a->repn = VStringRepn;
			    a->value = a->name + name_size + 1;
			    strcpy (a->value, str);
		    }
	    }else{
			ungetc (buf[0], f);
			a = VMalloc (sizeof (VAttrRec) + name_size +
				     strlen ("") + 1);
			a->repn = VStringRepn;
			a->value = a->name + name_size + 1;
			strcpy (a->value, "");
	    }
	    
	    /* Copy the attribute's name into the newly allocated node: */
	    strcpy (a->name, name_buf);
	    
	    /* Place the new node on the end of the growing attribute list: */
	    a->next = NULL;
	    if ((a->prev = list->prev))
		    a->prev->next = a;
	    else list->next = a;
	    list->prev = a;
    }

    /* Swallow the terminating "}": */
    if (fscanf (f, " %1s", buf) != 1 || buf[0] != '}') {
	VWarning ("VReadFile: Missing }");
Fail:	VDestroyAttrList (list);
	return NULL;
    }
    return list;
}


/*
 *  ReadString
 *
 *  Read an word, string, or number that may or may not be quoted.
 *  If ch is " then  parse a quoted string, otherwise parse to the next
 *  whitespace character. The string is returned in a reusable buffer.
 *
 *  The first several characters are read into a fixed-sized buffer
 *  large enough to hold most values. When it overflows, we begin
 *  allocating and reallocating larger buffers as we continue to read
 *  the value.
 */

#define StringAllocIncrement	100	/* each round of buffer resizing */

static char *ReadString (FILE *f, char ch, VStringConst name)
{
    VBoolean escaped = (ch == '"');
    size_t len = 0;
    char *cp;
    static size_t max_len;
    static char *buf;

    if (! buf) {
	buf = VMalloc (StringAllocIncrement);
	max_len = StringAllocIncrement;
    }

    if (! escaped)
	ungetc (ch, f);

    cp = buf;
    while (1) {
	ch = fgetc (f);

	/* Check for premature EOF: */
	if (ch == EOF) {
	    VWarning ("VReadFile: EOF encountered in %s attribute", name);
	    return NULL;
	}

	/* Check for closing " or escape sequence: */
	if (escaped) {
	    if (ch == '"')
		break;
	    if (ch == '\\')
		switch (ch = fgetc (f)) {

		case '\n':
		    continue;

		case 'n':
		    ch = '\n';
		}
	} else if (isspace (ch))
	    break;

	/* If the buffer in which we're accumulating the value is full,
	   allocate a larger one: */
	if (++len == max_len) {
	    buf = VRealloc (buf, max_len += StringAllocIncrement);
	    cp = buf + len - 1;
	}

	/* Store the character in the current buffer: */
	*cp++ = ch;
    }
    *cp = 0;

    /* Allocate a node of the correct size, or trim one already allocated so
       it is the correct size: */
    return buf;
}


/*
 *  ReadDelimiter
 *
 *  Read a Vista data file delimiter.
 */

static VBoolean ReadDelimiter (FILE *f)
{
    int ch;
    const char *cp;
    static char *msg = "VReadFile: Vista data file delimiter not found";

    /* Skip whitespace up to the first character of the delimeter: */
    while ((ch = fgetc (f)) != VFileDelimiter[0])
	if (ch == EOF || ! isspace (ch & 0x7F)) {
	    VWarning (msg);
	    return FALSE;
	}

    /* Swallow remaining characters of the delimiter: */
    for (cp = & VFileDelimiter[1]; *cp; cp++)
	if (*cp != fgetc (f)) {
	    VWarning (msg);
	    return FALSE;
	}

    return TRUE;
}


/*
 *  ReadData
 *
 *  Read the binary data accompanying attributes.
 */
static VBoolean ReadData (FILE *f, VAttrList list, VReadFileFilterProc *filter)
{
  VAttrListPosn posn, subposn;
  VAttrList sublist;
  VBundle b;
  VRepnKind repn;
  VBoolean read_data, data_found, length_found;
  VLong data, length;
  VTypeMethods *methods;
  VPointer value;

  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    switch (VGetAttrRepn (& posn)) {

    case VAttrListRepn:

      /* Recurse on nested attribute list: */
      VGetAttrValue (& posn, NULL, VAttrListRepn, & sublist);
      if (! ReadData (f, sublist, filter))
	return FALSE;
      break;

    case VBundleRepn:
      VGetAttrValue (& posn, NULL, VBundleRepn, & b);
      repn = VLookupType (b->type_name);

      /* If a filter routine was supplied, ask it whether to bother
	 with the binary data: */
      read_data = ! filter || (*filter) (b, repn);

      /* Extract any data and length attributes in the object's value: */
      if ((data_found = VLookupAttr (b->list, VDataAttr, & subposn))) {
	if (! VGetAttrValue (& subposn, NULL, VLongRepn, & data)) {
	  VWarning ("VReadFile: "
		    "%s attribute's data attribute incorrect",
		    VGetAttrName (& posn));
	  return FALSE;
	}
	/* fprintf(stderr," ReadData, data= %ld\n",data); */
	VDeleteAttr (& subposn);
      }
      if ((length_found = VLookupAttr (b->list, VLengthAttr, & subposn))) {
	if (! VGetAttrValue (& subposn, NULL, VLongRepn, & length)) {
	  VWarning ("VReadFile: "
		    "%s attribute's length attribute incorrect",
		    VGetAttrName (& posn));
	  return FALSE;
	}
	VDeleteAttr (& subposn);
      }

      /* None or both must be present: */
      if (data_found ^ length_found) {
	VWarning ("VReadFile: %s attribute has %s but not %s",
		  VGetAttrName (& posn),
		  data_found ? "data" : "length",
		  data_found ? "length" : "data");
	return FALSE;
      }

      /* Read the binary data associated with the object: */
      if (data_found) {
	if (data < offset) {
	  VWarning ("VReadFile: "
		    "%s attribute's data attribute incorrect",
		    VGetAttrName (& posn));
	  return FALSE;
	}

	if (! read_data)
	  data += length;

	/* To seek forward to the start of the data block we first
	   try fseek. That will fail on a pipe, in which case we
	   seek by reading. */
	if (data != offset &&
	    fseek (f, (long) data - offset, SEEK_CUR) == -1 &&
	    errno == ESPIPE &&
	    ! MySeek (f, data - offset)) {
	  VSystemWarning ("VReadFile: Seek within file failed");
	  return FALSE;				   
	}

	if (read_data) {
	  /* fprintf(stderr," ReadData, length= %ld\n",length); */
	  
	  b->data = VMalloc (b->length = length);
	  if (fread (b->data, 1, length, f) != length) {
	    VWarning ("VReadFile: Read from stream failed");
	    return FALSE;
	  }
	  offset = data + length;
	} else
	  /* bug: read error occured when bundle was not read
	     by a filter function. FK 24/03/98 */
	  offset = data;
      }

      /* Recurse to read binary data for sublist attributes: */
      if (! ReadData (f, b->list, filter))
	return FALSE;

      /* If the object's type is registered and has a decode method,
	 invoke it to decode the binary data: */
      if (read_data && repn != VUnknownRepn &&
	  (methods = VRepnMethods (repn)) && methods->decode) {
	if (! (value = (methods->decode) (VGetAttrName (& posn), b)))
	  return FALSE;

	/* Replace the old typed value with the newly decoded one: */
	VSetAttrValue (& posn, NULL, repn, value);
	VDestroyBundle (b);
      }
      break;

    default:
      break;
    }
  }
  return TRUE;
}


/*
 *  VWriteObjects
 *
 *  Write a list of objects, plus some other attributes, to a Vista data file.
 */

VBoolean VWriteObjects (FILE *file, VRepnKind repn, VAttrList attributes,
			int nobjects, VPointer objects[])
{
    VAttrList list;
    VAttrListPosn posn;
    int i;
    VBoolean result;

    /* Create an attribute list if none was supplied: */
    list = attributes ? attributes : VCreateAttrList ();

    /* Prepend to the attribute list an attribute for each object: */
    for (i = nobjects - 1; i >= 0; i--)
	VPrependAttr (list, VRepnName (repn), NULL, repn, objects[i]);

    /* Write the attribute list: */
    result = VWriteFile (file, list);

    /* Remove the attributes just prepended: */
    VFirstAttr (list, & posn);
    for (i = 0; i < nobjects; i++)
	VDeleteAttr (& posn);
    if (list != attributes)
	VDestroyAttrList (list);

    return result;
}


/*
 *  VWriteFile
 */

VBoolean VWriteFile (FILE *f, VAttrList list)
{
    DataBlock *db;
    VBundle b;
    VTypeMethods *methods;
    VRepnKind repn;
    VPointer value, ptr;
    VBoolean result, free_it;

    /* Write the Vista data file header, attribute list, and delimeter
       while queuing on data_list any binary data blocks to be written: */
    offset = 0;
    data_list = VListCreate ();
    FailTest (fprintf (f, "%s %d ", VFileHeader, VFileVersion));
    if (! WriteAttrList (f, list, 1)) {
	VListDestroy (data_list, VFree);
	return FALSE;
    }
    FailTest (fputs ("\n" VFileDelimiter, f));
    fflush (f);

    /* Traverse data_list to write the binary data blocks: */
    for (db = VListFirst (data_list); db; db = VListNext (data_list)) {
	repn = VGetAttrRepn (& db->posn);
	if (repn == VBundleRepn) {

	    /* A typed value includes its binary data block explicitly: */
	    VGetAttrValue (& db->posn, NULL, VBundleRepn, & b);
	    ptr = b->data;
	    free_it = FALSE;

	} else {

	    /* For any other representation, obtain the binary data block
	       from its encode_data method: */
	    VGetAttrValue (& db->posn, NULL, repn, & value);
	    methods = VRepnMethods (repn);
	    ptr = (methods->encode_data)
		(value, db->list, db->length, & free_it);
	    if (! ptr)
		goto Fail;
	}

	/* Write the binary data and free the buffer containing it if it was
	   allocated temporarily by an encode_data method: */
	if (db->length > 0) {
	    result = fwrite (ptr, 1, db->length, f) == db->length;
	    if (free_it)
		VFree (ptr);
	    if (! result)
		goto Fail;
	}
    }
    VListDestroy (data_list, VFree);
    return TRUE;

Fail:
    VWarning ("VWriteFile: Write to stream failed");
    VListDestroy (data_list, VFree);
    return FALSE;
}


/*
 *  WriteAttrList
 *
 *  Write a list of attributes to a file. Lines are indented by indent tabs.
 */

static VBoolean WriteAttrList (FILE *f, VAttrList list, int indent)
{
    VAttrListPosn posn;
    int i;

    /* Write the { marking the beginning of the attribute list: */
    FailTest (fputs ("{\n", f));

    /* Write each attribute in the list: */
    for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn))
	if (! WriteAttr (f, & posn, indent))
	    return FALSE;

    /* Write the } marking the end of the attribute list: */
    for (i = indent - 1; i > 0; i--)
	FailTest (fputc ('\t', f));
    FailTest (fputc ('}', f));
    return TRUE;

Fail:
    VWarning ("VWriteFile: Write to stream failed");
    return FALSE;
}


/*
 *  WriteAttr
 *
 *  Write a single attribute to a file. If the attribute's value is
 *  itself an attribute list, it is indented by indent+1 tabs.
 */

static VBoolean WriteAttr (FILE *f, VAttrListPosn *posn, int indent)
{
    int i;
    char *str;
    VRepnKind repn;
    VAttrList sublist;
    VBundle b;
    DataBlock *db;
    VTypeMethods *methods;
    size_t length;
    VPointer value;
    VBoolean result;
    VAttrListPosn subposn;

    /* Indent by the specified amount: */
    for (i = 0; i < indent; i++)
	FailTest (fputc ('\t', f));
    indent++;

    /* Output the attribute's name: */
    FailTest (fprintf (f, "%s: ", VGetAttrName (posn)));

    /* Ouput its value: */
    switch (repn = VGetAttrRepn (posn)) {

    case VAttrListRepn:
	VGetAttrValue (posn, NULL, VAttrListRepn, (VPointer) & sublist);
	result = WriteAttrList (f, sublist, indent);
	break;

    case VBundleRepn:
	VGetAttrValue (posn, NULL, VBundleRepn, (VBundle) & b);
	if (! WriteString (f, b->type_name))
	    return FALSE;
	FailTest (fputc (' ', f));

	/* If it's a typed value with binary data... */
	if (b->length > 0) {

	    /* Include "data" and "length" attributes in its attribute list: */
	    VPrependAttr (b->list, VLengthAttr, NULL, VLongRepn,
			  (VLong) b->length);
	    VPrependAttr (b->list, VDataAttr, NULL, VLongRepn,
			  (VLong) offset);

	    /* Add it to the queue of binary data blocks to be written: */
	    offset += b->length;
	    db = VNew (DataBlock);
	    db->posn = *posn;
	    db->list = b->list;
	    db->length = b->length;
	    VListAppend (data_list, db);
	}

	/* Write the typed value's attribute list: */
	result = WriteAttrList (f, b->list, indent);

	/* Remove the "data" and "length" attributes added earlier: */
	if (b->length > 0) {
	    VFirstAttr (b->list, & subposn);
	    VDeleteAttr (& subposn);
	    VDeleteAttr (& subposn);
	}
	break;

    case VStringRepn:
	VGetAttrValue (posn, NULL, VStringRepn, (VPointer) & str);
	result = WriteString (f, str);
	break;

    default:
	if (! (methods = VRepnMethods (repn)) ||
	    ! methods->encode_attr || ! methods->encode_data) {
	    VWarning ("VWriteFile: "
		      "%s attribute has unwriteable representation: %s",
		      VGetAttrName (posn), VRepnName (repn));
	    return FALSE;
	}

	/* Write the type name: */
	if (! WriteString (f, VRepnName (repn)))
	    return FALSE;
	FailTest (fputc (' ', f));

	/* Invoke the object type's encode_attr method to obtain an
	   attribute list: */
	VGetAttrValue (posn, NULL, repn, & value);
	sublist = (methods->encode_attr) (value, & length);

	/* If binary data is indicated... */
	if (length > 0) {

	    /* Include "data" and "length" attributes in the attr list: */
	    VPrependAttr (sublist, VLengthAttr, NULL, VLongRepn,
			  (VLong) length);
	    VPrependAttr (sublist, VDataAttr, NULL, VLongRepn,
			  (VLong) offset);

	    offset += length;
	}

	/* Add the object to the queue of binary data blocks to be written: */
	db = VNew (DataBlock);
	db->posn = *posn;
	db->list = sublist;
	db->length = length;
	VListAppend (data_list, db);

	/* Write the typed value's attribute list: */
	result = WriteAttrList (f, sublist, indent);

	/* Remove the "data" and "length" attributes added earlier: */
	if (length > 0) {
	    VFirstAttr (sublist, & subposn);
	    VDeleteAttr (& subposn);
	    VDeleteAttr (& subposn);
	}
    }

    /* Output a trailing newline: */
    if (result)
	FailTest (fputc ('\n', f));
    return result;

Fail:
    VWarning ("VWriteFile: Write to stream failed");
    return FALSE;
}


/*
 *  WriteString
 *
 *  Write a string, quoting it if necessary.
 */

static VBoolean WriteString (FILE *f, const char *str)
{
    const char *cp;
    int ch;

    /* Test for the presence of funny characters in the string value: */
    for (cp = str; (ch = *cp); cp++)
	if (! isalnum (ch) && ch != '.' && ch != '-' && ch != '+' && ch != '_')
	    break;

    /* If funny characters are present, output the string in quotes: */
    if (ch) {
	FailTest (fputc ('"', f));
	for (cp = (char *) str; (ch = *cp); cp++)
	    switch (ch) {

	    case '\n':
		FailTest (fputs ("\\n\\\n", f));
		break;

	    case '"':
		FailTest (fputs ("\\\"", f));
		break;

	    default:
		FailTest (fputc (ch, f));
	    }
	FailTest (fputc ('"', f));
    } else FailTest (fputs (str, f));
    return TRUE;

Fail:
    VWarning ("VWriteFile: Write to stream failed");
    return FALSE;
}


/*
 *  MySeek
 *
 *  Seek forward in a stream by reading it. Needed because fseek() won't
 *  work on pipes.
 */

static VBoolean MySeek (FILE *f, long bytes)
{
    size_t len;
    char buf[1000];

    while (bytes > 0) {
	len = VMin (bytes, sizeof (buf));
	if (fread (buf, 1, len, f) != len)
	    return FALSE;
	bytes -= len;
    }
    return TRUE;
}
