/*
 *  $Id: os.h 726 2004-03-08 13:12:45Z lohmann $
 *
 *  This header file tries to establish an environment containing the
 *  ANSI C, POSIX, and X/Open features that Vista relies on.
 *
 *  In a Standard C Conforming POSIX environment it is equivalent to
 *  including:  <float.h> <limits.h> <stdarg.h> <stdio.h> <stdlib.h>
 *		<string.h> <time.h> <unistd.h>
 *
 *  This file requires an ANSI C compiler.
 */

#ifndef V_os_h
#define V_os_h 1

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


/*
 *  Only supports ANSI C.
 */

#ifndef __STDC__
#error __FILE__ requires an ANSI C compiler
#endif


/*
 *  Header files that are on all platforms of interest.
 */

#define _XOPEN_SOURCE

#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>


/*
 *  Determine what platform we're on.
 */

#undef SunOS_4
#undef SunOS_5
#if defined (sun) || defined (__sun)
#ifdef SVR4
  #define SunOS_5
#else
  #define SunOS_4
#endif
#endif


/*
 *  Platform-specific definitions and header inclusions to get us
 *  up to ANSI C and POSIX.
 */

/* Define limits for floating point numbers, normally in <float.h>: */
#if defined (SunOS_4)					/* SunOS 4.1 */
  #define DBL_MAX (1.797693134862315708e+308)
  #define FLT_MAX ((float)3.40282346638528860e+38)
#else							/* ANSI C library */
  #include <float.h>
#endif

/* Define maximum number of bytes in pathname, normally in <limits.h>: */
#if defined (SunOS_4) && ! defined (PATH_MAX)		/* SunOS 4.1 */
  #define PATH_MAX _POSIX_PATH_MAX
#endif

/* Define function prototypes that are normally in <stdio.h>: */
#if defined (SunOS_4)
  extern int _filbuf (FILE *);
  extern int fclose (FILE *);
  extern int fflush (FILE *);
  extern int fgetc (FILE *);
  extern int fprintf (FILE *, const char *, ...);
  extern int fputc (int, FILE *);
  extern int fputs (const char *, FILE *);
  extern size_t fread (void *, size_t, size_t, FILE *);
  extern int fscanf (FILE *, const char *, ...);
  extern int fseek (FILE *, long int, int);
  extern size_t fwrite (const void *, size_t, size_t, FILE *);
  extern int printf (const char *, ...);
  extern void rewind (FILE *);
  extern int scanf (const char *, ...);
/* extern int sprintf (char *, const char *, ...); */
  extern int sscanf (const char *, const char *, ...);
  extern int ungetc (int, FILE *);
  extern int vsprintf (char *, const char *, va_list);
#endif

/* Define EXIT_SUCCESS, EXIT_FAILURE, which are normally in <stdlib.h>: */
#if defined (SunOS_4)					/* SunOS 4.1 */
  #define EXIT_SUCCESS 0
  #define EXIT_FAILURE 1
#endif

/* Get a definitions of strtox, which are normally in <stdlib.h>: */
#if defined (SunOS_4)					/* SunOS 4.1 */
  extern double strtod (const char *, char **);
  extern long strtol (const char *, char **, int);
#endif

/* Define strerror(), which is normally in <string.h>: */
#if defined (SunOS_4)					/* SunOS 4.1 */
  extern int sys_nerr;
  extern char *sys_errlist[];
  #define strerror(e)	((e) >= sys_nerr ? "Unknown" : sys_errlist[e])
#endif

/* Define memcpy(), etc., which are normally in <string.h>: */
#if defined (SunOS_4)					/* SunOS 4.1 */
  extern void *memchr (const void *, int, size_t);
  extern int memcmp (const void *, const void *, size_t);
  extern void *memcpy (void *, const void *, size_t);
  extern void *memmove (void *, const void *, size_t);
  extern void *memset (void *, int, size_t);
#endif

/* Get a definition of time(), which is normally in <time.h>: */
#if defined (SunOS_4)					/* SunOS 4.1 */
  extern time_t time (time_t *);
#endif


/*
 *  Definition of a random number generator.
 *
 *  On systems that have drand48(), use it; otherwise use rand().
 */

#if defined (sun) || defined (sgi) || defined (_XPG4)
  extern void srand48 (long);
  extern double drand48 (void);
  extern long mrand48 (void);
  #define VRandomSeed(seed)	srand48 ((long) seed)
  #define VRandomDouble()	((VDouble) drand48 ())
  #define VRandomLong()		((VLong) mrand48 ())
#else
  #define VRandomSeed(seed)	srand ((int) seed)
  #ifdef RAND_MAX
    #define VRandomDouble()	((unsigned int) rand () / (VDouble) RAND_MAX+1)
  #else
    #define VRandomDouble()	((unsigned int) rand () / (VDouble) 0x80000000)
  #endif
  #define VRandomLong()		((VLong) rand())
#endif


/*
 *  Other useful definitions.
 */

/* Include an identifying string in an object file without causing
   `unused variable' warnings: */
#define VRcsId(str) static char *rcsid = (0 ? (char *) & rcsid : (str))


#endif /* V_Vos_h */
