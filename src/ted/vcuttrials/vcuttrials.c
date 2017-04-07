/****************************************************************
 *
 * Copyright (C) Max Planck Institute
 * for Biological Cybernetics, Tuebingen, Germany
 *
 * Author: Gabriele Lohmann, 2015
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
 *****************************************************************/

/*
** trial average using spline interpolation
**
** G.Lohmann, May 2014
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define ABS(x) ((x) > 0 ? (x) : -(x))
#define SQR(x) ((x)*(x))

#define LEN     10000   /* buffer length        */
#define NTRIALS 10000   /* max number of trials */

typedef struct SpointStruct {
  VShort x;
  VShort y;
  VShort z;
} SPoint;

typedef struct TrialStruct {
  int   id;
  float onset;
  float duration;
  float height;
} Trial;



Trial trial[NTRIALS];
int ntrials=0;
int nevents=0;

int test_ascii(int val)
{
  if (val >= 'a' && val <= 'z') return 1;
  if (val >= 'A' && val <= 'Z') return 1;
  if (val >= '0' && val <= '9') return 1;
  if (val ==  ' ') return 1;
  if (val == '\0') return 1;
  if (val == '\n') return 1;
  if (val == '\r') return 1;
  if (val == '\t') return 1;
  if (val == '\v') return 1;
  return 0;
}



/* parse design file */
void ReadDesign(VString designfile)
{
  FILE *fp=NULL;
  int  i,j,k,id;
  char buf[LEN];
  float onset=0,duration=0,height=0;

  fp = fopen(designfile,"r");
  if (!fp) VError(" error opening design file %s",designfile);

  i = ntrials = nevents = 0;
  while (!feof(fp)) {
    for (j=0; j<LEN; j++) buf[j] = '\0';
    if (fgets(buf,LEN,fp) == NULL) break;
    if (strlen(buf) < 2) continue;
    if (buf[0] == '%' || buf[0] == '#') continue;
    if (! test_ascii((int)buf[0])) VError(" input file must be a text file");

    /* remove non-alphanumeric characters */
    for (j=0; j<strlen(buf); j++) {
      k = (int)buf[j];
      if (!isgraph(k) && buf[j] != '\n' && buf[j] != '\r' && buf[j] != '\0') {
	buf[j] = ' ';
      }
      if (buf[j] == '\v') buf[j] = ' '; /* remove tabs */
      if (buf[j] == '\t') buf[j] = ' ';
    }

    if (sscanf(buf,"%d %f %f %f",&id,&onset,&duration,&height) != 4)
      VError(" line %d: illegal input format",i+1);

    if (duration < 0.5 && duration >= -0.0001) duration = 0.5;
    trial[i].id       = id;
    trial[i].onset    = onset;
    trial[i].duration = duration;
    trial[i].height   = height;
    i++;
    if (i > NTRIALS) VError(" too many trials %d",i);

    if (id > nevents) nevents = id;
  }
  fclose(fp);

  ntrials = i;
}

int main (int argc,char *argv[])
{
  static VString  designfile = "";
  static VShort   cond_id = 1;
  static VShort   trial_id = 0;
  static VDouble  temporal_resolution = 1.0;
  static VDouble  start = 0;
  static VDouble  length = 20;
  static VOptionDescRec  options[] = {
    {"design", VStringRepn, 1, & designfile, VRequiredOpt, NULL,"Design file (ascii)" },
    {"cond", VShortRepn, 1, & cond_id, VOptionalOpt, NULL,"Id of experimental condition"},
    {"trial", VShortRepn, 1, & trial_id, VOptionalOpt, NULL,"Id of trial (starts at 0)"},
    {"resolution", VDoubleRepn, 1, & temporal_resolution, VOptionalOpt, NULL," output temporal resolution in secs"},
    {"start", VDoubleRepn, 1, & start, VOptionalOpt, NULL, "start relative to beginning of trial in secs"},
    {"length", VDoubleRepn, 1, & length, VOptionalOpt, NULL, "trial length in seconds"},
  };
  FILE *in_file=NULL,*out_file=NULL;
  VAttrList list=NULL;
  VAttrListPosn posn;
  int i,j,k,slice,row,col;
  int nslices=0,nrows=0,ncols=0,ntimesteps=0;
  double t=0;
  char *prg = GetLipsiaName("vcuttrials");
  fprintf (stderr, "%s\n", prg);


  /*  parse command line */
  VParseFilterCmd (VNumber (options),options,argc,argv,&in_file,&out_file);


  /* get image dimensions, read functional data */
  if (! (list = VReadFile (in_file, NULL))) exit (1);
  fclose(in_file);

  nslices = 0;
  VImage tmp=NULL;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & tmp);
    if (VPixelRepn(tmp) != VShortRepn) continue;
    nslices++;
  }
  if (nslices < 1) VError(" no slices");


  VAttrList geoinfo = VGetGeoInfo(list);
  VImage *src = (VImage *) VCalloc(nslices,sizeof(VImage));
  i = ntimesteps = nrows = ncols = 0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src[i]);
    if (VPixelRepn(src[i]) != VShortRepn) continue;
    if (VImageNBands(src[i]) > ntimesteps) ntimesteps = VImageNBands(src[i]);
    if (VImageNRows(src[i])  > nrows)  nrows = VImageNRows(src[i]);
    if (VImageNColumns(src[i]) > ncols) ncols = VImageNColumns(src[i]);
    i++;
  }
  fprintf(stderr," nslices: %d, nrows: %d, ncols: %d,  ntimesteps: %d\n",nslices,nrows,ncols,ntimesteps);


  /* read repetition time */
  double tr = 0;
  if (VGetAttr (VImageAttrList (src[0]), "repetition_time", NULL,
		VDoubleRepn, (VPointer) & tr) != VAttrFound) {
    VError(" attribute 'repetition_time' missing");
  }
  tr /= 1000.0;
  double experiment_duration = tr * (double)ntimesteps;


  /* read design file */
  ReadDesign(designfile);
  

  /* select trial id */
  int tid=0,j0=-1;
  for (j=0; j<ntrials; j++) {
    if (trial[j].id == cond_id && tid == trial_id) {
      j0 = j;
      break;
    }
    if (trial[j].id == cond_id) tid++;
  }
  if (j0 < 0) VError("trial not found");
  fprintf(stderr," trial onset: %f\n",trial[j0].onset);



  /* adjust start time */
  trial[j0].onset += start;
  if (trial[j0].onset < 0) {
    VWarning(" negative onset, reset to zero");
    trial[j0].onset = 0;
  }
  if (trial[j0].onset + length >= experiment_duration) {
    VWarning(" experiment_duration exceeded",experiment_duration);
    length = experiment_duration - trial[j0].onset-0.5;
    VWarning(" parameter '-length' set to %f\n",length);
  }
  
  int nt = (int)(length/temporal_resolution + 0.5);
  fprintf(stderr," number of output timesteps: %d\n",nt);  


  /* ini output data structs  */
  VImage *dest = (VImage *) VCalloc(nslices, sizeof(VImage));
  VAttrList out_list = VCreateAttrList();
  double *D=NULL;
  if (geoinfo != NULL) {   
    D = VGetGeoDim(geoinfo,NULL);
    D[0] = (double)4.0;
    D[1] = (double)ncols;
    D[2] = (double)nrows;
    D[3] = (double)nslices;
    D[4] = (double)(nt);
    VSetGeoInfo(geoinfo,out_list);
  }
  for (slice=0; slice<nslices; slice++) {
    dest[slice] = VCreateImage(nt,nrows,ncols,VShortRepn);
    VFillImage(dest[slice],VAllBands,0);
    VCopyImageAttrs (src[slice], dest[slice]);
    VAppendAttr(out_list,"image",NULL,VImageRepn,dest[slice]);
  }

  
  /* spline interpolation  */
  double tstep = temporal_resolution;
  double *xx = (double *) VCalloc(ntimesteps,sizeof(double));
  double *yy = (double *) VCalloc(ntimesteps,sizeof(double));
 
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, ntimesteps);


  for (slice=0; slice<nslices; slice++) {
    for (row=0; row<nrows; row++) {
      for (col=0; col<ncols; col++) {

	for (k=0; k<ntimesteps; k++) {
	  xx[k] = tr*((double)k);
	  yy[k] = (double)VPixel(src[slice],k,row,col,VShort);
	}
	gsl_spline_init (spline, xx, yy, ntimesteps);

	k=0;
	for (t=0; t<length; t += tstep) {
	  if (k >= nt) break;
	  double xi = trial[j].onset + t;
	  double yi = 0;
	  if (xi < experiment_duration && xi >= 0) {
	    if (gsl_spline_eval_e (spline,xi,acc,&yi) == GSL_EDOM) continue;
	  }
	  VPixel(dest[slice],k,row,col,VShort) = yi;
	  k++;
	}
      }
    } 
  }


  /* write to disk */
  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  exit(0);
}
