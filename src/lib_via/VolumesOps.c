/*! \file
  Operations on volume sets.

\par Reference:
G. Lohmann (1998). "Volumetric Image Analysis",
John Wiley & Sons, Chichester, England.

\par Author:
 Gabriele Lohmann, MPI-CBS
*/
/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>
#include <viaio/Volumes.h>

/* From the standard C libaray: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*!
\struct Volumes
\brief a linked list containing volumes

A Vista volume set is used to represent three-dimensional volumetric 
objects in 3D images in a three-dimensional run-length encoding.

A volume set in memory is referred to by the C pointer type Volumes. 
In data files it is  identified  by  the type name volumes, 
and in memory, by the VRepnKind code VolumesRepn. 
Since it is a standard object type with built-in support in the Vista library, 
volume sets can be read from  data  files,  written to data files, 
and manipulated as part of attribute lists by routines such as 
VReadFile(3Vi), VWriteFile(3Vi), and VGetAttr(3Vi).

A single volume, a member of a volume set, 
is referred to by the C pointer type Volume.

\par Author:
 Gabriele Lohmann, MPI-CBS
*/


/*!
\struct Volume
\brief a run-length encoding of a volumetric object
\param  short   label
\param  short   slice-dimension
\param  short   row-dimension
\param  short   column-dimension
\param  short   length of hash table (number of buckets) 
\param  short   total number of tracks in all buckets
\param  int     pointer to buckets
\param  VBucket pointer to next volume in linked list


The data struct 'volume' is used to represent  a single 3D object.
A linked list of single objects of type 'volume' make up a 'volumes'
data struct. The data struct uses 3D run-length encoding.

Use functions 'VAddVolume', 'VCreateVolume', 'VolumeExists',
'VFirstVolume', 'VNextVolume' to perform operations on volumes.
\par Author:
 Gabriele Lohmann, MPI-CBS
*/

/*!
\struct VTrack
\brief part of a run-length encoded volume

The data type 'VTrack' represents a single track of foreground voxels.
A linked list of tracks make up a volume. 

\par Author:
 Gabriele Lohmann, MPI-CBS
*/


#define IMAX(a,b) ((a) > (b) ? (a) : (b))

extern void VDestroyVolume(Volume);

/*!
  \fn Volumes VImage2Volumes(VImage src)
  \brief convert a raster image to volume repn. 
   Each grey value produces one volume.
  \param src   input image (ubyte or short repn)
*/
Volumes VImage2Volumes(VImage src)
{
  Volumes volumes;
  Volume *vol;
  VBit   *bin_pp;
  VUByte *ubyte_pp;
  VShort *short_pp;
  VInteger *int_pp;
  VLong *long_pp;
  VSByte *sbyte_pp;
  VTrack t;
  long label,maxlabel=0;
  long nbuckets,nbands,nrows,ncols;
  long *b0,*b1;
  long *r0,*r1;
  int i,b,r,c;
  VRepnKind repn;

  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);

  repn = VPixelRepn(src);
  if (repn == VFloatRepn || repn == VDoubleRepn) 
    VError("VImage2Volumes: illegal pixel repn");


  /* get number of labels */
  maxlabel = 0;
  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	label = VGetPixel(src,b,r,c);
	if (label > maxlabel) maxlabel = label;
      }
    }
  }
  maxlabel++;

  b0 = (long *) VCalloc(maxlabel,sizeof(long));
  b1 = (long *) VCalloc(maxlabel,sizeof(long));
  r0 = (long *) VCalloc(maxlabel,sizeof(long));
  r1 = (long *) VCalloc(maxlabel,sizeof(long));

  /*
  ** bounding boxes of all objects to determine optimal hash table length,
  ** and get largest label
  */
  for (i=0; i<maxlabel; i++) {
    b1[i] = 0;
    b0[i] = nbands;
    r1[i] = 0;
    r0[i] = nrows;
  }

  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	label = VGetPixel(src,b,r,c);
	if (b < b0[label]) b0[label] = b;
	if (b > b1[label]) b1[label] = b;
	if (r < r0[label]) r0[label] = r;
	if (r > r1[label]) r1[label] = r;
      }
    }
  }



  /*
  ** create volume set
  */
  vol = (Volume *) VMalloc(sizeof(VolumeRec) * (maxlabel));

  for (label = 1; label < maxlabel; label++) {
    vol[label] = NULL;
    nbuckets = (b1[label] - b0[label] + 1) * (r1[label] - r0[label] + 1);
    if (nbuckets < 1) continue;
    vol[label] = VCreateVolume(label,nbands,nrows,ncols,nbuckets);
    if (vol[label] == NULL) VError("error creating volume");
    vol[label]->ntracks = 0;
  }


  /*
  ** add tracks to hash table
  */
  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      c = 0;
      while (c < ncols) {
	label = VGetPixel(src,b,r,c);

	if (label != 0) {
	  t = (VTrack) VMalloc(sizeof(VTrackRec));
	  if (t == NULL)
	    VError("VImage2Volumes: error allocating track");

	  t->band = b;
	  t->row  = r;
	  t->col  = c;

	  if (repn == VBitRepn) {
	    bin_pp = (VBit *) VPixelPtr(src,b,r,c);
	    while ((c < ncols) && (*bin_pp++ > 0)) c++;
	  }
	  if (repn == VUByteRepn) {
	    ubyte_pp = (VUByte *) VPixelPtr(src,b,r,c);
	    while ((c < ncols) && (*ubyte_pp++ == label)) c++;
	  }
	  if (repn == VShortRepn) {
	    short_pp = (VShort *) VPixelPtr(src,b,r,c);
	    while ((c < ncols) && (*short_pp++ == label)) c++;
	  }
	  if (repn == VIntegerRepn) {
	    int_pp = (VInteger *) VPixelPtr(src,b,r,c);
	    while ((c < ncols) && (*int_pp++ == label)) c++;
	  }
	  else if (repn == VLongRepn) {
	    long_pp = (VLong *) VPixelPtr(src,b,r,c);
	    while ((c < ncols) && (*long_pp++ == label)) c++;
	  }
	  else {
	    sbyte_pp = (VSByte *) VPixelPtr(src,b,r,c);
	    while ((c < ncols) && (*sbyte_pp++ == label)) c++;
	  }
	  t->length = c - t->col;
	  AddTrack(vol[label],t);
	  c--;
	}
	c++;
      }
    }
  }
 
  /* add volumes to list of volumes */
  VFree(b0);
  VFree(b1);
  VFree(r0);
  VFree(r1);

  volumes = VCreateVolumes(nbands,nrows,ncols);

  for (label = 1; label < maxlabel; label++) {
    if (vol[label] == NULL) continue;
    if (vol[label]->ntracks > 0) {
      VAddVolume(volumes,vol[label]);
    }
  }
  VolumesAttrList(volumes) = VCopyAttrList (VImageAttrList (src));

  return volumes;
}



/*!
  \fn VImage Volumes2Image(Volumes src,VRepnKind repn)
  \brief convert a volume set to a raster image repn
  \param src   input volume set
  \param repn  output pixel repn (any repn).
*/
VImage
Volumes2Image(Volumes src,VRepnKind repn)
{
  Volume v;
  VTrack t;
  VImage dest;
  int c,c0,c1,i,nbands,nrows,ncols;
  VLong label;

  /*
  long maxlabel = (VLong) VRepnMaxValue(repn);
  long minlabel = (VLong) VRepnMinValue(repn);
  */

  nbands = src->nbands;
  nrows  = src->nrows;
  ncols  = src->ncolumns;

  dest = VCreateImage(src->nbands,src->nrows,src->ncolumns,repn);
  if (!dest) VError("Volumes2Image: Error creating output image");
  VFillImage(dest,VAllBands,0);

  for (v = src->first; v != NULL; v = v->next) {

    label = v->label;

    /*
    if (label < minlabel) {
      VWarning("clipping label values: %d set to %d",label,minlabel);
      label = minlabel;
    }
    else if (label > maxlabel) {
      VWarning("clipping label values: %d set to %d",label,maxlabel);
      label = maxlabel;
    }
    */

    for (i=0; i<v->nbuckets; i++) {
      for (t = v->bucket[i].first; t != NULL; t = t->next)  {
	if (t->band < 0 || t->band >= nbands) continue;
	if (t->row  < 0 || t->row  >= nrows) continue;

	c0 = t->col;
	c1 = t->col + t->length;
	if (c0 < 0) c0 = 0;
	if (c1 >= ncols) c1 = ncols - 1;

	for (c=c0; c <c1; c++) {
	  VSetPixel(dest,t->band,t->row,c,label);
	}
      } 
    }
  }
  VImageAttrList (dest) = VCopyAttrList (VolumesAttrList(src));
  return dest;
}



/*!
  \fn VImage Volumes2Bin(Volumes src)
  \brief convert a volume set to a binary raster image
  \param src   input volume set
*/
VImage
Volumes2Bin(Volumes src)
{
  Volume v;
  VTrack t;
  VImage dest;
  int c,c0,c1,i;

  dest = VCreateImage(src->nbands,src->nrows,src->ncolumns,VBitRepn);
  if (!dest) VError("Volumes2Bin: Error creating output image");
  VFillImage(dest,VAllBands,0);

  for (v = src->first; v != NULL; v = v->next) {
    for (i=0; i<v->nbuckets; i++) {
      for (t = v->bucket[i].first; t != NULL; t = t->next)  {
	c0 = t->col;
	c1 = t->col + t->length;

	for (c=c0; c <c1; c++) {
	  VPixel(dest,t->band,t->row,c,VBit) = 1;
	}
      } 
    }
  }

  VImageAttrList (dest) = VCopyAttrList (VolumesAttrList(src));
  return dest;
}


/*!
  \fn VImage Volume2Bin(Volume src)
  \brief convert a single volume to a binary raster image
  \param src   input volume
*/
VImage Volume2Bin(Volume v)
{
  VTrack t;
  VImage dest;
  long c,c0,c1,i;

  dest = VCreateImage(v->nbands,v->nrows,v->ncolumns,VBitRepn);
  if (!dest) VError("Volume2Bin: Error creating output image");
  VFillImage(dest,VAllBands,0);

  for (i=0; i<v->nbuckets; i++) {
    for (t = v->bucket[i].first; t != NULL; t = t->next)  {
      c0 = t->col;
      c1 = t->col + t->length;

      for (c=c0; c <c1; c++) {
	VPixel(dest,t->band,t->row,c,VBit) = 1;
      }
    } 
  }
  return dest;
}




/*!
  \fn Volume VBin2Volume(VImage src)
  \brief convert a binary raster image to a single volume representation
  \param src   input raster image (bit repn)
*/
Volume VBin2Volume(VImage src)
{
  Volume vol;
  VBit *srcpix;
  VTrack t;
  long len,nbands,nrows,ncols,label;
  long b0,b1,r0,r1,c0,c1;
  long b,r,c;

  /*
  ** get bounding box
  */
  nbands = VImageNBands(src);
  nrows = VImageNRows(src);
  ncols = VImageNColumns(src);
  if (VPixelRepn(src) != VBitRepn)
    VError(" input image must be of type bit.");

  b1 = 0;
  r1 = 0;
  c1 = 0;
  b0 = nbands;
  r0 = nrows;
  c0 = ncols;
  srcpix = (VBit *) VPixelPtr(src,0,0,0);
  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (*srcpix++ > 0) {
	  if (b < b0) b0 = b;
	  if (b > b1) b1 = b;
	  if (r < r0) r0 = r;
	  if (r > r1) r1 = r;
	  if (c < c0) c0 = c;
	  if (c > c1) c1 = c;
	}
      }
    }
  }

  /*
  ** create hash table
  */
  len = (b1 - b0 + 1) * (r1 - r0 + 1);
  label = 1;
  vol = VCreateVolume(label,nbands,nrows,ncols,len);
  if (!vol) VError("error creating volume");
  vol->label = 1;

  /*
  ** add tracks to hash table
  */
  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      c = 0;
      while (c < ncols) {
	if (VPixel(src,b,r,c,VBit) > 0) {
	  t = (VTrack) VMalloc(sizeof(VTrackRec));
	  t->band = b;
	  t->row  = r;
	  t->col  = c;
	  while ((c < ncols) && (VPixel(src,b,r,c,VBit)) > 0) c++;
	  t->length = c - t->col;
	  AddTrack(vol,t);
	}
	c++;
      }
    }
  }
  return vol;
}




/*!
\fn VBoolean VolumeInside(Volume v, long b, long r, long c)
\brief test if voxel [b,r,c] is inside volume v,
return TRUE if so, return FALSE, otherwise.
\param v input volume
\param b slice address of voxel
\param r row address of voxel
\param c column address of voxel
*/
VBoolean VolumeInside(Volume v, long b, long r, long c)
{
  long i;
  VTrack t;

  i = VolumeHash(v->nbands,b,r,v->nbuckets);

  for (t = v->bucket[i].first; t != NULL; t = t->next)  {
    if ((t->band == b) 
	&& (t->row == r)
	&& (t->col <= c) && (c < (t->col + t->length)))
      return TRUE;

    if ((t->band <= b) 
	&& (t->row <= r)
	&& (c >= (t->col + t->length)))
      return FALSE;
  }
  return FALSE;
}


/*
\fn VBoolean VolumeEqual(Volume v1, Volume v1) 
\brief test if volume <v1> equal volume <v2>,
return TRUE if so, return FALSE, otherwise.
\param v1 first volume
\param v2 second volume
*/
VBoolean VolumeEqual(Volume v1, Volume v2) 
{
  long i;
  VTrack t1,t2;

  if (VolumeNBuckets(v1) != VolumeNBuckets(v2)) return FALSE;
  if (VolumeNTracks(v1) != VolumeNTracks(v2)) return FALSE;

  for (i=0; i<VolumeNBuckets(v1); i++) {
    for (t1 = VFirstTrack(v1,i); VTrackExists(t1); t1 = VNextTrack(t1)) {
      for (t2 = VFirstTrack(v1,i); VTrackExists(t2); t2 = VNextTrack(t2)) {
	if (t1->band != t2->band) return FALSE;
	if (t1->row != t2->row) return FALSE;
	if (t1->col != t2->col) return FALSE;
	if (t1->length != t2->length) return FALSE;
      }
    }
  }
  return TRUE;
}



/*
** return ptr to track contraining voxel [b,r,c].
** return NULL, if no such track exists, i.e. if [b,r,c] 
** is not an interior point.
*/
VTrack VolumeGetTrack(Volume v, long b, long r, long c)
{
  long i;
  VTrack t;

  /* find track of given voxel */
  i = VolumeHash(v->nbands,b,r,v->nbuckets);
  for (t = v->bucket[i].first; t != NULL; t = t->next)  {
    if ((t->band == b) 
	&& (t->row == r)
	&& (t->col <= c) && (c < (t->col + t->length))) {
      return t;
    }

    if ((t->band <= b) 
	&& (t->row <= r)
	&& (c >= (t->col + t->length))) {
      return NULL;
    }
  }
  return NULL;
}


/*
**  test if voxel [b,r,c] is a border voxel using 
**  26-adjacency for interior voxels, and 6-adjacency "outside".
**  [b,r,c] is a border voxel, if there exists an exterior voxel
**  in its 6(?)-neighbourhood. 
*/
VBoolean VolumeBorder(Volume v, long b0, long r0, long c0)
{
  VTrack t;
  long b,r;

  t = VolumeGetTrack(v,b0,r0,c0);
  if (! t) return FALSE;           /* not an interior voxel */

  /* check northern voxel */
  if (t->col > c0-1) return TRUE;

  /* check southern voxel */
  if (c0+1 >= t->col + t->length) return TRUE;

  /* check top voxel */
  if (b0 > 0) {
    b = b0 - 1;
    if (! VolumeGetTrack(v,b,r0,c0)) return TRUE;
  }

  /* check bottom voxel */
  if (b0 < v->nbands - 1) {
    b = b0 + 1;
    if (! VolumeGetTrack(v,b,r0,c0)) return TRUE;
  }

  /* check eastern voxel */
  if (r0 > 0) {
    r = r0 - 1;
    if (! VolumeGetTrack(v,b0,r,c0)) return TRUE;
  }

  /* check western voxel */
  if (r0 < v->nrows - 1) {
    r = r0 + 1;
    if (! VolumeGetTrack(v,b0,r,c0)) return TRUE;
  }
  return FALSE;
}



/*
** count the number of border voxels in volume v
*/
double VolumeBorderSize(Volume v) 
{
  long b,r,c;
  double nborder;

  nborder = 0;
  for (b=0; b<v->nbands; b++) {
    for (r=0; r<v->nrows; r++) {
      for (c=0; c<v->ncolumns; c++)
	if (VolumeBorder(v,b,r,c)) nborder++;
    }
  }
  return nborder;
}



VTrack VCopyTrack(VTrack t)
{
  VTrack s;

  s = (VTrack) VMalloc(sizeof(VTrackRec));
  s->band   = t->band;
  s->row    = t->row;
  s->col    = t->col;
  s->length = t->length;
  return s;
}

/*
Volume
VCopyVolume(Volume src)
{
  Volume dest=NULL;
  int i;
  VTrack t,s;
  dest = VCreateVolume(src->label,src->nbands,src->nrows,src->ncolumns,src->nbuckets);

  for (i=0; i<src->nbuckets; i++) {
    for (t = VFirstTrack(src,i); VTrackExists(t); t = VNextTrack(t)) {
      s = VCopyTrack(t);
      AddTrack(dest,s);
    }
  }
  return dest;
}
*/


/*
** extract a slice from volume;
** add all tracks belonging to slice <band> to volume <dest>.
** <dest> must have been created prior to the call.
*/
void VolumeGetSlice(Volume src, Volume dest, int band)
{ 
  long r,i,len;
  VTrack t,s;

  if (!dest) return;

  len = src->nbuckets;

  for (r=0; r<src->nrows; r++) {
    i = VolumeHash(src->nbands,band,r,len);

    for (t = src->bucket[i].first; t != NULL; t = t->next) {
      if (t->band != band) continue;
      s = VCopyTrack(t);
      AddTrack(dest,s);
    }
  }
  return;
}


/*
** get radius**2 of smallest enclosing sphere
*/
double VolumeRadius(Volume vol,double *mean) 
{
  long i;
  double rmax,d;
  double u,v,w;
  VTrack t;

  rmax = 0;
  for (i=0; i<vol->nbuckets; i++) { 
    for (t = vol->bucket[i].first; t != NULL; t = t->next)  {

      u = t->band - mean[0];
      v = t->row  - mean[1];
      w = t->col  - mean[2];
      d = u * u + v * v + w * w;
      if (d > rmax) rmax = d;
      w = t->col + t->length  - mean[2];
      d = u * u + v * v + w * w;
      if (d > rmax) rmax = d;
    }
  }

  rmax = sqrt(rmax);
  return rmax;
}

