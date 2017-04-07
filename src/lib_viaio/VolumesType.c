/*
 *
 *  This file contains methods for the volume set (Volumes) type.
 *
 *  Author: Gabriele Lohmann, MPI/CNS, Nov. 1995
 */

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/os.h>
#include <viaio/Volumes.h>


/*
 *  Table of methods.
 */

/* Later in this file: */
static VDecodeMethod VolumesDecodeMethod;
static VEncodeAttrMethod VolumesEncodeAttrMethod;
static VEncodeDataMethod VolumesEncodeDataMethod;

/* Used in Type.c to register this type: */
VTypeMethods VolumesMethods = {
  (VCopyMethod *) VCopyVolumes,	        /* copy a Volumes */
  (VDestroyMethod *) VDestroyVolumes,   /* destroy a Volumes */
  VolumesDecodeMethod,		        /* decode a Volumes's value */
  VolumesEncodeAttrMethod,		/* encode a Volumes's attr list */
  VolumesEncodeDataMethod		/* encode a Volumes's binary data */
};


/*
 *  VolumesDecodeMethod
 *
 *  The "decode" method registered for the "Volumes" type.
 *  Convert an attribute list plus binary data to a Volumes object.
 */

static VPointer VolumesDecodeMethod (VStringConst name, VBundle b)
{
  Volumes volumes;
  VAttrList list;
  long nels,i,j;
  size_t length;
  VShort hashlen,nbands,nrows,ncols;
  long nt1,nt2;
  long label;
  VLong nvolumes;
  VLong  ntracks;
  VLong *p;
  VTrack t;
  VPointer data;
  Volume vol;

#define Extract(name, dict, locn, required)	\
  VExtractAttr (b->list, name, dict, VShortRepn, & locn, required)

    /* Extract the required attribute values for Volumes. */
    if (! Extract (VolNVolumesAttr, NULL, nvolumes, TRUE) ||
	! Extract (VolNBandsAttr, NULL, nbands, TRUE) ||
	! Extract (VolNRowsAttr, NULL, nrows, TRUE) ||
	! Extract (VolNColumnsAttr, NULL, ncols, TRUE))
      return NULL;
  if (nvolumes <= 0) {
    VWarning ("VolumesReadDataMethod: Bad Volumes file attributes");
    return NULL;
  }

  /* Create the Volumes data structure. */
  volumes = VCreateVolumes (nbands,nrows,ncols);
  if (! volumes)
    return NULL;

  /* Give it whatever attributes remain: */
  list = VolumesAttrList (volumes);
  VolumesAttrList (volumes) = b->list;
  b->list = list;

  /* check amount of binary data: */
  nels =  b->length / (VRepnPrecision (VShortRepn) / 8);
  length = nels * (VRepnPrecision (VShortRepn) / 8);

  /* Allocate storage for the Volumes binary data: */
  data = VMalloc(nels * sizeof(VShort));

  /* Unpack the binary data: */
  if (! VUnpackData (VShortRepn, nels, b->data, VMsbFirst, & length,
		     & data, NULL))
    VError("error unpacking data");
  
  p = (VLong *) data;
  for (i = 0; i < nvolumes; i++) {
    nt1 = *p++;
    nt2 = *p++;
    ntracks = (long) nt1 + (long) nt2 - (long) VRepnMinValue(VShortRepn);
    hashlen = *p++;
    nbands  = *p++;
    label   = *p++;
    vol     = VCreateVolume(label,nbands,nrows,ncols,hashlen);
    for (j=0; j<ntracks; j++) {
      t = VMalloc(sizeof(VTrackRec));
      t->band   = *p++;
      t->row    = *p++;
      t->col    = *p++;
      t->length = *p++;
      AddTrack(vol,t);
    }
    VAddVolume(volumes,vol);
  }
  return volumes;

#undef Extract
}


/*
 *  VolumesEncodeAttrMethod
 *
 *  The "encode_attrs" method registered for the "Volumes" type.
 *  Encode an attribute list value for a Volumes object.
 */

static VAttrList VolumesEncodeAttrMethod (VPointer value, size_t *lengthp)
{
  Volumes volumes = value;
  Volume v;
  VAttrList list;
  VLong ntracks;

  /* Temporarily prepend several attributes to the edge set's list: */
  if ((list = VolumesAttrList (volumes)) == NULL)
    list = VolumesAttrList (volumes) = VCreateAttrList ();
  VPrependAttr (list, VolNVolumesAttr, NULL, VShortRepn,
		(VShort) volumes->nvolumes);
  VPrependAttr (list, VolNBandsAttr, NULL, VShortRepn,
		(VShort) volumes->nbands);
  VPrependAttr (list, VolNRowsAttr, NULL, VShortRepn,
		(VShort) volumes->nrows);
  VPrependAttr (list, VolNColumnsAttr, NULL, VShortRepn,
		(VShort) volumes->ncolumns);

  /* Compute the file space needed for the Volumes's binary data: */

  ntracks = 0;
  for (v = volumes->first; v != NULL; v = v->next) 
    ntracks += v->ntracks;

  /*
  *lengthp = (VRepnPrecision (VShortRepn) / 8)
     * ((ntracks * 4) + (volumes->nvolumes * 2))
     + (2 * VRepnPrecision (VShortRepn) / 8) * volumes->nvolumes;
     */

  *lengthp = (VRepnPrecision (VShortRepn) / 8) * 
    (4 * ntracks + 5 * volumes->nvolumes);

  return list;
}


/*
 *  VolumesEncodeDataMethod
 *
 *  The "encode_data" method registered for the "Volumes" type.
 */

static VPointer 
VolumesEncodeDataMethod (VPointer value, VAttrList list,
			 size_t length, VBoolean *free_itp)
{
  Volumes volumes = value;
  Volume v;
  size_t len;
  VPointer ptr, p;
  VShort idata[4],hashlen,nbands,ntracks,label;
  long nt1,nt2;
  long i;
  VTrack t;

  /* Allocate a buffer for the encoded data: */
  p = ptr = VMalloc (length);

  /* Pack each edge: */
  for (v = volumes->first; v != NULL; v = v->next) {

    /* Pack the number of tracks in volume: */
    nt2 = 0;
    nt1 = v->ntracks + (long) VRepnMinValue(VShortRepn);
    if (nt1 > (long) VRepnMaxValue(VShortRepn)) {
      nt1 = (long) VRepnMaxValue(VShortRepn) - 1;
      nt2 = v->ntracks - nt1 + (long) VRepnMinValue(VShortRepn);
    }
    if (nt2 > (VLong) VRepnMaxValue(VShortRepn)) {
      VError("too many tracks in volume: ntracks:%d  %d %d",nt1,nt2,v->ntracks);
    }

    ntracks = (VShort) nt1;
    len =  VRepnPrecision (VShortRepn) / 8;
    if (! VPackData (VShortRepn, 1, & ntracks, VMsbFirst, & len, & p, NULL))
      return NULL;
    p = (char *) p + len;
    length -= len;

    ntracks = (VShort) nt2;
    len =  VRepnPrecision (VShortRepn) / 8;
    if (! VPackData (VShortRepn, 1, & ntracks, VMsbFirst, & len, & p, NULL))
      return NULL;
    p = (char *) p + len;
    length -= len;

    /* Pack the hashtable length: */
    hashlen = (VShort) v->nbuckets;
    len =  VRepnPrecision (VShortRepn) / 8;
    if (! VPackData (VShortRepn, 1, & hashlen, VMsbFirst, & len, & p, NULL))
      return NULL;
    p = (char *) p + len;
    length -= len;

    /* Pack the number of bands (used in hash function): */
    nbands = (VShort) v->nbands;
    len =  VRepnPrecision (VShortRepn) / 8;
    if (! VPackData (VShortRepn, 1, & nbands, VMsbFirst, & len, & p, NULL))
      return NULL;
    p = (char *) p + len;
    length -= len;

    /* Pack label id */
    label = (VShort) v->label;
    len =  VRepnPrecision (VShortRepn) / 8;
    if (! VPackData (VShortRepn, 1, & label, VMsbFirst, & len, & p, NULL))
      return NULL;
    p = (char *) p + len;
    length -= len;

    for (i=0; i<hashlen; i++) {
      for (t = v->bucket[i].first; t != NULL; t = t->next)  {
	idata[0] = t->band;
	idata[1] = t->row;
	idata[2] = t->col;
	idata[3] = t->length;	
	len = 4 * VRepnPrecision (VShortRepn) / 8;
	if (! VPackData (VShortRepn, 4, & idata[0], VMsbFirst, & len, & p, NULL))
	  return NULL;
	p = (char *) p + len;
	length -= len;
      }
    }
  }

  *free_itp = TRUE;

  return ptr;
}
