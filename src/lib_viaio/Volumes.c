/*
** Methods for the volume set (Volumes) data type.
**
** Author:
** G.Lohmann, MPI/CNS, <lohmann@cns.mpg.de>, Dez. 1995
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




/*
** create a new volume list and return a ptr to it.
*/
Volumes VCreateVolumes(long nbands, long nrows, long ncols)
{
  Volumes volumes;

  /* allocate memory: */
  volumes = VMalloc(sizeof(VolumesRec));
  if (!volumes) return NULL;

  /* initialize structure: */
  volumes->attributes = VCreateAttrList();
  volumes->first    = NULL;
  volumes->nvolumes = 0;
  volumes->nbands   = nbands;
  volumes->nrows    = nrows;
  volumes->ncolumns = ncols;

  return volumes;
}

/*
** create a single volume and initialize it, i,e. allocate space for its
** hashtable.
*/
Volume VCreateVolume(long label,long nbands,long nrows,long ncolumns,long nbuckets)
{
  Volume vol;
  long i;

  vol = VMalloc(sizeof(VolumeRec));
  vol->label    = label;
  vol->nbands   = nbands;
  vol->nrows    = nrows;
  vol->ncolumns = ncolumns;
  vol->nbuckets = nbuckets;
  vol->ntracks  = 0;
  vol->next     = NULL;

  if (nbuckets > MAXHASHLEN) {
    VWarning(" VCreateVolume, nbuckets exceeds hashlength, %ld %d",nbuckets,MAXHASHLEN);
    vol->nbuckets = MAXHASHLEN;
  }
  else vol->nbuckets  = (long) nbuckets;
  vol->bucket = (VBucket) VMalloc(sizeof(VBucketRec) * vol->nbuckets);
  if (!vol->bucket) return NULL;

  for (i=0; i<vol->nbuckets; i++) {
    vol->bucket[i].ntracks = 0;
    vol->bucket[i].first   = NULL;
    vol->bucket[i].last    = NULL;
  }
  return vol;
}


/*
** copy src volume to dest volume
*/
Volume VCopyVolume(Volume src)
{
  Volume dest;
  VTrack s,t,s_prev;
  long i,j,ntracks;

  dest = VCreateVolume(src->label,src->nbands,src->nrows,src->ncolumns,
		       src->nbuckets);
  dest->ntracks = src->ntracks;

  for (i=0; i<src->nbuckets; i++) {

    ntracks = src->bucket[i].ntracks;
    dest->bucket[i].ntracks = ntracks;
    dest->bucket[i].first = NULL;
    dest->bucket[i].last  = NULL;

    if (ntracks == 0) continue;

    s = VMalloc(sizeof(VTrackRec));
    dest->bucket[i].first = s;
    s_prev = NULL;
    t = src->bucket[i].first;

    for (j=0; j<ntracks-1; j++) {
      s->band   = t->band;
      s->row    = t->row;
      s->col    = t->col;
      s->length = t->length;

      s->previous = s_prev;
      s_prev = s;
      s->next = VMalloc(sizeof(VTrackRec));
      s = s->next;
      t = t->next;
    }
    s->band   = t->band;
    s->row    = t->row;
    s->col    = t->col;
    s->length = t->length;      
    s->previous = s_prev;
    s->next = NULL;
    dest->bucket[i].last = s;
  }
  return dest;
}



#define smaller(s,t)  \
     ((s)->band < (t)->band) \
  || (((s)->band == (t)->band) && ((s)->row < (t)->row)) \
  || (((s)->band == (t)->band) && ((s)->row == (t)->row) \
      && ((s)->col < (t)->col))

/*
** Add a new track to a hashtable while keeping the bucket sorted,
** The track must have been created prior to the call.
*/
void AddTrack(Volume v, VTrack t)
{
  long i;
  VTrack s,r;

  i = VolumeHash(v->nbands,t->band,t->row,v->nbuckets);

  if (v->bucket[i].ntracks < 1) {
    t->next = NULL;
    t->previous = NULL;
    v->bucket[i].first = t;
    v->bucket[i].last  = t;
  }
  else {

    s = v->bucket[i].last;
    while ((s != NULL) && (smaller(t,s))) {
      s = s->previous;
    }

    if ((s == v->bucket[i].last)) {  /* append to the end */
      t->previous = s;
      t->next = NULL;
      v->bucket[i].last = t;
      s->next = t;
    }

    else if ((s == NULL)) {         /* prepend to start */
      r = v->bucket[i].first;
      t->previous = NULL;
      t->next = r;
      r->previous = t;
      v->bucket[i].first = t;
    }

    else {            /* insert somewhere into the middle */
      r = s->next;
      r->previous = t;
      t->previous = s;
      t->next = r;
      s->next = t;
    }
  }
  v->bucket[i].ntracks++;
  v->ntracks++;
}

/*
** append a volume to the end of a volume list
*/
void VAddVolume(Volumes volumes, Volume vol)
{
  Volume u=NULL,v=NULL;

  /* find end of volume list: */
  v = volumes->first;
  while (v != NULL) {
    u = v;
    v = v->next;
  }

  /* append new volume to end of list: */
  if (u != NULL) u->next = vol;
  else volumes->first = vol;
  volumes->nvolumes++;
  return;
}


/*
** copy a list a volumes
*/
Volumes VCopyVolumes(Volumes src)
{
  Volumes dest;
  Volume v;

  dest = VCreateVolumes(src->nbands,src->nrows,src->ncolumns);
  if (!dest) return NULL;

  for (v=src->first; v != NULL; v = v->next) {
    VAddVolume(dest,VCopyVolume(v));
  }
  
  if (VolumesAttrList(src))
    VolumesAttrList(dest) = VCopyAttrList(VolumesAttrList(src));
  else
    VolumesAttrList(dest) = NULL;

  return dest;
}


/*
** destroy a single volume 
*/
void VDestroyVolume(Volume v)
{
  long i,n;
  VTrack t,s;

  for (i=0; i<v->nbuckets; i++) { 
    t = v->bucket[i].first;
    n = 0;
    while (t != NULL) { 
      s = t;
      t = t->next;
      VFree(s);
      n++;
    }
    if (n > v->bucket[i].ntracks) 
      VError("VDestroyVolume: free error");
  }
  
  v->nbuckets = 0;
  v->ntracks  = 0;
  v->next     = NULL;
  VFree(v->bucket);
  v->bucket = NULL;
  VFree(v);
}


/*
** destroy a set of volumes
*/
void VDestroyVolumes(Volumes volumes)
{
  Volume v,w;

  v = volumes->first;

  while (v != NULL) {
    w = v;
    v = v->next; 
    VDestroyVolume(w); 
  }
  volumes->first = NULL;
  VDestroyAttrList(volumes->attributes);
  VFree(volumes);
}


long VReadVolumes (FILE *file, VAttrList *attributes, Volumes **volumes)
{
  return VReadObjects (file, VolumesRepn, attributes,(VPointer **) volumes);
}

VBoolean VWriteVolumes (FILE *file, VAttrList attributes,long nvolumes,Volumes *volumes)
{
  return VWriteObjects (file, VolumesRepn, attributes, nvolumes,(VPointer *) volumes);
}

