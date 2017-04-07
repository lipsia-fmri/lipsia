/*
** Include file to be used in conjuntion with the volume set (Volumes)
** data type.
**
** Author:
**  G.Lohmann, <lohmann@cns.mpg.de>, Jan. 1996
*/

#include <stdio.h>

#define MAXHASHLEN 50000  /* max length of hash table */

typedef struct VTrackStruct {
  long band;
  long row; 
  long col;
  long length;
  struct VTrackStruct *next;
  struct VTrackStruct *previous;
} *VTrack, VTrackRec;


typedef struct VBucketStruct {
  long  ntracks;      /* number of tracks in one hashtable bucket */
  VTrack first;        /* ptr to first track in bucket             */
  VTrack last;         /* ptr to last track in bucket              */
} *VBucket, VBucketRec;

typedef struct VolumeStruct {
  long label;
  long nbands;
  long nrows;
  long ncolumns;
  long nbuckets;      /* length of hash table (number of buckets) */
  long   ntracks;       /* total number of tracks in all buckets   */
  VBucket bucket;      /* ptrs to buckets      */
  struct VolumeStruct *next;
} VolumeRec, *Volume;


typedef struct V_VolumesRec {
  VAttrList attributes;
  long nvolumes;        /* number of volumes in list       */
  long nbands;
  long nrows;
  long ncolumns;
  Volume first;          /* ptr to first volume in list     */
} VolumesRec;

#define VolumesAttrList(volumes) ((volumes)->attributes)
#define VolumesNum(volumes) ((volumes)->nvolumes)
#define VolumesNBands(volumes) ((volumes)->nbands)
#define VolumesNRows(volumes) ((volumes)->nrows)
#define VolumesNColumns(volumes) ((volumes)->ncolumns)
#define VolumesNVolumes(volumes) ((volumes)->ntracks)
/*
#define VolumesNTracks(volumes) ((volumes)->ntracks)
*/

#define VolumeNBands(volume) ((volume)->nbands)
#define VolumeNRows(volume) ((volume)->nrows)
#define VolumeNColumns(volume) ((volume)->ncolumns)
#define VolumeNBuckets(volume) ((volume)->nbuckets)
#define VolumeNTracks(volume) ((volume)->ntracks)
#define VolumeLabel(volume) ((volume)->label)
#define VFirstVolume(volumes) ((volumes)->first)
#define VNextVolume(volume) ((volume)->next)
#define VolumeExists(volume) ((volume) != NULL)

#define VTrackLength(track) ((track)->length)
#define VTrackExists(track) ((track) != NULL)
#define VFirstTrack(volume,i) ((volume)->bucket[(i)].first)
#define VNextTrack(track) ((track)->next)
#define VPreviousTrack(track) ((track)->previous)

#define VolumesAttr     "volumes"
#define VolNVolumesAttr "nvolumes"
#define VolNTracksAttr  "ntracks"
#define VolNBandsAttr   "nbands"
#define VolNRowsAttr    "nrows"
#define VolNColumnsAttr "ncolumns"

extern Volumes VCreateVolumes(long,long,long);
extern Volumes VCopyVolumes(Volumes);
extern void VDestroyVolumes(Volumes);
extern VBoolean VWriteVolumes(FILE *, VAttrList, long, Volumes *);

extern long VReadVolumes(FILE *, VAttrList *, Volumes **);

extern Volume VCreateVolume(long,long,long,long,long);
extern Volume VCopyVolume(Volume);
extern void VAddVolume(Volumes, Volume);
extern void AddTrack(Volume,VTrack);

extern double VolumeBorderSize(Volume);
extern VBoolean VolumeBorder(Volume,long,long,long);
extern VTrack VolumeGetTrack(Volume,long,long,long);
extern VBoolean VolumeInside(Volume,long,long,long);
extern double VolumeRadius(Volume,double *);

/*
** hash function
*/
#define VolumeHash(nbands, b, r, len) (((b) * (nbands) + (r)) % (len))
