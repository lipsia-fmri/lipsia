#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include <nifti/nifti1.h>

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352


/* slicetiming information */
void VGetSlicetimes(int slice_start,int slice_end,int slice_code,double slice_duration,double *slicetimes)
{
  int i;
  double tsum = 0;

  for (i=slice_start; i<=slice_end; i++) slicetimes[i] = 0;

  if (slice_end < slice_start) 
    VError(" inconsistent slicetime info, slice_start=%d, slice_end=%d",slice_start,slice_end);

  switch(slice_code) {

  case NIFTI_SLICE_SEQ_INC:
    tsum=0;
    for (i=slice_start; i<=slice_end; i++) {
      slicetimes[i] = tsum;
      tsum += slice_duration;
    }
    break;

  case NIFTI_SLICE_SEQ_DEC:
    tsum = 0;
    for (i=slice_end; i>=slice_start; i--) {
      slicetimes[i] = tsum;
      tsum += slice_duration;
    }
    break;

  case NIFTI_SLICE_ALT_INC:
    tsum=0;
    for (i=slice_start; i<=slice_end; i+=2) {
      slicetimes[i] = tsum;
      tsum += slice_duration;
    }
    for (i=slice_start+1; i<=slice_end; i+=2) {
      slicetimes[i] = tsum;
      tsum += slice_duration;
    }
    break;

  case NIFTI_SLICE_ALT_DEC:
    tsum = 0;
    for (i=slice_end; i>=slice_start; i-=2) {
      slicetimes[i] = tsum;
      tsum += slice_duration;
    }
    for (i=slice_end-1; i>=slice_start; i-=2) {
      slicetimes[i] = tsum;
      tsum += slice_duration;
    }
    break;

  case NIFTI_SLICE_ALT_INC2:
    tsum=0;
    for (i=slice_start+1; i<=slice_end; i+=2) {
      slicetimes[i] = tsum;
      tsum += slice_duration;
    }
    for (i=slice_start; i<=slice_end; i+=2) {
      slicetimes[i] = tsum;
      tsum += slice_duration;
    }
    break;

  case NIFTI_SLICE_ALT_DEC2:
    tsum = 0;
    for (i=slice_end-1; i>=slice_start; i-=2) {
      slicetimes[i] = tsum;
      tsum += slice_duration;
    }
    for (i=slice_end; i>=slice_start; i-=2) {
      slicetimes[i] = tsum;
      tsum += slice_duration;
    }
    break;

  default:
    VWarning(" unknown slicetiming code (%d)",slice_code);
  }
}
