/*
** read nifti-1 to vista format
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include <nifti/nifti2.h>
#include <nifti/nifti1_io.h>

#define ABS(x) ((x) > 0 ? (x) : -(x))


void VByteSwapData(char *data,size_t ndata,size_t nsize)
{
  size_t i;

  if (nsize < 2) return;

  if (nsize == 2) {
    for (i=0; i<ndata; i+=nsize) {
      nifti_swap_2bytes(1,&data[i]);
    }
  }
  else if (nsize == 4) {
    for (i=0; i<ndata; i+=nsize) {
      nifti_swap_4bytes(1,&data[i]);
    }
  }
  else if (nsize == 8) {
    for (i=0; i<ndata; i+=nsize) {
      nifti_swap_8bytes(1,&data[i]);
    }
  }
  else {
    VError(" no byte swap available for %d byte data",nsize);
  }
}

void VByteSwapNiftiHeader(nifti_1_header *hdr)
{
  int i;
  fprintf(stderr," byte swap...\n");
  nifti_swap_2bytes(1,&hdr->datatype);
  nifti_swap_2bytes(1,&hdr->dim_info);
  
  for (i=0; i<8; i++) {
    nifti_swap_2bytes(1,&hdr->dim[i]);
    nifti_swap_4bytes(1,&hdr->pixdim[i]);
  }

  nifti_swap_2bytes(1,&hdr->bitpix);
  nifti_swap_4bytes(1,&hdr->cal_max);
  nifti_swap_4bytes(1,&hdr->cal_min);

  nifti_swap_2bytes(1,&hdr->slice_start);
  nifti_swap_2bytes(1,&hdr->slice_end);
  nifti_swap_4bytes(1,&hdr->slice_duration);
  nifti_swap_4bytes(1,&hdr->toffset);

  nifti_swap_4bytes(1,&hdr->vox_offset);
  nifti_swap_2bytes(1,&hdr->intent_code);

  nifti_swap_4bytes(1,&hdr->scl_slope);
  nifti_swap_4bytes(1,&hdr->scl_inter);

  nifti_swap_4bytes(1,&hdr->cal_min);
  nifti_swap_4bytes(1,&hdr->cal_max);

  nifti_swap_2bytes(1,&hdr->qform_code);
  nifti_swap_2bytes(1,&hdr->sform_code);  
  nifti_swap_4bytes(1,&hdr->quatern_b);
  nifti_swap_4bytes(1,&hdr->quatern_c);
  nifti_swap_4bytes(1,&hdr->quatern_d);
  nifti_swap_4bytes(1,&hdr->qoffset_x);
  nifti_swap_4bytes(1,&hdr->qoffset_y);
  nifti_swap_4bytes(1,&hdr->qoffset_z);

  for (i=0; i<4; i++) {
    nifti_swap_4bytes(1,&hdr->srow_x[i]);
    nifti_swap_4bytes(1,&hdr->srow_y[i]);
    nifti_swap_4bytes(1,&hdr->srow_z[i]);
  }
}


void VByteSwapNifti2Header(nifti_2_header *hdr)
{
  int i;
  fprintf(stderr," byte swap...\n");
  nifti_swap_8bytes(1,&hdr->datatype);
  nifti_swap_8bytes(1,&hdr->dim_info);
  
  for (i=0; i<8; i++) {
    nifti_swap_8bytes(1,&hdr->dim[i]);
    nifti_swap_8bytes(1,&hdr->pixdim[i]);
  }

  nifti_swap_8bytes(1,&hdr->bitpix);
  nifti_swap_8bytes(1,&hdr->slice_start);
  nifti_swap_8bytes(1,&hdr->slice_end);

  nifti_swap_8bytes(1,&hdr->vox_offset);
  nifti_swap_8bytes(1,&hdr->scl_slope);
  nifti_swap_8bytes(1,&hdr->scl_inter);

  nifti_swap_8bytes(1,&hdr->cal_min);
  nifti_swap_8bytes(1,&hdr->cal_max);

  nifti_swap_8bytes(1,&hdr->slice_duration);
  nifti_swap_8bytes(1,&hdr->toffset);

  nifti_swap_8bytes(1,&hdr->slice_start);
  nifti_swap_8bytes(1,&hdr->slice_end);
  nifti_swap_4bytes(1,&hdr->slice_code);

  nifti_swap_4bytes(1,&hdr->xyzt_units);
  nifti_swap_4bytes(1,&hdr->intent_code);

  nifti_swap_4bytes(1,&hdr->qform_code);
  nifti_swap_4bytes(1,&hdr->sform_code);

  nifti_swap_8bytes(1,&hdr->quatern_b);
  nifti_swap_8bytes(1,&hdr->quatern_c);
  nifti_swap_8bytes(1,&hdr->quatern_d);
  nifti_swap_8bytes(1,&hdr->qoffset_x);
  nifti_swap_8bytes(1,&hdr->qoffset_y);
  nifti_swap_8bytes(1,&hdr->qoffset_z);

  for (i=0; i<4; i++) {
    nifti_swap_8bytes(1,&hdr->srow_x[i]);
    nifti_swap_8bytes(1,&hdr->srow_y[i]);
    nifti_swap_8bytes(1,&hdr->srow_z[i]);
  }
}
