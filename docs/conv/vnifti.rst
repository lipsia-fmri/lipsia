
Convert to/from nifti files
===========================

The program **vnifti** converts nifti1 files \*.nii to vista files \*.v and back. 
It recognizes the file type from the file name suffix.
If the input file has the suffix \*.nii it assumes that it is a nifti1 file
that should be converted to a vista file.
If the input file has a \*.v ending, it assumes that it is a vista file that
should be converted to a nifti1 file. Currently, only nifti1 is supported.
Nifti2 is not supported. 

In some nifti files, the repetition time may not be correctly given. In this
case the repetition time can be directly specified on the command line using the
parameter "-tr". If this parameter is not set, the program assumes that the repetition time
is correctly specified in the nii header.

The program is typically called as follows:

 ::

  vnifti -in data.v -out data.nii 
  vnifti -in data.nii -out data.v 
  vnifti -in data.nii.gz -out data.v




Parameters of *vnifti*


 -help         Usage information.
 -in           Input image. 
 -out          Output image.
 -tr           Repetition time in seconds. Default: 0
 -scale        Whether to scale 4D data to 16bit integer. Default: false
 -attrtype     Whether to output 4D data to lipsia 4D format. Default: true

.. index:: nifti
