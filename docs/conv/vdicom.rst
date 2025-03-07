Convert from dicom files
-------------------------

This program converts dicom files to vista files.


The program is typically called as follows:

 ::

   vdicom -in directory -p file


If the directory contains a single image,
this call will create one new file named "file.v"
If multiple images exists within the same directory, it will create
multiple output files named "file-1.v", "file-2.v", etc.


Parameters of 'vdicom'
~~~~~~~~~~~~~~~~~~~~~~~~~~~

 -help         Usage information.
 -in           Input directory containing raw DICOM files.
 -p            Prefix of output image name.
 -t            Target dirctory. Default: curent directory.
 -f            Force overwrite of existing files

.. index:: dicom
