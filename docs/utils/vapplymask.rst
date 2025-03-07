Masking an image
----------------------------


The program **vapplymask** applies a mask to an input image.
Here is an example:

 ::

   vapplymask -in zmap.v -mask mask.v -out result.v


The input into 'vapplymask' is a 3D image, for example a zmap resulting from a statistical analysis. 
The mask must have the same dimensions as the input image.
The output is a 3D image in which parts of the input image are zeroed out, namely those parts that are not covered by the mask.


Parameters of 'vapplymask'
``````````````````````````


    -help    Prints usage information.
    -in      Input file. Default: (none)
    -out     Output file. Default: (none)
    -mask    Mask.


.. index:: vapplymask
