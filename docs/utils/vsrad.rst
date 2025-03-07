Speckle filter
--------------


The program **vsrad** implements a speckle reducing filter.
It is used in conjuncton with **vdenoise**.
It used to reduce artefacts which may result from a prior call to 'vdenoise'.
'vdenoise' is a wavelet filter and may produce Gibbs ringing artefacts.
Denoising using 'vdenoise' and 'vsrad' has been tested on structural images 
(T1-weighted,Flash) and diffusion weighted images (DWI).
In principle, it should also be usable on fMRI data, but it has not yet been
tested in this context.


The program 'vsrad' should generally be called after 'vdenoise',
but it may also be used standalone, i.e. without first applying 'vdenoise'.
The number of iterations in 'vsrad' controls the filter strength.
Fewer iterations produce a lesser effect.
Here is a typical calling sequence.

 
 ::.

  vdenoise -in input.v -out tmp.v -level 0.5
  vsrad -in tmp.v -out output.v -iter 5



**Reference:**  *Lohmann et al (2010). Magn Res Med. 64(1), pp:15-22.*



Parameters of 'vsrad'

 -help     Usage information.

 -in       Input file.

 -out      Output file.

 -iter     Max number of iterations. Default: 7.

 -type     Type of diffusion filter (0 or 1). Default 0.

 -rho      Speckle parameter. Default: 0.1666.


.. index:: srad, speckle reduction

