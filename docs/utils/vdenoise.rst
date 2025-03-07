Denoising
---------

The program **vdenoise** removes noise from MRI data using wavelet
filtering. Wavelet filters may produce pseduo-Gibbs artefacts.
Such artefacts should be reduced using a second program called **vsrad**.
Denoising has been tested on structural images (T1-weighted,Flash) and diffusion weighted images.
In principle, it should also be usable on fMRI data, but it has not yet been
tested in this context.

The filter strength is controlled by the parameter '-level'. Reducing the
value of '-level' reduces the filter effect and also the possibility of artefacts.
The  type of wavelet may be specified by the parameter '-wavelet'. The 'haar'
wavelet is a good choice.
Here is a typical calling sequence.

 ::

   vdenoise -in input.v -out tmp.v -level 0.5
   vsrad -in tmp.v -out output.v -iter 5



**Reference:**
*Lohmann et al (2010). Magn Res Med. 64(1), pp:15-22.*



Parameters of 'vdenoise'

 -help     Usage information.
 -in       Input file.
 -out      Output file.
 -level    Filter strength. Default: 0.5.
 -wavelet  Wavelet function to be used. Default: haar.
           Options: [haar | daub4 | daub6 | daub8 | bspline103 | bspline202].
					 
					 
.. index:: denoise
