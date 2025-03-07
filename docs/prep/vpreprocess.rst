
Spatial smoothing and temporal filtering
-----------------------------------------

Spatial and temporal filtering can be applied to fMRI data using the program **vpreprocess**.
Spatial filtering is often used to improve the signal-to-noise-ratio (SNR), and to reduce
the inter-subject variability. Temporal filtering is used to remove low frequency drifts in an fMRI time
series. Temporal filtering is implemented in **vpreprocess** using a Fast Fourier Transform (FFT).
Specified frequencies are removed in the frequency domain. Finally, an
inverse FFT transforms the image back into the original space.

Spatial and temporal filtering can be done in one call. Example:

 ::

   vpreprocess -in input.v -out output.v -fwhm 6 -low 0 -high 100


This call will perform a spatial filter using a Gaussian kernel with a FWHM of 6 mm, and
a temporal highpass filter that allows frequencies faster than 1/100 Hz to pass.
Note that the cutoff is specified in *seconds* and not in *Hertz*. In here, the lowpass filter is set to 0, which effectively disables the lowpass filter. Spatial filtering can also be disabled by setting '-fwhm 0', the same is true for disabling the highpass filter '-high 0'.

**To summarize:** If '-high' is set to a positive number, 'vpreprocess' performs a high pass filter.
If '-low' is set to a positive number, 'vpreprocess' performs a low pass filter that
allows low frequencies to pass.
Highpass and lowpass filtering can be combined into a bandpass filter in which
frequencies *outside* the interval '-high' and '-low' are suppressed.
Example:

 ::

   vpreprocess  -in input.v -out output.v -low 3 -high 100


Instead of a band pass filter, a band stop filter can be applied using the option '-stop' (i.e. the opposite of a band pass filter).
In this case, frequencies between '-high' and '-low' are suppressed. Example:


 ::

   vpreprocess  -in input.v -out output.v -low 3 -high 100 -stop true

   
The parameter '-minval' can be used to exclude voxels (e.g. voxels outside the brain) from the analysis.
If a region of interest mask is specified, it takes precedence over the parameter '-minval'.



Parameters of 'vpreprocess'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 -help     Usage information.
 -in       Input file.
 -out      Output file.
 -fwhm     Spatial filter: FWHM in mm. Default: 0
 -lowpass  Suppress fast oscillations (in secs). Default: 0
 -highpass  Suppress slow oscillations (in secs). Default: 0
 -stop     Invert bandpass filter. Default: false
 -mask     Region of interest mask (optional).
 -minval   Signal threshold. Default: -1e+06.


.. index:: preprocess
