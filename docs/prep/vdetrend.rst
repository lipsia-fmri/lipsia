
Detrending
------------

The program **vdetrend** can be used for removing a baseline drift.
The drift can be removed either by subtracting a moving average or
by subtracting a polynomial model of the drift.
Available models are: linear, cubic or a polynomial of order 5.
The parameter '-del' specifies the number
of initial timepoints that will be replaced by a mean value.
The parameter '-linear' can be used to apply a linear detrending 
as a first step before other types of detrending are used.


Example:

 ::

   vdetrend -in func.v -out detrended_func.v -type demean -linear true -del 2



Parameters of 'vdetrend'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 -help     Usage information.
 -in       Input file.
 -out      Output file.
 -type     Type of detrending [demean | linear | cubic | order5]  Default: demean
 -linear   Whether to first subtract linear drift. Default: true
 -window   Window size in seconds (only for demeaning). Default: 100
 -del      Number of initial timepoints to ignore. Default: 0



.. index:: vdetrend
