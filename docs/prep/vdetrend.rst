
Detrending
============

The program **vdetrend** can be used for removing a baseline drift.
The baseline drift is modelled either as a linear drift, or as a cubic 
polynomial or as a polynomial of order 5. The parameter '-del' specifies the number
of initial timepoints that will be replaced by a mean value.


Example:

 ::

   vdetrend -in func.v -out detrended_func.v -type cubic -del 0



Parameters of 'vdetrend'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 -help     Usage information.
 -in       Input file.
 -out      Output file.
 -type     Type of detrending [linear | cubic | order5]  Default: linear
 -del      Number of initial timepoints to ignore. Default: 1
 -minval   Signal threshold. The filters are only applied above this value. Default: 0


.. index:: vdetrend
