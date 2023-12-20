Reading data for use in Semi-blind machine learning (SML-EL)
=============================================================

The program **vreadconnectome** prepares inputs for the program **vsml**.
The purpose of **vreadconnectome** is to convert the connectome data into the format expected by the
program  **vsml**.
It is assumed that the connectome data has been prepared using some other software package.

The input into **vreadconnectome** is a text-file in csv format containing the 
comma-separated connectome of one subject. Since the connectome is symmetric,
the connectome of one subject has K*(K-1)/2 entries, where K is the number of parcels/components of the atlas.
E.g. if the K=200, then there should be 19900 entries in the connectome.
Typically, such an entry is a number specifying the (tangent) correlation between two fMRI time series.

No comments ("#" or "%") are allowed in the input file. Also, there should not be any blank spaces.

The output filename must have the ending ".v". This is the convention for the vista-format.



Example:
``````````

 :: 


  vreadconnectome -in connectome_subj007.csv -ncomponents 200 -out connectome_subj007.v




Parameters of 'vreadconnectome':
````````````````````````````````

 -help          Prints usage information.
 -in            Input file (csv format)
 -out           Output file (vista-format)
 -ncomponents   Number of components/parcels in the atlas/decomposition



.. index:: vreadconnectome

