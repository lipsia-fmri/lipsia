Semi-blind machine learning - ensemble learning (SML-EL)
-----------------------------------------------------------

The program **vsml_statistics** implements a new approach for predictive modelling based on fMRI.
The main idea is to supplement fMRI data with readily available non-imaging information so 
that reliable predictive modeling becomes feasible even for smaller sample sizes.

The difference between *vsml* and *vsm_statistics* is that *vsm_statistics* investigates a number of
of randomly sampled training and test set pairs so that a statistic about the accuracy of predictions
can be made. The program *vsml* on the other hand applies SML to only one specific training and
test set pair.

The input into **vsml_statistics** is a collection of connectomes of N subjects together with a
textfile containing the target variable of interest, e.g. intelligence (one line per subject).
A second textfile containing non-imaging supplementary information is also needed (also one line per subject).
The number of rows in those two files must equal the number of subjects specified by the parameter '-nsubjects'.

The program randomly selects K non-overlapping training/test sets, where K is specified by the parameter '-numsamples'.
The number of elements in the training sets is specified by the parameter '-ntrain'.
Likewise, the number of elements in the training sets is specified by the parameter '-ntest'.
Those two numbers should not exceed the total numbers of subjects N.
The parameters '-dimX' and '-npls' are used to control the partial least squares regression (PLS),
where '-dimX' determines the number of features (edges of the connectome) that are inout into the PLS,
while '-npls' determines the number of latent components.
The parameter -nensembles' determines the number of ensembles in the ensemble learning process.


The output is a textfile with K rows (i.e. one row for each selection of training/test sets).
Each such row shows the resulting linear correlation and R^2 between the predicted and the
observed target variable (e.q. IQ). More precisely, each output row has 8 columns headed
'A','B','C','X','AX','BX',CX','alpha'.
The column headed 'A' contains the correlation between predicted and observed IQ without using
supplementary information.
The column headed 'B' contains the correlation between predicted and observed IQ with using
supplementary information, but without bias control.
The column headed 'C' contains the correlation between predicted and observed IQ with using
supplementary information, this time with bias control.
The column headed 'X' contains the correlation between predicted and observed IQ only using
supplementary information (without using fMRI data).

The column headed 'AX' contains information about the bias in column 'A',
i.e. the correlation of the prediction with the supplementary information. This correlation
should be approximately the same as the correlation of the observed IQ with the supplementary information.
Likewise, the column headed 'BX' contains information about the bias in column 'B'.
And the column headed 'CX' contains information about the bias in column 'C'.
The column headed 'alpha' shows the parameter alpha after adjustment for bias control.



Example:
^^^^^^^^^^^^

 :: 
 
   vsml_statistics -in func_*.v -ntrain 290 -ntest 100 -regressor IQ.txt -xx education.txt \
    -numsamples 50 -dimX 800 -npls 10 -nensembles 1000 -seed 12345 -out results.txt





Parameters of 'vsml':
^^^^^^^^^^^^^^^^^^^^^^^

 -help       Prints usage information.
 -in         Input files
 -out        Output textfile
 -ntrain     Number of subjects in the training set
 -ntest      Number of subjects in the test set
 -regressor  Textfile containing the target variable (e.g. IQ, one line per subject)
 -xxx        Supplementary information (e.g. education levels, one line per subject)
 -dimX       Number of features per ensemble
 -npls       Number of components for PLS
 -nensembles   Number of ensembles
 -numsamples Number of samples
 -seed       Seed for random number generator


.. index:: vsml_statistics



Reference:
^^^^^^^^^^^^^^^^^^^^^^

 Lohmann, G. et al (2023), bioRxiv, Improving the reliability of fMRI-based predictions of intelligence via semi-blind machine learning, https://doi.org/10.1101/2023.11.03.565485

