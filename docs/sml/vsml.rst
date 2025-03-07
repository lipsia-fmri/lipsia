Semi-blind machine learning - ensemble learning (SML-EL)
-----------------------------------------------------------

The program **vsml** implements a new approach for predictive modelling based on fMRI.
The main idea is to supplement fMRI data with readily available non-imaging information so 
that reliable predictive modeling becomes feasible even for smaller sample sizes.

The difference between *vsml* and *vsm_statistics* is that *vsm_statistics* investigates a number of
of randomly sampled training and test set pairs so that a statistic about the accuracy of predictions
can be made. The program *vsml* on the other hand applies SML to only one specific training and
test set pair.


The input into **vsml** are collections of connectomes together with a
textfiles containing the target variable of interest, e.g. intelligence,
and additional textfiles containing non-imaging supplementary information.

The program **vsml** expects as input a list of connectomes for training (parameter '-train'),
and a list of connectomes for testing (parameter '-test').
The connectomes must be in vista-format. The program **vreadconnectome** can be used to
convert those inputs into the required format.

Furthermore, the program **vsml** requires as input a text-file containing the target variable of
interest (e.g. IQ, parameter '-ytrain').
This file is used for training. It must contain one number per subject of the training set,
so that the number of rows in this file equals the number of training connectomes.

Optionally, a text-file containing the target variable for
interest for the test set can also be supplied (parameter '-ytest'). If available, this information can be used
to assess the accuracy of the prediction.

The order in which the connectomes are listed as input into the '-train' and '-test' parameters 
must coincide with the order of the rows in the respective text-files.

Likewise, **vsml** requires as input text-files containing information about the supplementary info.
There should be one file for the training set (parameter '-xtrain') and one file for the test set (parameter '-xtest').

The parameters '-dimX' and '-npls' are used to control the partial least squares regression (PLS),
where '-dimX' determines the number of features (edges of the connectome) that are input into the PLS,
while '-npls' determines the number of latent components.
The parameter -nensembles' determines the number of ensembles in the ensemble learning process.


The output is a text-file containing the predicted values of the target variable for the given test set.



Example:
```````````

 :: 


   vsml -train train_*.v -test test_*.v -ytrain IQ_train.txt -ytest IQ_test.txt \
    -xtrain Edu_train.txt -xtest Edu_test.txt -dimX 800 -npls 10 -nensembles 1000 -seed 12345  \
    -out results.txt




Parameters of 'vsml':
`````````````````````````

   -help       Prints usage information.
   -train      Input fMRI files, training set (Required).
   -test       Input fMRI files, test set (Optional).
   -out        Output textfile.
   -ytrain     Textfile containing the target variable of the training set.
   -ytest      Textfile containing the target variable of the test set.
   -xtrain     Textfile containing the supplementary info of the training set.
   -xtest      Textfile containing the supplementary info of the test set.
   -dimX       Number of features per ensemble.
   -npls       Number of components for PLS.
   -nensembles  Number of ensembles.
   -seed       Seed for random number generator.



.. index:: vsml




Reference:
^^^^^^^^^^^^^^^^^^^^^^

 Lohmann, G. et al (2023), bioRxiv, Improving the reliability of fMRI-based predictions of intelligence via semi-blind machine learning, https://doi.org/10.1101/2023.11.03.565485 
