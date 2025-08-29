Overview
======================================

Lipsia is a collection of tools for the analysis of functional magnetic resonance imaging (fMRI) data.
Its primary focus lies in implementing novel algorithms, including laminar-specific fMRI analysis (cylarim),
statistical inference (LISA), and Eigenvector centrality mapping (ECM).
Lipsia is designed with a focus on compactness and ease of installation,
making it readily accessible for researchers to incorporate these advanced analysis methods into their workflows.




Installation
```````````````````````
Lipsia supports Linux and other operating systems via Singularity and Docker, see the files
"singularity-recipe.txt" and/or "Dockerfile". Follow the instructions here: :doc:`install`.



Lipsia file format
```````````````````````````````````````
Lipsia uses its own data format, which is called vista (extension *.v).
Many lipsia programs also accept gzipped files or nifti-files as input (*.v.gz or *.nii.gz).
The output is always in unzipped vista-format.
You can easily convert your NIfTI data from and to lipsia with the program *vnifti*::

  vnifti -in data.nii -out data.v
  vnifti -in data.nii.gz -out data.v
  vnifti -in result.v -out result.nii




Getting started
```````````````````````

Lipsia generally uses the "-in" parameter to specify input images
and "-out" for output images. Many programs within Lipsia offer additional parameters,
such as brain masks or thresholds, to fine-tune the analysis.

As an example, consider Eigenvector Centrality Mapping (ECM).
To compute an ECM map, you would typically use the following command-line instruction::

  vecm -in functional.nii.gz -mask brainmask.nii.gz -out ecm.v

This command processes the input fMRI data (functional.nii.gz) within the specified brain mask
(brainmask.nii.gz) and generates an ECM map in the vista format (ecm.v).

To further process or visualize the ECM map with tools that may not support the vista format,
you can convert it to NIfTI format using the following command::

  vnifti -in ecm.v -out ecm.nii



*Reference:*
 
Lohmann, G. et al. (2010), "Eigenvector centrality mapping for analyzing connectivity patterns in fMRI data of the human brain"
 PLoS ONE, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0010232>





Laminar-specific fMRI analysis (Cylarim) 
````````````````````````````````````````````````````````````````

Cylarim is an innovative software for laminar fMRI analysis that
addresses a major limitation of many existing approaches: their reliance on pre-defined regions of interest (ROIs).
The definition of ROIs often involves arbitrary choices,
making them difficult to justify and prone to bias.
To overcome this problem, Cylarim utilizes small, overlapping cylinders that traverse the cortical ribbon,
enabling a comprehensive mapping of laminar-specific activity across extensive cortical regions.

The program **vcylarim** is the core tool of this package. It generates the cylinders and produces
various statistics based on the cylinder structure. Several other programs are used in preparatory steps,
or are used to further refine the output produced by "vcylarim".

For a more detailed description see (:doc:`cylarim/index`).



*Reference:*


Lohmann, G. et al (2025), bioRxiv, Beyond ROIs: Laminar fMRI Mapping with Cylarim, https://doi.org/10.1101/2025.03.25.645237






Statistical inference (LISA) in examples:
```````````````````````````````````````````````````

*Onesample test at the 2nd level* (:doc:`stats/vlisa_onesample`). 
Example: the input is a set of contrast maps called "data_*.nii.gz"::

  vlisa_onesample -in data_*.nii.gz -mask mask.nii -out result.v
  vnifti -in result.v -out result.nii


*Twosample test at the 2nd level* (:doc:`stats/vlisa_twosample`). 
Example: input are two sets of contrast maps called "data1_*.nii.gz" and "data2_*.nii.gz"::

  vlisa_twosample -in1 data1_*.nii.gz -in2 data2_*.nii.gz -mask mask.nii -out result.v
  vnifti -in result.v -out result.nii


*Single subject test (1st level)* (:doc:`stats/vlisa_prewhitening`). 
Example: input are two runs acquired in the same session called "run1.nii.gz" and "run2.nii.gz".


  vlisa_prewhitening -in run1.nii.gz run2.nii.gz -design des1.txt des2.txt -mask mask.nii -out result.v 
  vnifti -in result.v -out result.nii


Note: preprocessing should include a correction for baseline drifts::



*Reference:*
  
Lohmann et al (2018) "LISA improves statistical analysis for fMRI",
Nature Comm, https://www.nature.com/articles/s41467-018-06304-z






Semi-blind machine learning (SML):
```````````````````````````````````````````````````

SML is implemented in the program **vsml**. In the following, its usage is illustrated in an example.  
The program **vsml** expects as input connectome data for all subjects of the training and the test set.
It is assumed that these connectomes have been precomputed using some other software tool and exist as text-files in csv-format.
The first step of the processing chain is to convert these connectomes into the lipsia-format. 
This is done using the program **vreadconnectome** as shown in the example below. 
In our example, the training set consists of 400 subjects, the test set has 100 subjects. 

The information about the target variable of interest (e.g.IQ) must be supplied as a text-file for all subjects of the
training set ("IQ_train.txt"). If this information is also available for the test set, 
it can optionally be supplied ("IQ_test.txt") and will be used to report the accuracy of the prediction.

Likewise, information about supplementary non-imaging information 
(e.g. educational levels) must be supplied as text-files ("Edu_train.txt", "Edu_test.txt"). 
Each row in these files contains a numerical value corresponding to a subject's attribute (e.g. IQ or educational level). 
The order of rows in these text files should align with the sequence in which the connectomes are input into **vsml**.

A few more parameters can optionally be supplied to **vsml** to adjust the partial least squares regression and ensemble learning process,
but default settings of these parameters should usually work well enough. 

The output of *vsml* is a text file ("results.txt") showing the predictions of the target variable for the subjects of the test set.



*Example usage:*
  
  for i in {1...400}; do
    **vreadconnectome -in traindata_${i}.csv -out traindata_${i}.v -ncomponents 100; done**
  
  for i in {1...100}; do
    **vreadconnectome -in testdata_${i}.csv -out testdata_${i}.v -ncomponents 100; done**
  
  **vsml -train train_\*.v -test test_\*.v -ytrain IQ_train.txt -ytest IQ_test.txt -xtrain Edu_train.txt -xtest Edu_test.txt -out results.txt**




*Reference:*
  
Lohmann et al (2023) "Improving the reliability of fMRI-based predictions of intelligence via semi-blind machine learning", bioRxiv, https://doi.org/10.1101/2023.11.03.565485



* :ref:`genindex`
* :ref:`search`
