


Higher-level  level design files
================================

Higher-level designs are used to compute higher-level statistics such as two-sample or paired tests
so that groups of data sets can be statistically assessed.


For a higher level analysis,
a design file has to be created that can be used in 'vlisa_2ndlevel'.
For all higher level design files, the *number of rows* must coincide
with the number of input files to 'vlisa_2ndlevel'.

Note that convenience programs exist that can be used instead of the more complicated way of
specifying high-level statistics via 'vlisa_2ndlevel'. These programs are:

 - `vlisa_onesample`_
 - `vlisa_twosample`_


However, higher-level design files are needed in more complex situations.
Below, some simple design files are listed for tutorial purposes.




One-sample test
--------------------

 ::

        %
        % second level design file, one-sample test
        %
        1
        1
        1
        1
        1
        1
        1
        1

This example design can be used to perform a one-sample test for 8 scans.
Each row in this file corresponds to one scan.
The contrast vector in this case is '1' (to test for > 0) or '-1' (to test for < 0).
For such simple scenarios, an easier way to do a one-sample test is to use `vlisa_onesample`_.



Two-sample test
--------------------
 ::

        %
        % second level design file, two-sample test
        %
        1 0
        1 0
        1 0
        1 0
        1 0
        1 0
        1 0
        1 0
        0 1
        0 1
        0 1
        0 1
        0 1
        0 1
        0 1
        0 1



This example design can be used to perform a two-sample test
for comparing two groups of 8 subjects each (e.g. patients vs. healthy controls). The contrast vector comparing the two groups should be '1 -1'
(to test for patient > control) or '-1 1' (to test for patient < control).
**Important**: the order of the scans in this example is expected to be

::

       patient 1
       patient 2
       patient 3
       patient 4
       patient 5
       patient 6
       patient 7
       patient 8
       control 1
       control 2
       control 3
       control 4
       control 5
       control 6
       control 7
       control 8

For such simple scenarios, an easier way to the same test is to use  `vlisa_twosample`_.



Paired test
----------------

::

       %
       % second level design file, paired test
       %
       1 0 1 0 0 0 0 0 0 0
       0 1 1 0 0 0 0 0 0 0
       1 0 0 1 0 0 0 0 0 0
       0 1 0 1 0 0 0 0 0 0
       1 0 0 0 1 0 0 0 0 0
       0 1 0 0 1 0 0 0 0 0
       1 0 0 0 0 1 0 0 0 0
       0 1 0 0 0 1 0 0 0 0
       1 0 0 0 0 0 1 0 0 0
       0 1 0 0 0 0 1 0 0 0
       1 0 0 0 0 0 0 1 0 0
       0 1 0 0 0 0 0 1 0 0
       1 0 0 0 0 0 0 0 1 0
       0 1 0 0 0 0 0 0 1 0
       1 0 0 0 0 0 0 0 0 1
       0 1 0 0 0 0 0 0 0 1


This example design file can be used to perform a paired test. Here it specifies 8 subjects which have been scanned twice (e.g. the same subject scanned in the morning vs. scanned in the afternoon).
The design file much resembles the one from above for the twosample test. However, there are eight additional columns (one for each subject).
A value of "1" in these columns indicates that this scan belongs to the subject's column. In our example, the third column indicates the first subject (the two columns before indicate the experimental condition, as in the twosample test). Only the first and second scan (row in the design matrix) have a "1" here. This way of writing the design matrix is equivalent to a paired twosample test, as the additional subject-wise regressors take out what is commonly found in both scans of the same subject. Having more columns in the design matrix means that the contrast vector is longer (it should always correspond to the number of columns).
Furthermore, as we don't want to constrain the permutations within one subject (avoiding mixing of within- and across subjects), we need to supply the permutation group IDs with the **-grp** argument. We supply a text file with:


::

       1
       1
       2
       2
       3
       3
       4
       4
       5
       5
       6
       6
       7
       7
       8
       8




To estimate the effect, we could use as contrast vector '1 -1 0 0 0 0 0 0 0 0' (indicating morning > afternoon).
An easier way to do a paired test is to use  **vlisa_twosample -test paired**.

**Important**: the order of the scans in this example is expected to be


::

       morning subj 1
       afternoon subj 1
       morning subj 2
       afternoon subj 2
       morning subj 3
       afternoon subj 3
       morning subj 4
       afternoon subj 4
       morning subj 5
       afternoon subj 5
       morning subj 6
       afternoon subj 6
       morning subj 7
       afternoon subj 7
       morning subj 8
       afternoon subj 8


ANOVA (interactions and main effects, repeated measures)
-------------------------------------------------------------

::

       %
       % second level design file, ANOVA
       %
       1 0 0 0 1 0 0 0 0 0 0 0
       0 1 0 0 1 0 0 0 0 0 0 0
       1 0 0 0 0 1 0 0 0 0 0 0
       0 1 0 0 0 1 0 0 0 0 0 0
       1 0 0 0 0 0 1 0 0 0 0 0
       0 1 0 0 0 0 1 0 0 0 0 0
       1 0 0 0 0 0 0 1 0 0 0 0
       0 1 0 0 0 0 0 1 0 0 0 0
       0 0 1 0 0 0 0 0 1 0 0 0
       0 0 0 1 0 0 0 0 1 0 0 0
       0 0 1 0 0 0 0 0 0 1 0 0
       0 0 0 1 0 0 0 0 0 1 0 0
       0 0 1 0 0 0 0 0 0 0 1 0
       0 0 0 1 0 0 0 0 0 0 1 0
       0 0 1 0 0 0 0 0 0 0 0 1
       0 0 0 1 0 0 0 0 0 0 0 1



This example design can be used to compute an ANOVA (2x2, repeated measures). Here, we have two groups of subjects, with 4 subjects per group (hence, we have 2*2*4=16 scans = a design file with 16 rows). For instance, we have 4 control subjects scanned in the morning and afternoon, furthermore we have 4 patients scanned in the morning and afternoon (each subject was scanned twice, this is a repeated measures design).
Similar to the design file above for the paired test, we have two types of columns in the design matrix: the first four columns indicate the experimental condition for each of the two groups separately (morning control, afternoon control, morning patient, afternoon patient), the next 8 columns indicate the subject (there were 2*4 subjects). As with the paired test, this design effectively models  what is commonly found in both scans for any given subject.
Specifying the type of test (interaction or main effect) is done via the contrast vector (which has the same length as the number of columns in the design).
To study the interaction (morning control - afternoon control) - (morning patients - afternoon patients), i.e. which brain areas show a larger change morning minus afternoon or the controls, use this contrast vector: '1 -1 -1 1 0 0 0 0 0 0 0 0'. Note the signs of the third and fourth element (they differ from the first and second), as we are computing a difference of differences. To test for the main effect morning - afternoon, we specify the contrast vector: '1 -1 1 -1 0 0 0 0 0 0 0 0'.
The design files can either be generated using the tool `GLM <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/GLM>`_ of the FSL software package, or via manual specification.

Furthermore, as we don't want to constrain the permutations within one subject (avoiding mixing of within- and across subjects), we need to supply the permutation group IDs with the **-grp** argument. We supply a text file with:


::

       1
       1
       2
       2
       3
       3
       4
       4
       5
       5
       6
       6
       7
       7
       8
       8



**Important**: the order of the scans in this example is expected to be


::


       morning control 1
       afternoon control 1
       morning control 2
       afternoon control 2
       morning control 3
       afternoon control 3
       morning control 4
       afternoon control 4
       morning patient 1
       afternoon patient 1
       morning patient 2
       afternoon patient 2
       morning patient 3
       afternoon patient 3
       morning patient 4
       afternoon patient 4



.. _designformat_first: designformat_first.rst

.. _designformat_higher: designformat_higher.rst


.. _vlisa_2ndlevel: vlisa_2ndlevel.rst

.. _vlisa_onesample: vlisa_onesample.rst

.. _vlisa_twosample: vlisa_twosample.rst

.. _vlisa_precoloring: vlisa_precoloring.rst

.. _vlisa_prewhitening: vlisa_prewhitening.rst


.. index:: design files
