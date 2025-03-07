======================================
LIPSIA 3.1.1: fMRI analysis tools
======================================

Lipsia is a collection of tools for the analysis of functional magnetic resonance imaging (fMRI) data.
Its primary focus lies in implementing novel algorithms, including laminar-specific fMRI analysis (:doc:`cylarim <cylarim/index>`),
statistical inference (:doc:`LISA <stats/vlisa_generic>`), and Eigenvector centrality mapping (:doc:`ECM <nets/vecm>`).
Lipsia is designed with a focus on compactness and ease of installation,
making it readily accessible for researchers to incorporate these advanced analysis methods into their workflows.

The code is available on GitHub: `GitHub <https://github.com/lipsia-fmri/lipsia>`_.
If you like lipsia, please go to the repository and star it!

.. raw:: html

   <a class="github-button" href="https://github.com/lipsia-fmri/lipsia" data-icon="octicon-star" data-show-count="true" aria-label="Star lipsia-fmri/lipsia on GitHub">Star</a>
   <script async defer src="https://buttons.github.io/buttons.js"></script>



Installation and first steps: 
------------------------------------
:doc:`Installation <install>` and :doc:`Getting started <overview>`

Layer-specific fMRI analysis:
------------------------------------
:doc:`Cylarim <cylarim/index>`

Statistical inference:
------------------------------------
:doc:`LISA <stats/vlisa_generic>`

Semi-blind machine learning:
------------------------------------
:doc:`SML <sml/vsml>`

Task-related edge density (TED):
------------------------------------
:doc:`TED <ted/vted>`

Network tools:
------------------------------------
- :doc:`ECM <nets/vecm>`
- :doc:`CCM <nets/vccm>`
- :doc:`BCM <nets/vbcm>`


__________________________________________________________________


Credits
-------
If you use Lipsia in your research, please cite the relevant publications:





`Lohmann et al (2023) “Improving the reliability of fMRI-based predictions of intelligence via semi-blind machine learning”, bioRxiv <https://doi.org/10.1101/2023.11.03.565485>`_

.. code-block:: latex
  
  @article{Lohmann2023,
    title = {Improving the reliability of fMRI-based predictions of intelligence via semi-blind machine learning},
    url = {http://dx.doi.org/10.1101/2023.11.03.565485},
    DOI = {10.1101/2023.11.03.565485},
    publisher = {Cold Spring Harbor Laboratory},
    author = {Lohmann,  Gabriele and Heczko,  Samuel and Mahler,  Lucas and Wang,  Qi and Steiglechner,  Julius and Kumar,  Vinod J. and Roost,  Michelle and Jost,  J\"{u}rgen and Scheffler,  Klaus},
    year = {2023},
    month = nov 
  }

`Lohmann et al (2018) “LISA improves statistical analysis for fMRI”, Nature Comm <https://www.nature.com/articles/s41467-018-06304-z>`_

.. code-block:: latex
  
  @article{Lohmann2018,
    title = {LISA improves statistical analysis for fMRI},
    volume = {9},
    ISSN = {2041-1723},
    url = {http://dx.doi.org/10.1038/s41467-018-06304-z},
    DOI = {10.1038/s41467-018-06304-z},
    number = {1},
    journal = {Nature Communications},
    publisher = {Springer Science and Business Media LLC},
    author = {Lohmann,  Gabriele and Stelzer,  Johannes and Lacosse,  Eric and Kumar,  Vinod J. and Mueller,  Karsten and Kuehn,  Esther and Grodd,  Wolfgang and Scheffler,  Klaus},
    year = {2018},
    month = oct 
  }

`Lohmann, G. et al (2018), Eigenvector centrality mapping for ultrahigh resolution fMRI data of the human brain. bioRxiv <http://dx.doi.org/10.1101/494732>`_

.. code-block:: latex
  
  @article{Lohmann2018,
    title = {Eigenvector centrality mapping for ultrahigh resolution fMRI data of the human brain},
    url = {http://dx.doi.org/10.1101/494732},
    DOI = {10.1101/494732},
    publisher = {Cold Spring Harbor Laboratory},
    author = {Lohmann,  Gabriele and Loktyushin,  Alexander and Stelzer,  Johannes and Scheffler,  Klaus},
    year = {2018},
    month = dec 
  }

`Lohmann et al (2012), “Connectivity Concordance Mapping: A New Tool for Model-Free Analysis of fMRI Data of the Human Brain.” Frontiers in Systems Neuroscience. 2012;6:13. <http://dx.doi.org/10.3389/fnsys.2012.00013>`_

.. code-block:: latex
  
  @article{Lohmann2012,
    title = {Connectivity Concordance Mapping: A New Tool for Model-Free Analysis of fMRI Data of the Human Brain},
    volume = {6},
    ISSN = {1662-5137},
    url = {http://dx.doi.org/10.3389/fnsys.2012.00013},
    DOI = {10.3389/fnsys.2012.00013},
    journal = {Frontiers in Systems Neuroscience},
    publisher = {Frontiers Media SA},
    author = {Lohmann,  Gabriele and Ovadia-Caro,  Smadar and Jungeh\"{u}lsing,  Gerhard Jan and Margulies,  Daniel S. and Villringer,  Arno and Turner,  Robert},
    year = {2012}
  }

`Lohmann, G. et al (2010), Eigenvector centrality mapping for analyzing connectivity patterns in fMRI data of the human brain. PLoS ONE 5(4) <http://dx.doi.org/10.1371/journal.pone.0010232>`_

.. code-block:: latex

  @article{Lohmann2010,
    title = {Eigenvector Centrality Mapping for Analyzing Connectivity Patterns in fMRI Data of the Human Brain},
    volume = {5},
    ISSN = {1932-6203},
    url = {http://dx.doi.org/10.1371/journal.pone.0010232},
    DOI = {10.1371/journal.pone.0010232},
    number = {4},
    journal = {PLoS ONE},
    publisher = {Public Library of Science (PLoS)},
    author = {Lohmann,  Gabriele and Margulies,  Daniel S. and Horstmann,  Annette and Pleger,  Burkhard and Lepsien,  Joeran and Goldhahn,  Dirk and Schloegl,  Haiko and Stumvoll,  Michael and Villringer,  Arno and Turner,  Robert},
    editor = {Sporns,  Olaf},
    year = {2010},
    month = apr,
    pages = {e10232}
  }


.. toctree::
  :hidden:
  :caption: Getting started
  
  overview
  install

.. toctree::
  :hidden:
  :caption: Laminar-specific fMRI (Cylarim)

  cylarim/index
  cylarim/vcylarim
  cylarim/vcylarim_stats
  cylarim/vcylarim_plot
  cylarim/vcylarim_getmask
  cylarim/vmetric
  cylarim/vrim

.. toctree::
  :hidden:
  :caption: LISA statistical inference

  stats/vlisa_generic
  stats/vlisa_2ndlevel
  stats/vlisa_onesample
  stats/vlisa_twosample
  stats/vlisa_precoloring
  stats/vlisa_prewhitening

.. toctree::
  :hidden:
  :caption: Semi-blind machine learning (SML):

  sml/vsml
  sml/vreadconnectome
  sml/vsml_statistics

.. toctree::
  :hidden:
  :caption: Task-related edge density (TED)
    
  ted/vcuttrials
  ted/vted
  ted/vhubness
  ted/vtedfdr

.. toctree::
  :hidden:
  :caption: Network tools

  nets/vbcm
  nets/vccm
  nets/vecm

.. toctree::
  :hidden:
  :caption: Preprocessing

  prep/vdetrend
  prep/vpreprocess

.. toctree::
  :hidden:
  :caption: File format converters

  conv/vdicom
  conv/vnifti

.. toctree::
  :hidden:
  :caption: Utilities

  utils/vapplymask
  utils/vdenoise
  utils/vsrad
  stats/designformat
  stats/designformat_first
  stats/designformat_higher

