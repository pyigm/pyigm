.. highlight:: rest

***********
Metallicity
***********

.. index:: metallicity

Notebooks
=========

.. toctree::
   :maxdepth: 1

   PDF Examples <MetallicityPDF_examples>

Overview
========

pyigm contains (or will contain) a set of classes and
methods to enable analysis of the metallicity of gas
in the IGM/CGM.


MetallicityPDF
==============

The MetallicityPDF Class holds the PDF of metallicity
values for a system (which need not be specified).
By default, the analysis is performed in log10 space.
Here is a simple instantiation::

   ZH  = np.linspace(-5, 0., 25)
   pdf = np.exp(-(ZH+1.5)**2/0.2**2)
   # Class
   mpdf = MetallicityPDF(ZH, pdf)
   # Mean
   print('Mean of the PDF is {:g}'.format(mpdf.meanZH))



Attributes/Properties
---------------------

========   ============== ============================================
Variable   Type           Description
========   ============== ============================================
meanZH     float          Weighted Mean of [Z/H]
========   ============== ============================================


Methods
-------

One can calulate the bounds of a confidence interval with::

   mpdf.confidence_limits(0.68)

Plots
-----

One can generate a quick histogram plot of the PDF with::

   mpdf.hist_plot()

