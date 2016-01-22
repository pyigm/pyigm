""" Class for Metallicity PDF
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import warnings

import numpy as np
from scipy.interpolate import interp1d


import pdb


class MetallicityPDF(object):
    """A Class for Metallicity PDF

    Attributes
    ----------
    ZH : ndarray
      Values where the PDF has been evaluated
    phi_ZH : ndarray
      PDF of the metallicity
    """

    def __init__(self, ZH, pdf_ZH, parent=None, dZH=None, normalize=True):
        """
        Parameters
        ----------
        ZH : ndarray
          Values of ZH (log10, typically but not required) where the PDF has been evaluated
        pdf_ZH : ndarray
          PDF of the metallicity
        parent : object, optional
          Parent of the metalliicty PDF, e.g. an AbsSystem
        dZH : float or ndarray, optional
          Width of the ZH 'bins' (log 10)
        normalize : bool, optional
          Normalize the input PDF?

        Returns
        -------
        MetallicityPDF

        """
        # Tests
        if ZH.size != pdf_ZH.size:
            raise IOError("ZH and pdf_ZH must have the same size")
        # Set
        self.ZH = ZH
        self.pdf_ZH = pdf_ZH
        self.parent = parent

        # dZH
        if dZH is None:
            dZH = self.ZH - np.roll(self.ZH,1)
            dZH[0] = dZH[1]
            self.dZH = dZH
        else:
            self.dZH = dZH

        # Normalize pdf_ZH
        if normalize:
            self.normalize()
        else:
            warnings.warn("Not normalizing the PDF")

    @property
    def meanZH(self):
        """ Calculate and return average ZH from the PDF

        Weighted in log-space

        Returns
        -------
        meanZH : float
        """
        meanZH = np.sum(self.ZH*self.pdf_ZH*self.dZH)
        return meanZH

    @property
    def medianZH(self):
        """ Calculate and return median ZH from the PDF

        Weighted in log-space

        Returns
        -------
        medianZH : float
        """
        cumsum = np.cumsum(self.pdf_ZH*self.dZH)
        # Interpolate
        fint = interp1d(cumsum, self.ZH)
        medianZH = float(fint(0.5))
        return medianZH

    def confidence_limits(self, cl):
        """ Calculate the bounds of a given confidence interval

        Spline interpolation is used

        Parameters
        ----------
        cl : float
          Confidence interval, ranging from 0 to 1

        Returns
        -------
        ZH_min, ZH_max : float, float
          Bounds corresponding to the input confidence limit
        """
        from scipy.interpolate import interp1d

        if (cl <= 0.) or (cl >=1):
            raise IOError("cl must range from 0-1")
        # Spline the PDF cumulative sum vs ZH
        cumul = np.cumsum(self.dZH*self.pdf_ZH)
        f = interp1d(cumul, self.ZH)
        # Caculate
        frac = (1.- cl) / 2.
        ZH_min = float(f(frac))
        ZH_max = float(f(1-frac))
        # Return
        return ZH_min, ZH_max

    def normalize(self):
        """ Normalize the PDF

        """
        #
        norm = np.sum(self.dZH*self.pdf_ZH)
        self.pdf_ZH /= norm

    def hist_plot(self, bokeh=False):
        """ Simple histogram plot of the PDF

        Parameters
        ----------
        bokeh : bool, optional
          Generate a bokeh plot?

        Returns
        -------

        """
        if not bokeh:
            from matplotlib import pyplot as plt
            # imports
            try:
                import seaborn as sns; sns.set_style("white")
            except:
                pass
            # Giddy up
            plt.clf()
            plt.bar(self.ZH-self.dZH/2., self.pdf_ZH, width=self.dZH)
            plt.xlabel("[Z/H]")
            plt.ylabel("PDF([Z/H]")
            plt.show()
            plt.close()
        else:
            from bokeh.io import show
            from bokeh.plotting import figure
            p = figure(plot_width=400, plot_height=400, title='[Z/H] PDF')
            p.quad(top=self.pdf_ZH, bottom=0, left=self.ZH-self.dZH/2.,
                   right=self.ZH+self.dZH/2.)
            p.xaxis.axis_label = '[Z/H]'
            # Show
            show(p)

    def __add__(self, other):
        """ Combine the existing PDF with an input one

        Parameters
        ----------
        other : MetallicityPDF

        Returns
        -------
        MetallicityPDF
          The

        """
        # Check that the ZH arrays are identical
        if not np.allclose(self.ZH,other.ZH):
            raise IOError("ZH arrays need to be identical (for now)")
        # Add em'
        new = MetallicityPDF(self.ZH, self.pdf_ZH+other.pdf_ZH, normalize=False)
        # Return
        return new

    def __repr__(self):
        repr = '<{:s}: meanZH={:g}'.format(self.__class__.__name__, self.meanZH)
        try:
            repr = repr + ', ' + self.parent.name
        except:
            pass
        repr = repr + '>'
        return repr


