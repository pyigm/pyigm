""" Class for Metallicity PDF
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import warnings

import numpy as np
from scipy.interpolate import interp1d


import pdb

class GenericPDF(object):
    """A Class for any PDF
    """

    def __init__(self, x, pdf, parent=None, dx=None, normalize=True):
        """
        Parameters
        ----------
        x : ndarray
          x-values of the PDF
          Must be sorted!
        pdf : ndarray
          PDF
        parent : object, optional
          Parent of the PDF, e.g. an AbsSystem
        dx : float or ndarray, optional
          Width of the x 'bins'
        normalize : bool, optional
          Normalize the input PDF?

        Returns
        -------
        PDF
        """
        # Tests
        if x.size != pdf.size:
            raise IOError("x and pdf must have the same size")
        # Set
        self.x = x
        self.pdf = pdf
        self.parent = parent

        # dZH
        if dx is None:
            dx = self.x - np.roll(self.x,1)
            dx[0] = dx[1]
            self.dx = dx
        else:
            self.dx = dx

        # Normalize pdf_ZH
        if normalize:
            self.normalize()
        else:
            warnings.warn("Not normalizing the PDF")

    @property
    def mean(self):
        """ Calculate and return average from the PDF

        Weighted in log-space

        Returns
        -------
        mean : float
        """
        mean = np.sum(self.x*self.pdf*self.dx)
        return mean

    @property
    def gmean(self):
        """ Calculate and return the geometric mean from the PDF

        Weighted in log-space

        Returns
        -------
        mean : float
        """
        gmean = np.sum((10**self.x)*self.pdf*self.dx)
        return gmean

    @property
    def median(self):
        """ Calculate and return median ZH from the PDF

        Weighted in log-space

        Returns
        -------
        median : float
        """
        cumsum = np.cumsum(self.pdf*self.dx)
        # Interpolate
        fint = interp1d(cumsum, self.x, assume_sorted=True)  # Sorted is important
        median = float(fint(0.5))
        return median

    def confidence_limits(self, cl):
        """ Calculate the bounds of a given confidence interval

        Spline interpolation is used

        Parameters
        ----------
        cl : float
          Confidence interval, ranging from 0 to 1

        Returns
        -------
        x_min, x_max : float, float
          Bounds corresponding to the input confidence limit
        """
        if (cl <= 0.) or (cl >= 1):
            raise IOError("cl must range from 0-1")
        # Spline the PDF cumulative sum vs ZH
        cumul = np.cumsum(self.dx*self.pdf)
        f = interp1d(cumul, self.x, assume_sorted=True)
        # Caculate
        frac = (1.- cl) / 2.
        x_min = float(f(frac))
        x_max = float(f(1-frac))
        # Return
        return x_min, x_max

    def normalize(self):
        """ Normalize the PDF
        """
        #
        norm = np.sum(self.dx*self.pdf)
        self.pdf /= norm

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
            plt.bar(self.x-self.dx/2., self.pdf, width=self.dx)
            plt.xlabel("x")
            plt.ylabel("PDF(x)")
            plt.show()
            plt.close()
        else:
            from bokeh.io import show
            from bokeh.plotting import figure
            p = figure(plot_width=400, plot_height=400, title='x PDF')
            p.quad(top=self.pdf, bottom=0, left=self.x-self.dx/2.,
                   right=self.x+self.dx/2.)
            p.xaxis.axis_label = 'x'
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
        if not np.allclose(self.x,other.x):
            raise IOError("ZH arrays need to be identical (for now)")
        # Add em'
        new = GenericPDF(self.x, self.pdf+other.pdf, normalize=False)
        # Return
        return new

    def __repr__(self):
        repr = '<{:s}: mean={:g}'.format(self.__class__.__name__, self.mean)
        try:
            repr = repr + ', ' + self.parent.name
        except:
            pass
        repr = repr + '>'
        return repr


class MetallicityPDF(GenericPDF):
    """A Class for a Metallicity PDF
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
        # Generate
        GenericPDF.__init__(self, ZH, pdf_ZH, parent=parent, dx=dZH, normalize=normalize)
        # Copy a few names
        self.ZH = self.x
        self.pdf_ZH = self.pdf
        self.dZH = self.dx

    @property
    def meanZH(self):
        """ Calculate and return average ZH from the PDF

        Weighted in log-space

        Returns
        -------
        meanZH : float
        """
        return self.mean

    @property
    def medianZH(self):
        """ Calculate and return median ZH from the PDF

        Weighted in log-space

        Returns
        -------
        medianZH : float
        """
        return self.median

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

class DensityPDF(GenericPDF):
    """A Class for a Density PDF
    """

    def __init__(self, nH, pdf_nH, parent=None, dnH=None, normalize=True):
        """
        Parameters
        ----------
        nH : ndarray
          Values of nH (log10, typically but not required) where the PDF has been evaluated
        pdf_nH : ndarray
          PDF of the nH values
        parent : object, optional
          Parent of the metalliicty PDF, e.g. an AbsSystem
        dnH : float or ndarray, optional
          Width of the nH 'bins' (log 10)
        normalize : bool, optional
          Normalize the input PDF?

        Returns
        -------
        MetallicityPDF

        """
        # Tests
        if nH.size != pdf_nH.size:
            raise IOError("nH and pdf_nH must have the same size")
        # Generate
        GenericPDF.__init__(self, nH, pdf_nH, parent=parent, dx=dnH, normalize=normalize)
        # Copy a few names
        self.nH = self.x
        self.pdf_nH = self.pdf
        self.dnH = self.dx

    @property
    def meannH(self):
        """ Calculate and return average nH from the PDF

        Weighted in log-space

        Returns
        -------
        meannH : float
        """
        return self.mean

    @property
    def mediannH(self):
        """ Calculate and return median nH from the PDF

        Weighted in log-space

        Returns
        -------
        mediannH : float
        """
        return self.median

    def __add__(self, other):
        """ Combine the existing PDF with an input one

        Parameters
        ----------
        other : DensityPDF

        Returns
        -------
        DensityPDF

        """
        # Check that the nH arrays are identical
        if not np.allclose(self.nH,other.nH):
            raise IOError("nH arrays need to be identical (for now)")
        # Add em'
        new = DensityPDF(self.nH, self.pdf_nH+other.pdf_nH, normalize=False)
        # Return
        return new

    def __repr__(self):
        repr = '<{:s}: meannH={:g}'.format(self.__class__.__name__, self.meannH)
        try:
            repr = repr + ', ' + self.parent.name
        except:
            pass
        repr = repr + '>'
        return repr


