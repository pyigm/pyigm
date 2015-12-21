""" Class for Metallicity PDF
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
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

    def __init__(self, ZH, pdf_ZH, parent=None, dZH=None):
        """
        Parameters
        ----------
        ZH : ndarray
          Values of ZH (log10, typically but not required) where the PDF has been evaluated
        pdf_ZH : ndarray
          PDF of the metallicity
        parent : object, optional
          Parent of the metalliicty PDF, e.g. an AbsSystem
        dZH : float or ndarray
          Width of the ZH 'bins' (log 10)

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
        self.normalize()

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

    def normalize(self):
        """ Normalize the PDF

        First checks to see if it is normalized
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


    def __repr__(self):
        repr = '<{:s}: meanZH={:g}'.format(self.__class__.__name__, self.meanZH)
        try:
            repr = repr + ', ' + self.parent.name
        except:
            pass
        repr = repr + '>'
        return repr


