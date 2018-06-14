"""Contains a class that can used to generate random numbers from an
arbitrary (discrete) distribution, mostly taken from Neil Crighton's
old astro package.  Taken from pyntejos on GitHub

"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

npversion = map(int, np.__version__.split('.'))


class RanDist(object):
    """ Take random samples from an arbitrary discrete distribution."""

    def __init__(self, x, dist):
        """
        Inputs:
        x = values sampling the probability density
        dist = probability density distribution at each x value

        This finds the normalised probability density function
        (actually the mass density function because it's discrete) and
        the cumulative distribution function.

        Make sure the probability distribution is sampled densely
        enough (i.e. there are enough x values to properly define the
        shape of both the cdf and its inverse), because linear
        interpolation is used between the provided x values and cdf to
        infer new x values when generating the random numbers. A log
        sampling of x values is appropriate for distributions like
        inverse power laws, for example.
        """

        # normalise such that area under pdf is 1.
        self.pdf = dist / np.trapz(dist, x=x)
        # cumulative probability distribution
        self.cdf = dist.cumsum()
        self.cdf = self.cdf / float(self.cdf[-1])
        self.x = x

    def random(self, N=1, seed=None):
        """Return N random numbers with the requested distribution."""
        if seed is not None:  np.random.seed(seed)
        i = np.random.rand(N)
        y = np.interp(i, self.cdf, self.x)
        return y

    def plot_pdf(self):
        import matplotlib.pyplot as pl
        pl.plot(self.x, self.pdf)

    def self_test(self, N=1e4, log=False, seed=None, nbins=50):
        """ Make plots of the CDF, the PDF, and a histogram of N
        random samples from the distribution.
        """
        import matplotlib.pyplot as pl
        pl.figure()
        pl.subplots_adjust(hspace=0.001)
        # The cdf
        ax = pl.subplot(211)
        if log:
            ax.semilogx(self.x, self.cdf, 'b-')
        else:
            ax.plot(self.x, self.cdf, 'b-')
        ax.set_ylabel('cdf')
        ax.set_ylim(0, 1)
        # ax.set_xticklabels([])

        # The actual generated numbers
        ax = pl.subplot(212, sharex=ax)
        y = self.random(N, seed=seed)
        if log:
            bins = np.logspace(np.log10(min(y)), np.log10(max(y)), nbins)
        else:
            bins = nbins

        if npversion < [1, 3, 0]:
            vals, edges = np.histogram(y, bins=bins, normed=True, new=True)
        else:
            vals, edges = np.histogram(y, normed=True, bins=bins)
        if log:
            ax.loglog(self.x, np.where(self.pdf > 0, self.pdf, 1e-20), 'r-',
                      label='pdf')
            ax.loglog(edges[:-1], np.where(vals > 0, vals, 1e-20),
                      ls='steps-post', label='random values')
        else:
            ax.plot(self.x, self.pdf, 'r-', label='pdf')
            ax.plot(edges[:-1], vals, ls='steps-post', label='random values')
        ax.legend(frameon=False)
        ax.set_ylabel('pdf')



