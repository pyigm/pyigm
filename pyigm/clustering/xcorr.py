"""Module containing 2PCF Survey Class"""

class Survey2D2PCF(object):
    """
    Class for 2D-2PCF for absorber-galaxy two-point
    cross and auto-correlations. The Class loads several independent
    IgmGalaxyField objects together. It internally calculates the total numbers
    of cross- and auto-pairs of real data and produce randoms. It calculates
    the 2D2PCF using the Landy & Szalay estimator (or others), and it estimates
    the uncertainty with a bootstrap, jacknife or Landy & Szalay approximation."""

    def __init__(self, field):
        f = copy.deepcopy(field)
        self.fields = [f]

        self.DgDg = field.DgDg
        self.DgRg = field.DgRg
        self.RgRg = field.RgRg
        self.DaDg = field.DaDg
        self.DaRg = field.DaRg
        self.RaDg = field.RaDg
        self.RaRg = field.RaRg
        self.DaDa = field.DaDa
        self.DaRa = field.DaRa
        self.RaRa = field.RaRa

        self.Ngal = len(field.galreal)
        self.Nabs = len(field.absreal)

        self.rbinedges = field.rbinedges
        self.tbinedges = field.tbinedges
        self.rcbins = 0.5*(self.rbinedges[:-1] + self.rbinedges[1:])
        self.tcbins = 0.5*(self.tbinedges[:-1] + self.tbinedges[1:])
        self.rdiff  = self.rbinedges[1:] - self.rbinedges[:-1]
        self.tdiff  = self.tbinedges[1:] - self.tbinedges[:-1]
