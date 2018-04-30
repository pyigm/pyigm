""" Module containing the ClusteringField class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

from pyigm.field.igmfield import IgmGalaxyField

class ClusteringField(IgmGalaxyField):
    """The Field class is meant to contain information from a given
    galaxy field (with one or more QSO sightlines containing IGM
    information).

    Input parameters:
    ---
    absreal:   catalog of real absorption lines (numpy rec array). It has to have
               dtype.names RA,DEC,ZABS,LOGN,B
    galreal:   catalog of real galaxies (numpy rec array). It has to have
               dtype.names RA, DEC, ZGAL, MAG
    wa:        numpy array with QSO spectral coverage (wavelenght)
    fl:        numpy aray with the QSO spectrum flux it can be normalized or not)
    er:        numpy array with the QSO spectrum error (it can be normalized
               or not, but has to be consistent with fl)
    R:         resolution of the QSO spectrum spectrograph
    Ngal_rand: number of random galaxies per real one that will be created.
    Nabs_rand: number of random absorbers per real one that will be created.
    proper:    calculates everything in physical Mpc rather than co-moving (Boolean).


    Description of the Class:

    It first creates absrand and galrand using random_abs() and
    random_gal() functions. CRA and CDEC are the center position of the
    field, and so will define the coordinate system. Galaxies and
    absorbers (RA,DEC,z) are transformed to (X,Y,Z) co-moving
    coordinates (assuming proper=False). It will then calculate
    cross-pairs DaDg,DaRg,RaDg,RaRg and auto-pairs DgDg, DgRg, RgRg and
    DaDa, DaRa, RaRa where D means 'data' R means 'random', a means
    'absorber' and g means 'galaxy'. The class has also some
    implemented plots that are useful.

    """
    def __init__(self, R=20000, Ngal_rand=10, Nabs_rand=100, proper=False, field_name='', wrapped=True):
        self.absreal = absreal  # np rec array with absorber properties (single ion)
        self.galreal = galreal  # np rec array with galaxy properties
        self.Ngal_rand = Ngal_rand  # Ngal_rand x len(galreal) = NRANDOM (Gal)
        self.Nabs_rand = Nabs_rand  # Nabs_rand x len(absreal) = NRANDOM (Abs)
        self.absrand = random_abs(self.absreal, self.Nabs_rand, wa, fl, er, R=R)
        self.galrand = random_gal(self.galreal, self.Ngal_rand)
        self.CRA = np.mean(self.absreal.RA)
        self.CDEC = np.mean(self.absreal.DEC)
        self.name = field_name
        self.proper = proper
        self.wrapped = wrapped

        IGMSystem.__init__(self, radec, zabs, vlim, NHI=NHI, abs_type='DLA', **kwargs)
