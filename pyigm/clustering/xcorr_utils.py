"""Utilities for xcorr calculations"""
from __future__ import print_function, absolute_import, division, unicode_literals

import time
import numpy as np
import pdb
import numba as nb

from scipy.ndimage import gaussian_filter as gf
from scipy.interpolate import CubicSpline


from pyigm.clustering.randist import RanDist

Ckms = 299792.458  # speed of light in km/s

@nb.jit(nopython=True, cache=True)
def nb_auto_pairs_rt(x, y, z, rbinedges, tbinedges, wrap=True, track_time=False,
                  original=True):
    """
    [NT: give a nice description]

    Find the number of pairs in 2d and 3d for galaxies with
    coordinate x, y, z.

    x,y,z: comoving coordinates in Mpc.  x is (roughly) the redshift
           direction.

    rbinedges and tbinedges correspond to the arrays defining the bins
    edges """

    #start=time.clock()
    #x = np.array(X)
    #y = np.array(Y)
    #z = np.array(Z)
    #rbinedges = np.array(rbinedges)
    #tbinedges = np.array(tbinedges)

    npair_rt = np.zeros((len(rbinedges) - 1, len(tbinedges) - 1), dtype=nb.types.float64)#float)
    vals = np.zeros((len(rbinedges) - 1, len(tbinedges) - 1), dtype=nb.types.float64)#float)
    #npair_rt = np.zeros((len(rbinedges) - 1, len(tbinedges) - 1), float)
    #vals = np.zeros((len(rbinedges) - 1, len(tbinedges) - 1), float)

    # Ugly for loop
    for i in range(len(x) - 1):
        if i % 1000 == 0:
            print(i)
        # radial separation
        if wrap:
            rsep = np.abs(x[i+1:] - x[i])
        else:
            rsep = x[i+1:] - x[i]
        # transverse separation
        tsep = np.hypot(y[i+1:] - y[i], z[i+1:] - z[i])

        # Histogram me
        #vals, _ = np.histogramdd((rsep, tsep), (rbinedges, tbinedges),
        #                      range=None, normed=False, weights=None)
        # Loop me
        vals[:] = 0.
        #r_logic = []

        # Cut down to within the box!
        ok_sep = (rsep < rbinedges[-1]) & (tsep < tbinedges[-1])
        rsep = rsep[ok_sep]
        tsep = tsep[ok_sep]
        #
        r_logic = np.zeros((len(rsep), len(rbinedges) - 1), dtype=nb.types.boolean)#float)
        #r_logic = np.zeros((len(rsep), len(rbinedges) - 1), dtype=bool)
        for jj in range(len(rbinedges)-1):
            r_logic[:,jj]  = (rsep > rbinedges[jj]) & (rsep <= rbinedges[jj+1])
        #t_logic = np.zeros((len(tsep), len(tbinedges) - 1), dtype=bool)
        t_logic = np.zeros((len(tsep), len(tbinedges) - 1), dtype=nb.types.boolean)
        for jj in range(len(tbinedges)-1):
            t_logic[:,jj]  = (tsep > tbinedges[jj]) & (tsep <= tbinedges[jj+1])
        # Finish
        for jj in range(len(rbinedges)-1):
            for kk in range(len(tbinedges)-1):
                vals[jj,kk] = np.sum(r_logic[:,jj] & t_logic[:,kk])
        #pdb.set_trace()
        npair_rt += vals
        #npair_rt += vals2

    #end=time.clock()
    #print('\t Time elapsed = {} seconds.'.format(end-start))
    return npair_rt

def auto_pairs_rt(X, Y, Z, rbinedges, tbinedges, wrap=True, track_time=False,
                  original=True):
    """
    [NT: give a nice description]

    Find the number of pairs in 2d and 3d for galaxies with
    coordinate x, y, z.

    x,y,z: comoving coordinates in Mpc.  x is (roughly) the redshift
           direction.

    rbinedges and tbinedges correspond to the arrays defining the bins
    edges """

    start=time.clock()
    x = np.array(X)
    y = np.array(Y)
    z = np.array(Z)
    rbinedges = np.array(rbinedges)
    tbinedges = np.array(tbinedges)

    npair_rt = np.zeros((len(rbinedges) - 1, len(tbinedges) - 1), float)

    # Ugly for loop
    #   Could use map here..
    for i in range(len(x) - 1):
        # radial separation
        if wrap:
            rsep = np.abs(x[i+1:] - x[i])
        else:
            rsep = x[i+1:] - x[i]
        # transverse separation
        tsep = np.hypot(y[i+1:] - y[i], z[i+1:] - z[i])

        # Cut me (3x speed-up)
        ok_sep = (rsep <= rbinedges[-1]) & (tsep <= tbinedges[-1])
        rsep = rsep[ok_sep]
        tsep = tsep[ok_sep]

        # Histogram me (most expensive step, but faster than numba!)
        vals, _ = np.histogramdd((rsep, tsep), (rbinedges, tbinedges),
                              range=None, normed=False, weights=None)
        npair_rt += vals

    end=time.clock()
    print('\t Time elapsed = {} seconds.'.format(end-start))
    return npair_rt



def cross_pairs_rt(x1, y1, z1, x2, y2, z2, rbinedges, tbinedges,wrapped=True):
    """ Find the number of pairs in 2d and 3d for galaxies with
    coordinate (x1, y1, z1) and (x2,y2,z2).

    x,y,z: comoving coordinates in Mpc.  x is (roughly) the redshift
           direction
    """
    start=time.clock()
    x1 = np.array(x1)
    y1 = np.array(y1)
    z1 = np.array(z1)
    x2 = np.array(x2)
    y2 = np.array(y2)
    z2 = np.array(z2)
    rbinedges = np.array(rbinedges)
    tbinedges = np.array(tbinedges)

    npair_rt = np.zeros((len(rbinedges) - 1, len(tbinedges) - 1), float)

    if len(x1)>len(x2):
        auxx=x1
        auxy=y1
        auxz=z1
        x1=x2
        y1=y2
        z1=z2
        x2=auxx
        y2=auxy
        z2=auxz

    # Ugly for loop
    for i in range(len(x1) - 1):
        # radial separation
        if wrapped:
            rsep = np.abs(x1[i] - x2)
        else:
            rsep = x1[i] - x2
        # transverse separation
        tsep = np.hypot(y1[i]-y2, z1[i]-z2)

        # Cut me (3x speed-up)
        ok_sep = (rsep <= rbinedges[-1]) & (tsep <= tbinedges[-1])
        rsep = rsep[ok_sep]
        tsep = tsep[ok_sep]

        # Histogram
        vals, _ = np.histogramdd((rsep, tsep), (rbinedges, tbinedges),
                              range=None, normed=False, weights=None)
        npair_rt += vals

    end=time.clock()
    print('\t Time elapsed = {} seconds.'.format(end-start))
    return npair_rt


def W12(DD,RR,f=None):
    """ It computes,

    (1) the Natural/Peebles estimator: DD/(f*RR) - 1
    (2) its poissonian error

    for a given number of pairs. """

    Ndd = np.sum(DD) # different than Ngal*(Ngal-1)/2
    Nrr = np.sum(RR) # different than Nrdm*(Nrdm-1)/2

    if f==None:
        f = Ndd/Nrr
    RR = np.where(RR==0,1e-10,RR)

    W1 = DD / (f*RR)  - 1.
    W1 = np.where(RR==0,-2,W1)

    t1 = DD / ((f*RR)**2)
    t2 = (DD**2) / (RR**3)

    err_W1 = np.sqrt(t1 + t2)
    err_W1 = np.where(err_W1==0,1,err_W1)
    return W1, err_W1



def W3b(DD,RR,DR,RD,Ndd=None,Nrr=None,Ndr=None,Nrd=None):
    """ It computes,

    (1) the Landy&Szalay estimator: (DD - DR - RD + RR) / RR
    (2) its poissonian error

    for a given number of pairs. f1,f2,f3 corresponds to the
    normalization factors such that: W3 = DD/(f1*RR) - DR/(f2*RR) -
    RD/(f3*RR) + 1"""

    if Ndd is None:
        Ndd = np.sum(DD)
    if Nrr is None:
        Nrr = np.sum(RR)
    if Ndr is None:
        Ndr = np.sum(DR)
    if Nrd is None:
        Nrd = np.sum(RD)

    #normalize counts
    DD = DD / Ndd
    RR = RR / Nrr
    DR = DR / Ndr
    RD = RD / Nrd


    #assert np.sum(RR==0)==0, "RR has values equal to zero. Try again."
    RR = np.where(RR==0,1e-10,RR)

    W3 = DD/RR  - DR/RR - RD/RR + 1.
    W3 = np.where(DD==0,-1,W3)

    t1 = DD*Ndd / ((Nrr*RR)**2)
    t2 = DR*Ndr / ((Nrr*RR)**2)
    t3 = RD*Nrd / ((Nrr*RR)**2)
    t4 = ((DD*Ndd - DR*Ndr - RD*Nrd)**2 ) / ((Nrr*RR)**3)

    err_W3 = np.sqrt(t1 + t2 + t3 + t4)
    err_W3 = clean_matrix(err_W3,n=0)
    err_W3 = np.where(err_W3<=0,1.9,err_W3)

    return W3,err_W3


def W3(DD,RR,DR,RD,Ndd,Nrr,Ndr,Nrd):
    """Returns the Landy & Scalay estimator for the cross-correlation
        and its error. The pair counts DD, RR, DR, RD are not required to
    be normalized normalized. Ndd,Nrr,Ndr,Nrd are the normalizations
    factors for each pair count

    Parameters
    ----------
    DD : ndarray
    RR : ndarray
    DR : ndarray
    RD : ndarray
    Ndd : int
    Nrr : int
    Ndr : int
    Nrd : int

    Returns
    -------

    """

    #normalize counts
    DD = DD / Ndd
    RR = RR / Nrr
    DR = DR / Ndr
    RD = RD / Nrd

    RR = np.where(RR==0,1e-20,RR)
    W3 = DD/RR  - DR/RR - RD/RR + 1.

    #error from LS93, eqs. (48) (46) (43) and definition of Gp
    #Here Gp = RR when normalized
    err_W3 = np.sqrt( (1 + W3)**2 / (Ndd*RR) )
    #err_W3 = np.sqrt( (1 + W3) / (Ndd*DD) )

    #integral constrain correction, LS93 eq. (48) (46)
    C = np.sum(W3*RR)
    #print 'integral constraint: %s'%C
    W3 = (1 + C)*(1+W3) - 1
    if np.fabs(C)>1:
        print('integral constraint: {}'.format(C))


    #clean W3 and errors
    #W3     = clean_matrix(W3,n=0)
    W3     = np.where(RR==0,0,W3)
    W3     = np.where(DD==0,-1,W3)
    err_W3 = clean_matrix(err_W3,n=0)
    #err_W3 = np.where(err_W3<=0,1.9,err_W3)
    return W3, err_W3


def W_sig(W,err_W,s=1.):

    """ It return the values of W which are s times greater than
    W/err_W, the rest of the values are set to 0."""

    Ws = W
    Ws = np.where(W/err_W < s, -1, Ws)
    return Ws



def clean_matrix(W,n=-99):
    W = np.where(np.isnan(W),n,W)
    W = np.where(np.isinf(W),n,W)
    return W


def spline_sensitivity(galreal, Nmin=20, magmin=17., magmax=26., delta_mag=0.5, DZ=0.01,
                       smooth_scale=10., debug=False):
    """
    For a given galaxy with a
    given magnitude (and other properties), it calculates the redshift
    sensitivity function from galaxies in a magnitude band around the
    selected one (i.e., including slightly brighter and fainter
    galaxies),
    For extremely bright or faint galaxies (rare) the sensitivity function
    is calculated from at least Nmin (=20) galaxies (i.e. the magnitude
    band is increased).

    Parameters
    ----------
    galreal : np.recarray
    Nmin : int, optional
      Minimum number of galaxies to include in a bin
    DZ : float, optional
      delta z for the histogram for getting the
    smooth_scale : float, optional
       smoothing scale for the histogram (in number of bins, so depends on DZ)
    magmin : float, optional
       Minimum magnitude for the binning (bright end)
    magmax : float, optional
       Maximum magnitude for the binning (faint end)
    delta_mag : float, optional
       Step size for magnitude binning
    delta_mag : float (optional)
    debug : bool (optional)

    Returns
    -------
    VALS : dict
       dict containing the redshift histograms used in the CubicSpline analysis
    SPL : dict
      CubicSplines
    magbins : ndarray
      Mag bins used in the generation of the Splines
    """
    Nmin = int(Nmin)  # minimum number of galaxys for the fit
    zmin = np.max([np.min(galreal.ZGAL[galreal.ZGAL > 0]), 1e-9])
    zmax = np.max(galreal.ZGAL)

    # spline in z
    galreal.sort(order=str('MAG'))  # np.recarray.sort()

    bins = np.append(np.linspace(0, zmin, 20), np.arange(zmin + DZ, zmax + 10 * DZ, DZ))

    # Make subhistograms depending on magnitude. Use dictionary.
    VALS = dict()
    SPL = dict()
    aux_hist, _ = np.histogram(galreal.ZGAL, bins)
    VALS['all'] = gf(aux_hist.astype(float), smooth_scale)  # smooth the histogram
    SPL['all'] = CubicSpline(0.5 * (bins[:-1] + bins[1:]), VALS['all'].astype(float))

    # Generate magnitude bins
    magbins = np.arange(magmin, magmax, delta_mag)
    for mag in magbins:
        delta_mag2 = delta_mag
        q = 0
        # Insure there are Nmin galaxies in the bin by extending the magnitude bin if needed
        while True:
            cond = (galreal.MAG <= mag + delta_mag2 * 0.5) & (galreal.MAG > mag - delta_mag2 * 0.5)
            if np.sum(cond) >= Nmin:
                break
            else:
                delta_mag2 += 0.25
            q += 1
            assert q < 1000, 'Something wrong with the redshift distribution'

        aux_hist, _ = np.histogram(galreal.ZGAL[cond], bins)
        VALS['{}'.format(mag)] = gf(aux_hist.astype(float), smooth_scale)  # smooth the histogram
        SPL['{}'.format(mag)] = CubicSpline(0.5 * (bins[:-1] + bins[1:]), VALS['{}'.format(mag)].astype(float))
        spl = SPL['{}'.format(mag)]
        if debug:
            import matplotlib.pyplot as pl
            pl.plot(bins, spl(bins), '-', label='{}'.format(mag))
            pl.plot(bins[:-1], aux_hist, drawstyle='steps-mid')
            pl.xlim(0, 2)
            pl.legend()
            pl.show()

    if debug:
        import matplotlib.pyplot as pl
        pl.legend()
        pl.show()
    # Return
    return magbins, VALS, SPL


def random_abs_zmnx(absreal, Nrand, zmnx, wrest, dv_Galactic=100.):
    """From a real absorber catalog it creates a random catalog.
    Absorbers are simply placed randomly between redshift limits
    zmnx[0] and zmnx[1]

    Parameters
    ----------
    absreal
    Nrand
    zmnx

    Returns
    -------

    """

    absrand = absreal.repeat(Nrand)

    randz = np.random.uniform(low=zmnx[0], high=zmnx[1], size=2*len(absrand))  # Buffer to reject Galaxy

    # Avoid Galactic
    galactic = mask_galactic()
    Galz = galactic/wrest - 1.
    z_Gal = np.outer(np.ones_like(randz), Galz)
    # Diff
    diff = z_Gal - np.outer(randz, np.ones_like(Galz))
    mdiff = np.amin(np.abs(diff), axis=1)
    # Good ones
    gdz = mdiff > (dv_Galactic/Ckms) / (randz+1)
    gdi = np.where(gdz)[0]
    absrand.ZABS = randz[gdi[0:len(absrand)]]
    # Return
    return absrand


def random_abs_W(absreal, Nrand, wa, fl, er, sl=3., R=20000, FWHM=10., ion='HI'):
    """From a real absorber catalog it creates a random catalog.  For
    a given real absorber with (z_obs,logN_obs,b_obs) it places it at
    a new z_rand, defined by where the line could have been
    observed.

    Input parameters:
    ---
    absreal: numpy rec array with the absorber catalog.
    Nrand:   number of random lines per real one generated (integer).
    wa:      numpy array of wavelength covered by the spectrum.
    fl:      numpy array of normalized flux.
    er:      numpy array of error in the normalized flux of the spectrum for
             a given wavelength.
    sl:      significance level for the detection of the absorption line.
    R:       resolution of the spectrograph, assumed constant
    FWHM:    Full-width at half maximum in pixels (assumed constant). This
             parameter defines the smoothing scale for Wmin.
    ion:     Name of the ion. Function only valid for HI so far.

    From the error we calculate the Wmin = sl * wa * er / (1+z) / R,
    where z = wa/w0 - 1 (w0 is the rest frame wavelenght of the
    transition) and R is the resolution of the spectrograph. We then
    smooth Wmin with a boxcar along FWHM pixels.

    For the given absorber we transform (logN_obs,b_obs) to a W_obs assuming
    linear part of the curve-of-growth.

    We then compute the redshifts where W_obs could have been observed
    according to the given Wmin, and place Nrand new absorbers with
    the same properties as the given one accordingly.
    """
    from linetools.analysis.absline import Wr_from_N_b_transition
    from linetools.lists.linelist import LineList
    ism = LineList('ISM')

    Nrand = int(Nrand)
    absrand = absreal.repeat(Nrand)

    zmin = np.min(absreal.ZABS)
    zmax = np.max(absreal.ZABS)

    if ion == 'HI':
        z_Lya, Wmin_Lya = compute_Wmin(wa, fl, er, sl=sl, R=R, FWHM=FWHM, ion='HI')
        z_Lyb, Wmin_Lyb = compute_Wmin(wa, fl, er, sl=sl, R=R, FWHM=FWHM, ion='HILyb')

    for i in range(len(absreal)):

        if absreal.ZABS[i] > np.max(z_Lya):  # lines that were observed through Lyb
            Wr = Wr_from_N_b_transition(absreal.LOGN[i], absreal.B[i], 'HI 1025', linelist=ism)
            #Wr = logN_b_to_Wr(absreal.LOGN[i], absreal.B[i], ion='HILyb')
            z = z_Lyb
            z = np.where(z <= z_Lya, -1., z)  # mask out region with Lya coverage
            Wmin = Wmin_Lyb
        else:  # lines that were observed through Lya only
            Wr = Wr_from_N_b_transition(absreal.LOGN[i], absreal.B[i], 'HI 1215', linelist=ism)
            #Wr = logN_b_to_Wr(absreal.LOGN[i], absreal.B[i], ion='HI')
            z = z_Lya
            Wmin = Wmin_Lya

        zgood = (Wr > Wmin) & (z >= zmin) & (z < zmax)

        if np.sum(zgood) == 0:
            from matplotlib import pyplot as plt
            plt.plot(z, Wmin, drawstyle='steps-mid')
            plt.axis([z[0], z[-1], 0, 0.1])
            plt.show()
            print(Wmin)
        assert np.sum(zgood) > 0, \
            'There are not regions in the spectrum with Wmin<{} A. The minimum is {}. Adjust significance.'.format(
                Wr, np.min(Wmin))

        # Random time
        rand_z = RanDist(z, zgood * 1.)
        zrand = rand_z.random(Nrand)
        absrand.ZABS[i * Nrand:(i + 1) * Nrand] = zrand

    return absrand

def mask_galactic():
    # masked regions (potential Galactic absorption)
    galactic = np.array([1334.5323,  # CII
                         1335.7077, # CII*
                         1238.821,  # NV
                         1242.804,  # NV
                         1302.1685,  # OI
                         #1304.8576,  # OI*
                         #1306.0286,  # OI**
                         1304.3702,  # SiII
                         1260.4221,  # SiII
                         1526.707,  # SiII
                         1548.204,  # CIV
                         1550.781,  # CIV
                         1334.8132,  # PIII
                         1259.519,  # SII
                         1253.811,  # SII
                         1250.584,  # SII
                         1670.7886, # AlII
                         1608.4511])  # FeII
    return galactic


def compute_Wmin(wa, fl, er, sl=3., R=20000, FWHM=10, ion='HI', dv_mask=200.):
    """For a given spectrum and transition, it computes the minimun
    rest-frame equivalent width for that transition to be observed. It
    return a tuple of redshift and Wmin (z,Wmin)"""
    from scipy.ndimage import uniform_filter as uf
    #
    if ion == 'HI':
        w0 = 1215.67  # HI Lya w0 in angstroms
    if ion == 'HILyb':
        w0 = 1025.72  # HI Lyb w0 in angstroms
    # Mask galactic
    galactic = mask_galactic()
    zgal = galactic / w0 - 1.
    dzgal = (zgal + 1) * dv_mask / Ckms

    wa = np.array(wa)
    fl = np.array(fl)
    er = np.array(er)

    z = wa / w0 - 1.  # spectrum in z coordinates

    Wmin = sl * w0 * er / R / fl  # sl*wa / (1. + z) / R / (S/N)
    Wmin = np.where(Wmin <= 0, 1e10, Wmin)
    Wmin = np.where(np.isnan(Wmin), 1e10, Wmin)
    Wmin = np.where(np.isinf(Wmin), 1e10, Wmin)
    Wmin = uf(Wmin.astype(float), FWHM)  # smoothed version (uniform preferred over gaussian)
    for zi, dzi in zip(zgal, dzgal):
        cond = (z > zi - dzi) & (z < zi + dzi)
        Wmin = np.where(cond, 1e10, Wmin)
    return z, Wmin

def random_gal(galreal, Nrand, magbins, SPL):
    """
    Preferred random galaxy generator.

    Places Nrand new galaxies at a random redshift given
    by a smoothed version of the observed sensitivity function.

    Parameters
    ----------
    galreal
    Nrand : int
      Number of random galaxies generated per real galaxy
    magbins : ndarray
    SPL : dict

    Returns
    -------
    galrand : np rec array
      Copy of galreal, Nrand times
    """

    # Init
    galrand = galreal.repeat(Nrand)
    zmin = np.max([np.min(galreal.ZGAL[galreal.ZGAL > 0]), 1e-9])
    zmax = np.max(galreal.ZGAL)
    rvals = np.linspace(0, zmax, 1e4)

    # TODO -- Use better masking than +/- 90 mag
    # Awful for loop, but it appears to run fast
    for i in range(len(galreal)):
        if (galreal.MAG[i] > 90) or (galreal.MAG[i] < -90):  # no magnitude, use the whole distribution
            #vals = VALS['all']
            spl = SPL['all']
        else:
            # And this is ugly expensive
            ind_mag = np.where(np.fabs(galreal.MAG[i] - magbins) == np.min(np.fabs(galreal.MAG[i] - magbins)))[0][0]
            mag = magbins[ind_mag]
            #vals = VALS['{}'.format(mag)]
            spl = SPL['{}'.format(mag)]
        if i % 1000 == 0:
            print('{}/{}'.format(i + 1, len(galreal)))

        dist = np.array(spl(rvals))
        dist = np.where((rvals < zmin) | (rvals > zmax), 0, dist)  # get rid of redshifts beyond observed
        rand_z = RanDist(rvals, dist)
        zrand = rand_z.random(Nrand)
        galrand.ZGAL[i * Nrand:(i + 1) * Nrand] = zrand
    # Return
    return galrand


def collapse_along_LOS(DD, nbins=None, s=0):
    """Sums pair counts over the first nbins along the sightline
    dimension. Returns an array with the values for transverse bins. If
    nbins is None then collapses the whole array.

    Parameters
    ----------
    DD : ndarray
      Pair counts to sum
    nbins :  int, optional
      Number of bins in the radial dimension to collapse along
      If None, take them all
    s : float, optional
      For Gaussian filtering

    Returns
    -------
    DD_1D : ndarray
      Collapsed pair counts
    """
    # Gaussian filter
    if s > 0:
        sDD = gf(DD, s)
    else:
        sDD = DD
    if nbins is None:
        nbins = sDD.shape[0]
    # Avoidable loop?
    #old_DD_1D = np.array([np.sum(sDD.T[i][:nbins]) for i in range(sDD.shape[1])])
    DD_1D = np.sum(sDD[:nbins,:], axis=0)
    # Return
    return DD_1D
