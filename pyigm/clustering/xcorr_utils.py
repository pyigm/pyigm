"""Utilities for xcorr calculations"""
import time
import numpy as np

def auto_pairs_rt(X, Y, Z, rbinedges, tbinedges, wrap=True, track_time=False):
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

    for i in range(len(x) - 1):
        # radial separation
        if wrap:
            rsep = np.abs(x[i+1:] - x[i])
        else:
            rsep = x[i+1:] - x[i]
        # transverse separation
        tsep = np.hypot(y[i+1:] - y[i], z[i+1:] - z[i])

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

    for i in range(len(x1) - 1):
        # radial separation
        if wrapped:
            rsep = np.abs(x1[i] - x2)
        else:
            rsep = x1[i] - x2
        # transverse separation
        tsep = np.hypot(y1[i]-y2, z1[i]-z2)

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

def W3(DD,RR,DR,RD,Ndd=None,Nrr=None,Ndr=None,Nrd=None):
    """Returns the Landy & Scalay estimator for the cross-correlation
    and its error. The pair counts DD, RR, DR, RD are not required to
    be normalized normalized. Ndd,Nrr,Ndr,Nrd are the normalizations
    factors for each pair count (if None, the normalizations are taken
    from the sum over the whole array."""

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


def random_gal(galreal, Nrand, Nmin=20, DZ=0.01, smooth_scale=10.):
    """
    Prefered random galaxy generator. For a given galaxy with a
    given magnitude (and other properties), it calculates the redshift
    sensitivity function from galaxies in a magnitude band around the
    selected one (i.e., including slightly brighter and fainter
    galaxies), and places Nrand new galaxies at a random redshift given
    by a smoothed version of the observed sensitivity function. For
    extremely bright or faint galaxies (rare) the sensitivity function
    is calculated from at least Nmin (=20) galaxies (i.e. the magnitude
    band is increased).

    Parameters
    ----------
    galreal
    Nrand
    Nmin : int, optional
    DZ : float, optional
      delta z for the histogram for getting the
    smooth_scale : float, optional
       smoothing scale for the histogram (in number of bins, so depends on DZ)

    Returns
    -------
    galrand : np rec array
      Copy of galreal, Nrand times
    """
    from pyigm.clustering.randist import RanDist
    #from pyntejos.fit import InterpCubicSpline
    from scipy.interpolate import CubicSpline
    from scipy.ndimage import gaussian_filter as gf

    debug = 0
    Ckms = 299792.458
    Nmin = int(Nmin)  # minimum number of galaxys for the fit
    zmin = np.max([np.min(galreal.ZGAL[galreal.ZGAL > 0]), 1e-9])
    zmax = np.max(galreal.ZGAL)
    # spline in z
    galreal.sort(order='MAG')  # np.recarray.sort()
    galrand = galreal.repeat(Nrand)
    delta_mag = 0.5

    bins = np.append(np.linspace(0, zmin, 20), np.arange(zmin + DZ, zmax + 10 * DZ, DZ))
    # bins = np.arange(0, zmax+DZ, DZ)

    # Make subhistograms depending on magnitude. Use dictionary.
    VALS = dict()
    SPL = dict()
    aux_hist, _ = np.histogram(galreal.ZGAL, bins)
    VALS['all'] = gf(aux_hist.astype(float), smooth_scale)  # smooth the histogram
    #SPL['all'] = InterpCubicSpline(0.5 * (bins[:-1] + bins[1:]), VALS['all'].astype(float))
    SPL['all'] = CubicSpline(0.5 * (bins[:-1] + bins[1:]), VALS['all'].astype(float))

    delta_mag2 = delta_mag
    magbins = np.arange(17, 26, delta_mag2)
    # magbins = [15,19,21,22,26]
    rvals = np.linspace(0, zmax, 1e4)
    for mag in magbins:
        q = 0
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
        vals = VALS['{}'.format(mag)]
        spl = SPL['{}'.format(mag)]
        rand_z = RanDist(rvals, spl(rvals))
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

    for i in range(len(galreal)):
        if (galreal.MAG[i] > 90) or (galreal.MAG[i] < -90):  # no magnitude, use the whole distribution
            vals = VALS['all']
            spl = SPL['all']
        else:
            ind_mag = np.where(np.fabs(galreal.MAG[i] - magbins) == np.min(np.fabs(galreal.MAG[i] - magbins)))[0][0]
            mag = magbins[ind_mag]
            vals = VALS['{}'.format(mag)]
            spl = SPL['{}'.format(mag)]
        if i % 1000 == 0:
            print('{}/{}'.format(i + 1, len(galreal)))

        dist = np.array(spl(rvals))
        dist = np.where((rvals < zmin) | (rvals > zmax), 0, dist)  # get rid of redshifts beyond observed
        rand_z = RanDist(rvals, dist)
        zrand = rand_z.random(Nrand)
        galrand.ZGAL[i * Nrand:(i + 1) * Nrand] = zrand
    return galrand
