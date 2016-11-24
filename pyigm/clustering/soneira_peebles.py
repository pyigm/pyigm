""" Module for Soneira Peebles algorithm
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb
import yaml
import imp

import numpy as np
import matplotlib.pyplot as plt


def get_xy(ncenters, ndim, rstate):
    """ Generate 'xyz' positions in a n-dimension 'sphere' with
    brute force.

    Parameters
    ----------
    ncenters : int
      Number of 'spheres' to generate
    ndim : int
      Number of dimensions
    rstate : RandomState
      Random number state

    Returns
    -------
    'xyz' centers:  (ncenters, ndim)

    """
    flag = True
    while flag:
        buff = 10 * ncenters
        args = [buff,ndim]
        randxy = rstate.rand(*args)*2 - 1.  # Multi-dimensional array of -1 to 1 values
        # Generate Radius**2
        radius2 = np.sum(randxy**2,1)
        # Grab the good ones
        gdr2 = np.where(radius2 < 1.)[0]
        # Check we have enough
        if len(gdr2) > ncenters:
            flag = False
    # Return the first ncenters
    return randxy[gdr2[0:ncenters],:]


def sp_cluster(L, R, eta, l, ndim, seed=None):
    """
    Parameters
    ----------
    L : int
      Number of levels to generate
    R : float
      Radius of 0th sphere
    eta : int
      Sphere multiplier
    l : float
      Sphere radius shrinking factor
    ndim

    Returns
    -------

    """
    # Random state
    if seed is None:
        seed = 12345
    rstate = np.random.RandomState(seed)

    #
    fdict = dict(pos=np.zeros(ndim), parent=0, level=0, radius=R)
    idict = dict(L=L, R=R, eta=eta, l=l, ndim=ndim)
    sp_dict = {0:[fdict], 'input':idict}

    #
    radius = R
    for kk in range(L):
        level = kk+1
        sp_dict[level] = []
        nparents = len(sp_dict[kk])
        # Generate sphere positions
        xyval = get_xy(nparents*eta, ndim, rstate)
        # Offset and save
        cnt = 0
        for jj in range(nparents):
            for ii in range(eta):
                # Make the dict
                xpos = sp_dict[kk][jj]['pos'] + xyval[cnt,:] * radius
                new_dict = dict(pos=xpos, radius=radius/l, parent=jj, level=level, idx=cnt)
                sp_dict[level].append(new_dict)
                cnt += 1
                #pdb.set_trace()
        # Shrink the radius
        radius = radius / l
    # Return
    return sp_dict


def plot_1d(sp_dict):
    """ Check 1D output
    Parameters
    ----------
    sp_dict
    """
    if sp_dict['input']['ndim'] != 1:
        raise IOError("Only lines")

    clrs = ['r','b','g','y','k']
    plt.clf()
    ax = plt.gca()
    nlevel = len(sp_dict.keys()) - 1
    for level in range(nlevel):
        off = level*0.1
        for item in sp_dict[level]:
            xy = item['pos']
            radius = item['radius']
            ax.fill_between([xy[0]-radius, xy[0]+radius], [-1+off,-1+off], [1-off,1-off], alpha=0.25, color='b') #color=clrs[level]

    ax.set_xlim(sp_dict[0][0]['radius']*-2, sp_dict[0][0]['radius']*2)
    ax.set_ylim(-2., 2)
    plt.show()
    plt.close()


def plot_circles(sp_dict):

    #
    if sp_dict['input']['ndim'] != 2:
        raise IOError("Only circles")

    clrs = ['r','b','g','y','k']
    plt.clf()
    nlevel = len(sp_dict.keys()) - 1
    for level in range(nlevel):
        for item in sp_dict[level]:
            xy = item['pos']
            radius = item['radius']
            circle=plt.Circle((xy[0],xy[1]),radius,color=clrs[level], fill=False)
            plt.gca().add_artist(circle)
    ax = plt.gca()
    ax.set_xlim(sp_dict[0][0]['radius']*-2, sp_dict[0][0]['radius']*2)
    ax.set_ylim(sp_dict[0][0]['radius']*-2, sp_dict[0][0]['radius']*2)
    plt.show()
    plt.close()


# MAIN
if __name__ == "__main__":

    flg_tst = 0
    flg_tst += 2**0  # Test 1D
    #flg_tst += 2**1  # Test 2D

    #
    if (flg_tst % 2**1) >= 2**0:  # Test 2D
        sp_dict = sp_cluster(2, 1., 4, 4., ndim=1)
        plot_1d(sp_dict)

    if (flg_tst % 2**2) >= 2**1:  # Test 2D
        sp_dict = sp_cluster(2, 1., 4, 4., 2)
        plot_circles(sp_dict)
