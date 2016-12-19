""" Class for f(N) constraints
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import imp
import pdb

from astropy.io import fits
from astropy import cosmology

# Path for pyigm
pyigm_path = imp.find_module('pyigm')[1]


class FNConstraint(object):
    """A Class for fN constraints

    Parameters
    ----------
    fN_dtype : str
    zeval : float
        Redshift where the constraint is evaluated
    ref : str
        Reference
    flavor : str
        Specific type of constraint

    Attributes
    ----------
    fN_dtype : str
        Constraint type for the fN
           'fN' -- Standard f(N) evaluation
           'MFP' -- MFP
           'LLS' -- LLS incidence 
           'teff' -- tau effective
           'beta' -- slope constraint
    flavor : str
        Specific type of constraint
    comment : str
    ref : str
        Reference
    cosm : astropy.cosmology, optional
        Cosmology
    zeval : float
        Redshift where the constraint is evaluated
    data : dict
        Dictionary containing the constraints
    """

    @classmethod
    def from_row(cls, ftype, row):
        """  Read from binary FITS table
        Parameters
        ----------
        row : astropy.table.row
        """
        slf = cls(ftype)
        # Common items
        common = ['REF','COSM','TYPE','COMMENT']
        slf.ref = row['REF']
        if row['COSM'] in ['VANILLA', '', 'h=0.7, Om=0.3, OL=0.7']:
            cosmo = cosmology.core.FlatLambdaCDM(70., 0.3, Ob0=0.0455)
            #sigma_8 = 0.80
        else:
            pdb.set_trace()
            raise ValueError('Not ready for this type of cosmology')
        slf.cosmo = cosmo
        slf.flavor = row['TYPE']
        slf.comment = row['COMMENT']

        # zeval
        if 'ZEVAL' in row.array.names:
            slf.zeval = row['ZEVAL']
        elif 'Z_LLS' in row.array.names:
            slf.zeval = row['Z_LLS']
        elif 'Z_MFP' in row.array.names:
            slf.zeval = row['Z_MFP']
        elif 'Z_TEFF' in row.array.names:
            slf.zeval = row['Z_TEFF']
        else:
            raise ValueError('fN.data: No redshift info!')

        # zip the rest
        slf.data = dict(zip(row.array.names,row))
        for item in common:
            slf.data.pop(item) # No need to duplicate
        # Return
        return slf

    @classmethod
    def from_fitsfile(cls, fits_file):
        """ Build up a list of fN constraints from a multi-extension FITS file

        Parameters
        ----------
        fits_file : str or list
           Name of FITS file

        Returns
        -------
        fN_list : list
           List of FNConstraint objects

        """
        # List of constraints
        fN_cs = []

        # Read
        if isinstance(fits_file, list):
            for ifile in fits_file:
                tmp_cs = cls.from_fitsfile(ifile)
                for cs in tmp_cs:
                    fN_cs.append(cs)
        else:
            hdus = fits.open(fits_file)
            if len(hdus) == 1:
                raise ValueError('Expecting a multi-extension fits file -- {:s}'.format(
                                fits_file))
            # Loop through hdu
            for hdu in hdus[1:]:
                data = hdu.data
                # Get ftype
                if 'FN' in data.dtype.names:
                    ftype = 'fN' # Standard f(N) data
                elif 'TAU_LIM' in data.dtype.names:
                    ftype = 'LLS' # LLS survey
                elif 'MFP' in data.dtype.names:
                    ftype = 'MFP' # MFP measurement
                elif 'TEFF' in data.dtype.names:
                    ftype = 'teff' # tau effective (Lya)
                else:
                    raise ValueError('Cannot figure out ftype')

                # Loop on the Table
                for row in data:
                    fNc = cls.from_row(ftype, row)
                    fN_cs.append(fNc)

        # Return
        return fN_cs

    @classmethod
    def load_defaults(cls):
        """ Load default constraints on f(N)

        Returns
        -------
        all_fN_cs : list
          list of FNConstraint objects
        """
        fn_file = pyigm_path+'/data/fN/fN_constraints_z2.5_vanilla.fits'
        k13r13_file = pyigm_path+'/data/fN/fN_constraints_K13R13_vanilla.fits'
        n12_file = pyigm_path+'/data/fN/fN_constraints_N12_vanilla.fits'
        # Load
        all_fN_cs = cls.from_fitsfile([fn_file,k13r13_file, n12_file])
        # Return
        return all_fN_cs

    # Initialize with type
    def __init__(self, fN_dtype, zeval=0., ref='', flavor='', cosmo=None):
        if fN_dtype not in ['fN', 'MFP', 'LLS', 'teff', 'beta']:
            raise IOError('Bad f(N) constraint')
        self.fN_dtype = fN_dtype
        self.zeval = zeval
        self.ref = ref
        self.flavor = flavor
        self.cosmo = cosmo


    def __repr__(self):
        return ('<{:s}: {:s}_{:s} zeval={:g}, ref={:s}>'.format(
                self.__class__.__name__, self.fN_dtype, self.flavor,
                self.zeval, self.ref))


def fN_data_from_ascii_file(infile):
    """
    Parameters
    ----------
    infile

    Returns
    -------

    """

    assert False # Bad code
    #makes new fN constraint with data type fN
    fNc = FNConstraint('fN')
    ftype = fNc.fN_dtype.encode('ascii')
    fNc.fN_dtype = ftype
    fNc.ref=infile.encode('ascii')
    
    # Open file
    f = open(infile, 'r')

    # Read and ignore header lines
    firstline = f.readline()
    # get rid of newline /n symbol
    firstline =firstline.strip()
    #get zeval and DX from first line
    values = firstline.split()
    fNc.zeval = float(values[0])
    ZEVAL = float(values[0])
    DX = float(values[1])

    #declaration of variables
    BINS1 =[]
    BINS2 = []
    fn = []
    SIG_FN1 = []
    SIG_FN2 = []
    count = 0
    numlines=0

    # Loop over lines and extract info
    for line in f:
        line = line.strip()
        columns = line.split()
        BINS1.append(float(columns[0]))
        BINS2.append(float(columns[1]))
    fn.append(float(columns[2]))
    SIG_FN1.append(float(columns[3]))
    SIG_FN2.append(float(columns[3]))
    numlines +=1
    if (float(columns[0])!=0) or (float(columns[1])!=0) or (float(columns[2])!=0) or (float(columns[3])!=0):
        count +=1
    f.close()

    NPT = int(count)
    bins = []
    bins.append(BINS1)
    bins.append(BINS2)
    sig_fn = []
    sig_fn.append(SIG_FN1)
    sig_fn.append(SIG_FN2)
    
    BINS = np.ndarray(shape=(2, numlines), dtype=float, buffer=np.array(bins))
    SIG_FN = np.ndarray(shape=(2, numlines), dtype=float, buffer=np.array(sig_fn))
    FN = np.ndarray(shape=(numlines,), dtype=float, buffer=np.array(fn))
    
    #makes array with names in ASCII not unicode
    arrayofnames = ['BINS','FN','SIG_FN','DX','NPT','ZEVAL']
    names = []
    for name in arrayofnames:
        newname = name.encode('ascii')
        names.append(newname)

    values = [BINS,FN,SIG_FN,DX,NPT,ZEVAL]

    fNc.data = dict(zip(names, values))
    
    return fNc


