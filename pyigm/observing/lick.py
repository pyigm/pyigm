'''
#;+ 
#; NAME:
#; lick
#;    Version 1.1
#;
#; PURPOSE:  Simple scripts for Lick obs
#;   2015 Written by JXP
#;-
#;------------------------------------------------------------------------------
'''

# Import libraries
import numpy as np
from numpy.ma.core import MaskedConstant
import os, subprocess

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

from xastropy.xutils import xdebug as xdb
from xastropy.obs import finder as x_finder
from xastropy.obs import radec as x_rad

# def wiki :: Generate a Wiki Table
# def starlist :: Generate a starlist file

#### ###############################
def wiki(targs, keys, fndr_pth=None, dbx_pth=None, outfil=None, skip_finder=False):
    """
    Generate a Wiki table for Lick observing.
    Should work for any of the Wiki pages

    Parameters:
    ----------
    targs: Table (RA, DEC keys required)
    keys: List
      List of keys to include in the Table + order
    fndr_pth: string
      Folder for finder charts
    dbx_pth: string
      Dropbox path for the finders
    skip_finder: False
      Skip making the finders

    Writes a file to disk that can be pasted into the Wiki
    """
    reload(x_finder)
    # Outfil
    if outfil is None:
        outfil = 'tmp_wiki.txt'
    f = open(outfil, 'w')

    # Finders?
    if not fndr_pth is None:
        if dbx_pth is None:
            dbx_pth = './'
            dbx_folder = './'
        else: # Expecting Public
            ifind = dbx_pth.find('Observing/')
            if ifind == -1:
                xdb.set_trace()
            else:
                dbx_folder = os.getenv('DROPBOX_DIR')+'/Public/'+dbx_pth[ifind:]
        #
        print('lick.wiki: Will copy finders to {:s}'.format(dbx_folder))
        # Get name tag
        name_tag = get_name_tag(targs.dtype.names)
        # Type
        #if isinstance(targs['RA'][0], basestring):
        #    radec = 1 # : separated strings
        #else:
        #    radec = 2 # decimal degrees
        # Finders
        fndr_files = []
        for targ in targs:
            # Finder
            #xdb.set_trace()
            if not skip_finder:
                x_finder.main([targ[name_tag], targ['RA'], targ['DEC']], fpath=fndr_pth)
            # Copy? + Save
            fil1 = fndr_pth+targ[name_tag]+'.pdf'
            fil2 = dbx_folder
            if not skip_finder:
                subprocess.call(["cp", fil1, dbx_folder])
            fndr_files.append(dbx_pth+targ[name_tag]+'.pdf')
        
    # Header
    lin = '||' 
    for key in keys:
        lin = lin+str(key)+'||'
    if 'fndr_files' in locals():
        lin=lin+'finder||comment||'
    f.write(str(lin+'\n'))
    
    # Targets
    for ii,targ in enumerate(targs):
        lin = '||' 
        for key in keys:
            lin = lin+str(targ[key])+'||'
        # Finder chart
        if 'fndr_files' in locals():
            lin = lin+'['+fndr_files[ii]+' pdf_finder]|| ||' # Lick formatting is different
        # Write
        f.write(str(lin+'\n'))

    # Close
    print('lick.wiki: Wrote {:s}'.format(outfil))
    f.close()

#### ###############################
def starlist(targs, outfil=None):
    '''
    Generates a Lick approved starlist
    FORMAT:
      1-25: Name 
      Then: RA, DEC, EPOCH [WITH COLONS!]

    Parameters:
    ----------
    targs: Table (targets with RA, DEC)
    outfil: string (None)
    '''
    reload(x_rad)
    # Init
    if outfil is None:
        outfil = 'starlist.txt'
    # Open
    f = open(outfil, 'w')

    # Name tag
    name_tag = get_name_tag(targs.dtype.names)

    # Loop
    for jj,targ in enumerate(targs):
        ras,decs = x_rad.dtos1((targ['RA'], targ['DEC']))
        #decs = targ['DEC'].replace(':',' ')
        if not decs[0] in ['+','-']:
            decs = '+'+decs
        # Name
        mask = []
        while type(targ[name_tag]) is MaskedConstant:
            mask.append(name_tag)
            name_tag = get_name_tag(targs.dtype.names,mask=mask)
        lin = targ[name_tag][0:3]+'_J'+ras.replace(':','')[0:4]+decs.replace(':','')[0:5]
        #xdb.set_trace()
        # RA
        lin = lin + '   ' + ras
        # DEC
        lin = lin + '  ' + decs
        # EPOCH
        lin = lin + '  2000.'
        # Write
        f.write(str(lin+'\n'))

    # Close
    print('lick.wiki: Wrote starlist -- {:s}'.format(outfil))
    f.close()

#
def get_name_tag(tbl_nms,mask=[None]):
    '''
    Find the name tag + return

    Parmaters:
    --------
    tbl_nms: List of names
    '''
    # Target name key
    name_tag = None
    for tag in ['Target', 'Name', 'NAME', 'QSO']:
        if tag in tbl_nms: 
            name_tag = tag
            if not name_tag in mask:
                break

    # Return
    return name_tag
