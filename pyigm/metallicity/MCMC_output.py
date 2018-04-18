import numpy as np
import random

from six import string_types

##See also "README_python_get_pdf_cdf.py" for a longer example!!
##Note there are two methods: one where you send an individual file
##  and one where you send a LIST of files

##Example:
# import sys
# sys.path.append('../../../IDL/')
# from MCMC_output import MCMC_output
# mcmc = MCMC_output()
# mcmc.infiles = 'J0910+1014_242_34_z0.26412_emcee.pkl'
# mcmc.NHIlow = 15.0
# mcmc.quantity = 'met'
# walkers, nsys, limit_code = mcmc.get_output()
# print("{} {} {}".format(len(walkers), nsys, limit_code))



################################################
################################################
################################################


##You can call this function independently
##  (e.g., for use in plotting scripts) like:
##  >>> from MCMC_output import NHI_defs
##  >>> NHIdict = NHI_defs()


def NHI_defs():
    
    NHIdict = {}
    
    NHIdict['names']    = ['pLLS', 'LLS', 'SLLS', 'DLA']
    # NHIdict['names']    = ['pLLS', 'LLS', 'pLLS+LLS', 'SLLS', 'DLA']
    
    NHIdict['absorber'] = [15.0, 19.0]
    NHIdict['SLFS']     = [15.0, 16.2]
    NHIdict['pLLS']     = [16.2, 17.2]
    NHIdict['LLS']      = [17.2, 19.0]
    NHIdict['pLLS+LLS'] = [NHIdict['pLLS'][0], NHIdict['LLS'][1]]
    NHIdict['SLLS']     = [19.0, 20.3]
    NHIdict['DLA']      = [20.3, 22.0]
    
    return NHIdict



################################################
################################################
################################################



class MCMC_output(object):
    """
    Get data from MCMC output files (*.pkl or *.hdf5)
    
    Parameters - Inputs
    ----------
    infiles : str, or list of strings [str, str, ...]
        The filename or list of filenames (*.pkl or *.hd5) you want to read
    NHItype : str
        The NHI range by "type" ('pLLS', 'LLS', ...) as defined
        in NHI_defs() above, if you don't use NHIlow/NHIhigh
    NHIlow : float
        The NHI range low end, if you don't use NHItype
    NHIhigh : float
        The NHI range high end, if you don't use NHItype
    quantity : str
        By default, it return metallicities. You can also return
        any of the other quantities looped over by the MCMC (as of
        this writing, these are (col, red, met, dens, carbalpha)).
        Or, if you want it to return EVERYTHING, set self.quantity=None
    err_ci_detections : int/decimal, 0 to 100
        The confidence interval (C.I.) that should be returned for the
        errors on the median for detections.
        Defaults to 50 (interquartile range).
    err_ci_limits : int/decimal, 0 to 100
        The confidence interval (C.I.) that should be returned for the
        errors on the median for upper/lower limits.
        Defaults to 0.80.
    limit_threshold : float -- between 0 and 1
        The "threshold" that defines an upper/lower limit. If the
        most-extreme (left/right) bin is:
        >>> hist[0] >= limit_threshold*max(hist)
        --> upper limit, -1
        or
        >>> hist[-1] >= limit_threshold*max(hist)
        --> lower limit, -2
    smash : bool -- True or False
        Whether to smash the resulting PDFs (i.e., return
        a single PDF that combines all others, or return
        a PDF for each input file).
        Requires ndraw to be set.
    binsize : float
        The binsize for the histogram
    base_hist : np.array()
        The EDGES of the bins for the histogram. There will
        be (n_bins + 1) values in this numpy array. Can generate
        like: base_hist = np.arange(min, max+0.00001, binsize)
    ndraw : int
        Number of walkers to randomly draw, if you want to get the
        same number for each sightline (e.g., creating an MDF)
        The smash option requires this to be set.
    
    
    Outputs
    -----------------
    You might find them useful at some point, anyway.
    
    infiles_good : The list of files for which nsys=1 (not 0), indicating
                   they are within the desired NHI range
    info        : Output of get_info(); mcmcout['info']
    median      : Output of get_median(); median of mcmcout['tags']
    err_low     : Output of get_median(); lower bound of median
    err_high    : Output of get_median(); upper bound of median
    cdf         : Output of get_pdf_cdf(); cumulative distribution of
                  mcmcout['tags'] (may be smashed if smash==True)
    pdf         : Output of get_pdf_cdf(); probability distribution of
                  mcmcout['tags'] (may be smashed if smash==True)
    bins        : The centroids of the histogram bins, whereas base_hist are the edges
    walkers     : Output of get_output(); the mcmcout['tags'] walkers' positions
    limit_code  : Output of test_if_limit(); whether the PDF is an
                  upper limit (-1), lower limit (-2), or detection (0).
                  Only gets set if get_output() is called.
    nsys        : The number of systems within the defined NHI range
    nsys_mask   : An array of 0's and 1's for whether each input file
                  is within the defined NHI range (nsys=1) or not (nsys=0)
    
    
    Functions - To be called by user
    --------------------------
    load()          : Loads *.pkl or *.hd5 file. Only useful
                      if you're doing stuff manually.
                      Calls load_pickle() or load_h5py().
                      Returns: mcmcout (*.pkl file object or *.hd5 file object)
    get_output()    : Get all walkers in MCMC output for mcmcout['tags']
                      (may be randomly-drawn --- e.g., if smash==True).
                      Also creates medians, error bars on medians, etc.
                      This is the primary function.
                      Calls get_info(), get_walkers(), get_median(), and
                      get_pdf_cdf() --- in that order.
                      Returns: walkers, nsys, limit_code
    
    Functions - Behind-the-scenes, but MIGHT want to be called by user
    --------------------------
    get_info()      : Get the mcmcout['info'] dictionary,
                      including NHI, eNHI, z, errz, initial guesses, etc.
                      Calls get_info_pickle() or get_info_h5py().
                      Called by get_output(), but it's useful separately if you
                      don't want to spend time getting all of the walkers
                      first. Otherwise, just use get_output() and then get the
                      value of self.info
                      Returns: info, nsys
    get_walkers()   : Get all walkers in MCMC output for mcmcout['tags']
                      (may be randomly-drawn --- e.g., if smash==True).
                      Calls get_walkers_pickle() or get_walkers_h5py().
                      Returns: walkers, nsys, limit_code
    get_median()    : Get median value (and errors) of a single mcmcout['tags'],
                      including col, met, dens, red, carbalpha
                      Only useful if you want to smash the medians to get
                      an OVERALL median!! Otherwise, just use get_output() and then
                      get the value of self.median
                      Returns: median, err_low, err_high, nsys
    get_pdf_cdf()   : Get PDF and CDF of a single mcmcout['tags'],
                      including col, met, dens, red, carbalpha
                      Only useful if you want to smash the medians to get
                      an OVERALL median!! Otherwise, just use get_output() and then
                      get the value of self.pdf and self.cdf
                      NOT UPDATED for when self.quantity=None
                      Returns: pdf, cdf, nsys, base_hist, bins

    Functions - Behind-the-scenes, error checking/testing
    --------------------------
    ndraw_error()   : Checks if we NEED to randomly-draw walkers based
                      on smash==True (you might want to in other cases,
                      but this forces it for smash==True)
                      Returns: (Nothing)
    test_if_limit() : Tests if the PDF is an upper/lower limit
                      Returns: (Nothing)
    test_if_infiles_is_list() : Tests if self.infiles is a list (we have
                      to loop over) or not (a single string/file)
                      Returns: Bool (True or False)
    
    Functions - Behind-the-scenes, utilities to load files, etc.
    --------------------------
    load_pickle()   : Loads *.pkl file
                      Returns: mcmcout (*.pkl file object)
    load_h5py()     : Loads *.hd5 file
                      Returns: mcmcout (*.hd5 file object)
    get_info_pickle() : Gets mcmcout['info'] for *.pkl file
                      Returns: (Nothing)
    get_info_h5py() : Gets mcmcout['info'] for *.hd5 file
                      Returns: (Nothing)
    get_walkers_h5py() : Gets walkers for *.hd5 file
                      Returns: walkers, nsys
    get_walkers_pickle() : Gets walkers for *.pkl file
                      Returns: walkers, nsys
    
    """
    
    
    def __init__(self, infiles=False, NHItype=None, NHIlow=0.0, NHIhigh=99.0, quantity='met', err_ci_detections=50, err_ci_limits=80, limit_threshold=0.70, smash=False, binsize=0.20, base_hist=None, ndraw=None):
        
        ##When this class is FIRST called, these all get run.
        ##It basically sets defaults.
        ##You can overwrite these by hand, too.
        
        ############
        ##Things you may want to re-define each time
        ############
        self.infiles = infiles
        self.NHItype = NHItype
        self.NHIlow = float(NHIlow)
        self.NHIhigh = float(NHIhigh)
        self.quantity = quantity
        self.limit_threshold = limit_threshold
        ##err_ci is the Confidence Interval default
        self.err_ci_detections = float(err_ci_detections)
        self.err_ci_limits = float(err_ci_limits)
        
        self.smash = smash
        if (smash) and (ndraw is None):
            self.ndraw_error()
        else:
            self.ndraw = ndraw
        
        self.binsize=binsize
        if base_hist is not None:
            self.base_hist=base_hist
        else:
            # self.base_hist=np.arange(-3.6, 1.0001, self.binsize)
            self.base_hist=np.arange(-5.0, 2.5001, self.binsize)
        
        
        ############
        ##The things below will be populated later, as outputs
        ############
        self.infiles_good = []
        self.walkers = np.array([])
        self.info = []
        self.nsys = 0
        self.nsys_mask = []
        self.limit_code = 999
        
        ##For histogram
        self.bins=self.base_hist[0:-1] + (self.binsize/2.0)
        self.pdf = []
        self.cdf = []
        
        ##For median
        self.median = -999
        self.err_low = -999
        self.err_high = -999
    
    
    
    
    ################################################
    ################################################
    ################################################
    ################################################
    
    
    
    
    def inputs(self):
        print("Inputs:")
        print("    self.infiles            (str or list of str)")
        print("    self.NHItype            (str, like 'pLLS')")
        print("    self.NHIlow             (float, like 17.2)")
        print("    self.NHIhigh            (float, like 19.0)")
        print("    self.quantity           (str, like 'met')")
        print("    self.err_ci_detections  (float/int, 0 to 100)")
        print("    self.err_ci_limits      (float/int, 0 to 100)")
        print("    self.limit_threshold    (float, 0 to 1)")
        print("    self.smash              (bool)")
        print("    self.binsize            (float)")
        print("    self.base_hist          (np.array)")
        print("    self.ndraw              (int)")
        
        return
    
        
    ################################################
    ################################################
    
    
    def outputs(self):
        print("Outputs:")
        print("    self.infiles_good  (list of str)")
        print("    self.info          (list of dict)")
        print("    self.walkers       (list of np.array())")
        print("    self.pdf           (list of np.array())")
        print("    self.cdf           (list of np.array())")
        print("    self.bins          (np.array())")
        print("    self.median        (list of (lists of) int)")
        print("    self.err_low       (list of float)")
        print("    self.err_high      (list of float)")
        print("    self.nsys          (int)")
        print("    self.nsys_mask     (list of int: 0, 1)")
        print("    self.limit_code    (list of int: 0, -1, -2)")
        
        return
    
        
    ################################################
    ################################################
    
    
    def get_output(self):
        self.get_info()
        self.get_walkers()
        self.get_median()
        self.get_pdf_cdf()
        return self.walkers, self.nsys, self.limit_code
    
    
    ################################################
    ################################################
    
    
    def load(self):
        """
        Loads a *.pkl or *.hd5 file. Only useful if
        doing something manually.
        """
        
        ##Check if the "infiles" is a single string (otherwise,
        ##  it's a list we want to loop over)
        infiles_is_list = self.test_if_infiles_is_list()
        
        all_mcmcout = []
        if infiles_is_list is True:
            ##Then we have a LIST of input files to loop over and get
            
            old_infiles = self.infiles
            
            for fil in old_infiles:
                self.infiles = fil
                all_mcmcout.append(self.load())
            
            ############
            ##Done looping over files
            ############
            
            mcmcout = all_mcmcout[:]
            
            ##Reset the infiles
            self.infiles = old_infiles
            
        else:
            ##Get file extension for automatically detecting file format (*.pkl or *.hd5)
            fileext = self.infiles.split('.')[-1]
            
            if fileext == 'pkl':
                mcmcout = self.load_pickle()
            else:
                mcmcout = self.load_h5py()
        
        
        return mcmcout
    
    
    ################################################
    ################################################
    
    
    def get_info(self):
        
        ##Check to see if user gave (NHItype) or (NHIlow and NHIhigh).
        ##i.e., if NHItype was set, then that's what we'll go with.
        if self.NHItype:
            NHIdict = NHI_defs()
            self.NHIlow = float(NHIdict[self.NHItype][0])
            self.NHIhigh = float(NHIdict[self.NHItype][1])
        
        ##Check if the "infiles" is a single string (otherwise,
        ##  it's a list we want to loop over)
        infiles_is_list = self.test_if_infiles_is_list()
        
        if infiles_is_list is True:
            ##Then we have a LIST of input files to loop over and get
            
            ##Check if we need to reset ndraw.
            ##  That is, if (smash is set) and (ndraw is not set), we have a problem.
            if self.smash:
                # self.smash_error()
                if self.ndraw is None:
                    self.ndraw_error()
            
            ##Save the old infiles so we can get individual walkers
            ##  by overwriting self.infiles temporarily
            old_infiles = self.infiles
            
            all_files_good = []
            all_info = []
            all_nsys = 0
            all_nsys_mask = []
            for fil in old_infiles:
                self.infiles = fil
                
                ##Get file extension for automatically detecting file format (*.pkl or *.hd5)
                fileext = self.infiles.split('.')[-1]
                
                if fileext == 'pkl':
                    self.get_info_pickle()
                else:
                    self.get_info_h5py()
                
                
                ##Check to make sure we actually want this sightline
                ##  i.e., that it is within our NHI range
                if self.nsys > 0:
                    all_files_good.append(self.infiles)
                    all_info.append(self.info)
                    all_nsys += self.nsys
                    all_nsys_mask.append(1)
                else:
                    all_nsys_mask.append(0)
            
            ############
            ##Done looping over files
            ############
            
            self.infiles_good = all_files_good
            self.nsys = all_nsys
            self.nsys_mask = all_nsys_mask
            self.info = all_info
            
            ##Reset the infiles
            self.infiles = old_infiles
        
        else:
            ##Then we have a single file to get
            
            ##Get file extension for automatically detecting file format (*.pkl or *.hd5)
            fileext = self.infiles.split('.')[-1]
            
            if fileext == 'pkl':
                self.get_info_pickle()
            else:
                self.get_info_h5py()
            
            if self.nsys > 0:
                self.infiles_good = self.infiles
                self.nsys_mask = [1]
            else:
                self.infiles_good = []
                self.info = []
                self.nsys_mask = [0]
        
        
        
        ########################
        
        return self.info, self.nsys
    
    
    ################################################
    ################################################
    
    
    def get_walkers(self):
        
        ##Check to see if user gave (NHItype) or (NHIlow and NHIhigh).
        ##i.e., if NHItype was set, then that's what we'll go with.
        if self.NHItype:
            NHIdict = NHI_defs()
            self.NHIlow = float(NHIdict[self.NHItype][0])
            self.NHIhigh = float(NHIdict[self.NHItype][1])
        
        ##Check if the "infiles" is a single string (otherwise,
        ##  it's a list we want to loop over)
        infiles_is_list = self.test_if_infiles_is_list()
        
        if infiles_is_list is True:
            ##Then we have a LIST of input files to loop over and get
            
            ##Save the old infiles so we can get individual walkers
            ##  by overwriting self.infiles temporarily
            ##  i.e.,
            ##
            ##>>> WE'RE USING RECURSION HERE <<<
            ##
            old_infiles = self.infiles
            
            all_files_good = []
            all_walkers = []
            all_nsys = 0
            all_nsys_mask = []
            all_limit_code = []
            for fil in old_infiles:
                self.infiles = fil
                self.get_walkers()
                if self.nsys > 0:
                    all_files_good.append(self.infiles)
                    all_walkers.append(self.walkers)
                    all_nsys += self.nsys
                    all_nsys_mask.append(1)
                    all_limit_code.append(self.limit_code)
                
                else:
                    all_nsys_mask.append(0)
            
            ############
            ##Done looping over files
            ############
            
            ##Reset the infiles
            self.infiles = old_infiles
            ##Save the output
            self.infiles_good = all_files_good
            self.walkers = all_walkers
            self.nsys = all_nsys
            self.nsys_mask = all_nsys_mask
            self.limit_code = all_limit_code
        
        else:
            ########################
            ##Just a single file
            ########################
            
            ##This NEEDS to be reset each time
            self.walkers = np.array([])
            self.nsys = 0
            self.nsys_mask = []
            self.limit_code = 999
            
            ########################
            
            ##Get file extension for automatically detecting file format (*.pkl or *.hd5)
            fileext = self.infiles.split('.')[-1]
            
            if fileext == 'pkl':
                walkers, nsys = self.get_walkers_pickle()
            else:
                walkers, nsys = self.get_walkers_h5py()
            
            
            ########################
            
            if nsys > 0:
                ##Check if we need to reset ndraw.
                ##  That is, if (smash is set) and (ndraw is not set), we have a problem.
                if (self.smash) and (self.ndraw is None):
                    self.ndraw_error()
                
                ##Randomly draw "ndraw" walkers rather than
                ##  returning all of them, if user desires
                if self.ndraw is None:
                    self.walkers = walkers
                else:
                    
                    ##Initialize an empty list
                    random_walkers = []
                    
                    ##Now randomly get them
                    for idraw in range(self.ndraw):
                        random_walkers.append(random.choice(walkers))
                    
                    self.walkers = random_walkers
                
                
                self.infiles_good = self.infiles
                self.nsys = nsys
                self.nsys_mask.append(nsys)
            
            else:
                ##These will be [], [], and 0, respectively
                self.infiles_good = []
                self.walkers = walkers
                self.nsys = nsys
                self.nsys_mask.append(nsys)
            
            ########################
            
            ##Get limit code --- make sure this is done AFTER setting
            ##  self.nsys (and self.walkers)
            self.test_if_limit()
        
        
        ########################
        
        return self.walkers, self.nsys, self.limit_code
    
    
    ################################################
    ################################################
    
    
    def get_median(self):
        ##Check to see if get_output() has been called yet
        if len(self.walkers) == 0:
            print("You need to call get_walkers() first! Returning...")
            return [[], [], [], []]
        
        ##Check if the "infiles" is a single string (otherwise,
        ##  it's a list we want to loop over)
        infiles_is_list = self.test_if_infiles_is_list()
        
        if infiles_is_list is True:
            ##Then we have a LIST of input files to loop over and get
            
            ##Check if we need to reset ndraw.
            ##  That is, if (smash is set) and (ndraw is not set), we have a problem.
            if self.smash:
                # self.smash_error()
                if self.ndraw is None:
                    self.ndraw_error()
            
            all_medians = []
            all_err_low = []
            all_err_high = []
            for idx in range(len(self.walkers)):
                # ##Check to make sure we actually want this sightline
                # ##  i.e., that it is within our NHI range
                # if self.nsys[idx] > 0:
                
                ##ASSUME IT IS WITHIN OUR NHI RANGE!! COULD BE DANGEROUS!!
                ##Check if this is a detection or a limit
                if self.limit_code[idx] == 0:
                    ci_low = 50 - (self.err_ci_detections/2.0)
                    ci_high = 50 + (self.err_ci_detections/2.0)
                else:
                    ci_low = 50 - (self.err_ci_limits/2.0)
                    ci_high = 50 + (self.err_ci_limits/2.0)
                
                temp_all_walkers = np.array(self.walkers[idx])
                ##Need to check if we only have ONE walker axis, or if we
                ##  have EVERYTHING (in which case we get a median for each axis)
                if self.quantity:
                    all_medians.append(np.median(temp_all_walkers))
                    all_err_low.append(np.percentile(temp_all_walkers,ci_low,axis=0))
                    all_err_high.append(np.percentile(temp_all_walkers,ci_high,axis=0))
                else:
                    all_medians.append([np.median(temp_all_walkers[:,i]) for i in range(len(temp_all_walkers[0]))])
                    all_err_low.append([np.percentile(temp_all_walkers[:,i],ci_low,axis=0) for i in range(len(temp_all_walkers[0]))])
                    all_err_high.append([np.percentile(temp_all_walkers[:,i],ci_high,axis=0) for i in range(len(temp_all_walkers[0]))])
            
            ############
            ##Done looping over files
            ############
            
            ##Check if we want to return a single median of ALL walkers,
            ##  or a list of medians, one for each sightline.
            if self.smash:
                ##Assume detections
                ci_low = 50 - (self.err_ci_detections/2.0)
                ci_high = 50 + (self.err_ci_detections/2.0)
                
                temp_all_walkers = np.array(self.walkers)
                if self.quantity:
                    self.median = np.array(np.median(temp_all_walkers))
                    self.err_low = np.array(np.percentile(temp_all_walkers,ci_low,axis=0))
                    self.err_high = np.array(np.percentile(temp_all_walkers,ci_high,axis=0))
                else:
                    ##Get the median for EACH quantity
                    self.median = [np.median(temp_all_walkers[:,:,i].flatten()) for i in range(len(self.walkers[0][0]))]
                    self.err_low = [np.percentile(temp_all_walkers[:,:,i].flatten(),ci_low,axis=0) for i in range(len(self.walkers[0][0]))]
                    self.err_high = [np.percentile(temp_all_walkers[:,:,i].flatten(),ci_high,axis=0) for i in range(len(self.walkers[0][0]))]
            else:
                self.median = all_medians[:]
                self.err_low = all_err_low[:]
                self.err_high = all_err_high[:]
        
        else:
            ##Then we have a single file to get
            if self.nsys > 0:
                ##Check if this is a detection or a limit
                if self.limit_code == 0:
                    ci_low = 50 - (self.err_ci_detections/2.0)
                    ci_high = 50 + (self.err_ci_detections/2.0)
                else:
                    ci_low = 50 - (self.err_ci_limits/2.0)
                    ci_high = 50 + (self.err_ci_limits/2.0)
                
                temp_all_walkers = np.array(self.walkers)
                ##Need to check if we only have ONE walker axis, or if we
                ##  have EVERYTHING (in which case we get a median for each axis)
                if self.quantity:
                    self.median = np.median(temp_all_walkers)
                    self.err_low = np.percentile(temp_all_walkers,ci_low,axis=0)
                    self.err_high = np.percentile(temp_all_walkers,ci_high,axis=0)
                else:
                    self.median = [np.median(temp_all_walkers[:,i]) for i in range(len(temp_all_walkers[0]))]
                    self.err_low = [np.percentile(temp_all_walkers[:,i],ci_low,axis=0) for i in range(len(temp_all_walkers[0]))]
                    self.err_high = [np.percentile(temp_all_walkers[:,i],ci_high,axis=0) for i in range(len(temp_all_walkers[0]))]
                
                self.nsys_mask = [1]
            else:
                self.median = -999
                self.err_low = -999
                self.err_high = -999
                self.nsys_mask = [0]
        
        
        
        ########################
        
        return self.median, self.err_low, self.err_high, self.nsys
    
    
    ################################################
    ################################################
    
    
    def get_pdf_cdf(self):
        
        ##Check to see if get_output() has been called yet
        if len(self.walkers) == 0:
            print("You need to call get_walkers() first! Returning...")
            return [[], [], [], [], []]
        
        ##Recalculate these in case the user changed base_hist.
        ##  They're not used below, only outputted for convenience.
        self.bins = self.base_hist[0:-1] + (self.binsize/2.0)
        
        
        ##If we don't check now and self.quantity == None, then the
        ##  next check -- if it's a list -- will always pass because
        ##  there will of COURSE be a list of walkers. But, that does NOT
        ##  mean there are NOT multiple input files, so it will get confused.
        ##  We COULD adjust the code below to handle this, but... that's difficult.
        if not self.quantity:
            try:
                from get_bash_color import bash_color
                errtxt = bash_color('ERROR:', 'red', bold=True)
            except:
                errtxt = "ERROR:"
            
            print("\n{} get_pdf_cdf() is not set up to create PDFs/CDFs for EVERY quantity at once!! Returning...\n".format(errtxt))
            return [[], [], [], [], []]
        
        
        ##Check if the "walkers" is a single list or
        ##  a 2D list (of numpy arrays)
        infiles_is_list = False
        if not isinstance(self.walkers[0], float):
            infiles_is_list = True
        
        
        if infiles_is_list is True:
            ##Then we have a LIST of input files to loop over and get
            
            ##Check if we need to reset ndraw.
            ##  That is, if (smash is set) and (ndraw is not set), we have a problem.
            if self.smash:
                # self.smash_error()
                if self.ndraw is None:
                    self.ndraw_error()
            
            all_pdf = []
            all_walkers = []
            all_nsys = 0
            for idx in range(len(self.walkers)):
                # ##Check to make sure we actually want this sightline
                # ##  i.e., that it is within our NHI range
                # if self.nsys[idx] > 0:
                
                ##ASSUME IT IS WITHIN OUR NHI RANGE!! COULD BE DANGEROUS!!
                all_nsys += 1
                all_walkers.append(self.walkers[idx])
                pdf_temp, edge_temp = np.histogram(np.sort(self.walkers[idx]), bins=self.base_hist, density=True)
                
                ##Normalize to the number of walkers
                all_pdf.append(pdf_temp/len(self.walkers[idx]))
            
            ############
            ##Done looping over files
            ############
            
            self.nsys = all_nsys
            
            ##Normalize the area
            smashed_pdf = np.sum(all_pdf, axis=0)/all_nsys
            area = np.sum(smashed_pdf, axis=0)*self.binsize
            ##
            smashed_pdf = smashed_pdf/area
            all_pdf = np.array([ap/area for ap in all_pdf])
            
            ##Check if we want to return a single (smashed) PDF/CDF or individual ones
            if self.smash:
                ##Normalize
                self.pdf = smashed_pdf
                self.cdf = np.cumsum(smashed_pdf)*self.binsize
            else:
                self.pdf = all_pdf
                # self.cdf = np.array([np.cumsum(my_pdf) for my_pdf in self.pdf])
                ##Actually, let's ALWAYS assume we want the OVERALL CDF (i.e., smashed)!!
                self.cdf = np.cumsum(smashed_pdf)*self.binsize
        
        else:
            ##Then we have a single file to get
            pdf_temp, edge_temp = np.histogram(np.sort(self.walkers), bins=self.base_hist, density=True)
            
            pdf_temp = pdf_temp/len(self.walkers)
            
            # ##Normalize
            # smashed_pdf = np.sum(pdf_temp, axis=0)/self.nsys
            # area = np.sum(smashed_pdf, axis=0)*self.binsize
            # ##
            # smashed_pdf = smashed_pdf/area
            # pdf_temp = pdf_temp/area
            
            ##Normalize to the number of walkers
            self.pdf = pdf_temp
            self.cdf = np.cumsum(self.pdf)*self.binsize
            self.nsys = 1
            self.nsys_mask = [self.nsys]
        
        
        ########################
        
        return self.pdf, self.cdf, self.nsys, self.base_hist, self.bins
    
    
    
    
    ################################################
    ################################################
    ################################################
    ################################################
    
    
    
    
    # def smash_error(self):
    #     
    #     print("(Remember because you set smash=True, you'll want to call get_pdf_cdf() and get_median() separately, not just get_output()...)")
    #     
    #     return
    
    
    ################################################
    ################################################
    
    


    def ndraw_error(self):
        
        self.ndraw = 10000
        
        print("Smashing requires the same number of walkers from each sightline. Setting ndraw={}".format(self.ndraw))
        
        return
    
    
    ################################################
    ################################################
    
    
    def test_if_limit(self):
        
        ##limit_threshold:
        ##If the height at the extreme edge is, e.g., 70%
        ##  the max height of the histogram, it is probably
        ##  an upper or lower limit
        
        if self.nsys <= 0:
            self.limit_code = 999
        
        ##Need to check if we only want ONE walker axis (e.g., self.quantity == 'met'),
        ##  or if we want to get EVERYTHING (i.e., self.quantity == None)
        if self.quantity:
            ##Then we have a single axis to get
            hist, edges = np.histogram(self.walkers, bins=30)
            max_height = hist.max()
            
            if hist[0] >= max_height*self.limit_threshold:
                ##Upper limit
                self.limit_code = -1
            elif hist[-1] >= max_height*self.limit_threshold:
                ##Lower limit
                self.limit_code = -2
            else:
                self.limit_code = 0
        else:
            temp_limit_code = []
            for idx in range(len(self.walkers[0])):
                ##Get all the walkers for THIS axis
                hist, edges = np.histogram(np.array(self.walkers)[:,idx], bins=30)
                max_height = hist.max()
                
                if hist[0] >= max_height*self.limit_threshold:
                    ##Upper limit
                    temp_limit_code.append(-1)
                elif hist[-1] >= max_height*self.limit_threshold:
                    ##Lower limit
                    temp_limit_code.append(-2)
                else:
                    temp_limit_code.append(0)
                
            self.limit_code = temp_limit_code
        
        # return self.limit_code
        return
    
    
    
    
    ################################################
    ################################################
    
    
    def test_if_infiles_is_list(self):
        
        ##Check if the "infiles" is a single string (otherwise,
        ##  it's a list we want to loop over)
        infiles_is_list = False
        try:
            ##Python2
            if not isinstance(self.infiles, basestring):
                infiles_is_list = True
        except:
            ##Python3
            if not isinstance(self.infiles, string_types):
                infiles_is_list = True
        
        return infiles_is_list
    
    
    
    
    ################################################
    ################################################
    ################################################
    ################################################
    
    
    
    
    def load_pickle(self):
        import pickle
        try:
            ##Python2
            fil=open(self.infiles)
            mcmcout=pickle.load(fil)
        except:
            ##Python3
            fil=open(self.infiles,'rb')
            mcmcout=pickle.load(fil, encoding='latin1')
        
        fil.close()
        
        return mcmcout
    
    
    ################################################
    ################################################
    
    
    def load_h5py(self):
        import h5py
        mcmcout=h5py.File(self.infiles, 'r')
        
        return mcmcout
    
    
    ################################################
    ################################################
    
    
    def get_info_pickle(self):
        mcmcout = self.load_pickle()
        
        ##Check that it's in the NHI range
        mcmctags = np.array(mcmcout['tags'])
        nhiloc = np.where(mcmctags == 'col')[0][0]
        if float(self.NHIlow) <= mcmcout['guess'][nhiloc] < float(self.NHIhigh):
            self.nsys = 1
            self.info = mcmcout['info'].copy()
        else:
            self.nsys = 0
            self.info = []
        
        return
                    
    
    ################################################
    ################################################
    
    
    def get_info_h5py(self):
        # ##We don't save the info in the *.hd5 files...
        # ##  ...so we have to look to the *.pkl files for this
        # old_infiles = self.infiles
        # self.infiles = old_infiles.replace('.hd5', '.pkl')
        # 
        # self.get_info_pickle()
        # 
        # ##Return back to *.hd5 file
        # self.infiles = old_infiles
        
        mcmcout = self.load_h5py()
        
        ##Get the keys for the attributes dictionary
        hd5_attrs = [str(at) for at in mcmcout['inputs'].attrs]
        
        ##Build a temp info dictionary
        this_info = {}
        for ha in hd5_attrs:
            this_info[ha] = mcmcout['inputs'].attrs[ha]
        
        if float(self.NHIlow) <= this_info['NHI'] < float(self.NHIhigh):
            self.nsys = 1
            self.info = this_info.copy()
        else:
            self.nsys = 0
            self.info = []

        
        return
    
    
    ################################################
    ################################################
    
    
    def get_walkers_pickle(self):
        
        mcmcout = self.load_pickle()
        
        ##Check that it's in the NHI range
        mcmctags = np.array(mcmcout['tags'])
        nhiloc = np.where(mcmctags == 'col')[0][0]
        if float(self.NHIlow) <= mcmcout['guess'][nhiloc] < float(self.NHIhigh):
            ##Need to check if we only want ONE walker axis (e.g., self.quantity == 'met'),
            ##  or if we want to get EVERYTHING (i.e., self.quantity == None)
            if self.quantity:
                qnt= np.where(mcmctags == self.quantity)[0]
                ##Check if the axis even exists for this absorber. If it does not,
                ##  return an empty array
                if len(qnt) == 0:
                    ##Print an error to the screen
                    errtxt = "Warning:"
                    try:
                        from get_bash_color import bash_color
                        errtxt = bash_color(errtxt, 'red', bold=True)
                    except:
                        pass
                    
                    print("{} 'self.quantity={}' not found in {}".format(errtxt, self.quantity, self.infiles))
                    walkers = []
                    nsys = 0
                else:
                    qnt = qnt[0]
                    ##We need the pdf to be a numpy array before we try to do [:,qnt] on it
                    ##  it should already be, but let's be sure if that.
                    walkers = np.array(np.array(mcmcout['pdfs'])[:,qnt])
                    nsys = 1
            else:
                walkers = np.array(np.array(mcmcout['pdfs'])[:,:])
                nsys = 1
        else:
            walkers = []
            nsys = 0
        
        return walkers, nsys
    
    
    ################################################
    ################################################
    
    
    def get_walkers_h5py(self):
        
        ## Load
        mcmcout = self.load_h5py()
        
        ##Check that it's in the NHI range
        try:
            mcmctags = np.array([str(i,'utf-8') for i in mcmcout['outputs']['tags'].value])
        except:
            mcmctags = np.array([str(i) for i in mcmcout['outputs']['tags'].value])
        
        if float(self.NHIlow) <= mcmcout['inputs']['guess'][0] < float(self.NHIhigh):
            ##Need to check if we only want ONE walker axis (e.g., self.quantity == 'met'),
            ##  or if we want to get EVERYTHING (i.e., self.quantity == None)
            if self.quantity:
                qnt= np.where(mcmctags == self.quantity)[0]
                ##Check if the axis even exists for this absorber. If it does not,
                ##  return an empty array
                if len(qnt) == 0:
                    ##Print an error to the screen
                    errtxt = "Warning:"
                    try:
                        from get_bash_color import bash_color
                        errtxt = bash_color(errtxt, 'red', bold=True)
                    except:
                        pass
                    
                    print("{} 'self.quantity={}' not found in {}".format(errtxt, self.quantity, self.infiles))
                    walkers = []
                    nsys = 0
                else:
                    qnt = qnt[0]
                    ##We need the pdf to be a numpy array before we try to do [:,qnt] on it
                    ##  it should already be, but let's be sure if that.
                    walkers = np.array(mcmcout['outputs']['pdfs'][:,qnt])
                    nsys = 1
            else:
                walkers = np.array(mcmcout['outputs']['pdfs'][:,:])
                nsys = 1
        else:
            walkers = []
            nsys = 0
        
        return walkers, nsys



