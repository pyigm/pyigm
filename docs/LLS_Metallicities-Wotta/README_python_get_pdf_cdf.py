##Example: How to get walkers or create a PDF and CDF

import numpy as np
from matplotlib import pyplot as plt


############
##Import CBW's module/class that gets MCMC output data
############
from MCMC_output import MCMC_output



################################################
################################################
################################################



def make_pdf_cdf_alltogether(infile=False,
    NHItype=None, NHIlow=None, NHIhigh=None,
    base_hist=None, quantity=None,
    limit_threshold=None
    ):
    
    """
    This is the all-in-one method. That is, it gets
    the data for a whole LIST of *.pkl/*.hd5 files.
    See make_pdf_cdf_onebyone for a one-by-one method.
    """
    
    
    ########################
    ##ONE WAY TO CALL THE CLASS AND SET OPTIONS
    ########################
    
    ##Call the class with defaults and save it to a variable
    mcmc=MCMC_output(
        NHItype = 'LLS',
        binsize = 0.20,
        base_hist=np.arange(-3.6, 1.0001, 0.20),
        smash = True,
        ndraw = 10000
        )
    
    
    # ########################
    # ##ANOTHER WAY TO CALL THE CLASS AND SET OPTIONS
    # ##THIS IS USEFUL IF YOU WANT TO CALL THESE MULTIPLE
    # ##  TIMES (e.g., ON MULTIPLE SAMPLES)
    # ########################
    
    # ##Call the class and save it to a variable
    # mcmc=MCMC_output()
    
    # ##Set some of its values/variables
    # # mcmc.NHIlow = '16.0'
    # # mcmc.NHIhigh = '17.2'
    # mcmc.NHItype = 'pLLS'
    # mcmc.binsize = 0.20
    # mcmc.base_hist=np.arange(-3.6, 1.0001, mcmc.binsize)
    # mcmc.smash = True
    # mcmc.ndraw = 10000
    
    
    ########################
    
    
    listfile_files = []
    with open('./pkl_files.list', 'r') as infile:
        for line in infile:
            listfile_files.append(line.strip())
    
    mcmc.infiles = listfile_files
    
    all_walkers, all_nsys, all_limit_code = mcmc.get_walkers()
    MDF, cumul_frac, nsys, base_hist, bins = mcmc.get_pdf_cdf()
    
    ########################
    # print(len(bins), len(MDF), len(cumul_frac))
    
    plt.plot(bins, MDF, drawstyle='steps-mid')
    plt.show()
    ##
    plt.plot(bins, cumul_frac)
    plt.show()
    
    print("Done!")
    
    # return MDF, cumul_frac, nsys, base_hist, bins
    return



################################################
################################################



def make_pdf_cdf_onebyone(infile=False,
    NHItype=None, NHIlow=None, NHIhigh=None,
    base_hist=None, quantity=None,
    NHI_defs_path=None,
    limit_threshold=None
    ):
    
    """
    This is the one-by-one method. That is, it gets
    the data for a single *.pkl/*.hd5 file, and then
    combines them after.
    See make_pdf_cdf_alltogether for an all-in-one method.
    """
    
    ##Call the class and save it to a variable
    mcmc=MCMC_output()
    
    ##Set some of its values/variables
    # mcmc.NHIlow = '16.0'
    # mcmc.NHIhigh = '17.2'
    mcmc.NHItype = 'pLLS'
    mcmc.NHI_defs_path = './'
    mcmc.binsize = 0.20
    mcmc.base_hist=np.arange(-3.6, 1.0001, mcmc.binsize)
    bins = mcmc.bins
    
    
    ########################
    
    
    ##Now loop over your list of files to [get and coadd]
    files = []
    with open('./pkl_files.list', 'r') as infile:
        for line in infile:
            files.append(line.strip())
    
    
    ##Now actually get the data
    walkers = []
    nsys = 0
    limit_code = []
    pdf = []
    cdf = []
    for file in files:
        ##Set the INDIVIDUAL (*.pkl/*.hdf5) file that
        ##  we will get the walkers from.
        mcmc.infiles = file
        
        ##Now get.
        this_walkers, this_nsys, this_limit_code = mcmc.get_walkers()
        
        ##Make sure we only add them IF this thing exists. Which includes
        ##  if it's a 'pLLS' (or whatever our NHI range is).
        if this_nsys > 0:
            walkers.append(this_walkers)
            nsys += this_nsys
            limit_code.append(this_limit_code)
            this_pdf, this_cdf, this_nsys, this_base_hist, this_bins = mcmc.get_pdf_cdf()
            
            ##We already know this_nsys > 0 from above.
            ##  No need to check again.
            pdf.append(this_pdf)
    
    
    ########################
    
    
    ##Convert to numpy arrays for easier math
    pdf = np.array(pdf)
    
    
    ##Smash all PDFs into a single one
    ##Normalize
    smashed_pdf = np.sum(pdf, axis=0)/nsys
    area = np.sum(smashed_pdf, axis=0)*mcmc.binsize
    MDF = smashed_pdf/area
    cumul_frac = np.cumsum(MDF)*mcmc.binsize
    
    
    ########################
    
    # print(len(bins),len(MDF),len(cumul_frac))
    
    plt.plot(bins, MDF)
    plt.show()
    ##
    plt.plot(bins, cumul_frac)
    plt.show()
    
    print("Done!")
    
    # return cumul_frac, MDF, nsys, mcmc.base_hist, mcmc.bins
    return




