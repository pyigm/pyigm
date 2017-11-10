# Understanding MCMC output (*.pkl format files)

    import pickle
    filename = 'J1619+3342_z0.2694_emcee.pkl'


## Read in the file
Python2
    file = open(filename)
    mcmcout = pickle.load(file)
    file.close()
    
Python3
    file = open(filename,'rb')
    mcmcout = pickle.load(file,encoding='latin1')
    file.close()


## Get all possible tags of `mcmcout`

    mcmcout.keys()


# INPUTS

## Sightline info

Get info about the sightline (name, redshift, etc.) and initial guesses

    mcmcout['info']

Get observational data (metal ion column densities)

    mcmcout['data']


## Input parameter tag order

Get input parameter tag order

    mcmcout['tags']

Get initial guesses on the input parameters

    mcmcout['guess']



# OUTPUTS

## Just the answers

Get the list of "which column" means "what percent"

    mcmcout['percent']

The actual results

    mcmcout['results']

The order is:

    mcmcout['results'][mcmcout_percent_index][mcmcout_tags_index]

E.g. for `percent=75` and for `tags='met'` then

    mcmcout['results'][4][2]


## Residuals

The order of the ions (the observational input)

    mcmcout['ions']

The ion column densities (logN) of the single, best-fit model. The list order is: `mcmcout['ions']`

    mcmcout['best_fit']


## The details

The grid position of every walker. Each value is a list of input parameters whose order is: `mcmcout['tags']`

    mcmcout['pdfs']

The "acceptance fraction" of the MCMC walkers (higher is better). Should be ~50%.

    mcmcout['acceptance']

The "effective logNHI" that the MCMC found (rather than the input logNHI). Mostly, as long as this is roughly the same as the input value, you can probably just use the input logNHI.

    mcmcout['effNHI']

