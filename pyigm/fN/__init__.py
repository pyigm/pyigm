# Dependent modules
try:
    import pymc
except ImportError:
    print('-----------------------------------------------------------')
    print('-----------------------------------------------------------')
    print('WARNING: Not loading mcmc in pyigm.fN   \n Install pymc if you want it')
    print('-----------------------------------------------------------')
else:
    import mcmc