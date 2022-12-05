""" Simple script to generate EW splines"""

from pyigm import utils as pyigm_utils
from pyigm.fN import tau_eff
from pyigm.fN import fnmodel

from IPython import embed

if __name__ == "__main__":
    # b=24
    #pyigm_utils.mk_ew_lyman_spline(
    #    24., ew_fil='EW_SPLINE_b24_2022.yml', chk=False)

    # Test
    fN_default = fnmodel.FNModel.default_model()
    zem = 3.
    ilambda = 1200.*(1+2.9)
    teff = tau_eff.lyman_ew(ilambda, zem, fN_default)
