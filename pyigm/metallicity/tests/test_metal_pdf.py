# Module to run tests on FNModel

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
from astropy import units as u

from ..pdf import MetallicityPDF


#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_init_pdf():
    # Dummy data
    ZH  = np.linspace(-5, 0., 25)
    pdf = np.exp(-(ZH+1.5)**2/0.2**2)
    # Class
    mpdf = MetallicityPDF(ZH, pdf)
    # Test
    np.testing.assert_allclose(mpdf.meanZH, -1.4998713559597918)

def test_stats():
    # Dummy data
    ZH  = np.linspace(-5, 0., 25)
    pdf = np.exp(-(ZH+1.5)**2/0.2**2)
    # Class
    mpdf = MetallicityPDF(ZH, pdf)
    cl_lim = mpdf.confidence_limits(0.68)

    np.testing.assert_allclose(cl_lim, (-1.7738954485146285, -1.4708227781158598))
