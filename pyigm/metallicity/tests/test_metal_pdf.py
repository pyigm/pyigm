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


def test_stats():
    # Dummy data
    ZH  = np.linspace(-5, 0., 25)
    pdf = np.exp(-(ZH+1.5)**2/0.2**2)
    # Class
    mpdf = MetallicityPDF(ZH, pdf)
    # Test mean
    np.testing.assert_allclose(mpdf.meanZH, -1.4998713559597918)
    #
    np.testing.assert_allclose(mpdf.medianZH, -1.5967047485449448)


def test_cl():
    # Dummy data
    ZH  = np.linspace(-5, 0., 25)
    pdf = np.exp(-(ZH+1.5)**2/0.2**2)
    # Class
    mpdf = MetallicityPDF(ZH, pdf)
    cl_lim = mpdf.confidence_limits(0.68)

    np.testing.assert_allclose(cl_lim, (-1.7738954485146285, -1.4708227781158598))


def test_add():
    # Dummy data
    ZH  = np.linspace(-5, 0., 25)
    pdf = np.exp(-(ZH+1.5)**2/0.2**2)
    mpdf = MetallicityPDF(ZH, pdf)
    pdf2 = np.exp(-(ZH+1.7)**2/0.3**2)
    mpdf2 = MetallicityPDF(ZH, pdf2)
    # Sum
    sum_pdf = mpdf + mpdf2
    sum_pdf.normalize()
    # Test
    np.testing.assert_allclose(sum_pdf.meanZH,-1.5999356764972328)
