# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
import pytest
import datetime
import numpy as np
import numpy.testing as nptest
from geospacepy.terrestrial_ellipsoidal import ecef_cart2geodetic

#Test on Example 3-3 of Vallado
@pytest.mark.parametrize('n_vecs',[1,3])
def test_ecef_cart2geodetic_on_vallado_example(n_vecs):
    R_ECEF_single = np.array([65248340,68628750,64482960]).reshape(1,-1)
    R_ECEF = np.repeat(R_ECEF_single,n_vecs,axis=0)

    gdlat_expected = 34.352496
    gclon_expected = 46.4464
    h_ellp_expected = 50852200
    geodetic_iteration_tol = 1e-7

    gdlats,gclons,h_ellps = ecef_cart2geodetic(R_ECEF,tol=geodetic_iteration_tol)

    test_tol_gdlat_degrees = 1e-6
    test_tol_gclon_degrees = 1e-4 #only 4 digits given in textbook
    test_tol_meters = 15

    nptest.assert_allclose(gdlat_expected,gdlats,rtol=0.,atol=test_tol_gdlat_degrees)
    nptest.assert_allclose(gclon_expected,gclons,rtol=0.,atol=test_tol_gclon_degrees)
    nptest.assert_allclose(h_ellp_expected,h_ellps,rtol=0.,atol=test_tol_meters)
    
