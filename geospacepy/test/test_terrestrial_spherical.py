# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
import pytest
import datetime
import numpy as np
import numpy.testing as nptest
from geospacepy.special_datetime import datetime2jd
from geospacepy.terrestrial_spherical import (eci2ecef,ecef2eci,
                                            ecef_cart2spherical,ecef_spherical2cart,
                                            ecef2enu,enu2ecef)

example_jd = datetime2jd(datetime.datetime(2010,5,29,9,13,30))

def _example_vector(x,y,z,n_vectors):
    """Makes a shape=(n_vector,3) array out of floats x,y,z"""
    R = np.array([x,y,z]).reshape(1,3)    
    return np.broadcast_to(R,(n_vectors,3))

@pytest.mark.parametrize('n_vectors,jds',[(1,example_jd),
                                          (4,example_jd),
                                          (4,np.full((4,),example_jd)),
                                          (4,np.full((1,4),example_jd))])
def test_eci_to_ecef_round_trip(n_vectors,jds):
    R_ECI_in = _example_vector(1.,0.,0.,n_vectors)
    R_ECEF = eci2ecef(R_ECI_in,jds)
    R_ECI_out = ecef2eci(R_ECEF,jds)
    nptest.assert_allclose(R_ECI_in,R_ECI_out,atol=1e-8,rtol=0.)

@pytest.mark.parametrize('n_vectors',[1,4])
def test_ecef_cart_spherical_round_trip(n_vectors):
    R_ECEF_in = _example_vector(.5,.5,.5,n_vectors)
    lats,lons,rs = ecef_cart2spherical(R_ECEF_in)
    R_ECEF_out = ecef_spherical2cart(lats,lons,rs)
    nptest.assert_allclose(R_ECEF_in,R_ECEF_out,atol=1e-8,rtol=0.)

@pytest.mark.parametrize('n_vectors,latlon_shape',[(1,1),
                                                   (4,1),
                                                   (4,(4,)),
                                                   (4,(1,4))])
def test_ecef_enu_round_trip(n_vectors,latlon_shape):
    lat,lon = 40.0150,-105.2705 #Boulder, Colorado, USA
    if latlon_shape == 1:
      lats,lons = lat,lon
    else:
      lats,lons = np.ones(latlon_shape)*lat,np.ones(latlon_shape)*lon

    R_ECEF_in = _example_vector(0.,1.,0.,n_vectors)
    R_ENU = ecef2enu(R_ECEF_in,lats,lons)
    R_ECEF_out = enu2ecef(R_ENU,lats,lons)
    nptest.assert_allclose(R_ECEF_in,R_ECEF_out,atol=1e-8,rtol=0.)

