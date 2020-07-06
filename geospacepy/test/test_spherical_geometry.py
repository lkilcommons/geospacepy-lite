# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
import pytest
import numpy as np
import numpy.testing as nptest

from geospacepy.spherical_geometry import (angle_difference,
                                            grid_surface_integral,
                                            great_circle_distance)

def test_angle_difference_across_zero():
    ang2 = 10.
    ang1 = 350.
    #naiively 10. - 350. = -340, but adding 360 = 20
    expected_difference = 20.
    difference = angle_difference(ang1,ang2,'deg')
    assert(np.abs(expected_difference-difference)<.001)

@pytest.mark.parametrize('lat1,lat2,lon,algorithm',[(-85.,85.,10,'lawofcosines'),
                                                    (-85.,85.,10,'haversine'),
                                                    (50.,60.,10,'lawofcosines'),
                                                    (50.,60.,10,'haversine'),
                                                    (50.,50.001,10,'lawofcosines'),
                                                    (50.,50.001,10,'haversine')])
def test_great_circle_distance_on_meridian(lat1,lat2,lon,algorithm):
    dist_expected = np.radians(np.abs(lat2-lat1))
    dist = great_circle_distance(lat1,lon,lat2,lon,'deg',algorithm=algorithm)
    pertol = .2
    percent_error = np.abs(dist-dist_expected)/dist_expected*100.
    assert percent_error<pertol

@pytest.mark.parametrize('lon1,lon2,algorithm',[(-30.,140.,'lawofcosines'),
                                                (-30.,140.,'haversine'),
                                                (50.,80.,'lawofcosines'),
                                                (50.,80.,'haversine'),
                                                (50.,50.001,'lawofcosines'),
                                                (50.,50.001,'haversine')])
def test_great_circle_distance_on_equator(lon1,lon2,algorithm):
    dist_expected = np.radians(np.abs(lon2-lon1))
    dist = great_circle_distance(1.,lon1,1.,lon2,'deg',algorithm=algorithm)
    pertol = .2
    percent_error = np.abs(dist-dist_expected)/dist_expected*100. 
    assert percent_error<pertol

# @pytest.mark.parametrize('lat1,lat2,lon',[(50.,60.,10)])
# def test_great_circle_midpoint_on_meridian(lat1,lat2,lon):
#     lat_mid_expected = lat1+angle_difference(lat1,lat2,'deg')
#     lon_mid_expected = lon

#     lat_mid,lon_mid = great_circle_midpoint(lat1,lon,lat2,lon,'deg')
#     tol = 1e-6
#     assert np.abs(lat_mid-lat_mid_expected)<tol
#     assert np.abs(lon_mid-lon_mid_expected)<tol

@pytest.mark.parametrize('lat1,lon1,lat2,lon2',[(50.,10.,60.,20),
                                               (50.,10.,50.001,10.001),
                                               (50.,10.,50.0001,10.0001),])
def test_haversine_and_lawofcosines_produce_same_great_circle_distance(lat1,lon1,lat2,lon2):    
    dist_haversine = great_circle_distance(lat1,lon1,lat2,lon2,'deg',
                                            algorithm='haversine')
    dist_lawofcosines = great_circle_distance(lat1,lon1,lat2,lon2,'deg',
                                                algorithm='lawofcosines')
    pertol = .003
    percent_difference = np.abs(1.-dist_haversine/dist_lawofcosines)*100.
    assert percent_difference < pertol

def test_grid_surface_integral_of_half_unit_sphere_is_2pi():
    lats = [15.,45.,75.]
    lons = [45.,135.,225.,315.]
    grid_lats,grid_lons = np.meshgrid(lats,lons,indexing='ij')
    #The area of a half-sphere is 2*pi*R**2
    grid_values = np.ones_like(grid_lats)
    expected_integrated_values = 2*np.pi
    integrated_values = grid_surface_integral(grid_lats,
                                                grid_lons,
                                                grid_values,
                                                1.,'deg')
    tol = .001
    assert(np.abs(integrated_values-expected_integrated_values)<tol)



