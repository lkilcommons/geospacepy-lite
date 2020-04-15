import pytest
import numpy as np
import numpy.testing as nptest

from geospacepy.spherical_geometry import angle_difference,grid_surface_integral

def test_angle_difference_across_zero():
    ang2 = 10.
    ang1 = 350.
    #naiively 10. - 350. = -340, but adding 360 = 20
    expected_difference = 20.
    difference = angle_difference(ang1,ang2,'deg')
    assert(np.abs(expected_difference-difference)<.001)

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