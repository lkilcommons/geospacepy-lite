from geospacepy import astrodynamics2
import pytest
import numpy as np
from numpy import testing as nptest
import datetime,os

def test_julian_date_matches_vallado():
	"""
	Test that the julian date function in astrodyanmics matches 
	an example from Vallado: Fundamentals of Astrodynamics and Applications
	"""
	#Test using example in vallado pp. 194
	y,mo,d,h,m,s = 1992,8,20,12,14,0
	jd_ut1 = astrodynamics2.ymdhms2jd(y,mo,d,h,m,s)
	expected_jd = 2448855.009722
	#should be 2448855.009722" % (jd_ut1)
	assert abs(jd_ut1-expected_jd) < .000001

def test_greenwich_siderial_time_matches_vallado():
	y,mo,d,h,m,s = 1992,8,20,12,14,0
	jd_ut1 = astrodynamics2.ymdhms2jd(y,mo,d,h,m,s)
	theta_gst = astrodynamics2.jd2gst(jd_ut1,deg=True)
	expected_gst = 152.578787886
	assert (theta_gst-expected_gst) < .0001

if __name__ == '__main__':
	pytest.main()
