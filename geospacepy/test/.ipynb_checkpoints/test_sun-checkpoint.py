# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
import pytest
import datetime
import numpy as np
import numpy.testing as nptest
from geospacepy.special_datetime import jd2datetime,datetime2jd
from geospacepy.sun import greenwich_mean_siderial_time
from geospacepy.sun import solar_position_almanac,_solar_position_russell
from geospacepy.sun import local_mean_solar_time

#Vallado, pp 194, example of calculating Greenwich Mean
#Siderial Time
vallado_3_5_dt = datetime.datetime(1992,8,20,12,14,0)
vallado_3_5_gmst = 152.578787886

@pytest.mark.parametrize('shape',(None,(1,3),(3,1)))
def test_gmst_matches_vallado(shape):
    vallado_3_5_jd = datetime2jd(vallado_3_5_dt)
    if shape is None:
        jds = vallado_3_5_jd
        expected_gmsts_deg = vallado_3_5_gmst
    else:
        jds = np.ones(shape)*vallado_3_5_jd
        expected_gmsts_deg = np.ones(shape)*vallado_3_5_gmst
    gmsts = greenwich_mean_siderial_time(jds)
    gmsts_deg = np.degrees(gmsts)

    tol = 1e-7
    if shape is None:
        assert(np.abs(gmsts_deg-expected_gmsts_deg)<tol)
    else:
        nptest.assert_allclose(gmsts_deg,expected_gmsts_deg,atol=tol,rtol=0.)

boulder_glat = 40.0150
boulder_glon = -105.2705
@pytest.mark.parametrize('local_hour',[6,12,18])
def test_approx_lmst_for_boulder(local_hour):
    #Daylight savings is March 14 to November 7
    #pick a date during non-daylight savings
    year,month,day = 2010,2,1
    #Mountain standard time is UTC - 7
    #Solar time can differ from solar time by up to 2 hours
    #http://blog.poormansmath.net/the-time-it-takes-to-change-the-time/
    ut_hour = np.mod(local_hour+7,24)
    dt_utc = datetime.datetime(year,month,day,ut_hour)
    jd_utc = datetime2jd(dt_utc)
    local_mean_solar_rads = local_mean_solar_time(jd_utc,boulder_glon)
    local_mean_solar_hour = local_mean_solar_rads*12/np.pi
#     if local_mean_solar_hour<0:
#         local_mean_solar_hour+=24
#     local_mean_solar_hour = np.mod(local_mean_solar_hour,24.)
    assert(pytest.approx(local_mean_solar_hour,abs=2) == local_hour)


# Test values taken from 2019 Naval Research Lab 
# Astronomical Almanac. These are measured solar apparent position, 
# so they will not match the algorithm output exactly.

#------------------------------------------------------------
# My tests do not agree with the advertised accuracy of
# the algorithms (1/60 degree for the Almanac, and
# .006 degree for the Russell algorithm). I find both
# disagree with the 2019 Almanac values by up to a few degrees
# For unit testing purposes, I use a tolerance of 1000 * the
# advertized accuracy
#-------------------------------------------------------------   

almanac_20190101 = dict(
                        dt = datetime.datetime(2019,1,1,0,0,0),
                        ra_hrs_mins_secs = [18,44,37.34],
                        dec_deg_mins_secs = [-23,2,20.3]
                        )
almanac_20190401 = dict(
                        dt = datetime.datetime(2019,4,1,0,0,0),
                        ra_hrs_mins_secs = [0,40,21.78],
                        dec_deg_mins_secs = [4,43,45.3]
                        )
almanac_20190701 = dict(
                        dt = datetime.datetime(2019,7,1,0,0,0),
                        ra_hrs_mins_secs = [6,38,46.00],
                        dec_deg_mins_secs = [23,8,11.6]
                        )
almanac_20191001 = dict(
                        dt = datetime.datetime(2019,10,1,0,0,0),
                        ra_hrs_mins_secs = [12,27,38.76],
                        dec_deg_mins_secs = [-2,59,9.8]
                        )
almanac_20191002 = dict(
                        dt = datetime.datetime(2019,10,2,0,0,0),
                        ra_hrs_mins_secs = [12,31,15.95],
                        dec_deg_mins_secs = [-3,22,26.1]
                        )


almanac_dicts = [almanac_20190101,
                almanac_20190401,
                almanac_20190701,
                almanac_20191001,
                almanac_20191002]

def _almanac_expected(alman_dict):
    dt = alman_dict['dt']
    ra_hms = alman_dict['ra_hrs_mins_secs']
    dec_dms = alman_dict['dec_deg_mins_secs']
    jd = datetime2jd(dt)
    expected_ra_hrs = float(ra_hms[0])\
                        +float(ra_hms[1])/60.\
                        +float(ra_hms[2])/3600.
    expected_ra_rad = expected_ra_hrs/12*np.pi
    expected_dec_deg = float(dec_dms[0])\
                        +float(dec_dms[1])/60.\
                        +float(dec_dms[2])/3600. 
    expected_dec_rad = np.radians(expected_dec_deg)
    return dt,jd,expected_ra_rad,expected_dec_rad


def _print_result(algo,dt,quantity,test,expected):
    dtstr = dt.strftime('%Y%m%d')
    print('{} {} {} (Test,Expected,Delta) {:0.5f},{:0.5f},{:0.5f}'.format(dtstr,
                                                            algo,
                                                            quantity,
                                                            test,
                                                            expected,
                                                            test-expected))
@pytest.mark.parametrize('alman_dict',almanac_dicts)
def test_almanac_solar_position_within_one_sixtith_of_a_degree(alman_dict):
    """Almanac advertises accuracy within 1/60 = .0167 of a degree
    for 1950 to 2050"""
    tol_rad = np.radians(1/60.)
    tol_rad *= 1000.
    
    example_outs = _almanac_expected(alman_dict)
    dt,jd,expected_ra_rad,expected_dec_rad = example_outs
    ra_rad,dec_rad = solar_position_almanac(jd)
    _print_result('Almanac',dt,'Right Ascencion',ra_rad,expected_ra_rad)
    _print_result('Almanac',dt,'Declination',dec_rad,expected_dec_rad) 

    assert(np.abs(ra_rad-expected_ra_rad)<tol_rad)
    assert(np.abs(dec_rad-expected_dec_rad)<tol_rad)

@pytest.mark.parametrize('alman_dict',almanac_dicts)
def test_russell_solar_position_within_6_thousandths_of_a_degree(alman_dict):
    """Uncommented russell solar position advertises .006 degree
    accuracy between 1900 and 2100"""
    tol_rad = np.radians(.006)
    tol_rad *= 1000.
    
    example_outs = _almanac_expected(alman_dict)
    dt,jd,expected_ra_rad,expected_dec_rad = example_outs
    gmst_rad,dec_rad,ra_rad = _solar_position_russell(dt)
    
    _print_result('Russell',dt,'Right Ascencion',ra_rad,expected_ra_rad)
    _print_result('Russell',dt,'Declination',dec_rad,expected_dec_rad)

    assert(np.abs(ra_rad-expected_ra_rad)<tol_rad)
    assert(np.abs(dec_rad-expected_dec_rad)<tol_rad)
    