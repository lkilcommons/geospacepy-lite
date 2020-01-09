from geospacepy import special_datetime
import pytest
import numpy as np
from numpy import testing as nptest
import datetime,os,pkgutil

#Example 3-4: Converting from Gregorian date (YMD HMS) to Julian date
#Vallado, Fundamentals of Astrodynamics and Applications 3rd Ed., pp. 190
vallado_3_4_dt = datetime.datetime(1996,10,26,14,20,0)
vallado_3_4_jd = 2450383.09722222
#Example 3-13 : Converting from Julian date to Gregorian (YMD HMS)
#Vallado, Fundamentals of Astrodynamics and Applications 3rd Ed., pp. 209
vallado_3_13_jd = 2449877.3458762
vallado_3_13_dt = datetime.datetime(1995,6,8,20,18,3,703700)
#Example 3-10 : Converting HMS to Time of Day (SOD)
#Vallado, Fundamentals of Astrodynamics and Applications 3rd Ed., pp. 206
vallado_3_10_dt = datetime.datetime(2000,1,1,13,22,45,98000)
vallado_3_10_sod = 48165.98
#Example 3-12 : Converting Year-Month-Day Hour-Minute-Second to Days of Year
#Vallado, Fundamentals of Astrodynamics and Applications 3rd Ed., pp. 20
vallado_3_12_dt = datetime.datetime(2001,3,18,12,14,0)
vallado_3_12_year = 2001
vallado_3_12_doy = 77.5097222
#Test output of matlab 'datevec' to matlab 'datenum'
# dn = datenum(1995,6,8,20,18,3.7037); fprintf('%10.15f\n',dn);
# 728818.845876200241037
matlab_dt = datetime.datetime(1995,6,8,20,18,3,703700)
matlab_datenum = 728818.845876200241037

def ex_as_arr(example,shape):
    """Make arrays of an example input and output for one of the
    scalar time conversion functions for use in testing the vectorized
    version of that conversion function"""
    if shape is None:
        return example
    elif isinstance(shape,int):
        return [example for i in range(shape)]
    elif isinstance(shape,tuple):
        array_type = 'object' if isinstance(example,datetime.datetime) else float
        example_arr = np.full(shape,example,dtype=array_type)
        return example_arr
    else:
        raise ValueError('Invalid inputs')

def test_julian_date_matches_vallado():
    """
    Test that the julian date function matches
    an example from Vallado: Fundamentals of Astrodynamics and Applications
    """
    dt = vallado_3_4_dt
    jd_ut = special_datetime.datetime2jd(dt)
    expected_jd = vallado_3_4_jd
    assert abs(jd_ut-expected_jd) < .000001

def test_gregorian_date_matches_vallado():
    """
    Test that the gregorian-julian date function in astrodyanmics matches
    an example from Vallado: Fundamentals of Astrodynamics and Applications
    """
    #Test using example in vallado pp. 209
    given_jd = vallado_3_13_jd
    dt = special_datetime.jd2datetime(given_jd)
    expected_dt = vallado_3_13_dt
    delta_t = (dt-expected_dt).total_seconds()
    #Better than millisecond accuracy, at least for this date!
    assert abs(delta_t) < 1.0e-4

@pytest.mark.parametrize('example_dt,example_jd',[(vallado_3_4_dt,
                                                    vallado_3_4_jd),
                                                    (vallado_3_13_dt,
                                                    vallado_3_13_jd)])
@pytest.mark.parametrize('shape',[3,(3,),(3,1),(1,3)])
def test_vectorized_datetime_to_julian_date(example_dt,example_jd,shape):
    """
    Test the vectorization function factory works to create
    datetimearr2jd, and that this function works with scalars,
    1D lists and 1D numpy arrays
    """
    dts = ex_as_arr(example_dt,shape)
    expected_jds = ex_as_arr(example_jd,shape)
    jds = special_datetime.datetimearr2jd(dts)
    #Comparison is in units of days (use tolerance of a milisecond)
    atol = .001/86400.
    nptest.assert_allclose(jds,expected_jds,rtol=0.,atol=atol)

@pytest.mark.parametrize('example_jd,example_dt',[(vallado_3_4_jd,
                                                    vallado_3_4_dt),
                                                    (vallado_3_13_jd,
                                                    vallado_3_13_dt)])
@pytest.mark.parametrize('shape',[3,(3,),(3,1),(1,3)])
def test_vectorized_julian_date_to_datetime(example_jd,example_dt,shape):
    """
    Test the vectorization function factory works to create
    datetimearr2jd, and that this function works with scalars,
    1D lists and 1D numpy arrays
    """
    jds = ex_as_arr(example_jd,shape)
    expected_dts = ex_as_arr(example_dt,shape)
    dts = special_datetime.jdarr2datetime(jds)
    timedeltas = (expected_dts-dts).flatten()
    delta_ts = np.array([abs(tdelta.total_seconds()) for tdelta in timedeltas])
    #Comparison in units of seconds use tolerance of a milisecond
    nptest.assert_allclose(delta_ts,np.zeros_like(delta_ts),rtol=0.,atol=1e-3)

def test_datetime_to_day_of_year():
    dt = vallado_3_12_dt
    expected_doy = vallado_3_12_doy
    doy = special_datetime.datetime2doy(dt)
    tol = .01/86400. #Accurate to 10 miliseconds (about 9 miliseconds different)
    assert abs(doy-expected_doy)<tol

@pytest.mark.parametrize('example_doy,example_dt,example_year',[(vallado_3_12_doy,
                                                                vallado_3_12_dt,
                                                                vallado_3_12_year)])
@pytest.mark.parametrize('shape',[3,(3,),(3,1),(1,3)])
def test_vectorized_day_of_year_to_datetime(example_doy,example_dt,example_year,
                                            shape):
    """Test that vectorization works correctly in the presence of additional
    input arguments (the year, for the day of year to datetime conversion)"""
    doy = ex_as_arr(example_doy,shape)
    year = ex_as_arr(example_year,shape)
    expected_dts = ex_as_arr(example_dt,shape)
    dts = special_datetime.doyarr2datetime(doy,year)
    timedeltas = (expected_dts-dts).flatten()
    delta_ts = np.array([abs(tdelta.total_seconds()) for tdelta in timedeltas])
    #Comparison in units of seconds use tolerance of ten miliseconds
    nptest.assert_allclose(delta_ts,np.zeros_like(delta_ts),rtol=0.,atol=1e-2)

def test_day_of_year_to_datetime():
    doy = vallado_3_12_doy
    year = vallado_3_12_year
    expected_dt = vallado_3_12_dt
    dt = special_datetime.doy2datetime(doy,year)
    delta_t = (dt-expected_dt).total_seconds()
    #This is accurate to about 9 miliseconds, not sure if this
    #is the tabulated values fault or what
    assert abs(delta_t) < .01

def test_datenum_to_datetime():
    dn = matlab_datenum
    expected_dt = matlab_dt
    dt = special_datetime.datenum2datetime(dn)
    delta_t = (dt-expected_dt).total_seconds()
    #Accurate to within 2 microseconds
    #(likely b/c more than microsecond precision in tabulated datenum)
    assert abs(delta_t) < .000002

def test_datetime_to_datenum():
    dt = matlab_dt
    expected_dn = matlab_datenum
    dn = special_datetime.datetime2datenum(dt)
    delta_t = (dn-expected_dn)/86400. #days to seconds
    #Accurate to a microsecond
    assert abs(delta_t) < .000001
