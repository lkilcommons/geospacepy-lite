# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
import numpy as np
import datetime
from geospacepy.special_datetime import datetime2jd,dt_j2000
from geospacepy.array_management import BroadcastLenOneInputsToMatchArrayInputs

def solar_position_almanac(jds):
    """
    Finds the apparent solar right ascension and 
    declination for any number of julian dates.

    This algorithm is entitled 'Low precision formulas for the Sun'
    and can be found in section C of the Naval Research Laboratory
    Astronomical Almanac (2019). This formula should also be in
    other editions of the Almanac.

    The Almanac describes this formula as yeilding a precision
    better than 1' (1/60 degrees) for the years 1950 to 2050

    Parameters
    ----------

    jds : np.ndarray or float
        Julian dates for which to calculate solar positions

    Returns
    -------

    alpha : np.ndarray or float (matches input)
        Solar apparent right ascension (angle in equatorial
        plane measured clockwise from the vernal equinox direction),
        in radians

    delta : np.ndarray or float (matches input)
        Solar declination (equiv. to subsolar latitude), in radians

    """

    jd_j2000_epoch = datetime2jd(dt_j2000)
    jd2000 = jds - jd_j2000_epoch #J2000 epoch = 2451545.0 

    #Solar mean longitude (degrees)
    L = 280.460 + .9856474*jd2000
    L = np.mod(L,360.)

    #Solar mean anomaly (degrees)
    g = 357.528+0.9856003*jd2000
    g = np.mod(g,360.)

    #Solar ecliptic longitude (degrees)
    g_rad = np.radians(g)
    lam = L + 1.915*np.sin(g_rad)+.020*np.sin(2*g_rad)

    #Solar ecliptic latitude (degrees)
    beta = 0.

    #Obliquity of the ecliptic (degrees)
    epsilon = 23.439 - .0000004*jd2000

    #Right Ascension (there are 2 formulae...prefer the second b/c arctan)
    #Similar to subsolar longitude but measured counterclockwise from
    #vernal equinox direction instead of counterclockwise from
    #prime meridian
    epsilon_r = np.radians(epsilon)
    lam_r = np.radians(lam)
    t = np.tan(epsilon_r/2)**2.
    f = 180./np.pi
    lam_r = np.radians(lam)
    alpha = lam - f*t*np.sin(2.*lam_r) + (f/2.)*t**2*np.sin(4*lam_r)
    #alpha = np.arctan2(np.cos(epsilon_r)*np.sin(lam_r),np.cos(lam_r))
    alpha_r = np.radians(alpha)

    #Declination (equiv. to subsolar latitude, always < 90.)
    delta = np.degrees(np.arcsin(np.sin(epsilon_r)*np.sin(lam_r)))
    delta_r = np.radians(delta)
    return alpha_r,delta_r

def _solar_position_russell(dt):
    """This function is DEPRECATED use solar_position_almanac instead.
    There is reliable documentation for the solar_position_almanac
    algorithm, whereas the Russell algorithm only references
    'private communication'
    
    The following is the fortran code from which this code was translated:

    From C.T. Russell, (1971) "Geophysical Coordinate Transformations",
    Cosmic. Electrodyn. 2, 184-196
    ...
    G.D. Mead (private communication) has written a simple subroutine to\
    calculate the position of the Sun in GEI coordinates. It is accurate
    for years 1901 through 2099, to within 0.006 deg. The input is the
    year, day of year and seconds of the day in UT. The output is
    Greenwich Mean Sideral Time in degrees, the ecliptic longitude,
    apparent right ascension and declination of the Sun in degrees.
    The listing of this program follows. We note that the cartesian
    coordinates of the vector from the Earth to the Sun are:

      X = cos(SRASN) cos(SDEC)
      Y = sin(SRASN) cos(SDEC)
      Z = sin(SDEC)

      SUBROUTINE SUN(IYR, IDAY, SECS, GST, SLONG, SRASN, SDEC)
    C PROGRAM TO CALCULATE SIDEREAL TIME AND POSITION OF THE SUN.
    C GOOD FOR YEARS 1901 THROUGH 2099. ACCURACY 0.006 DEGREE.
    C INPUT IS IYR, IDAY (INTEGERS), AND SECS, DEFINING UN. TIME.
    C OUTPUT IS GREENWICH MEAN SIDEREAL TIME (GST) IN DEGREES,
    C LONGITUDE ALONG ECLIPTIC (SLONG), AND APPARENT RIGHT ASCENSION
    C AND DECLINATION (SRASN, SDEC) OF THE SUN, ALL IN DEGREES
      DATA RAD /57.29578/
      DOUBLE PRECISION DJ, FDAY
      IF(IYR. LT. 1901. OR. IYR. GT. 2099) RETURN
      FDAY = SECS/86400
      DJ = 365* (IYR-1900) + (IYR-1901)/4 + IDAY + FDAY -0.5D0
      T = DJ / 36525
      VL = DMOD (279.696678 + 0.9856473354*DJ, 360.D0)
      GST = DMOD (279.690983 + 0.9856473354*DJ + 360.*FDAY + 180., 360.D0)
      G = DMOD (358.475845 + 0.985600267*DJ, 360.D0) / RAD
      SLONG = VL + (1.91946 -0.004789*T)*SIN(G) + 0.020094*SIN (2.*G)
      OBLIQ = (23.45229 -0.0130125*T) / RAD
      SLP = (SLONG -0.005686) / RAD
      SIND = SIN (OBLIQ)*SIN (SLP)
      COSD = SQRT(1.-SIND**2)
      SDEC = RAD * ATAN (SIND/COSD)
      SRASN = 180. -RAD*ATAN2
      (COTAN (OBLIQ)*SIND/COSD, -COS (SLP)/COSD)
      RETURN
      END
    """
    iyear = dt.year
    iday = dt.timetuple().tm_yday
    secs = dt.hour*3600.+dt.minute*60.+dt.second
    fday = secs/86400.
    dj = 365*(iyear-1900)+(iyear-1901)/4 + iday + fday - .5
    t = dj/36525.
    vl = np.mod(279.696678 + 0.9856473354*dj, 360)
    gst = np.mod(279.690983 + 0.9856473354*dj + 360.*fday + 180., 360.)
    g = np.mod(358.475845 + 0.985600267*dj, 360.) * np.pi/180.
    slong = vl + (1.91946 -0.004789*t)*np.sin(g) + 0.020094*np.sin(2.*g)
    obliq = (23.45229 -0.0130125*t) * np.pi/180.

    slp = (slong - 0.005686) * np.pi/180.
    sin_d = np.sin(obliq)*np.sin(slp)
    cos_d = np.sqrt(1-sin_d**2)
    sdec = np.arctan(sin_d/cos_d)
    sransn = np.pi - np.arctan2(1/np.tan(obliq)*sin_d/cos_d,
                                -1*np.cos(slp)/cos_d)
    #GST is in degrees    
    gst = np.radians(gst)
    return gst,sdec,sransn

def greenwich_mean_siderial_time(jds):
    """Calculate the angle in the plane of the equator
    between the vernal equinox direction and the prime meridian (the 
    line of longitude through Greenwich, England).
    
    Parameters
    ----------

    jds : float or np.ndarray
        The julian date(s) of the times for which the GMST should
        be calculated

    Returns
    -------

    theta_GST : float or np.ndarray
        The Greenwich Mean Siderial Time in radians

    Notes
    -----

    Because this calculation depends on the actual exact number of earth
    rotations since the J2000 epoch, the time (julian date) strictly speaking
    should be in the UT1 system (the time system determined from observations
    of distant stars), because this system takes into account the small changes
    in earth's rotation speed.
    
    Generally though, UTC is available instead of UT1. UTC is determined
    from atomic clocks, and is kept within +- 1 second of UT1 
    by the periodic insertion of leap seconds.
    
    """
    jd_j2000 = datetime2jd(dt_j2000)
    t_ut1 = (jds-jd_j2000)/36525. #Get Julian centuries since the j2000.0 epoch
    #Note that this formula can be broken up into a two part (hours and seconds) version using a two part
    #T_UT1. Where 876600 is multiplied by 3600., and in the exponentiation, the accuracy can be increased
    #by breaking up the T_UT1
    theta_GST_s = 67310.54841+(876600.*3600.+8640184.812866)*t_ut1+.093104*t_ut1**2-6.2e-6*t_ut1**3
    
    # NOTE: In Python (and Numpy), the output of modulus is 
    # always the same sign as the divisor (360. in the case of angles)
    # so we can treat negative out of bounds the same
    # as positive
    
    #Make sure abs(theta_GST) <= 86400 seconds
    theta_GST_s = np.mod(theta_GST_s,86400.)

    #Convert theta_GST to degrees from seconds
    theta_GST = theta_GST_s/240.

    # Ensure in 0 to 360.
    theta_GST = np.mod(theta_GST,360.)

    # Radians
    theta_GST = theta_GST * np.pi / 180.
    return theta_GST

@BroadcastLenOneInputsToMatchArrayInputs
def local_hour_angle(jds,glons):
    """
    Finds local hour angle in radians. The sign convention is that of astronomy 
    (positive to the west, meaning angle increases opposite the 
    rotation direction of earth).

    Under this sign convention the hour angle (in units of hours) is the
    time it will be before for the sun is directly overhead at the specified
    location and time.

    Parameters
    ----------

    jds : np.ndarray or float
        Time (as julian date)

    glons : np.ndarray or float
        Geographic longitude

    Returns
    -------

    lhas : np.ndarray or float
        Local hour angles for locations at specified times (in radians)
    
    .. note::

        See also Vallado Figure 3-9 (pp.157)

    """
    sra,dec = solar_position_almanac(jds)
    gmst = greenwich_mean_siderial_time(jds) 
    
    #Greenwich Mean Sideral Time, Right Ascension, and longitude are
    #all measured with positive angles counterclockwise about the north pole
    #(with the earth's rotation direction).
    
    #But the hour angle is defined as positive opposite the earth's rotation
    phi = np.radians(glons)
    lhas = (gmst+phi) - sra
    return lhas

@BroadcastLenOneInputsToMatchArrayInputs
def local_mean_solar_time(jds,glons):
    """
    Find the local solar time (using the mean equinox)

    Parameters
    ----------

    jds : np.ndarray or float
        Time (as julian date)

    glons : np.ndarray or float
        Geographic longitude

    Returns
    -------

    lmsts : np.ndarray or float
            Local mean solar time at locations at specified times (in radians)

    Notes
    -----

    As an angle, the solar local time increases in the same direction
    as the ISO 6709 longitude convention (eastward positive).

    The differences between hour angle and local time is that hour angle
    is:
    1.) Measured with positive in the westward direction
    2.) Zero when the sun is directly overhead (instead of 12 for solar time)

    See also Vallado pp. 184

    """
    lhas = local_hour_angle(jds,glons)
    #lhas *= -1. #Convert to positive in the eastward direction
    lmsts = lhas + np.pi  #Equiv to + 12 in hours units
    return lmsts

@BroadcastLenOneInputsToMatchArrayInputs
def solar_zenith_angle(jds,glats,glons):
    """
    Finds solar zenith angle using Astronomical Almanac low-accuracy
    solar position.

    Parameters
    ----------

    jds : np.ndarray or float
        Time (as julian date)

    glats : np.ndarray or float
        Geographic (geocentric-spherical) latitude of the location

    glons : np.ndarray or float
        Geographic longitude of the location

    Returns
    -------

    szas : np.ndarray or float
        Solar Zenith angles for the time/location combinations
        specified (in radians)

    """
    lam = np.radians(glats)
    phi = np.radians(glons)
    
    sra,sdec = solar_position_almanac(jds)
    sha = local_hour_angle(jds,glons)

    cossza = np.sin(lam)*np.sin(sdec) + np.cos(lam)*np.cos(sdec)*np.cos(sha)
    szas = np.arccos(cossza)
    return szas