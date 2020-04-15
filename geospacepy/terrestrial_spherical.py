import numpy as np
from numpy import (cos,sin,tan,arccos,arcsin,arctan2)
from rotations import rot1,rot2,rot3
from scipy.spatial.transform import Rotation

R_EARTH = 6371200 #m

def ecef2eci(R_ECEF,theta_GST):
    """Rotate a vector from cartesian Earth Centered Earth Fixed (ECEF) to
    Cartesian Earth Centered Inertial (ECI)

    PARAMETERS
    ----------
        
        Vec_ECEF - np.ndarray
            Array of three component vectors [n x 3]
        theta_GST - float or np.ndarray
            Greenwich Sidereal Time angle 
    """
    R_ECI = rot3(-1*theta_GST,R_ECEF)
    return R_ECI

def eci2ecef(R_ECI,theta_GST):
    """Rotate a vector from cartesian Earth Centered Inertial (ECI) to 
    cartesian Earth Centered Earth Fixed (ECEF)"""
    R_ECEF = rot3(theta_GST,R_ECI)
    return R_ECEF

def ecef_cart2spherical(R_ECEF,deg=True):
    #R_ECEF is the cartesian Earth Centered Earth Fixed vector in any units
    #For clarity, function is not vectorized
    R_ECEF = R_ECEF.flatten() #Make sure the vector is 1-d
    r = sqrt(R_ECEF[0]**2+R_ECEF[1]**2+R_ECEF[2]**2)
    x = R_ECEF[0]
    y = R_ECEF[1]
    z = R_ECEF[2]
    longitude = arctan2(y,x) #Longitude is angle in x,y plane

    latitude =  arcsin(z/r) #Latitude is angle z-ward from x,y plane
    #Convert to degrees for return if deg switch on
    if deg:
        longitude = longitude*180./pi
        latitude = latitude*180./pi
    return array([r,latitude,longitude])

def ecef2enu(R_ECEF,lat,lon):
    #Rotate a vector from ecef to local east, north, up
    #coordinates centered at lat,lon

    lonrot = 90.+lon
    #lonrot = lon
    #lonrot[lonrot > 360.] = lonrot[lonrot>360.]-360.
    if lonrot > 360.:
        lonrot = lonrot-360.

    colat = 90.-lat
    R_ENU = rot1(colat,
                rot3(lonrot,R_ECEF,deg=True)
                ,deg=True)

    return R_ENU

def enu2ecef(R_ENU,lat,lon):
    #Rotate a vector from ecef to local east, north, up
    #coordinates centered at lat,lon

    lonrot = 90.+lon
    #lonrot = lon
    #lonrot[lonrot > 360.] = lonrot[lonrot>360.]-360.
    if lonrot > 360.:
        lonrot = lonrot-360.

    colat = 90.-lat
    R_ECEF = rot1(-1*colat*pi/180., rot3(-1*lonrot*pi/180.,R_ENU) )
    return R_ECEF