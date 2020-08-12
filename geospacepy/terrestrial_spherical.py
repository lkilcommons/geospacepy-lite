# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
import numpy as np
from numpy import (cos,sin,tan,arccos,arcsin,arctan2)
from geospacepy.array_management import (CheckInputsAreThreeComponentVectors,
                                        BroadcastLenOneInputsToMatchArrayInputs)
from geospacepy.rotations import rot1,rot2,rot3
from geospacepy.sun import greenwich_mean_siderial_time
#from scipy.spatial.transform import Rotation

R_EARTH = 6371200 #m

@CheckInputsAreThreeComponentVectors('R_ECEF')
def ecef2eci(R_ECEF,jds):
    """Rotate a vector from cartesian Earth Centered Earth Fixed (ECEF) to
    Cartesian Earth Centered Inertial (ECI)

    PARAMETERS
    ----------
        
    R_ECEF : np.ndarray
        Array of n three component vectors ( shape=(n,3) ) in ECEF
    jds : float or np.ndarray
        Time stamp (as julian date) for all (if float) or each (if array)
        of the vectors in R_ECEF

    RETURNS
    -------

    R_ECI : np.ndarray
        n three component vectors ( shape=(n,3) ) in ECI

    """
    theta_GST = greenwich_mean_siderial_time(jds)
    R_ECI = rot3(-1*theta_GST,R_ECEF)
    return R_ECI

@CheckInputsAreThreeComponentVectors('R_ECI')
def eci2ecef(R_ECI,jds):
    """Rotate a vector from cartesian Earth Centered Inertial (ECI) to 
    cartesian Earth Centered Earth Fixed (ECEF)

    PARAMETERS
    ----------
        
    R_ECI : np.ndarray
        Array of n three component vectors ( shape=(n,3) ) in ECI
    jds : float or np.ndarray
        Time stamp (as julian date) for all (if float) or each (if array)
        of the vectors in R_ECEF

    RETURNS
    -------

    R_ECEF : np.ndarray
        n three component vectors ( shape=(n,3) ) in ECEF

    """
    theta_GST = greenwich_mean_siderial_time(jds)
    R_ECEF = rot3(theta_GST,R_ECI)
    return R_ECEF

@CheckInputsAreThreeComponentVectors('R_ECEF')
def ecef_cart2spherical(R_ECEF):
    """Transform a vector in cartesian Earth Centered Earth Fixed (ECEF)
    to geocentric spherical (a.k.a geographic coordinates)

    PARAMETERS
    ----------

    R_ECEF : np.ndarray
        Array of n three component vectors ( shape=(n,3) ) in cartesian
        ECEF

    RETURNS
    -------

    lats : np.ndarray
        Array of n latitude values (shape=(n,))
    lons : np.ndarray
        Array of n longitude values (shape=(n,)). Sign convention
        is ISO 6709 (west is negative)
    rs : np.ndarray
        Array of n radii (distance from center of the earth) 
        (shape=(n,))
        (same units as R_ECEF)

    """
    xs = R_ECEF[:,0].flatten()
    ys = R_ECEF[:,1].flatten()
    zs = R_ECEF[:,2].flatten()
    rs = np.sqrt(xs**2.+ys**2.+zs**2.)
    lons = np.degrees(np.arctan2(ys,xs)) #Longitude is angle in x,y plane
    lats = np.degrees(np.arcsin(zs/rs)) #Latitude is angle z-ward from x,y plane
    return lats,lons,rs

@BroadcastLenOneInputsToMatchArrayInputs
def ecef_spherical2cart(lats,lons,rs):
    """Transform n positions in geocentric spherical, a.k.a. spherical Earth
    Centered Earth Fixed (ECEF), a.k.a geographic to cartesian Earth Centered Earth
    Fixed (ECEF)

    PARAMETERS
    ----------

    lats : float or np.ndarray
        Array of n latitude values (shape=(n,))
    lons : float or np.ndarray
        Array of n longitude values (shape=(n,)). Sign convention
        is ISO 6709 (west is negative)
    rs : float or np.ndarray
        Array of n radii (distance from center of the earth) 
        (shape=(n,))

    RETURNS
    -------

    R_ECEF : np.ndarray
        Array of n three component vectors ( shape=(n,3) ) in cartesian
        ECEF

    .. note::

        Within the code the ISO 31-11 convention for theta (polar, zero
        at positive z axis), and phi (azimuthal, zero at positive x axis) 
        is used
    """
    thetas = np.radians(90.-lats)
    phis = np.radians(lons)
    x = rs*np.sin(thetas)*np.cos(phis)
    y = rs*np.sin(thetas)*np.sin(phis)
    z = rs*np.cos(thetas)
    R_ECEF = np.concatenate([x.reshape(-1,1),
                             y.reshape(-1,1),
                             z.reshape(-1,1)],axis=1)
    return R_ECEF
    
@CheckInputsAreThreeComponentVectors('R_ECEF')
def ecef2enu(R_ECEF,lats,lons):
    """Rotate n vectors from Earth Centered Earth Fixed (ECEF)
    to local east, north, up coordinates centered at one location if only a
    single lat/lon pair is passed, or n locations if lats and lons are
    arrays of length n

    PARAMETERS
    ----------

    R_ECEF : np.ndarray
        Array of n 3-component vectors (shape=(n,3)) in cartesian ECEF

    lats : float or np.ndarray
        Latitude(s) of ENU coordinates
    
    lons : float or np.ndarray
        Longitude(s) of ENU coordinates
    
    RETURNS
    -------

    R_ENU : np.ndarray

        Array of n 3-component vecotrs (shape=(n,3)) in cartesian ENU
        (R_ENU[:,0] is eastward component, R_ENU[:,1] is northward, etc.)

    """
    #This works even if lon is < -360, as numpy modulus returns
    #the same sign as the denominator (360.) 
    lonrot = np.mod(90.+lons,360.)
    colat = 90.-lats

    R_ENU = rot1(np.radians(colat),rot3(np.radians(lonrot),R_ECEF))
    return R_ENU


@CheckInputsAreThreeComponentVectors('R_ENU')
def enu2ecef(R_ENU,lats,lons):
    """Rotate n vectors from East North Up (ENU) (relative to location(s)
    specified by lats and lons) to Earth Centered Earth Fixed (ECEF) coordinates
    
    PARAMETERS
    ----------

    R_ENU : np.ndarray
        Array of n 3-component vectors (shape=(n,3)) in cartesian
        ENU

    lats : float or np.ndarray
        Latitude(s) which define ENU directions
    
    lons : float or np.ndarray
        Longitude(s) which define ENU direction
    
    RETURNS
    -------

    R_ECEF : np.ndarray

        Array of n 3-component vectors (shape=(n,3)) in cartesian ECEF
    """

    #This works even if lon is < -360, as numpy modulus returns
    #the same sign as the denominator (360.) 
    lonrot = np.mod(90.+lons,360.)
    colat = 90.-lats
    R_ECEF = rot3(np.radians(-1*lonrot), rot1(np.radians(-1*colat),R_ENU) )
    return R_ECEF