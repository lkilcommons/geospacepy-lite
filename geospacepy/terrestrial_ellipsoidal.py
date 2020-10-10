import numpy as np
from warnings import warn
from geospacepy.array_management import (CheckInputsAreThreeComponentVectors,
                                        BroadcastLenOneInputsToMatchArrayInputs)

ECC_EARTH_SQUARED = .006694385
R_EARTH_MEAN_EQ = 6378137 #in m

@CheckInputsAreThreeComponentVectors('R_ECEF')
def ecef_cart2geodetic(R_ECEF,tol=1e-7,maxiters=100):
    """This implements Algorithm 12 (pp.179) of Fundamentals of Astrodynamics 
    and Applications (3rd Edition) by David Vallado. The algorithm's purpose
    is to transform the position vector of a spacecraft in earth-centered
    earth-fixed cartesian coordinates to geodetic latitude, geocentric
    longitude, and height above the surface of the earth.
    
    Parameters
    ----------

    R_ECEF : np.ndarray, shape=(n,3)
        n position vectors in cartesian ECEF, units must be meters
    tol : float
        Tolerance for iterative determination of geodetic latitude, in radians.
        The algorithm stops when the change in geodetic latitude from one
        iteration to the next is less than the tolerance. Defaults to 1e-7, 
        which is the value used in an example from the Vallado textbook 
        (Example 3-3)
    maxiters : int
        Maximum number of times the algorithm will attempt to refine the
        geodetic latitude before raising RuntimeError
    
    Returns
    -------

    gdlats : np.ndarray, shape=(n,)
        Geodetic latitudes of each position in R_ECEF
    gclons : np.ndarray, shape=(n,)
        Geocentric longitudes of each position in R_ECEF
    h_ellps : np.ndarray, shape=(n,)
        Ellipsoidal height of each position in meters 
        (the distance between each position and the ground, 
        measured perpendicular to the surface of an ellipsoidal 
        approximation of the earth)

    Warns
    -----

    UserWarning
        If ECEF positions are below earth's surface (because it may be
        the position was in kilometers instead of meters)

    """
    if np.any(np.logical_not(np.isfinite(R_ECEF))):
        raise ValueError(('Non-finite earth-centered-earth-fixed positions')
                          +' can not be converted to geodetic')

    X = R_ECEF[:,0].flatten()
    Y = R_ECEF[:,1].flatten()
    Z = R_ECEF[:,2].flatten()

    R = np.sqrt(X**2.+Y**2.+Z**2.)
    if np.any(R<R_EARTH_MEAN_EQ):
        warnstr = ('Apparently underground positions '
                  +'were passed to ecef_cart2geodetic '
                  +'( || R_ECEF || < {} m) '.format(R_EARTH_MEAN_EQ)
                  +'if you did not mean to calculate underground '
                  +'positions, check R_ECEF is in units of meters.')
        warn(warnstr,UserWarning)

    # Projection of spacecraft position on equatorial plane
    R_eq = np.sqrt(X**2+Y**2)

    #Determine the right ascension of the spacecraft
    sin_alpha = Y/R_eq
    cos_alpha = X/R_eq
    alpha = np.arctan2(Y,X)

    #Longitude is just right ascension, even for an ellipsoidal earth
    glons = alpha

    #Determine the declination of the spacecraft
    delta = np.arctan2(Z,R_eq)

    def C_earth(gdlats):
        """Ellipsoidal parameter used in formula for the equatorial
        projection of a position vector in terms of ellipsoid
        (geodetic) parameters.
            R_eq = (C_Earth+h_ellp)*cos(gdlat) (Vallado Eq 3-7)
        """
        return R_EARTH_MEAN_EQ/np.sqrt(1.-ECC_EARTH_SQUARED*np.sin(gdlats)**2)
    
    def refine_geodetic_latitude_estimate(gdlats):
        return np.arctan2(Z+C_earth(gdlats)*ECC_EARTH_SQUARED*np.sin(gdlats),R_eq)

    tolerance_reached = False
    prev_estimated = None
    gdlats_estimated = delta
    gdlats = None
    for i in range(maxiters):
        #print(i,gdlats_estimated)
        prev_estimated = gdlats_estimated.copy()
        gdlats_estimated = refine_geodetic_latitude_estimate(gdlats_estimated)
        if np.all(np.abs(gdlats_estimated-prev_estimated)<tol):
            #print("Converged after {} iterations".format(i))
            gdlats = gdlats_estimated
            break

    if gdlats is None:
        warnstr=('Geodetic latitude estimation failed to converge'
                +'to iteration-to-iteration tolerance {}'.format(tol)
                +'after {} iterations'.format(maxiters))
        raise RuntimeError(warnstr)

    #Find height above the surface of the ellipsoid
    h_ellps = R_eq/np.cos(gdlats)-C_earth(gdlats)
    return np.degrees(gdlats),np.degrees(glons),h_ellps
