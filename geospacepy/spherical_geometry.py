# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
import numpy as np
from numpy import (sin,cos,tan,arcsin,arccos,arctan2)

def _azifac(aziunit):
    if aziunit=='hour':
        return np.pi/12
    elif aziunit=='deg':
        return np.pi/180.
    elif aziunit=='rad':
        return 1.
    else:
        raise ValueError(('Invalid aziunit {}'.format(aziunit)
                          +' valid values are deg or hour or rad'))

def angle_difference(ang1,ang2,aziunit):
    """Difference between two angles in degrees or hours (ang2-ang1),
    taking into account wrapping

    PARAMETERS
    ----------
    
    ang1 : float or np.ndarray
        Angle(s) to subtract
    ang2 : float or np.ndarray
        Angle(s) to subtract from
    degorhour : str, optional
        'deg' for input angles and result in degrees
        'hour' for input angles in hours

    RETURNS
    -------

    diff : float or np.ndarray
        Difference (ang2-ang1)

    """
    ang2rad = _azifac(aziunit)
    y = np.sin(ang2*ang2rad-ang1*ang2rad)
    x = np.cos(ang2*ang2rad-ang1*ang2rad)
    diff = np.arctan2(y,x)/ang2rad
    return diff

def angle_midpoint(ang1,ang2,aziunit):
    """
    Midpoint between two angles in degrees or hours
    """
    return ang1 + angle_difference(ang1,ang2,aziunit)/2.

def _great_circle_distance(location1,location2,lonorlt='lt'):
    """Return angular distance in radians between n-by-2 numpy arrays
    location1, location2 (calculated row-wise so diff between
    location1[0,] and location2[0,]
    assuming that these arrays have the columns lat[deg],localtime[hours]
    and that they are points on a sphere of constant radius
    (the points are at the same altitude)
    """
    azi2rad = np.pi/12. if lonorlt=='lt' else np.pi/180
    wrappt = 24. if lonorlt=='lt' else 360.
    #Bounds check
    over = location1[:,1] > wrappt
    under = location1[:,1] < 0.
    location1[over,1]=location1[over,1]-wrappt
    location1[under,1]=location1[under,1]+wrappt

    if location1.ndim == 1 or location2.ndim == 1:
        dphi = np.abs(location2[1]-location1[1])*azi2rad
        a = (90-location1[0])/360*2*np.pi #get the colatitude in radians
        b = (90-location2[0])/360*2*np.pi
        C =  np.pi - np.abs(dphi - np.pi)#get the angular distance in longitude in radians
    else:
        dphi = np.abs(location2[:,1]-location1[:,1])*azi2rad
        a = (90-location1[:,0])/360*2*np.pi #get the colatitude in radians
        b = (90-location2[:,0])/360*2*np.pi
        C =  np.pi - np.abs(dphi - np.pi)#get the angular distance in longitude in radians
    return np.arccos(np.cos(a)*np.cos(b)+np.sin(a)*np.sin(b)*np.cos(C))

def _great_circle_distance_law_of_cosines(theta1,phi1,theta2,phi2):
    """Computes great circle distance between any number of paired 
    locations using the law of cosines. Warning:
    can be inaccurate for short distances. All angles in radians.
    Greek letter convention is theta for colatitude, phi for azimuth/longitude.
    """
    a = theta1
    b = theta2
    dphi = np.abs(phi2-phi1)
    C =  np.pi - np.abs(dphi - np.pi)
    dist = np.arccos(np.cos(a)*np.cos(b)+np.sin(a)*np.sin(b)*np.cos(C))
    return dist

def _great_circle_distance_haversine_formula(theta1,phi1,theta2,phi2):
    """Computes great circle distance between any number of paired
    locations using the haversine formula, which is known to be accurate
    for short distances. All angles in radians.
    Greek letter convention is theta for colatitude, phi for azimuth/longitude
    """
    lambda1 = np.pi/2-theta1
    lambda2 = np.pi/2-theta2
    delta_lambda = lambda2-lambda1
    delta_phi = phi2-phi1
    def hav(x):
        """Haversine function"""
        return np.sin(x*.5)**2.
    def archav(theta):
        """Inverse haversine"""
        return 2.*np.arcsin(np.sqrt(theta))
    def archav2(theta):
        """Inverse haversine using arctan2"""
        return 2.*np.arctan2(np.sqrt(theta),np.sqrt(1.-theta))

    hav_dist = hav(delta_lambda)+np.cos(lambda1)*np.cos(lambda2)*hav(delta_phi)
    
    dist = archav2(hav_dist)
    return dist

def great_circle_distance(lat1,azi1,lat2,azi2,aziunit,algorithm='lawofcosines'):
    """Computes great circle distance between any number of paired
    locations

    PARAMETERS
    ----------

    lat1 : float or np.ndarray
        Latitude(s) of first location(s)
    azi1 : float or np.ndarray
        Azimuth(s) (longitude,localtime) of first location(s) in pair(s)
    theta2 : float or np.ndarray
        Latitude(s) of second location(s)
    azi2 : float or np.ndarray
        Azimuth(s) (longitude,localtime) of second location(s) in pair(s) 
    aziunit : str
        Unit of passed azimuth(s), valid values: hour, deg or rad
    algorithm : str,optional
        Which formula to use in computing the great circle distance.
        Valid values are 'lawofcosines' or 'haversine'.

    RETURNS
    -------

    dist : float or np.ndarray
        Distance (in aziunit) between location(s) 1 and 2

    """
    algorithms = {
                    'lawofcosines':_great_circle_distance_law_of_cosines,
                    'haversine':_great_circle_distance_haversine_formula
                }
    if algorithm not in algorithms:
        raise ValueError('Invalid algorithm {}'.format(algorithm)
                        +' valid values: {}'.format([key for key in algorithms]))

    azi2rad = _azifac(aziunit)
    theta1,theta2 = np.radians(90.-lat1),np.radians(90.-lat2)
    phi1,phi2 = azi1*azi2rad,azi2*azi2rad
    dist = algorithms[algorithm](theta1,phi1,theta2,phi2)
    return dist

# There is currently a problem with this algorithm
# def great_circle_midpoint(lat1,azi1,lat2,azi2,aziunit):
#     """
#     Finds the midpoint latitude and azimuth (longitude or localtime) between 
#     two locations along a great circle arc.

#     PARAMETERS
#     ----------

#     lat1 : float or np.ndarray
#         Latitude(s) of first location(s)
#     azi1 : float or np.ndarray
#         Azimuth(s) (longitude,localtime) of first location(s) in pair(s)
#     theta2 : float or np.ndarray
#         Latitude(s) of second location(s)
#     azi2 : float or np.ndarray
#         Azimuth(s) (longitude,localtime) of second location(s) in pair(s) 
#     aziunit : str
#         Unit of passed azimuth(s), valid values: hour, deg or rad

#     RETURNS
#     -------
    
#     lat_mid : float or np.ndarray
#         Latitude(s) of midpoint(s)
#     azi_mid : float or np.ndarray
#         Longitude(s) (or localtimes) of midpoint(s) (unit determined by aziunit)

#     """
#     azi2rad = _azifac(aziunit)

#     a = np.radians(90.-lat1)
#     b = np.radians(90.-lat2)
#     c = great_circle_distance(lat1,azi1,lat2,azi2,aziunit,
#                                 algorithm='lawofcosines')
#     dphi = (azi2-azi1)*azi2rad
#     C =  np.pi - np.abs(dphi - np.pi)
    
#     #original: g = arccos((cos(b)-cos(c)*cos(a)*sin(c/2))*sin(c/2)/sin(c)+cos(a)*cos(c/2))
#     cos_g = cos(a)*cos(c/2)+((cos(b)-cos(a)*cos(c))/sin(c))*sin(c/2)
#     g = arccos(cos_g)
#     sin_I = (sin(c/2)*sin(b)*sin(C)/(sin(c)*sin(g)))
#     I = arcsin(sin_I)
#     lat_mid = 90.-np.degrees(g)
#     azi_mid = azi1+I/azi2rad
#     return lat_mid, azi_mid

def great_circle_rectangle_area(lats_bottom,lats_top,azis_left,azis_right,r,aziunit):
    """Calculate the surface area of any number of 'great-circle rectangles'
    (spherical quadrilaterals, whose 4 sides are all great circle arcs) 
    on the surface of a sphere of radius r.

    PARAMETERS
    ----------
        lats_bottom - float or np.ndarray 
            Latitude(s) specifying the bottom arc 
            of the great-circle rectangle(s)
        lats_top - float or np.ndarray
            Latitude(s) specifying the top arc 
            of the great-circle rectangle(s)
        azis_left - float or np.ndarray
            Azimuth/Longitude(s) specifying the left side arc 
            of the great-circle rectangle(s)
        azis_right - float or np.ndarray
            Azimuth/Longitude(s) specifying the right side arc 
            of the great-circle rectangle(s)
        r - float or np.ndarray
            Radius of the sphere on which to calculate the area
        aziunit - str
            Units of the azimuth/longitude ('deg' for degrees,
            'hour' for hours, 'rad' for radians)

    RETURNS
    -------
        
        areas - float or np.ndarray
            Surface areas of the great circle rectangles

    """
    azi2rad = _azifac(aziunit)

    if np.any((lats_top-lats_bottom) < 0):
        raise ValueError(('Latitudes for bottom edge cannot be > '
                          +'latitudes for top edge'))
    
    dazi = angle_difference(azis_left,azis_right,aziunit)
    dazi = np.mod(dazi,2*np.pi/azi2rad)

    if np.any(dazi*azi2rad>np.pi):
        raise ValueError(('Angles > 180 between longitudes for the right edge'
                          +' and longitudes for the left edge were found.'
                          +' usually this means the right and left edge inputs'
                          +' got switched.')) 

    theta_top = (90.-lats_top)*np.pi/180. #theta / lat - converted to radians
    theta_bottom = (90.-lats_bottom)*np.pi/180. #theta / lat - converted to radians
    dphi = np.abs(dazi)*azi2rad #delta phi / lon - converted to radians
    areas = np.abs(r**2*dphi*(np.cos(theta_top)-np.cos(theta_bottom)))
    return areas

def grid_surface_integral(grid_lats,grid_azis,grid_values,sphere_radius,aziunit):
    """Calculate the approximate surface integral of a field of values
    specified on a grid of m latitudes and n longitudes

    PARAMETERS
    ----------

        grid_lats - np.ndarray, shape=(m,n)
            The latitude values of the grid points. Latitudes must
            change along dimension 0, and stay constant along dimension 1
        grid_azis - np.ndarray, shape=(m,n)
            The azimuth/longitude values of the grid points. Longitudes must
            change along dimension 1, and stay constant along dimension 0
        grid_values - np.ndarray, shape=(m,n)
            The values of the field be integrated at the grid points
        sphere_radius - float
            The radius of the sphere on which the surface integral will
            be computed
        aziunit - str
            The unit of the azimuthal angles ('deg' for longitude, 'hour' for
            localtimes)

    RETURNS
    -------

        integrated_values - float
            The result of the surface integral, 
            units: units of grid_values * units of sphere_radius ** 2

    """

    if np.any(np.not_equal(grid_lats[:,0],grid_lats[:,1])):
        raise ValueError(('Latitudes are not the same in columns 0 and 1'
                          +'this function expects the grid arrays to vary'
                          +'in latitude along dimension 0'))
    if np.any(np.not_equal(grid_azis[0,:],grid_azis[1,:])):
        raise ValueError(('Azimuths/longitudes are not the same in rows 0 and 1'
                          +'this function expects the grid arrays to vary'
                          +'in longitude along dimension 1'))

    lats = grid_lats[:,0]
    dlats = np.diff(lats)
    dlat = np.abs(np.nanmedian(dlats))

    azis = grid_azis[0,:]
    dazis = angle_difference(azis[:-1],azis[1:],aziunit)
    dazis = np.mod(dazis,2*np.pi/_azifac(aziunit)) #Ensure > 0
    dazi = np.nanmedian(dazis)

    gridcell_bottom_lats = grid_lats-dlat/2.
    gridcell_top_lats = grid_lats+dlat/2.
    gridcell_left_azis = angle_difference(dazi/2.,grid_azis,aziunit)
    gridcell_right_azis = grid_azis+dazi/2.
    
    gridcell_areas = great_circle_rectangle_area(gridcell_bottom_lats,
                                                    gridcell_top_lats,
                                                    gridcell_left_azis,
                                                    gridcell_right_azis,
                                                    sphere_radius,
                                                    aziunit)

    integrated_values = np.nansum(gridcell_areas*grid_values)
    return integrated_values