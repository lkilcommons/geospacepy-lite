import numpy as np
from numpy import (sin,cos,tan,arcsin,arccos,arctan2)

def _azifac(degorhour):
    if degorhour=='hour':
        return np.pi/12
    elif degorhour=='deg':
        return np.pi/180.
    else:
        raise ValueError(('Invalid degorhour {}'.format(degorhour)
                          +' valid values are deg or hour'))

def angle_difference(ang1,ang2,degorhour='hour'):
    """Difference between two angles in degrees or hours (ang2-ang1),
    taking into account wrapping
    """
    ang2rad = _azifac(degorhour)
    y = np.sin(ang2*ang2rad-ang1*ang2rad)
    x = np.cos(ang2*ang2rad-ang1*ang2rad)
    return np.arctan2(y,x)/ang2rad

def angle_midpoint(ang1,ang2,degorhour='hour'):
    """
    Midpoint between two angles in degrees or hours
    """
    return ang1 + angle_difference(ang1,ang2,degorhour=degorhour)/2.

def great_circle_distance(location1,location2,lonorlt='lt'):
    """Returns n angular distances in radians between n-by-2 numpy arrays
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

def great_circle_midpoint(location1,location2,angDist='compute',lonorlt='lt'):
    """
    Finds the midpoint lat and lt or lon between two locations 
    along a great circle arc. Can pass angDist as an array to speed up process 
    if already computed, otherwise computes as needed using great_circle_distance
    """
    azi2rad = np.pi/12. if lonorlt=='lt' else np.pi/180
    wrappt = 24. if lonorlt=='lt' else 360.
    #Bounds check
    over = location1[:,1] > wrappt
    under = location1[:,1] < 0.
    location1[over,1]=location1[over,1]-wrappt
    location1[under,1]=location1[under,1]+wrappt

    if location1.ndim == 1 or location2.ndim == 1:
        a = (90-location1[0])*azi2rad
        b = (90-location2[0])*azi2rad
        if angDist is 'compute':
            c = great_circle_distance(location1,location2,lonorlt=lonorlt)
        else:
            c = angDist
        C = (location2[1]-location1[1])/24*2*np.pi
        azi1 = location1[1]
    else:
        a = (90-location1[:,0])*azi2rad
        b = (90-location2[:,0])*azi2rad
        if angDist is 'compute':
            c = great_circle_distance(location1,location2,lonorlt=lonorlt)
        else:
            c = angDist
        C = (location2[:,1]-location1[:,1])/24*2*np.pi
        azi1 = location1[:,1]

    #original: g = arccos((cos(b)-cos(c)*cos(a)*sin(c/2))*sin(c/2)/sin(c)+cos(a)*cos(c/2))
    cos_g = cos(a)*cos(c/2)+((cos(b)-cos(a)*cos(c))/sin(c))*sin(c/2)
    g = arccos(cos_g)
    sin_I = (sin(c/2)*sin(b)*sin(C)/(sin(c)*sin(g)))
    I = arcsin(sin_I)
    lat_mid = 90-g/(2*np.pi)*360
    azi_mid = azi1+I/(2*np.pi)*24
    return lat_mid, azi_mid

def great_circle_rectangle_area(lats_bottom,lats_top,azis_left,azis_right,r,degorhour):
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
        degorhour - str
            Units of the azimuth/longitude ('deg' for degrees,
            'hour' for hours)

    RETURNS
    -------
        
        areas - float or np.ndarray
            Surface areas of the great circle rectangles

    """
    ang2rad = _azifac(degorhour)

    if np.any((lats_top-lats_bottom) < 0):
        raise ValueError(('Latitudes for bottom edge cannot be > '
                          +'latitudes for top edge'))
    
    dazi = angle_difference(azis_left,azis_right,degorhour=degorhour)
    dazi = np.mod(dazi,2*np.pi/ang2rad)

    if np.any(dazi*ang2rad>np.pi):
        raise ValueError(('Angles > 180 between longitudes for the right edge'
                          +' and longitudes for the left edge were found.'
                          +' usually this means the right and left edge inputs'
                          +' got switched.')) 

    theta_top = (90.-lats_top)*np.pi/180. #theta / lat - converted to radians
    theta_bottom = (90.-lats_bottom)*np.pi/180. #theta / lat - converted to radians
    dphi = np.abs(dazi)*ang2rad #delta phi / lon - converted to radians
    areas = np.abs(r**2*dphi*(np.cos(theta_top)-np.cos(theta_bottom)))
    return areas

def grid_surface_integral(grid_lats,grid_azis,grid_values,sphere_radius,degorhour):
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
    dazis = angle_difference(azis[:-1],azis[1:],degorhour=degorhour)
    dazis = np.mod(dazis,2*np.pi/_azifac(degorhour)) #Ensure > 0
    dazi = np.nanmedian(dazis)

    gridcell_bottom_lats = grid_lats-dlat/2.
    gridcell_top_lats = grid_lats+dlat/2.
    gridcell_left_azis = angle_difference(dazi/2.,grid_azis,degorhour=degorhour)
    gridcell_right_azis = grid_azis+dazi/2.
    
    gridcell_areas = great_circle_rectangle_area(gridcell_bottom_lats,
                                                    gridcell_top_lats,
                                                    gridcell_left_azis,
                                                    gridcell_right_azis,
                                                    sphere_radius,
                                                    degorhour)

    integrated_values = np.nansum(gridcell_areas*grid_values)
    return integrated_values