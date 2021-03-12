# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
from numpy import *
import numpy as np
import bisect
import datetime
import logging
import matplotlib
from matplotlib.colors import Normalize, LogNorm
from geospacepy.spherical_geometry import (angle_difference,
                                            angle_midpoint,
                                            great_circle_distance)

def dipole_tilt_angle(dts):
    """
    Computes the dipole tilt angle (in degrees) given a single datetime or array of datetimes
    Approximation from:
    M. Nowada, J.-H. Shue, C.T. Russell, Effects of dipole tilt angle on geomagnetic activity,
    Planetary and Space Science, Volume 57, Issue 11, September 2009, Pages 1254-1259,
    ISSN 0032-0633, http://dx.doi.org/10.1016/j.pss.2009.04.007.
    """
    if isinstance(dts,np.ndarray):
        original_shape = np.shape(dts)
        dts = dts.flatten().tolist()
    else:
        dts = [dts] if not isinstance(dts,list) else dts
        original_shape = (len(dts),)

    phis=[]
    for dt in dts:
        doy = dt.timetuple().tm_yday
        ut_hr = dt.hour+dt.minute/60.+dt.second/60./60.
        #dipole tilt angle due to the time of year
        phi_year = 23.4*np.cos((doy-172.)*2*np.pi/365.25)
        phi_uthr = 11.2*np.cos((ut_hr-16.72)*2*np.pi/24.)
        phi = phi_year+phi_uthr
        phis.append(phi)

    if len(phis)==1:
        return phis[0]
    else:
        return np.array(phis).reshape(original_shape)

def sathat(ut,pos,secondCoord='Longitude',lattype='geocentric',up_is_geodetic=False):
    """
    Estimates the direction of the along and cross track unit vectors of a spacecraft
    in an east, north, up coordinate system.

    PARAMETERS
    ----------
    ut - numpy.ndarray (m rowns, 1 column)
        Spacecraft Timestamps in UT second of day
    pos - numpy.ndarray( m rows, 2 columns )
        Spacecraft locations in Lat,Lon or Lat,Localtime depending on secondCoord
    secondCoord - string, {'Localtime','Longitude'}
        Second coordinate of position arrays
    lattype - string, {'geocentric','geodetic'}
        Which style of latitude to use
    up_is_geodetic = bool
        If True, applies an additional transformation, which is appropriate if the satellite defines 'up'
        as geodetic normal, i.e. normal to the ellipse, and you would like the resulting unit vectors
        to transform the data into Z being radial, i.e. geocentric for you

    RETURNS
    -------
    s_along - numpy.ndarray (m-1 rows, 3 columns )
        Along track unit vector in ENU
    s_across - numpy.ndarray (m-1 rows, 3 columns )
        Across track left unit vector in ENU
    s_up - numpy.ndarray (m-1 rows, 3 columns)
        Geocentric upward unit vector in ENU

    NOTES
    -----
    Assumes no altitude change
    The algorithm is based on spherical trig
    The returned vectors are estimated for time ut2
    I'm aware that Localtime and longitude are not
    always simply related by a factor of 12/180, but
    since we only work with the difference between
    sequential longitudes or local times, this conversion
    is accurate enough.
    """
    import random
    #Define some useful constants and lambdas
    ecc_earth=.081819221456

    #Geocentric to Geodetic
    gc2gd = lambda gclat: np.arctan2(np.tan(gclat/180*np.pi),(1.-ecc_earth**2.))/np.pi*180

    #Geodetic to Geocentric
    gd2gc = lambda gdlat: np.arctan2(np.tan(gdlat/180*np.pi),1./(1.-ecc_earth**2.))/np.pi*180

    if lattype=='geocentric':
        gclat = pos[:,0]
        gdlat = gc2gd(gclat)
    elif lattype=='geodetic':
        gdlat = pos[:,0]
        gclat = gd2gc(gdlat)
    else:
        raise ValueError('%s is not a valid selection for lattype' % (lattype))

    ut1 = ut.flatten()[:-1]
    ut2 = ut.flatten()[1:]
    dt = ut2-ut1

    #Not sure why this problem occurs, but the algorithm
    #gives strange values when the two adjacent points are very close
    #together (i.e. delta_t < 1 second)

    colat1 = 90-gclat[:-1]
    colat2 = 90-gclat[1:]
    if secondCoord.lower()=='localtime':
        lon1 = pos[:-1,1]*180./12.
        lon2 = pos[1:,1]*180./12.
    else:
        lon1 = pos[:-1,1]
        lon2 = pos[1:,1]

    #Correct for earth rotation
    #0.0041780741 degrees per seconds from Wolfram Alpha
    lon1 = lon1 - (ut2-ut1)*0.0041780741  #Diff is forward difference

    arc_a = colat1*pi/180. #colat of point 1 in rad
    arc_b = colat2*pi/180. #colat of point 2 in rad
    ang_C = (lon2-lon1)*pi/180. #angle subtended by arc connecting points 1 and 2

    #Do some spherical trig
    cos_arc_c = cos(arc_a)*cos(arc_b) + sin(arc_a)*sin(arc_b)*cos(ang_C)
    sin_arc_c = sqrt(1.-cos_arc_c**2)
    cos_ang_A = (cos(arc_a)*sin(arc_b) - sin(arc_a)*cos(arc_b)*cos(ang_C))/sin_arc_c
    sin_ang_A = sin(arc_a)*(sin(ang_C)/sin_arc_c)
    cos_ang_B = (cos(arc_b)*sin(arc_a)-sin(arc_b)*cos(arc_a)*cos(ang_C))/sin_arc_c
    sin_ang_B = sin(arc_b)/sin_arc_c*sin(ang_C)

    #pdb.set_trace()
    # We assume a circular cross section and no alitude change, this is now ENU in whatever the original
    # geo(detic,centric) coordinate system was.
    all_zeros = np.zeros_like(sin_ang_A)
    all_ones = np.ones_like(sin_ang_A)

    #At point 1
    s_along_1 = column_stack((sin_ang_B,cos_ang_B,all_zeros)) #Along track
    s_cross_1 = column_stack((-1*cos_ang_B,sin_ang_B,all_zeros)) #Cross track (+ to the left)
    #s_up_1 = column_stack((all_zeros,all_zeros,all_ones))

    #At point 2
    s_along_2 = column_stack((sin_ang_A,-1*cos_ang_A,all_zeros)) #Along track
    s_cross_2 = column_stack((cos_ang_A,sin_ang_A,all_zeros)) #Cross track (+ to the left)
    #s_up_2 = column_stack((all_zeros,all_zeros,all_ones))

    #average
    s_along = (s_along_1+s_along_2)/2
    s_cross = (s_cross_1+s_cross_2)/2
    s_up = column_stack((all_zeros,all_zeros,all_ones))

    #Deal with the missing value from finite difference
    #(arrays will have same length as input time, lat and lon arrays)
    #s_along = np.row_stack((s_along[0,:],s_along))
    #s_cross = np.row_stack((s_cross[0,:],s_cross))
    #s_up = np.row_stack((s_up[0,:],s_up))

    s_along = np.row_stack((s_along,s_along[-1,:]))
    s_cross = np.row_stack((s_cross,s_cross[-1,:]))
    s_up = np.row_stack((s_up,s_up[-1,:]))

    #s_along[:,0] = moving_average(s_along[:,0],4)
    #s_along[:,1] = moving_average(s_along[:,1],4)
    #s_cross[:,0] = moving_average(s_cross[:,0],4)
    #s_cross[:,1] = moving_average(s_cross[:,1],4)

    if up_is_geodetic:
        # Now we apply the additional rotation that will apply the geodetic to geocentric transform
        # This is a rotation about the eastward direction by the difference between the geocentric and geodetic latitudes
        for r in np.arange(len(s_along[:,0])):
            dth = (gdlat[r]-gclat[r])/180*pi
            c = np.cos(dth)
            s = np.sin(dth)
            if r<10:
                print(gdlat[r] - gclat[r])
            rotmat = np.array([[1,   0,  0], [0,   c,  -1*s],[0,s,  c]])
            s = random.randint(0,86400/3)
            if s==1:
                print("Index %d: gdlat=%.3f gclat=%.3f" % (r,gdlat[r],gclat[r]))
                print((str(rotmat)))

            # Now we do the matrix product of the rotation matrix and each row of the along, across, and upward unit vectors

            s_along[r,:] = np.dot(rotmat,s_along[r,:])
            s_cross[r,:] = np.dot(rotmat,s_cross[r,:])
            s_up[r:,] = np.dot(rotmat,s_up[r,:])

    return s_along,s_cross,s_up
    
def simple_passes(latitude,half_or_full_orbit='half'):

    npts = len(latitude.flatten())
    entered_north = []
    entered_south = []

    for k in range(1,npts):
        #poleward crossing
        if latitude[k-1] < 0. and latitude[k] >= 0.:
            entered_north.append(k)
            #print "Entered Northern Hemisphere: ind:%d,lat:%.3f" % (k,latitude[k])
        elif latitude[k-1] > 0. and latitude[k] <= 0.:
            entered_south.append(k)
            #print "Entered Southern Hemisphere: ind:%d,lat:%.3f" % (k,latitude[k])

    if half_or_full_orbit is 'half':
        xings = entered_north+entered_south
    elif half_or_full_orbit is 'full':
        if entered_north[0] < entered_south[0]:
            xings = entered_north
        elif entered_south[0] < entered_north[0]:
            xings = entered_south

    xings.sort()

    return xings

def timepos_ticklabels(ax,t,lat,ltlon,fs=10):
    """
    Make Multi-Line Tick Labels for Spacecraft Data with time, lat and localtime/longitude
    """
    lat,ltlon,t = lat.flatten(),ltlon.flatten(),t.flatten()
    using_datetimes = isinstance(t.tolist()[0],datetime.datetime)
    if using_datetimes:
        #Need to convert datetimes to matplotlib dates if we are listing with datetimes
        mplt = np.array(matplotlib.dates.date2num(t.tolist()))
    else:
        mplt = t
    ticks = ax.get_xticks()
    #print ticks with multiple lines
    newlab=[]
    for tick in ticks:
        ind = (np.abs(mplt-tick)).argmin()
        if using_datetimes:
            newlab.append('%s\n%.1f\n%.1f' % (t[ind].strftime('%H:%M'),lat[ind],ltlon[ind]))
        else:
            newlab.append('%.1f\n%.1f\n%.1f' % (t[ind],lat[ind],ltlon[ind]))
    ax.set_xticklabels(newlab)
    matplotlib.artist.setp(ax.get_xmajorticklabels(),size=fs,rotation=0)

def multiline_timelabels(ax,tdata,xdata,strffmt='%H:%M',xfmt=['%.1f']):
    """Adds additional lines to the labels of an existing axes.
    
    INPUTS
    ------

    tdata : numpy.ndarray
        data that was passed to the plot function, must be an array of datetime objects
    xdata : numpy.ndarray
        any number of columns of additional data to be added to labels. Must have same number of rows as tdata.
    strffmt : str
        format specification to datetime.datetime.strftime for time labels
    xfmt : list
        list of formats for each column of xdata

    RETURNS
    -------

    """
    #Manually create the tick labels
    #There is probably a better way to do this with FuncFormatter, but I couldn't
    #figure out how to get all of the relavent lat and LT information into it
    from matplotlib import dates as mpldates

    #Get the tick marks
    xticks = ax.get_xticks()
    xticks_datetime = array(mpldates.num2date(xticks))
    xlabels = []
    for l in range(len(xticks)):
        tick = xticks_datetime[l]
        tick = tick.replace(tzinfo=None) #Remove the timezone so we can compare the two types
        ind = None
        for k in range(len(tdata)): #Can't get nonzero to work on this??
            if tdata[k] == tick:
                ind = k
        if ind is not None: #Sometimes tick is not found if it wants to tickmark outside of data range
            tickstr = tick.strftime(strffmt)
            if len(xdata.shape)>1: #If more than one column of additional data
                for c in range(len(xdata[0,:])):
                    tickstr+="\n"
                    tickstr+=xfmt[c] % (xdata[ind,c])
            else:
                tickstr+="\n"
                tickstr+=xfmt[0] % (xdata[ind])
            xlabels.append(tickstr)
        else:
            xlabels.append(tick.strftime(strffmt))

    ax.set_xticklabels(xlabels)

    return ax

def draw_dialplot(ax,minlat=50,padding=3,fslt=10,fslat=12,southern_hemi=False,
    draw_circles=True,latlabels=True):
    """Draws the dialplot and labels the latitudes
    
    PARAMETRS
    ---------
    
    ax : matplotlib.axes.Axes
        Axes object to plot on
    minlat : {60,50,40}, optional
        Latitude of largest ring of dialplot
    padding : int, optional
        Amount of extra space to put around the plot ( in plot units )
    fslt : int,optional
        Font size for hour labels
    fslat : int,optional
        Font size for latitude labels
    southern_hemi : bool,optional
        Defaults to False, put negative signs on latitude labels
    draw_circles : bool,optional
        Draw the circles at the labeled latitudes (default: True)
    latlabels : bool,optional
        Draw latitude labels on plot (default: True)

    """
    phi = linspace(0,2*pi,3000)

    ax.figure.set_facecolor('white')
    thecolor = 'grey'
    thelinestyle='solid'
    thezorder = -100 #make sure the lines are in the background of the plot

    if draw_circles:
        #Circles
        if minlat == 60:
            ax.plot(30*cos(phi),30*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
            ax.plot(20*cos(phi),20*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
            ax.plot(10*cos(phi),10*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
        elif minlat == 50:
            ax.plot(40*cos(phi),40*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
            ax.plot(30*cos(phi),30*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
            ax.plot(20*cos(phi),20*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
            ax.plot(10*cos(phi),10*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
        elif minlat == 40:
            ax.plot(50*cos(phi),50*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
            ax.plot(40*cos(phi),40*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
            ax.plot(30*cos(phi),30*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
            ax.plot(20*cos(phi),20*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
            ax.plot(10*cos(phi),10*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)
    else:
        #Just the border
        ax.plot((minlat-10)*cos(phi),(minlat-10)*sin(phi),color=thecolor,linestyle=thelinestyle,zorder=thezorder)

    #Labels
    tcolor = 'red'

    r_text = 90-minlat+1; th_text = 3*pi/2
    x = r_text*cos(th_text)
    y = r_text*sin(th_text)
    ax.text(x,y,'0',fontsize=fslt,color=tcolor,horizontalalignment='center',verticalalignment='top')

    r_text = 90-minlat+4; th_text = 7*pi/4
    x = r_text*cos(th_text)
    y = r_text*sin(th_text)
    ax.text(x,y,'3',fontsize=fslt,color=tcolor,horizontalalignment='center',verticalalignment='center')

    r_text = 90-minlat+1; th_text = 0*pi/2;
    x = r_text*cos(th_text)
    y = r_text*sin(th_text)
    ax.text(x,y,'6',fontsize=fslt,color=tcolor,horizontalalignment='left',verticalalignment='center')

    r_text = 90-minlat+4; th_text = pi/4;
    x = r_text*cos(th_text)
    y = r_text*sin(th_text)
    ax.text(x,y,'9',fontsize=fslt,color=tcolor,horizontalalignment='center',verticalalignment='center')

    r_text = 90-minlat+1; th_text = 1*pi/2;
    x = r_text*cos(th_text)
    y = r_text*sin(th_text)
    ax.text(x,y,'12',fontsize=fslt,color=tcolor,horizontalalignment='center',verticalalignment='bottom')

    r_text = 90-minlat+4; th_text = 3*pi/4;
    x = r_text*cos(th_text)
    y = r_text*sin(th_text)
    ax.text(x,y,'15',fontsize=fslt,color=tcolor,horizontalalignment='center',verticalalignment='center')

    r_text = 90-minlat+1; th_text = 2*pi/2;
    x = r_text*cos(th_text)
    y = r_text*sin(th_text)
    ax.text(x,y,'18',fontsize=fslt,color=tcolor,horizontalalignment='right',verticalalignment='center')

    r_text = 90-minlat+4; th_text = 5*pi/4;
    x = r_text*cos(th_text)
    y = r_text*sin(th_text)
    ax.text(x,y,'21',fontsize=fslt,color=tcolor,horizontalalignment='center',verticalalignment='center')

    # line([0 0],[-50 50],'Color',[0.7 0.7 0.7],'LineWidth',1.5);
    # line([-50 50],[0 0],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    for i in range(1,25):
        th = (i-1)*pi/12;
        r_min = 3;
        r_max = 90-minlat;
        ax.plot([r_min*cos(th), r_max*cos(th)],[r_min*sin(th), r_max*sin(th)],color=thecolor,linestyle=thelinestyle,zorder=thezorder,
             linewidth=1)

    if latlabels:
        sh = r'-' if southern_hemi else ''
        ax.text( 6,-5,sh+r'$80^o$',fontsize=fslat,color=tcolor);
        ax.text(16,-5,sh+r'$70^o$',fontsize=fslat,color=tcolor);
        ax.text(26,-5,sh+r'$60^o$',fontsize=fslat,color=tcolor);

        if minlat < 60:
            ax.text(36,-5,r'$50^o$',fontsize=fslat,color=tcolor);

        if minlat < 50:
            ax.text(46,-5,r'$40^o$',fontsize=fslat,color=tcolor);

    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axis('tight')
    ax.set_xlim([-1*(90.-minlat+padding),(90.-minlat+padding)])
    ax.set_ylim([-1*(90.-minlat+padding),(90.-minlat+padding)])
    return ax

def latlt2polar(lat,lt,hemisphere):
    """
    Converts an array of latitude and lt points to polar for a top-down dialplot (latitude in degrees, LT in hours)
    i.e. makes latitude the radial quantity and MLT the azimuthal

    get the radial displacement (referenced to down from northern pole if we want to do a top down on the north,
        or up from south pole if visa-versa)
    """
    if hemisphere=='N':
        r = 90.-lat
    elif hemisphere=='S':
        r = 90.-(-1*lat)
    else:
        raise ValueError('%s is not a valid hemisphere, N or S, please!' % (hemisphere))
    #convert lt to theta (azimuthal angle) in radians
    theta = lt/24. * 2*pi - pi/2

    #the pi/2 rotates the coordinate system from
    #theta=0 at negative y-axis (local time) to
    #theta=0 at positive x axis (traditional polar coordinates)
    return r,theta

def latlon2polar(lat,lon,hemisphere):
    """
    Converts an array of latitude and lt points to polar for a top-down dialplot (latitude in degrees, LT in hours)
    i.e. makes latitude the radial quantity and MLT the azimuthal
    """
    #Get the radial displacement (referenced to down from northern pole if we want to do a top down on the north,
    #   or up from south pole if visa-versa)
    if hemisphere=='N':
        r = 90.-lat
    elif hemisphere=='S':
        r = 90.-(-1*lat)
    else:
        raise ValueError('%s is not a valid hemisphere, N or S, please!' % (hemisphere))
    #convert lt to theta (azimuthal angle) in radians
    theta = lon/360. * 2*pi - pi/2

    #make sure theta is positive
    theta[theta<0.] = theta[theta<0.]+2*pi

    #the pi/2 rotates the coordinate system from
    #theta=0 at negative y-axis (local time) to
    #theta=0 at positive x axis (traditional polar coordinates)
    return r,theta

def latlt2cart(lat,lt,hemisphere):
    """
    Latitude and local time to cartesian for a top-down dialplot
    """
    r,theta = latlt2polar(lat,lt,hemisphere)
    return r*cos(theta),r*sin(theta)

def latlon2cart(lat,lon,hemisphere):
    """
    Latitude and longitude to cartesian for a top-down dialplot
    """
    r,theta = latlon2polar(lat,lon,hemisphere)
    return r*cos(theta),r*sin(theta)

def hairplot(ax,lat,lt,C,hemisphere,max_size=10,max_val=None,vmin=None,vmax=None,ref_units=None,dialplot=True,horizontal=False,min_displayed=0.,**kwargs):
    """
    Makes top-down polar plots with vertical lines to indicate intensity and color along spacecraft track.
    Can handle either an array of colors or a variable of colors.
    """
    
    #Draw the background dialplot on
    if dialplot:
        draw_dialplot(ax)

    if max_val is None:
        max_val = np.nanmax(np.abs(C))

    if vmax is None:
        vmax = nanmax(C)
    if vmin is None:
        vmin = nanmin(C)

    X,Y = latlt2cart(lat,lt,hemisphere)
    if not horizontal:
        X1 = zeros_like(X)
        Y1 = C/max_val*max_size
    else:
        Y1 = zeros_like(Y)
        X1 = C/max_val*max_size

    #Implement filtering very small values out
    above_min = np.abs(C) > min_displayed

    if 'alpha' not in kwargs:
        kwargs['alpha']=.75

    norm = Normalize(vmin=vmin,vmax=vmax)
    Q = ax.quiver(X[above_min],Y[above_min],X1[above_min],Y1[above_min],C[above_min],
        angles='xy',units='xy',width=.4,scale_units='xy',scale=1,headwidth=0,headlength=0,norm=norm,**kwargs)

    #Q appears to not actually create a mappable??
    mappable = matplotlib.cm.ScalarMappable(norm=Q.norm,cmap=Q.cmap)
    mappable.set_array(C[above_min])

    key_label = str(max_val)+'[%s]' % (ref_units) if ref_units is not None else ''
    if ref_units is not None:
        ax.quiverkey(Q, 0.45, 0, max_size, key_label,color=mappable.to_rgba(vmax),labelpos='E')

    return mappable

def crosstrackplot(ax,lat,lt,vcross,hemisphere,max_size=10,
        max_val=None,vmin=None,vmax=None,ref_units=None,key_label=None,
        dialplot=True,horizontal=False,min_displayed=0.,
        leftorright='left',label_start_end=False,no_quiverkey=False,
        alpha=.75,single_color=None,key_pos=(.05,.05),**kwargs):
    """
    Makes top-down polar plots with lines perpendicular to the
    spacecraft track to indicate intensity and color along spacecraft track.
    Can handle either an array of colors or a variable of colors.

    If single_color is not None, then will not make a colorbar
    """
    #Draw the background dialplot on
    if dialplot:
        draw_dialplot(ax)

    if max_val is None:
        max_val = np.nanmax(np.abs(vcross))

    if vmax is None:
        vmax = nanmax(vcross)
    if vmin is None:
        vmin = nanmin(vcross)

    X,Y = latlt2cart(lat,lt,hemisphere)
    dX,dY = np.diff(X),np.diff(Y)
    #Fix length
    dX,dY=np.concatenate(([dX[0]],dX)),np.concatenate(([dY[0]],dY))

    #Unit Vector Along Track
    alongX = dX/np.sqrt(dX**2+dY**2)
    alongY = dY/np.sqrt(dX**2+dY**2)

    #Cross Track Unit Vector
    acrossX = alongY if hemisphere == 'S' else -1*alongY
    acrossY = -1*alongX if hemisphere == 'S' else alongX

    if leftorright == 'right':
        acrossY=-1*acrossY
        acrossX=-1*acrossX

    Y1 = vcross/max_val*max_size*acrossY
    X1 = vcross/max_val*max_size*acrossX

    #Implement filtering very small values out
    above_min = np.abs(vcross) > min_displayed

    if not single_color:
        norm = Normalize(vmin=vmin,vmax=vmax)
        Q = ax.quiver(X[above_min],Y[above_min],X1[above_min],Y1[above_min],vcross[above_min],
            angles='xy',units='xy',width=.4,scale_units='xy',scale=1,alpha=alpha,headwidth=0,headlength=0,norm=norm,**kwargs)

        #Q appears to not actually create a mappable??
        mappable = matplotlib.cm.ScalarMappable(norm=Q.norm,cmap=Q.cmap)
        mappable.set_array(vcross[above_min])

    else:
        Q = ax.quiver(X[above_min],Y[above_min],X1[above_min],Y1[above_min],
            color=single_color,angles='xy',units='xy',width=.4,scale_units='xy',scale=1,alpha=alpha,headwidth=0,headlength=0,**kwargs)

    if ref_units is not None:
        key_label_pre = str(max_val)+'[%s]' % (ref_units)
        if key_label is not None:
            key_label = key_label_pre+' (%s)' % (key_label)
        else:
            key_label = key_label
        keycolor = mappable.to_rgba(vmax) if single_color is None else single_color
        if not no_quiverkey:
            ax.quiverkey(Q, key_pos[0], key_pos[1], max_size, key_label,color=keycolor,labelpos='E',coordinates='figure')

    return mappable if single_color is None else Q


def vector_plot(ax,data,satname='dmsp',color='blue',latlim=-50.,max_vec_len=12.,max_magnitude=1000.,min_displayed=0.,
    reference_vector_len=500.,reference_vector_label="500nT",labeljustify='left',labeltrack=True,fontsize=8,skip=2,plottime=False,
    timejustify='left',ntimes=1,col5isnorth=True,spacecraft_coords=False,alpha=.55,width=.2,secondCoord='localtime'):
    """
    Makes top-down dialplots of spacecraft tracks and vector data.

    PARAMETERS
    ----------
    
    ax : matplotlib.axes
        Thing we're going to plot on
    data : numpy.ndarray
        n x 5 array of spacecraft data
        column 1 = time (UT sec of day)
        column 2 = latitude (magnetic or otherwise)
        column 3 = localtime (magnetic or local solar)
        column 4 = eastward component of vector data
        column 5 = northward component of vector data
    satname : str, optional
        Text to place at end of spacecraft track
    latlim : float, optional
        Largest ring of dial plot (set to negative to indicate
            that data is for southern hemisphere)
    secondCoord : str, optional
        localtime or longitude
    max_vec_len : float, optional
        Maximum length of vector in plot coordinates (i.e. degrees latitude)
    max_magnitude : float, optional
        Value of sqrt(Vec_east**2+Vec_north**2) that will be associated with
        a vector of length max_vec_len on plot (scaling factor)
    min_displayed : float, optional
        Value of sqrt(Vec_east**2+Vec_north**2) that represents the threshold for
        'noise' level measurements. The code will not display vectors below this
        value to reduce visual clutter.
    reference_vector_len : float, optional
        The size of the reference vector in the lower left of the plot
    reference_vector_label : str, optional
        Label for reference vector
    labeljustify : {'left','right','auto'}, optional
        Where to position the satname at the end of the track
    labeltrack : bool, optional
        Draw the label at the end of track if True
    fontsize : int, optional
        Size of fonts for labels
    skip : int, optional
        Cadence of vectors to plot, i.e. skip=2 plots every other vector, except for the ten largest
    plottime : boolean, optional
        Plot the start time of the pass at the first point
    col2isnorth : boolean, optional
        Assume that the 5th column is northward, if false, assumes it's radial (i.e. equatorward, i.e. Apex d2)

    """
    #if labeltext=='default':
    labeltext=satname

    if latlim > 0:
        inlatlim = data[:,1] > latlim
    elif latlim < 0:
        inlatlim = data[:,1] < latlim
    plot_data = data[inlatlim,:]
    if len(plot_data) == 0:
        print("Warning: vector_plot called with no data in display region")
        return
    #convert to colat
    plot_data[:,1] = 90-abs(plot_data[:,1])
    #convert mlt to radians
    if secondCoord.lower() == 'localtime':
        plot_data[:,2] = plot_data[:,2]/24. * 2*pi - pi/2
    elif secondCoord.lower() == 'longitude':
        plot_data[:,2] = plot_data[:,2]/180. * pi

    #the pi/2 rotates the coordinate system from
    #theta=0 at negative y-axis (local time) to
    #theta=0 at positive x axis (traditional polar coordinates)

    #calculate the vector magnitudes
    magnitudes = sqrt(plot_data[:,3]**2+plot_data[:,4]**2)

    #calculate the scaling factors (the percentage of the maximum length each vector will be
    #maximum length corresponds to a magnitude equal to data_range(2)

    sfactors = (magnitudes)/(max_magnitude)

    #normalize so each vector has unit magnitude
    plot_data[:,3] = plot_data[:,3]/magnitudes;
    plot_data[:,4] = plot_data[:,4]/magnitudes;

    #stretch all vectors to the maximum vector length and then scale them
    #by each's individual scale factor
    plot_data[:,3] = (plot_data[:,3]*sfactors)*max_vec_len;
    plot_data[:,4] = (plot_data[:,4]*sfactors)*max_vec_len;

    #finally rotate the coordinate system for the datavar from E_hat,N_hat to r_hat,theta_hat

    #first if in the northern hemisphere, and we're using a coordinate system
    #that eastward northward (i.e. GEO), instead of eastward equatorward (i.e. Apex)
    #flip the sign of the N_hat component to be radially outward (away from the pole)

    if(latlim > 0) and col5isnorth:
        plot_data[:,4] = -1*plot_data[:,4]

    X = plot_data[:,1]*cos(plot_data[:,2])
    Y = plot_data[:,1]*sin(plot_data[:,2])
    r_hat = column_stack((cos(plot_data[:,2]),sin(plot_data[:,2])))
    th_hat = column_stack((-1*sin(plot_data[:,2]),cos(plot_data[:,2])))
    X1 = plot_data[:,4]*r_hat[:,0]+plot_data[:,3]*th_hat[:,0]
    Y1 = plot_data[:,4]*r_hat[:,1]+plot_data[:,3]*th_hat[:,1]

    #Make a mask that will remove any values that are below the minimum
    #magnitude to display
    g = magnitudes>min_displayed

    #Set up an order so the largest values get plotted first
    C = argsort(magnitudes[g])[::-1]

    #Keep the largest 5% of data and then use the skipping
    #to thin the vectors for speed and reduced clutter
    #n_largest_to_keep = ceil(.05*len(magnitudes))
    #Just keep the 10 largest magntiude vectors
    C = concatenate((C[:10],C[list(range(10,len(C),skip))]))

    ax.hold(True)
    ax.quiver(X[C],Y[C],X1[C],Y1[C],angles='xy',units='xy',width=width,color=color,scale_units='xy',scale=1,label=labeltext,alpha=alpha)
    if labeltrack:
        if labeljustify=='right':
            ax.text(X[-1],Y[-1],satname,color=color,va='top',ha='right',fontsize=fontsize)
        elif labeljustify=='left':
            ax.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)
        elif labeljustify=='auto':
            if len(X) > 1:
                if (X[-1] > 0 and X[-2] < X[-1]) or (X[-1] < 0 and X[-2] > X[-1]):
                    ax.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)
                elif (X[-1] > 0 and X[-2] > X[-1]) or (X[-1] < 0 and X[-2] < X[-1]):
                    ax.text(X[-1],Y[-1],satname,color=color,va='top',ha='right',fontsize=fontsize)
                else:
                    ax.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)
            else:
                ax.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)
    #plot the reference arrow
    ref_sfactor = reference_vector_len/max_magnitude
    ref_X1 = (sqrt(reference_vector_len**2/2)/reference_vector_len*ref_sfactor)*max_vec_len
    ref_Y1 = (sqrt(reference_vector_len**2/2)/reference_vector_len*ref_sfactor)*max_vec_len
    ax.quiver(-abs(latlim)+15,-abs(latlim)+5,ref_X1,ref_Y1,angles='xy',units='xy',width=.2,color='black',scale_units='xy',scale=1,label=reference_vector_label)
    ax.text(-abs(latlim)+15,-abs(latlim)+5,reference_vector_label,color='black',va='top',size=8)
    if plottime:
        for d in linspace(0,len(plot_data)-1,ntimes).tolist():
            t = datetime.datetime(2000,1,1)+datetime.timedelta(seconds=plot_data[d,0])
            ax.text(X[d],Y[d],t.strftime('%X'),color=color,va='top',ha=timejustify,fontsize=fontsize,alpha=.75)
    return ax

def vector_component_plot(ax_e,ax_n,data,satname='dmsp',color='blue',latlim=-50.,max_vec_len=12.,max_magnitude=1000.,reference_vector_len=500.,
    reference_vector_label="500nT",labeljustify='left',labeltrack=True,fontsize=8,skip=2,cmap='bwr',plottime=False,timejustify='left',ntimes=1):
    """
    Makes top-down dialplots of spacecraft tracks and vector data.

    PARAMETERS
    ----------
        ax_e : matplotlib.axes
            Thing we're going to plot eastward component of vector on
        ax_n : matplotlib.axes
            Thing we're going to plot eastward component of vector on
        data : numpy.ndarray
            n x 5 array of spacecraft data
            column 1 = time (UT sec of day)
            column 2 = latitude (magnetic or otherwise)
            column 3 = localtime (magnetic or local solar)
            column 4 = eastward component of vector data
            column 5 = northward component of vector data
        satname : str, optional
            Text to place at end of spacecraft track
        latlim : float, optional
            Largest ring of dial plot (set to negative to indicate
                that data is for southern hemisphere)
        max_vec_len : float, optional
            Maximum length of vector in plot coordinates (i.e. degrees latitude)
        max_magnitude : float, optional
            Value of sqrt(Vec_east**2+Vec_north**2) that will be associated with
            a vector of length max_vec_len on plot (scaling factor)
        reference_vector_len : float, optional
            The size of the reference vector in the lower left of the plot
        reference_vector_label : str, optional
            Label for reference vector
        labeljustify : {'left','right','auto'}, optional
            Where to position the satname at the end of the track
        labeltrack : bool, optional
            Draw the label at the end of track if True
        fontsize : int, optional
            Size of fonts for labels
        skip : int, optional
            Cadence of vectors to plot, i.e. skip=2 plots every other vector, except for the ten largest
        plottime : boolean, optional
            Plot the start time of the pass at the first point

    """
    #if labeltext=='default':
    labeltext=satname

    if latlim > 0:
        inlatlim = data[:,1] > latlim
    elif latlim < 0:
        inlatlim = data[:,1] < latlim
    plot_data = data[inlatlim,:]
    if len(plot_data) == 0:
        print("Warning: vector_plot called with no data in display region")
        return
    #convert to colat
    plot_data[:,1] = 90-abs(plot_data[:,1])
    #convert mlt to radians
    plot_data[:,2] = plot_data[:,2]/24. * 2*pi - pi/2
    #the pi/2 rotates the coordinate system from
    #theta=0 at negative y-axis (local time) to
    #theta=0 at positive x axis (traditional polar coordinates)

    #calculate the vector magnitudes
    magnitudes = sqrt(plot_data[:,3]**2+plot_data[:,4]**2);

    #calculate the scaling factors (the percentage of the maximum length each vector will be
    #maximum length corresponds to a magnitude equal to data_range(2)

    sfactors = (magnitudes)/(max_magnitude)

    #normalize so each vector has unit magnitude
    scaled_data_e = plot_data[:,3]/magnitudes;
    scaled_data_n = plot_data[:,4]/magnitudes;

    #stretch all vectors to the maximum vector length and then scale them
    #by each's individual scale factor
    scaled_data_e = (scaled_data_e*sfactors)*max_vec_len;
    scaled_data_n = (scaled_data_n*sfactors)*max_vec_len;

    #finally rotate the coordinate system for the datavar from E_hat,N_hat to r_hat,theta_hat

    #first if in the northern hemisphere, flip the sign of the N_hat component to be radially outward (away from the pole)

    if(latlim > 0):
        scaled_data_n = -1*scaled_data_n

    X = plot_data[:,1]*cos(plot_data[:,2])
    Y = plot_data[:,1]*sin(plot_data[:,2])
    r_hat = column_stack((cos(plot_data[:,2]),sin(plot_data[:,2])))
    th_hat = column_stack((-1*sin(plot_data[:,2]),cos(plot_data[:,2])))


    Y_E = ones_like(X)*scaled_data_e
    X_E = zeros_like(X)

    Y_N = ones_like(X)*scaled_data_n
    X_N = zeros_like(X)

    #Set up an order so the largest values get plotted first
    C_N = argsort(scaled_data_n)[::-1]
    C_E = argsort(scaled_data_e)[::-1]

    C_N = concatenate((C_N[:10],C_N[list(range(10,len(C_N),skip))]))
    C_E = concatenate((C_E[:10],C_E[list(range(10,len(C_E),skip))]))

    #Use the original component data as the color
    ax_e.quiver(X[C_E],Y[C_E],X_E[C_E],Y_E[C_E],plot_data[C_E,3],angles='xy',units='xy',width=.2,
        clim=[-.5*max_magnitude,.5*max_magnitude],scale_units='xy',scale=1,label=labeltext,alpha=.5,cmap=cmap)
    ax_n.quiver(X[C_N],Y[C_N],X_N[C_N],Y_N[C_N],plot_data[C_N,4],angles='xy',units='xy',width=.2,
        clim=[-.5*max_magnitude,.5*max_magnitude],scale_units='xy',scale=1,label=labeltext,alpha=.5,cmap=cmap)

    if labeltrack:
        if labeljustify=='right':
            ax_e.text(X[-1],Y[-1],satname,color=color,va='top',ha='right',fontsize=fontsize)
            ax_n.text(X[-1],Y[-1],satname,color=color,va='top',ha='right',fontsize=fontsize)

        elif labeljustify=='left':
            ax_e.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)
            ax_n.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)

        elif labeljustify=='auto':
            if len(X) > 1:
                if (X[-1] > 0 and X[-2] < X[-1]) or (X[-1] < 0 and X[-2] > X[-1]):
                    ax_e.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)
                    ax_n.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)

                elif (X[-1] > 0 and X[-2] > X[-1]) or (X[-1] < 0 and X[-2] < X[-1]):
                    ax_e.text(X[-1],Y[-1],satname,color=color,va='top',ha='right',fontsize=fontsize)
                    ax_n.text(X[-1],Y[-1],satname,color=color,va='top',ha='right',fontsize=fontsize)
                else:
                    ax_e.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)
                    ax_n.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)

            else:
                ax_e.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)
                ax_n.text(X[-1],Y[-1],satname,color=color,va='top',ha='left',fontsize=fontsize)

    #plot the reference arrow
    ref_sfactor = reference_vector_len/max_magnitude
    ref_X1 = (sqrt(reference_vector_len**2/2)/reference_vector_len*ref_sfactor)*max_vec_len
    ref_Y1 = (sqrt(reference_vector_len**2/2)/reference_vector_len*ref_sfactor)*max_vec_len
    ax_e.quiver(-abs(latlim)+15,-abs(latlim)+5,ref_X1,ref_Y1,angles='xy',units='xy',width=.2,color='black',scale_units='xy',scale=1,label=reference_vector_label)
    ax_e.text(-abs(latlim)+15,-abs(latlim)+5,reference_vector_label,color='black',va='top',size=8)
    if plottime:
        for d in linspace(0,len(plot_data)-1,ntimes).tolist():
            t = datetime.datetime(2000,1,1)+datetime.timedelta(seconds=plot_data[d,0])
            ax_e.text(X[d],Y[d],t.strftime('%X'),color=color,va='top',ha=timejustify,fontsize=fontsize,alpha=.75)
            ax_n.text(X[d],Y[d],t.strftime('%X'),color=color,va='top',ha=timejustify,fontsize=fontsize,alpha=.75)
    return ax_e,ax_n

def cubic_bez_arc(lat,azi1,azi2,lonorlt='lt'):
    """Returns the control point locations for a cubic bezier curve approximation of the arc between
        azi1 and azi2 at a radius of 90-abs(lat)"""
    #From: http://hansmuller-flex.blogspot.com/2011/04/approximating-circular-arc-with-cubic.html
    #The derivation of the control points for an arc of less than 90 degrees is a little more complicated.
    #If the arc is centered around the X axis, then the length of the tangent line is r * tan(a/2),
    #instead of just r.   The magnitude of the vector from each arc endpoint to its control point is k * r * tan(a/2).
    k = 4/3*(np.sqrt(2)-1) #Magic number for bezier curve arc approximations

    azi2rad = np.pi/12. if lonorlt=='lt' else np.pi/180.
    maxazi = 24. if lonorlt=='lt' else 360.

    r = 90.-np.abs(lat)

    if lonorlt=='lt':
        aziunit = 'hour' 
    elif lonorlt=='lon':
        aziunit = 'deg'
    else:
        raise ValueError('Invalid lonorlt {}'.format(lonorlt))
        
    theta = angle_difference(azi1,azi2,aziunit)
    midpoint = azi1+theta/2.
    tangent_len = r*np.tan(theta*azi2rad/2) # Length of tangent line
    r_cp = np.sqrt(r**2+(tangent_len/2)**2)
    th_cp = theta*k/2

    cp_lat = (90.-r_cp)*np.sign(lat)
    cp1_lonorlt = midpoint-th_cp
    cp2_lonorlt = midpoint+th_cp
    return cp_lat,cp1_lonorlt,cp2_lonorlt

def polarbinplot(ax,bin_edges,bin_values,hemisphere='N',lonorlt='lt',**kwargs):
    """
        Plots a collection of bins in polar coordinates on a dialplot
        
    INPUTS
    ------

    bin_edges : numpy.ndarray
        Must have shape [n x 4] with columns:
        bin_lat_start,bin_lat_end,bin_lonlt_start,bin_lonlt_end
    lonorlt : 'lon' or 'lt'
        Use longitude or localtime as azimuthal coordinate
    
    RETURNS
    -------

    mappable : matplotlib.collections.PatchCollection
        The collection of bin shapes (which can be fed to colorbar)

    """
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    from matplotlib.collections import PatchCollection

    #Compute number of bins
    nbins = len(bin_edges[:,0])

    azifun = latlon2cart if lonorlt=='lon' else latlt2cart
    azi2rad = np.pi/12. if lonorlt=='lt' else np.pi/180.

    #Get the color scale extents
    vmax = nanmax(bin_values) if 'vmax' not in kwargs else kwargs['vmax']
    vmin = nanmin(bin_values) if 'vmin' not in kwargs else kwargs['vmin']

    logscale = False
    if not logscale:
        norm = Normalize(vmin=vmin,vmax=vmax)
    else:
        norm = LogNorm(vmin=vmin,vmax=vmax)

    #mappable = matplotlib.cm.ScalarMappable(norm=norm,cmap=)
    #mappable.set_array(bin_values)
    hemifac = -1. if hemisphere == 'S' else 1.
    #Control point for 4 point bezier curves
    #Will put control points at the edge of the plot and at mlt1 and mlt2,
    #so that the curve approximates a circular arc

    cp12_lat,cp12_lonorlt1,cp12_lonorlt2 = cubic_bez_arc(bin_edges[:,0],bin_edges[:,2],bin_edges[:,3],lonorlt=lonorlt)
    cp34_lat,cp34_lonorlt1,cp34_lonorlt2 = cubic_bez_arc(bin_edges[:,1],bin_edges[:,3],bin_edges[:,2],lonorlt=lonorlt)

    #for n in range(5):
    #   print "lat: "+str(bin_edges[n,0])+" lt1: "+str(bin_edges[n,2])+" lt2: "+str(bin_edges[n,3])
    #   print "lat_cp: "+str(cp12_lat[n])+" lt1_cp: "+str(cp12_lonorlt1[n])+" lt2_cp: "+str(cp12_lonorlt2[n])

    X1,Y1 = azifun(bin_edges[:,0],bin_edges[:,2],hemisphere)
    X12c1,Y12c1 = azifun(cp12_lat,cp12_lonorlt1,hemisphere)
    X12c2,Y12c2 = azifun(cp12_lat,cp12_lonorlt2,hemisphere)
    X2,Y2 = azifun(bin_edges[:,0],bin_edges[:,3],hemisphere)
    X3,Y3 = azifun(bin_edges[:,1],bin_edges[:,3],hemisphere)
    X34c1,Y34c1 = azifun(cp34_lat,cp34_lonorlt1,hemisphere)
    X34c2,Y34c2 = azifun(cp34_lat,cp34_lonorlt2,hemisphere)
    X4,Y4 = azifun(bin_edges[:,1],bin_edges[:,2],hemisphere)

    control_codes = [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,
                    Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO]

    patches = [] #Path patches which make up the plot
    color_arr = []
    for ib in range(nbins):
        if np.isfinite(bin_values[ib]):
            verticies = [(X1[ib],Y1[ib]),(X12c1[ib],Y12c1[ib]),(X12c2[ib],Y12c2[ib]),(X2[ib],Y2[ib]),
                            (X3[ib],Y3[ib]),(X34c1[ib],Y34c1[ib]),(X34c2[ib],Y34c2[ib]),(X4[ib],Y4[ib]),(X1[ib],Y1[ib])]

            patches.append(PathPatch(Path(verticies,control_codes))) #Bool argument is explicitly closed shape
            color_arr.append(bin_values[ib])
    #Now make it into a Collection (which is a subclass of ScalarMappable)
    mappable = PatchCollection(patches,
        cmap=matplotlib.cm.jet if 'cmap' not in kwargs else kwargs['cmap'],
        alpha=.9 if 'alpha' not in kwargs else kwargs['alpha'],
        edgecolor="None" if 'edgecolor' not in kwargs else kwargs['edgecolor'],
        norm=norm)

    mappable.set_array(np.array(color_arr).flatten())

    ax.add_collection(mappable)

    return mappable

def polarbin_vectorplot(ax,bin_edges,bin_values_E,bin_values_N,
                            hemisphere='N',lonorlt='lt',color='black',
                            max_vec_len= 12.,max_magnitude= 1000.,
                            reference_vector_len= 500.,
                            reference_vector_label= "500nT",
                            alpha=.55,width=.35,zorder=10):

    #Find bin center lats and localtimes
    lats = angle_midpoint(bin_edges[:,0],bin_edges[:,1],'deg')
    latlim = 50.

    if lonorlt == 'lt':
        degorhour_azi = 'hour'
    elif lonorlt == 'lon':
        degorhour_azi = 'deg'
    else:
        raise ValueError('Invalid lonorlt %s' % (lonorlt))

    azis = angle_midpoint(bin_edges[:,2],bin_edges[:,3],degorhour_azi)

    azi2rad = np.pi/12. if lonorlt=='lt' else np.pi/180.

    r = 90.-np.abs(lats).flatten()
    theta = azis.flatten()*azi2rad-np.pi/2

    #calculate the vector magnitudes
    magnitudes = sqrt(bin_values_E**2+bin_values_N**2).flatten()

    #calculate the scaling factors (the percentage of the maximum length each vector will be
    #maximum length corresponds to a magnitude equal to data_range(2)

    sfactors = (magnitudes)/(max_magnitude)

    #normalize so each vector has unit magnitude
    scaled_data_e = bin_values_E.flatten()/magnitudes;
    scaled_data_n = bin_values_N.flatten()/magnitudes;

    #stretch all vectors to the maximum vector length and then scale them
    #by each's individual scale factor
    scaled_data_e = (scaled_data_e*sfactors)*max_vec_len;
    scaled_data_n = (scaled_data_n*sfactors)*max_vec_len;

    if hemisphere=='N':
        scaled_data_n = -1*scaled_data_n
    elif hemisphere=='S':
        scaled_data_n = scaled_data_n
    else:
        raise ValueError('Invalid Hemisphere %s (N or S)' % (hemisphere))

    X = r*np.cos(theta)
    Y = r*np.sin(theta)

    r_hat = column_stack((np.cos(theta),sin(theta)))
    th_hat = column_stack((-1*np.sin(theta),np.cos(theta)))

    X1 = scaled_data_n*r_hat[:,0]+scaled_data_e*th_hat[:,0]
    Y1 = scaled_data_n*r_hat[:,1]+scaled_data_e*th_hat[:,1]

    print('X',X.shape,'Y',Y.shape,'X1',X1.shape,'Y1',Y1.shape)

    ax.quiver(X,Y,X1,Y1,angles='xy',units='xy',
                width=width,color=color,
                scale_units='xy',scale=1,alpha=alpha,
                zorder=zorder)

    #plot the reference arrow
    ref_sfactor = float(reference_vector_len)/float(max_magnitude)
    ref_X1 = (sqrt(reference_vector_len**2/2)/reference_vector_len*ref_sfactor)*max_vec_len
    ref_Y1 = (sqrt(reference_vector_len**2/2)/reference_vector_len*ref_sfactor)*max_vec_len
    ax.quiver(-abs(latlim)+15,-abs(latlim)+10,ref_X1,ref_Y1,
                angles='xy',units='xy',width=width,color=color,
                scale_units='xy',scale=1,label=reference_vector_label,
                zorder=0.)
    ax.text(-abs(latlim)+15,-abs(latlim)+10,reference_vector_label,
                color='black',va='top',size=8,bbox={'alpha':0.})


def rolling_window(a, window):
    """Make for the lack of a decent moving average"""
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def moving_average(x,window_size):
    """Creates a weighted average smoothed version of x using the weights in window"""
    return np.nanmean(rolling_window(np.concatenate((x[:window_size//2],x,x[-window_size//2+1:])),window_size),-1)

def moving_median(x,window_size):
    """Creates a weighted average smoothed version of x using the weights in window"""
    return np.nanmedian(rolling_window(np.concatenate((x[:window_size//2],x,x[-window_size//2+1:])),window_size),-1)


