"""
astrodynamics2.py - python library of astrodynamical functions for ASEN 5050
Author - Liam Kilcommons
"""
from numpy import *
import numpy as np
import matplotlib.pyplot as pp
import itertools, datetime

G = 6.67e-11 #N m^2/s^2
m_earth = 5.9742e24 #kg
r_earth = 6371200 #m
mu = G*m_earth
#Cartesian Unit Vectors
I = array([1.,0.,0.])
J = array([0.,1.,0.])
K = array([0.,0.,1.])

def rot1(angle,vec,deg=False):
    #Angle in radians unless deg=True
    if deg:
        angle = angle*pi/180.
    c = cos(angle)
    s = sin(angle)
    rotmat = array([[1,   0,  0],
                   [0,   c,  s],
                   [0,-1*s,  c]])
    rotvec = dot(rotmat,vec.reshape((-1,1)))
    return rotvec.reshape(vec.shape)

def rot2(angle,vec,deg=False):
    #Angle in radians unless deg=True
    if deg:
        angle = angle*pi/180.
    c = cos(angle)
    s = sin(angle)
    rotmat = array([[c,   0,-1*s],
                   [0,   1,   0],
                   [s,   0,   c]])
    rotvec = dot(rotmat,vec.reshape((-1,1)))
    return rotvec.reshape(vec.shape)

def rot3(angle,vec,deg=False):
    #Angle in radians unless deg=True
    if deg:
        angle = angle*pi/180.
    c = cos(angle)
    s = sin(angle)
    rotmat = array([[   c,   s,  0],
                   [-1*s,   c,  0],
                   [   0,  0,  1]])
    rotvec = dot(rotmat,vec.reshape((-1,1)))
    return rotvec.reshape(vec.shape)

def rot_tests():
    for func in [rot1,rot2,rot3]:
        for vec in [array([1,0,0]),array([[1],[0],[0]]),array([1,0,0]).flatten()]:
            print("Applying %s with angle=pi/2 to %s" % (func.__name__,str(vec)))
            print("Result: %s" % (func(pi/2,vec)))

def convert(input,type,inputUnits,outputUnits):
    #Create a dictionary of conversion factors
    length_systems = ['earth radii','km','m']
    time_systems = ['']

def orbitPlotter(ecc,p=nan,a=nan,inputUnits='earth radii',step=.01,planetaryRadius=1.):
    #ecc = eccentricity
    #p = semiparameter
    #a = semimajor axis
    #nu = true anomoly

    #Parse input (nu)
    if isnan(p) and isnan(a):
        raise ValueError('Please specifiy either: p, the semiparameter, or a, the semimajor axis')
    elif isnan(p) and not isnan(a):
        p = a*(1-ecc**2)
    elif isnan(a) and not isnan(p):
        a = p*(1-ecc**2)

    nu = arange(0,2*pi,step)
    r = p/(1+ecc*cos(nu))

    #convert to cartesian
    x = r*cos(nu)
    y = r*sin(nu)

    planet_x = planetaryRadius*cos(nu)
    planet_y = planetaryRadius*sin(nu)

    fig = pp.figure()
    ax = pp.axes(aspect='equal')
    ax.plot(x,y,'b-')
    ax.hold(True)
    ax.plot(planet_x,planet_y,'g-')
    ax.set_xlabel(inputUnits)
    ax.set_title('Trajectory Plot: eccentricity=%.2f, semiparameter=%.2f, semimajor=%.2f [%s]' % (ecc,p,a,inputUnits))

    return fig

def truetoeccentric(nu,ecc,a=nan,b=nan,tolerence=.00001):
    #Convert true anomally in degrees to eccentric anomally in degress
    #a and b are unit independent
    nu = nu*pi/180.
    if ~isnan(a) and isnan(b):
        b = a*sqrt(1-ecc**2)
    elif ~isnan(b) and isnan(a):
        a = b/sqrt(1-ecc**2)

    p = b**2/a
    r = p/(1+ecc*cos(nu))

    Efroma = arccos((r*cos(nu)+a*ecc)/a)
    Efromb = arcsin(r*sin(nu)/b)
    Efromecc = 2*arctan(sqrt((1-ecc)/(1+ecc))*tan(nu/2))

    if abs(Efroma-Efromb) > tolerence:
        print("Warning: Eccentric anomally from semimajor (cosine) is not within %f rad of to Eccentric anomally from semiminor (sine)" %(tolerence))
    if abs(Efroma-Efromecc) > tolerence:
        print("Warning: Eccentric anomally from semimajor (cosine) is not within %f rad of to Eccentric anomally from eccentricity (tangent)" %(tolerence))
    if abs(Efromb-Efromecc) > tolerence:
        print("Warning: Eccentric anomally from semiminor (cosine) is not within %f rad of to Eccentric anomally from eccentricity (tangent)" %(tolerence))

    return Efroma*180/pi, Efromb*180/pi, Efromecc*180/pi

def eccentrictotrue(E,ecc,a=nan,b=nan,tolerence=.00001):
    #Convert eccentric anomally in degrees to true anomally in degrees
    #takes semimajor and semiminor axes in earth radii
    #a and b are unit independent
    E = E*pi/180.
    if ~isnan(a) and isnan(b):
        b = a*sqrt(1-ecc**2)
    elif ~isnan(b) and isnan(a):
        a = b/sqrt(1-ecc**2)

    r = a*(1-ecc*cos(E))
    nufroma = arccos((a*cos(E)-a*ecc)/r)
    nufromb = arcsin(b*sin(E)/r)
    nufromecc = 2*arctan(sqrt((1+ecc)/(1-ecc))*tan(E/2))

    if abs(nufroma-nufromb) > tolerence:
        print("Warning: True anomally from semimajor (cosine) is not within %f rad \n of Eccentric anomally from semiminor (sine)" %(tolerence))
    if abs(nufroma-nufromecc) > tolerence:
        print("Warning: True anomally from semimajor (cosine) is not within %f rad \n of Eccentric anomally from eccentricity (tangent)" %(tolerence))
    if abs(nufromb-nufromecc) > tolerence:
        print("Warning: True anomally from semiminor (cosine) is not within %f rad \n of Eccentric anomally from eccentricity (tangent)" %(tolerence))

    return nufroma*180/pi, nufromb*180/pi, nufromecc*180/pi

def kepler(ecc,a,E=nan,M=nan,tminustp=nan,tolerence=.001,dist_units="ER"):
    #ecc is eccentricity
    #a is semi-major axis in earth radii
    #nu is true anomally in degrees
    #E is eccentric anomally in degrees
    #M is mean anomally in degrees
    #tminustp is time since periapse in seconds
    #Returns (E,M,tminustp)

    #Convert Units
    if dist_units == "ER":
        a = a*r_earth #ER to meters
    elif dist_units == "m":
        a = a
    elif dist_units == "km":
        a = a*1000.
    else:
        raise ValueError("Invalid dist_units value: %s, valid options are ER,m or km" % (dist_units))

    if ~isnan(E):
        E = E*pi/180. #Radians
    if ~isnan(M):
        M = M*pi/180. #Radians

    #Compute mean motion
    n = sqrt(mu/a**3)

    if any(~isnan([E,M,tminustp])):

        if isnan(E):
            if isnan(M) and not isnan(tminustp):
            #Solve for M using tminustp via M = n(t-t_p)
                M = n*tminustp

            elif isnan(tminustp) and not isnan(M):
                tminustp = M/n

            #Now we have M and tminustp so we can solve for E using newton-raphson

            #Use Algorithm 2 in Vallado to guess for E


            if (M > -1*pi and M < 0) or M > pi:
                guessE = M-ecc
            else:
                guessE=M+ecc

            E = newtonraphsonkepler(ecc,M,guessE)
        else:
            M = E-ecc*sin(E)
            tminustp = M/n
        return E*180/pi,M*180/pi,tminustp
    else:
        raise ValueError('Must specify either M, E, or tminustp to solve keplers equation')
        return (nan,nan,nan)

def between_minus_pi_and_pi(angle,inunit='radians'):
    if inunit in ['Degrees','deg','degrees']:
        angle = angle*pi/180.
    if angle > 2*pi:
        angle = mod(angle,2*pi)
    if angle > pi:
        angle = angle-2*pi
    return angle

def newtonraphsonkepler(ecc,M,guess,tolerence=.001):
    delta=1000.
    num=1.
    Eprev = guess
    while delta>tolerence:
        Enext = Eprev + (M-Eprev+ecc*sin(Eprev))/(1-ecc*cos(Eprev))
        delta = abs(Eprev-Enext)
        print("Iteration %d: E=%.10f, delta=%.10f" % (num,Enext,delta))
        num+=1
        Eprev = Enext
    return Enext

def rv2coe(Rijk,Vijk,debug=True):
    #Units of R and V are km and km/sec respectively
    mu_km =  mu/(1000**3)
    r = linalg.norm(Rijk)
    v = linalg.norm(Vijk)
    if debug:
        print("|R|: %f" % (r))
        print("|V|: %f" % (v))

    a = (2./r-v**2/mu_km)**-1.
    ecc_vec = ((v**2-mu_km/r)*Rijk - dot(Rijk,Vijk)*Vijk)/mu_km
    ecc = linalg.norm(ecc_vec)
    if debug:
        print("mu_km: %f" % (mu_km))
        print("semimajor: %f" %(a))
        print("ecc: %f" % (ecc))

    #Angular Momentum
    h_vec = cross(Rijk,Vijk)
    h = linalg.norm(h_vec)
    if debug:
        print("angular mom: %f" % (h))
        print("angular mom vec: [%f,%f,%f]" % (h_vec[0],h_vec[1],h_vec[2]))

    #Inclination
    inc = arccos(dot(K,h_vec)/(linalg.norm(K)*h))
    if debug:
        print("inclination: %f" % (inc))

    #Right Ascention of Ascending Node
    n_vec = cross(K,h_vec) #node vector
    n = linalg.norm(n_vec)
    Omega = arccos(dot(I,h_vec)/(linalg.norm(I)*h))
    if n_vec[1] < 0.:
        Omega = 2*pi-Omega
    if debug:
        print("n_vec [%f,%f,%f]" % (n_vec[0],n_vec[1],n_vec[2]))
        print("n: %f" % (n))
        print("Omega: %f" %(Omega))

    #Argument of periapse
    w = arccos(dot(n_vec,ecc_vec)/(n*ecc))
    if ecc_vec[2] < 0.:
        w = 2*pi-w

    #True Anomaly
    nu = arccos(dot(ecc_vec,Rijk)/(linalg.norm(ecc_vec)*linalg.norm(Rijk)))
    if dot(Rijk,Vijk) < 0.:
        nu = 2*pi - nu

    #convert all angle to degrees
    inc = inc*180/pi
    Omega = Omega*180/pi
    w = w*180/pi
    nu = nu*180/pi

    return a,ecc,inc,Omega,w,nu

def readTLE(line1,line2,convertMeanMotion=True):
    card1 = int(line1[0])
    #1 blank
    satnum_1 = int(line1[2:6])
    satclass = line1[7]
    #8 blank
    international_designator = line1[9:16].strip()
    id_yr = int(line1[9:10])
    id_lchan_num = int(line1[11:13])
    id_piece = line1[14:16].strip()
    #17 blank
    epoch = float(line1[18:31])
    epoch_yr = int(line1[18:19])
    if epoch_yr < 50:
        epoch_yr = epoch_yr+2000.
    else:
        epoch_yr = epoch_yr+1900.
    epoch_day = float(line1[20:31])

    satnum_2 = int(line2[2:6])
    if satnum_1 != satnum_2:
        raise ValueError("Satellite Numbers not agree between TLE line 1 (%d) and TLE line 2 (%d)!" % (satnum_1,satnum_2))
    i = float(line2[8:15])
    RAAN = float(line2[17:24]) # Right Ascension of Ascending Node [deg]
    ecc = float("0."+line2[26:32]) # Eccentricity
    w = float(line2[34:41]) #Argument of Perigee [deg]
    M = float(line2[43:50]) #Mean Anomally [deg]
    n = float(line2[52:62]) #Mean Motion [rev/day]
    if convertMeanMotion:
        n = n*2*pi/86400. #Rev per day to rad per second
    revnum = float(line2[63:67]) #Revolution number at epoch [revs]

    return i,ecc,RAAN,w,M,n,epoch_yr,epoch_day

def coe2rv(a,ecc,i,Omega,w,nu,debug=True):
    #All distances in km, all angles in degrees
    #Follows Vallado 4th ed. pg. 125
    mu_km =  mu/(1000**3)

    #All angles to radians
    i = i*pi/180
    Omega = Omega*pi/180
    w = w*pi/180
    nu = nu*pi/180

    #Compute semiparameter
    p = a*(1-ecc**2)

    #Vectors in Perifocal frame
    Rpqw = array([p*cos(nu)/(1+ecc*cos(nu)),
                    p*sin(nu)/(1+ecc*cos(nu)),
                    0.])
    alpha = sqrt(mu_km/p)
    Vpqw = array([-1*alpha*sin(nu),alpha*(ecc+cos(nu)),0.])

    if debug:
        print("Perifocal R (R_pqw): [%f,%f,%f]" % (Rpqw[0],Rpqw[1],Rpqw[2]))
        print("Perifocal V (V_pqw): [%f,%f,%f]" % (Vpqw[0],Vpqw[1],Vpqw[2]))
    Rijk = rot3(-1*Omega,rot1(-1*i,rot3(-1*w,Rpqw)))
    Vijk = rot3(-1*Omega,rot1(-1*i,rot3(-1*w,Vpqw)))

    return Rijk,Vijk

#FUNCTIONS FROM HOMEWORK 4
#--------------------------
#Define a basic eci to ecef function
#I'll have it return the cartesian ECEF vector
def eci2ecef(R_ECI,theta_GST,deg=True):
    #R_ECI is ECI cartesian vector
    #Unit agnostic, use the deg switch to decide whether the angle will be degrees or radians
    R_ECEF = rot3(theta_GST,R_ECI,deg=deg) #Keeping it simple, pipe the deg argument through to rot3
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

#Define an ECEF spherical to cartesian transform
def ecef_spherical2cart(lat,lon,r,re=6371.2,deg=True):
    #Offer option to use custom earth radius in km
    #Convert to radians if nessecary
    if deg:
        lat = lat*pi/180.
        lon = lon*pi/180.
    #Height in km
    #r = height+re
    x = r*cos(lat)*cos(lon)
    y = r*cos(lat)*sin(lon)
    z = r*sin(lat)
    return array([x,y,z])

#Cartesian ECEF to Cartesian ECI
def ecef2eci(R_ECEF,theta_GST,deg=True):
    R_ECI = rot3(-1*theta_GST,R_ECEF,deg=deg)
    return R_ECI


#Define some conversions between geocentric (spherical) and geodetic latitude
def spherical2geodetic(gclat,deg=True,ecc_earth=.081819221456):
    #eccentricity from Vallado back cover
    #Valid for points on earth surface only
    if deg:
        gclat = gclat*pi/180.
    gdlat = arctan2(tan(gclat),(1.-ecc_earth**2.))
    if deg:
        gdlat = gdlat*pi/180.
    return gdlat

def geodetic2spherical(gdlat,deg=True,ecc_earth=.081819221456):
    #Valid for points on earth surface only
    if deg:
        gdlat = gdlat*pi/180.
    gclat = arctan2(tan(gdlat),1./(1.-ecc_earth**2.))
    if deg:
        gclat = gdlat*180./pi
    return gclat

#Main function
def ecef2topo(R_ECEF, gdlat, lon, height, deg=True,):
    #Assume input lat is geodetic, if it were geocentric/spherical, use above spherical2geodetic
    #gdlat = geodetic2spherical(gdlat)
    R_site_ecef = ecef_spherical2cart(gdlat,lon,height,deg=deg)
    #Find the ECEF vector of the site
    rho_ecef = R_ECEF-R_site_ecef
    #Compute ECEF range vector
    rho_topo = rot3(lon,rho_ecef,deg=True)
    #Rotate range vector about Z-axis by longitude
    rho_topo = rot2((90.-lat),rho_topo,deg=True)
    #Rotate range vector about y-axis by colatitude
    el = arcsin(rho_topo[2]/linalg.norm(rho_topo))
    #elevation is acos(rho_Z/|rho|), angle up from the SE plane
    beta = pi-arctan2(rho_topo[1],rho_topo[0])
    #Like theta for spherical coords, the azimuth is the angle of rho IN the SE plan
    #But since it's referenced to local north instead of south, it's pi - atan(y/x)
    betasin = arcsin(rho_topo[1]/sqrt(rho_topo[0]**2+rho_topo[1]**2))
    betacos = arccos(-1*rho_topo[0]/sqrt(rho_topo[0]**2+rho_topo[1]**2))
    rng = linalg.norm(rho_topo)
    #The range is just the distance to the spacecraft from the site, so it's just the length of rho vector
    #Convert to degrees for return
    el = el*180./pi
    beta = beta*180./pi
    print("Beta from sin: %.5f" % (betasin*180./pi))
    print("180-Beta from sin: %.5f" % (180.-betasin*180./pi))
    print("Beta from cos: %.5f" % (betacos*180./pi))
    print("-Beta from cos: %.5f" % (-1*betacos*180./pi))
    print("Beta from tan: %.5f" % (beta))
    return array([el,beta,rng]),rho_topo

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

#Compute the theta_GST from Year and fractional day of year
def doy2ymdhms(year,doy):
    #Not vectorized
    #January - 31 Days
    #February - 28 Days (Unless leap year - 29 Days)
    #March - 31 Days
    #April - 30 Days
    #May - 31 Days
    #June - 30 Days
    #July - 31 Days
    #August - 31 Days
    #September - 30 Days
    #October - 31 Days
    #November - 30 Days
    #December - 31 Days
    if len(doy)>1:
        raise ValueError('Not Vectorized!')
    decimaldoy = doy-floor(doy)
    doy = floor(doy)

    mons = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Nov','Dec']
    ndays = array([31,28,31,30,31,30,31,31,30,31,30,31])
    if mod(year,4) == 0:
        ndays[2] = 29
    doys = [sum(ndays[0:k]) for k in arange(len(ndays))]
    doys[0] += 1 #Add the first day in january since we're not zero based
    diffdoys = diffdoys-doy
    for (j,diffdoy) in enumerate(diffdoys):
        if diffdoy < ndays[j] and diffdoy > 0:
            mo = j
            d = diffdoy
            break
    #Hour, Minute, Second parse
    h = floor(decimaldoy*24)
    mn = floor(decimaldoy*24*60)
    s = floor(decimaldoy*24*60*60)
    return y,mo,d,h,mn,s

def ymdhms2jd(year,mon,day,hr,mn,sc):
    #Takes UTC ymdhms time and returns julian date
    #FIXME: Add leap second support
    leapsecond = False
    if year < 1900:
        raise ValueError('Year must be 4 digit year')
    t1 = 367.*year

    t2 = int(7.*(year+int((mon+9.)/12.))/4.)
    t3 = int(275.*mon/9.)
    t4 = day + 1721013.5
    if not leapsecond:
        t5 = ((sc/60.+mn)/60+hr)/24
    else:
        t5 = ((sc/61.+mn)/60+hr)/24
    #print t1,t2,t3,t4,t5
    return t1-t2+t3+t4+t5

def jd2ymdhms(jd):
    dt = jd2datetime(jd)
    return dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second

def jd2datetime(jd):
    #Takes julian date and returns datetime.datetime in UTC

    #The inverse of above, from Vallado pp 208. (algorithm 22)
    T1900 = (jd-2415019.5)/365.25
    year = 1900+int(T1900)
    leapyrs = int((year-1900-1)*.25)
    days = (jd-2415019.5)-((year-1900)*(365.0) + leapyrs)
    if days < 1.0:
        year-=1
        leapyrs = int((year-1900-1)*.25)
        days = (jd-2415019.5)-((year-1900)*(365.0) + leapyrs)
    #doy = int(days)
    return datetime.datetime(year,1,1,0,0,0)+datetime.timedelta(days=days-1)

def jd2gst(JD_UT1,deg=True):
    #Following Vallado pg. 194, gets Greenwich Mean Sideral Time (GMST) in degrees if deg=True or radians otherwise
    #Note that input time is in UT1 NOT UTC. If have UTC and want very accurate theta_gst, need to do UT1 = UTC + Delta_UT1
    #Delta_UT1 is obtainable from the Earth Orientation Parameters (EOP)
    T_UT1 = (JD_UT1-2451545.)/36525 #Get Julian centuries
    #Note that this formula can be broken up into a two part (hours and seconds) version using a two part
    #T_UT1. Where 876600 is multiplied by 3600., and in the exponentiation, the accuracy can be increased
    #by breaking up the T_UT1
    theta_GST_s = 67310.54841+(876600.*3600.+8640184.812866)*T_UT1+.093104*T_UT1**2-6.2e-6*T_UT1**3
    #Make sure abs(theta_GST) <= 86400 seconds
    if abs(theta_GST_s) > 86400.:
        theta_GST_s = mod(theta_GST_s,86400.)
    #Convert theta_GST to degrees from seconds
    theta_GST = theta_GST_s/240.
    if theta_GST < 0.:
        theta_GST = 360.-theta_GST
    if theta_GST > 360.:
        theta_GST = mod(theta_GST,360.)
    if not deg:
        theta_GST = theta_GST * pi / 180.
    return theta_GST


def groundtrack(year,decimaldoy,a,e,w,Omega,M0,n,timestep=60.,timelen=3*3600.,w_e=7.2921158553e-5):
    #year and decimaldoy are the UT1 timestamp/epoch for the orbital elements
    #w_e is earth rotation rate in rad/s
    #n is mean motion in rad/s
    #a is semimajor in km
    #timelen is length of time to propagate orbit for in seconds
    #timestep is the length of each time step in seconds

    ndeg = n * 180/pi #Convert n to degrees per second
    nsteps = floor(timelen/timestep)
    #Compute Julian Date
    yr,mo,dy,hr,mn,sc = doy2ymdhms(year,decimaldoy)
    jd = ymdhms2jd(yr,mo,dy,hr,mn,sc)

    #Init output arrays

    lat_arr = zeros(nsteps+1,1)
    lon_arr = zeros(nsteps+1,1)

    #Set initial values
    M = M0
    theta_GST = jd2gst(jd)

    for k in arange(nsteps):
        E, M_out, tminustp = kepler(ecc,a,M=M)
        nu_sin,nu_cos,nu_tan = eccentrictotrue(E,ecc,a=a)
        nu = quadrant_check(nu_sin,nu_cos)
        #def coe2rv(a,ecc,i,Omega,w,nu,debug=True):
        R_ECI,V_ECI = coe2rv(a,ecc,i,Omega,w,nu)
        R_ECEF = eci2ecef(R_ECI)
        r,lat_arr[k],lon_arr[k] = ecef_cart2spherical(R_ECEF)
        #Convert Spherical Latitude to Geodetic
        lat_arr[k] = spherical2geodetic(lat_arr[k],deg=True)
        #Increment theta_GST and M
        theta_GST = theta_GST + w_e*timestep
        M = M+n_deg*timestep

    return lat_arr,lon_arr

def hour_angle_approx(dt,lons):
    """
    Returns hour angle in degrees
    """
    lons[lons<0.] = lons[lons<0.]+360.
    gamma = 2 * pi / 365 * (dt.timetuple().tm_yday - 1 + float(dt.hour - 12) / 24)
    eqtime = 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) \
                - 0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma))
    decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) \
            - 0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) \
            - 0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma)
    time_offset = eqtime + 4 * lons
    tst = dt.hour * 60 + dt.minute + dt.second / 60 + time_offset
    ha = tst / 4 - 180.
    return ha


def lon2lt(dt,lons):
    """
    Converts and array of longitudes into solar local times
    """
    phi = np.radians(lons)
    #Returns in radians
    gst,sdec,sra = solar_position_approx(dt)
    #Calculate hour angle
    sha = sra - (gst+phi)
    #Convert to hours
    lts = sha*12./np.pi+12.
    return lts

def solar_position_approx(dt,degrees=False):
    """
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
    #Since GST is in degrees convert declination and right ascension
    if degrees:
        sdec = sdec * 180./np.pi
        sransn = sransn * 180./np.pi
        return gst,sdec,sransn
    else:
        gst = np.radians(gst)
        return gst,sdec,sransn

def solar_zenith_angle(dt,lats,lons,degrees=True):
    """
    Finds solar zenith angle using Russel solar position
    """
    lam = np.radians(lats)
    phi = np.radians(lons)
    gst,sdec,sra = solar_position_approx(dt)
    #Calculate hour angle
    sha = sra - (gst+phi)
    cossza = np.sin(lam)*np.sin(sdec) + np.cos(lam)*np.cos(sdec)*np.cos(sha)
    if degrees:
        return np.degrees(np.arccos(cossza))
    else:
        return np.arccos(cossza)

def quadrant_check(from_sin,from_cos,deg=True):
    if not deg:
        #Convert to deg
        from_sin = from_sin*180/pi
        from_cos = from_cos*180/pi
    #Compute other possibilities
    other_sin = 180-from_sin
    other_cos = -1*from_cos
    #Make sure positive
    if other_sin < 0:
        other_sin = other_sin+360.
    if other_cos < 0:
        other_cos = other_cos+360.
    if from_sin < 0:
        from_sin = from_sin+360.
    if from_cos < 0:
        from_cos = from_cos+360.

    #Iterate through all possible combinations
    if from_sin == from_cos:
        a = from_sin
    elif other_sin == from_cos:
        a = from_cos
    elif from_sin == other_cos:
        a = from_sin
    elif other_sin == other_cos:
        a = other_sin
    else:
        print("From sin: %f\n" % (from_sin))
        print("From cos: %f\n" % (from_cos))
        print("180-from_sin: %f\n" % (other_sin))
        print("-1*from_cos: %f\n" % (other_cos))
        raise RuntimeError("Unable to determine quadrant!")
    if not deg:
        a = a*pi/180. #Convert to radians

    return a
