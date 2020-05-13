# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
from numpy import *
import numpy as np
import matplotlib.pyplot as pp
import itertools, datetime
from geospacepy.special_datetime import jd2datetime,datetime2jd

G = 6.67e-11 #N m^2/s^2
m_earth = 5.9742e24 #kg
r_earth = 6371200 #m
mu = G*m_earth
#Cartesian Unit Vectors
I = array([1.,0.,0.])
J = array([0.,1.,0.])
K = array([0.,0.,1.])

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
