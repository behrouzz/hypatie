"""
Module transform
================
This module supplies classes and functions for dealing with coordinates
"""
import numpy as np
from datetime import datetime, timedelta
from collections.abc import Iterable
import re
from .time import datetime_to_jd
from .iau2000A import gcrs2tete, tete_rotmat

d2r = np.pi/180
r2d = 180/np.pi

def _rev(x):
    if x >= 360:
        while x >= 360:
            x = x - 360
    elif x < 0:
        while x < 0:
            x = 360 + x
    return x

rev = np.vectorize(_rev)


class RAhms:
    def __init__(self, h=0, m=0, s=0):
        self.h = h
        self.m = m
        self.s = s
        self.hours = h + m/60 + s/3600
        self.minutes = h*60 + m + s/60
        self.seconds = h*3600 + m*60 + s
        self.deg = 15 * self.hours
        self.rad = self.deg * (3.141592653589793/180)
        

class DECdms:
    def __init__(self, sign='+', d=0, m=0, s=0):
        self.sign = sign
        sgn = -1 if sign=='-' else 1
        self.d = d
        self.m = m
        self.s = s
        self.degrees = (d + m/60 + s/3600) * sgn
        self.minutes = (d*60 + m + s/60) * sgn
        self.seconds = (d*3600 + m*60 + s) * sgn
        self.deg = self.degrees
        self.rad = self.deg * (3.141592653589793/180)


class RAdeg:
    def __init__(self, ra):
        self.rad = ra * (3.141592653589793/180)
        self.h = int(ra/15)
        self.m = int(((ra/15)-self.h)*60)
        self.s = ((((ra/15)-self.h)*60)-self.m)*60
        self.str = str(self.h)+'h' + str(self.m)+'m' + str(self.s)+'s'
        

class DECdeg:
    def __init__(self, dec):
        self.rad = dec * (3.141592653589793/180)
        self.sign = '-' if dec<0 else '+'
        dec = abs(dec)
        self.d = int(dec)
        self.m = abs(int((dec-self.d)*60))
        self.s = (abs((dec-self.d)*60)-self.m)*60
        self.str = self.sign + str(self.d)+'d' + \
                   str(self.m)+'m' + str(self.s)+'s'
        

def to_epoch(ra, dec, epoch):
    """
    Convert coordinates to another epoch

    Note: you have to multiply d_ra and d_dec by (epoch-2000) and
    add the results to the initial ra and dec to get the new ra and
    dec in the desired epoch.

    Example:
    >>> ra2000, dec2000 = 11.795879583333333, -26.3824475
    >>> epoch = 1950
    >>> d_ra, d_dec = to_epoch(ra2000, dec2000, epoch=epoch)
    >>> ra = ra2000 + (epoch-2000)*d_ra
    >>> dec = dec2000 + (epoch-2000)*d_dec
    >>> print(ra,dec)
    11.183690048723463 -26.65500294330134
    
    Arguments
    ---------
        ra  (deg)    : ra in J2000
        dec (deg)    : dec in J2000
        epoch (year) : epoch to convert to
        
    Returns
    -------
        d_ra  : delta RA per year (deg/yr)
        d_dec : delta DEC per year (deg/yr)
    """
    T = (epoch - 2000) / 100
    m = 3.07496 + 0.00186*T
    n = 1.33621 - 0.00057*T
    nn = 20.0431 - 0.0085*T

    d_ra = m + n*np.sin(ra*d2r)*np.tan(dec*d2r)
    d_dec = nn*np.cos(ra*d2r)

    d_ra = RAhms(s=d_ra).deg
    if d_dec<0:
        d_dec = DECdms(sign='-', s=d_dec).deg
    else:
        d_dec = DECdms(s=d_dec).deg
    return d_ra, d_dec


def tete(t):
    """
    Get rotation matrix in order to convert GCRS J2000 coordinates to
    True Equator True Equinox (of date).

    You can convert a GCRS position (pos) to True Equator True Equinox
    using this formula:
    
    rot_mat = tete(t)
    tete_coord = np.matmul(rot_mat, pos)
    
    Arguments
    ---------
        t (datetime) : Observation time (UTC)
        
    Returns
    -------
        Rotation matrix
    """
    return tete_rotmat(t)


def to_tete(pos, t):
    """
    Convert GCRS J2000 coordinates to True Equator True Equinox (of date)
    
    Arguments
    ---------
        pos  (np.array) : [x, y, z] coordinates in GCRS J2000
        t (datetime)    : Observation time (UTC)
        
    Returns
    -------
        position at True Equator True Equinox
    """
    return gcrs2tete(pos, t)


def obliq(t):
    obl_cf = [23.439279444444445, -0.013010213611111111,
              -5.0861111111111115e-08, 5.565e-07,
              -1.6e-10, -1.2055555555555555e-11]
    T = (datetime_to_jd(t) - 2451545) / 36525 # centuries from J2000
    cf = np.array(obl_cf[::-1])
    f = np.poly1d(cf)
    return f(T)


def _time(t):
    if isinstance(t, datetime):
        return t
    elif isinstance(t, str) and bool(re.match("\d{4}-\d\d-\d\d \d\d:\d\d:\d\d", t)):
        return datetime.strptime(t, '%Y-%m-%d %H:%M:%S')
    elif t is None:
        return datetime.utcnow()
    else:
        raise Exception("only datetime or str: '%Y-%m-%d %H:%M:%S'")
   

def mag(x):
    """Returns magnitude of a vector"""
    return np.linalg.norm(np.array(x))

def unit(x):
    """Returns unit vector of a vector"""
    return x / mag(x)

def to_xy_plane(pos):
    """Returns new positions transformed in the xy plane"""
    mags = np.array([mag(i) for i in pos])
    pos_z0 = pos * np.array([1,1,0])
    u = np.array([unit(i) for i in pos_z0])
    new_pos = np.array([u[i] * mags[i] for i in range(len(mags))])
    return new_pos

def rotating_coords(pos, period, times):
    """
    Transforms plane position coordinates (with z=0) to rotating frame
    
    Arguments
    ---------
        pos (np.array): 3d position array
        perios (float): period in days
        times (list): times list
        
    Returns
    -------
        rotating frame position array
    """
    # calculating theta
    t = np.array(times) - times[0]
    t = np.array([i.total_seconds() for i in t])
    period = period * 86400
    omega = 2*np.pi / period
    theta = omega * t

    # calculating rotating coordinates
    x, y = pos[:, 0], pos[:, 1]
    rot_x = x*np.cos(-theta) - y*np.sin(-theta)
    rot_y = x*np.sin(-theta) + y*np.cos(-theta)
    rot_pos = np.c_[rot_x, rot_y, np.zeros(len(theta))]
    return rot_pos


def radec_to_altaz(ra, dec, obs_loc, t):
    """
    Convert ra/dec coordinates to az/alt coordinates

    Arguments
    ---------
        lon (float): longtitude of observer location
        lat (float): latitude of observer location
        ra (iter of float): right ascension value(s)
        dec (iter of float): declination value(s)
        t (datetime or str): time of observation in UTC. default now.

    Returns
    -------
        altitude(s), azimuth(s)
    """
    lon, lat = obs_loc
    t = _time(t)
    if isinstance(ra, Iterable):
        ra = np.array(ra)
        dec = np.array(dec)

    J2000 = datetime(2000,1,1,12)
    d = (t - J2000).total_seconds() / 86400 #day offset

    UT = t.hour + t.minute/60 + t.second/3600
    LST = (100.46 + 0.985647 * d + lon + 15*UT + 360) % 360
    ha = (LST - ra + 360) % 360
    
    x = np.cos(ha*d2r) * np.cos(dec*d2r)
    y = np.sin(ha*d2r) * np.cos(dec*d2r)
    z = np.sin(dec*d2r)
    xhor = x*np.cos((90-lat)*d2r) - z*np.sin((90-lat)*d2r)
    yhor = y
    zhor = x*np.sin((90-lat)*d2r) + z*np.cos((90-lat)*d2r)
    az = np.arctan2(yhor, xhor)*r2d + 180
    alt = np.arcsin(zhor)*r2d
    return az, alt


def altaz_to_radec(lon, lat, az, alt, t=None):
    """
    Convert az/alt coordinates to ra/dec coordinates

    Arguments
    ---------
        lon (float): longtitude of observer location
        lat (float): latitude of observer location
        az (iter of float): azimuth value(s)
        alt (iter of float): altitude value(s)
        t (datetime or str): time of observation in UTC. default now.

    Returns
    -------
        right ascension value(s), declination value(s)
    """

    t = _time(t)
    if isinstance(az, Iterable):
        az = np.array(az)
        alt = np.array(alt)

    ALT = alt * d2r
    AZ  = az  * d2r
    LAT = lat * d2r

    sin_DEC = np.sin(LAT)*np.sin(ALT) + np.cos(LAT)*np.cos(ALT)*np.cos(AZ)
    DEC = np.arcsin(sin_DEC)

    cos_HA = (np.sin(ALT) - np.sin(LAT)*np.sin(DEC)) / \
             (np.cos(LAT)*np.cos(DEC))
    
    HA = np.where(cos_HA>1, 0, np.arccos(cos_HA))

    dec = DEC * r2d
    ha  = HA  * r2d
    
    J2000 = datetime(2000,1,1,12)
    d = (t - J2000).total_seconds() / 86400
    UT = t.hour + t.minute/60 + t.second/3600
    LST = (100.46 + 0.985647 * d + lon + 15*UT + 360) % 360

    ra = np.where(az>=180, (LST-ha)%360, (LST+ha)%360)

    return ra, dec


def _in_vec(inp):
    if len(inp.shape)>1:
        if inp.shape[-1]==3:
            vec_out = inp[:,0], inp[:,1], inp[:,2]
        elif inp.shape[-1]==2:
            vec_out = inp[:,0], inp[:,1]
        else:
            raise Exception('Input imension not sopported!')
    elif len(inp.shape)==1:
        vec_out = inp
    else:
        raise Exception('Dimensionality problem! 2')
    return vec_out


def _out_vec(vec):
    if len(vec.shape)>1:
        vec = vec.T
        if vec.shape[-1]==3:
            x, y, z = vec[:,0], vec[:,1], vec[:,2]
            vec_out = np.vstack((x, y, z)).T
        elif vec.shape[-1]==2:
            x, y = vec[:,0], vec[:,1]
            vec_out = np.vstack((x, y)).T
        else:
            raise Exception('Dimensionality problem! 3')
    elif len(vec.shape)==1:
        vec_out = vec
    else:
        raise Exception('Dimensionality problem! 4')
    return vec_out
        


def sph2car(sph):
    lon, lat, r = _in_vec(sph)
    lon, lat = lon*d2r, lat*d2r
    x = r * np.cos(lat) * np.cos(lon)
    y = r * np.cos(lat) * np.sin(lon)
    z = r * np.sin(lat)
    car = _out_vec(np.array([x, y, z]))
    return car



def car2sph(car):
    x, y, z = _in_vec(car)
    lon = rev(np.arctan2(y, x)*r2d)
    lat = np.arctan2(z, np.sqrt(x**2 + y**2))*r2d
    r = np.sqrt(x**2 + y**2 + z**2)
    sph = _out_vec(np.array([lon, lat, r]))
    return sph


def _equsph2eclsph(radec):
    ra, dec = _in_vec(radec)
    ra, dec = ra*d2r, dec*d2r
    epsilon = (23 + 26/60 + 21.448/3600)*d2r
    sb = np.sin(dec)*np.cos(epsilon) - np.cos(dec)*np.sin(epsilon)*np.sin(ra)
    cbcl = np.cos(dec)*np.cos(ra)
    cbsl = np.sin(dec)*np.sin(epsilon) + np.cos(dec)*np.cos(epsilon)*np.sin(ra)
    lon = rev(np.arctan2(cbsl,cbcl)*r2d)
    r = np.sqrt(cbsl**2 + cbcl**2)
    lat = np.arctan2(sb,r)*r2d
    lonlat = _out_vec(np.array([lon, lat]))
    return lonlat


def equ_sph2ecl_sph(equsph):
    if len(equsph.shape)>1:
        eclsph = np.zeros((len(equsph), 3))
        eclsph[:,:2] = _equsph2eclsph(equsph[:,:2])
        eclsph[:, 2] = equsph[:,2]
    else:
        lon, lat = _equsph2eclsph(equsph[:2])
        eclsph = np.array([lon, lat, equsph[2]])
    return eclsph


def equ_car2ecl_car(equ_car):
    equ_sph = car2sph(equ_car)
    ecl_sph = equ_sph2ecl_sph(equ_sph)
    ecl_car = sph2car(ecl_sph)
    return ecl_car



def angular_separation(r1, d1, r2, d2):
    """
    Calculate angular separation between two point

    Arguments
    ---------
        r1 (float): right ascension of the first point in degrees
        d1 (float): declination of the first point in degrees
        r2 (float): right ascension of the second point in degrees
        d2 (float): declination of the second point in degrees

    Returns
    -------
        angular sepration in degrees
    """
    r1 = (np.pi/180) * r1
    d1 = (np.pi/180) * d1
    r2 = (np.pi/180) * r2
    d2 = (np.pi/180) * d2

    radical = np.sqrt(
        (np.cos(d2)**2)*(np.sin(r2-r1)**2) + ((np.cos(d1)*np.sin(d2) - np.sin(d1)*np.cos(d2)*np.cos(r2-r1))**2)
        )
    
    kasr = radical / ( (np.sin(d1)*np.sin(d2)) + (np.cos(d1)*np.cos(d2)*np.cos(r2-r1)) )

    sep = (180/np.pi) * np.arctan(kasr)
    
    return sep


def posvel(ra, dec, pmra, pmdec, distance, radvel):
    """
    Cartesian position and velocity vectors.
    
    Note:
    -----
    1) If you have parallax instead of distance, use this formulla:
    distance = (1/(plx/1000)) * 30856775814913.67 # km
    
    2) If you have redshift instead of radial velocity, use this formulla:
    radvel = z * 299792.458 # km/s
    
    Arguments
    ---------
        ra (iter of float): RA of objects
        dec (iter of float): DEC of objects
        pmra (iter of float): pmRA of objects
        pmdec (iter of float): pmDEC of objects
        distance (iter of float): Distance to objects in km
        radvel (iter of float): Radial Velocity of objects in km/s
        
    Returns
    -------
        position vectors (km), velocity vectors (km/s)
    """
    
    # Changes of RA and DEC in one second
    d_ra = ((pmra/1000)/3600)/31557600
    d_dec = ((pmdec/1000)/3600)/31557600
    
    r1, d1 = ra, dec
    r2 = r1 + d_ra
    d2 = d1 + d_dec
    
    # Changes in distance
    dist1 = distance
    dist2 = dist1 + radvel
    
    # Convert to cartesian
    x1, y1, z1 = radec_to_cartesian(r1, d1, dist1)
    x2, y2, z2 = radec_to_cartesian(r2, d2, dist2)
    pos1 = np.array([x1,y1,z1])
    pos2 = np.array([x2,y2,z2])
    
    # Velocity
    d_pos = pos2 - pos1
    vel = d_pos * 1 # because it's in one second
    return pos1, vel
