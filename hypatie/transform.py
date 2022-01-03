"""
Module transform
=============
This module supplies two functions (radec_to_altaz and altaz_to_radec)
to transform between coordinates systems.
"""
import numpy as np
from datetime import datetime
from collections.abc import Iterable
import re

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

def radec_to_altaz(lon, lat, ra, dec, t=None):
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
    t = _time(t)
    d2r = np.pi/180
    r2d = 180/np.pi

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
    return alt, az

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
    d2r = np.pi/180
    r2d = 180/np.pi

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

def radec_to_cartesian(ra, dec, r):
    """
    Convert ra/dec/distance coordinates to cartesian coordinates

    Arguments
    ---------
        ra (iter of float): right ascension value(s) in degrees
        dec (iter of float): declination value(s) in degrees
        r (iter of float): distance value(s)

    Returns
    -------
        x,y,z
    """
    d2r = np.pi/180
    x = r * np.cos(dec*d2r) * np.cos(ra*d2r)
    y = r * np.cos(dec*d2r) * np.sin(ra*d2r)
    z = r * np.sin(dec*d2r)
    return x,y,z

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

def hmsdms_to_deg(hmsdms):
    """
    Convert HMS (hours, minutes, seconds) and DMS (degrees, minutes, seconds) to
    RA, DEC in decimal degrees.
    Example:
        hmsdms_to_deg('06 45 08.91728 -16 42 58.0171')
    Return:
        (101.28715533333333, -15.28388413888889)
    """
    ls = hmsdms.split(' ')
    ra_h = int(ls[0])
    ra_m = int(ls[1])
    ra_s = float(ls[2])
    dec_d = int(ls[3])
    dec_m = int(ls[4])
    dec_s = float(ls[5])

    ra = 15*ra_h + 15*ra_m/60 + 15*ra_s/3600
    dec = dec_d + dec_m/60 + dec_s/3600

    return ra, dec

def posvel(ra, dec, pmra, pmdec, distance, radvel):
    """
    Cartesian position and velocity vectors.
    
    Note:
    -----
    1) If you have parallax instead of distance, use this formulla:
    distance = (1/(plx/1000)) * 30856775814913.67
    
    2) If you have redshift instead of radial velocity, use this formulla:
    radvel = z * 299792.458
    
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
        position vectors, velocity vectors
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
