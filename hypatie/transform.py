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
