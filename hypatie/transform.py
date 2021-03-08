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

def radec_to_altaz(lon, lat, ra, dec, t=None):
    """Convert RaDec to AltAz"""
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
    """Convert AltAz to RaDec"""

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

