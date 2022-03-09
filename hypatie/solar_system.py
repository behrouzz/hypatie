import numpy as np
from .time import utc2tt, datetime_to_jd
from datetime import datetime

d2r = np.pi/180
r2d = 180/np.pi

def rev(x):
    if x >= 360:
        while x >= 360:
            x = x - 360
    elif x < 0:
        while x < 0:
            x = 360 + x
    return x


def get_sun(t):
    """
    Get geocentric position of the Sun

    Approximation based on method presented by The Astronomical Almanac,
    to calculate the apparent coordinates of the Sun, mean equinox and
    ecliptic of date, to a precision of about 0°.01 (36″), for dates 
    between 1950 and 2050.
    Ref: https://en.wikipedia.org/wiki/Position_of_the_Sun

    Parameters
    ----------
        t (datetime): time of observation (UTC)

    Returns
    -------
        ra  : RA of the Sun (deg)
        dec : DEC of the Sun (deg)
        r   : Distance of the Sun (AU)
    """
    tt = utc2tt(t)
    jd = datetime_to_jd(tt)
    n = jd - 2451545
    L = rev(280.460 + 0.985647*n) # mean longitude
    g = rev(357.528 + 0.9856003*n) # mean anomaly
    # Ecliptic coordinates
    ecl_lon = rev(L + 1.915 * np.sin(g*d2r) + 0.020 * np.sin(2*g*d2r))
    ecl_lat = 0
    r = 1.00014 - 0.01671*np.cos(g*d2r) - 0.00014*np.cos(2*g*d2r)
    # Equatorial coordinates
    obl = 23.439 - 0.0000004*n
    ra = rev(np.arctan2(np.cos(obl*d2r)*np.sin(ecl_lon*d2r), np.cos(ecl_lon*d2r))*r2d)
    dec = np.arcsin(np.sin(obl*d2r)*np.sin(ecl_lon*d2r))*r2d
    return ra, dec, r


def dec_sun(t):
    t = utc2tt(t)
    t0 = datetime(t.year, 1, 1)
    dt = t - t0
    N = dt.days + dt.seconds/86400 + dt.microseconds/(1_000_000*86400)
    in_cos = 0.98565*d2r*(N+10) + 1.914*d2r * np.sin(0.98565*d2r*(N-2))
    dec = -np.arcsin(0.39779 * np.cos(in_cos)) * r2d
    return dec
