import pickle, io
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



def jd_to_sec(jd):
    """jd to seconds since J2000"""
    return (jd - 2451545.0) * 86400.0


def num2txt(arr):
    output = io.BytesIO()
    np.savetxt(output, arr)
    return output.getvalue().decode('utf-8')


class Segment:
    """
    Segment of SPK data created by numeph python package
    """
    def __init__(self, cet_tar, domain, coef):
        self.cet_tar = cet_tar
        self.center = cet_tar[0]
        self.target = cet_tar[1]
        self.domain = domain
        self.coef = coef

    def to_str(self):
        cfx = self.coef[0,:,:]
        cfy = self.coef[1,:,:]
        cfz = self.coef[2,:,:]
        cf_xyz = np.concatenate((cfx,cfy,cfz))
        cf_str = num2txt(cf_xyz)
        n_cols = cf_xyz.shape[1]
        ini_dom = self.domain[0,0]
        dt_rec = self.domain[0,1] - self.domain[0,0]
        dt_dom = self.domain[-1,-1] - self.domain[0,0]
        n_recs = cfx.shape[0]

        first_row = [self.center, self.target, n_cols, n_recs, dt_rec, ini_dom, dt_dom]
        first_row = ['segment'] + [str(int(i)) for i in first_row]
        first_row = ','.join(first_row)+'\n'
        return first_row + cf_str

    def get_pos(self, t):
        """
        Get position of the object at time t
        
        Parameters
        ----------
            t (datetime)    : time for which the position is requested
        
        Returns
        ----------
            pos (np.array): position of the object at t
        """
        jd = datetime_to_jd(t)
        t = jd_to_sec(jd)
        mask = np.logical_and(t>=self.domain[:,0], t<self.domain[:,1])
        rec = np.where(mask)[0][0] # record index
        cfx = self.coef[0,rec,:]
        cfy = self.coef[1,rec,:]
        cfz = self.coef[2,rec,:]
        fx = np.polynomial.chebyshev.Chebyshev(coef=cfx, domain=self.domain[rec])
        fy = np.polynomial.chebyshev.Chebyshev(coef=cfy, domain=self.domain[rec])
        fz = np.polynomial.chebyshev.Chebyshev(coef=cfz, domain=self.domain[rec])
        pos = np.vstack((fx(t),fy(t),fz(t))).T[0]
        return pos


def load_pickle(fname):
    """
    Load an ephemeris pickle file (created by numeph python package)
    
    Parameters
    ----------
        fname (str) : path and name of the pickle file
    
    Returns
    ----------
        Dictionary of Segments
    """
    f = open(fname, 'rb')
    data = pickle.load(f)
    f.close()
    dc = {}
    for k,v in data.items():
        dc[k] = Segment(k, v[0], v[1])
    return dc
