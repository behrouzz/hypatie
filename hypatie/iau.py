"""
Module iau
===============
Transform coordinates from GCRS to True Equator True Equinox and many more
IAU 2000A models
Ref: https://arxiv.org/abs/astro-ph/0602086
"""

import csv
import os.path
from urllib.request import urlretrieve
import numpy as np
from datetime import datetime, timedelta
from .time import utc2tdb, datetime_to_jd
from .utils import is_recent

NUT_URL = 'https://raw.githubusercontent.com/behrouzz/astrodatascience/main/data/iau2000A_nut.csv'
FINALS2000A_URL = 'https://maia.usno.navy.mil/ser7/finals2000A.all'
FINALS2000A_FNAME = 'finals2000A.all'

s2r = (1/3600) * (np.pi/180)
d2r = np.pi / 180

def get_T(utc):
    t = utc2tdb(utc)
    T = (datetime_to_jd(t) - 2451545) / 36525
    return T

def precession(T):
    e0 = 84381.406
    psi_A_cf = [0, 5038.481507, -1.0790069, -0.00114045, 0.000132851, -0.0000000951]
    w_A_cf = [e0, -0.025754, 0.0512623, -0.00772503, -0.000000467, 0.0000003337]
    chi_A_cf = [0, 10.556403, -2.3814292, -0.00121197, 0.000170663, -0.0000000560]

    psi_A = np.poly1d(psi_A_cf[::-1])(T)
    w_A = np.poly1d(w_A_cf[::-1])(T)
    chi_A = np.poly1d(chi_A_cf[::-1])(T)

    # convert to radian
    e0 = e0 * s2r
    psi_A = psi_A * s2r
    w_A = w_A * s2r
    chi_A = chi_A * s2r

    s1 = np.sin(e0)
    s2 = np.sin(-psi_A)
    s3 = np.sin(-w_A)
    s4 = np.sin(chi_A)

    c1 = np.cos(e0)
    c2 = np.cos(-psi_A)
    c3 = np.cos(-w_A)
    c4 = np.cos(chi_A)

    row1 = [c4*c2-s2*s4*c3, c4*s2*c1+s4*c3*c2*c1-s1*s4*s3, c4*s2*s1+s4*c3*c2*s1+c1*s4*s3]
    row2 = [-s4*c2-s2*c4*c3, -s4*s2*c1+c4*c3*c2*c1-s1*c4*s3, -s4*s2*s1+c4*c3*c2*s1+c1*c4*s3]
    row3 = [s2*s3, -s3*c2*c1-s1*c3, -s3*c2*s1+c3*c1]

    P = np.array([row1, row2, row3])
    return P




def create_phi(T): #new
    fname = 'iau2000A_nut.csv'
    if not os.path.exists(fname):
        print('Downloading '+fname+'...')
        urlretrieve(NUT_URL, fname)

    data = []
    with open(fname, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            data.append(row)
    data = data[1:]
    
    M = np.array([[int(i) for i in row[1:15]] for row in data])
    S = np.array([float(i[15]) for i in data]) * s2r
    Sd = np.array([float(i[16]) for i in data]) * s2r
    Cp = np.array([float(i[17]) for i in data]) * s2r
    C = np.array([float(i[18]) for i in data]) * s2r
    Cd = np.array([float(i[19]) for i in data]) * s2r
    Sp = np.array([float(i[20]) for i in data]) * s2r

    phi = {
    1 : 908103.259872 + 538101628.688982 * T,
    2 : 655127.283060 + 210664136.433548 * T,
    3 : 361679.244588 + 129597742.283429 * T,
    4 : 1279558.798488 + 68905077.493988 * T,
    5 : 123665.467464 + 10925660.377991 * T,
    6 : 180278.799480 + 4399609.855732 * T,
    7 : 1130598.018396 + 1542481.193933 * T,
    8 : 1095655.195728 + 786550.320744 * T,
    9 : 5028.8200 * T + 1.112022 * T**2,
    10 : 485868.249036 + 1717915923.2178 * T + 31.8792 * T**2 + 0.051635 * T**3 - 0.00024470 * T**4,
    11 : 1287104.79305 + 129596581.0481 * T - 0.5532 * T**2 + 0.000136 * T**3 - 0.00001149 * T**4,
    12 : 335779.526232 + 1739527262.8478 * T - 12.7512 * T**2 - 0.001037 * T**3 + 0.00000417 * T**4,
    13 : 1072260.70369 + 1602961601.2090 * T - 6.3706 * T**2 + 0.006593 * T**3 - 0.00003169 * T**4,
    14 : 450160.398036 - 6962890.5431 * T + 7.4722 * T**2 + 0.007702 * T**3 - 0.00005939 * T**4,
    }


    d_psi = 0
    d_e = 0
    for i in range(1,1366):
        PHI = 0
        for j in range(1,15):
            PHI = PHI + M[i-1,j-1] * phi[j]
        PHI = PHI * s2r
        d_psi = d_psi + ((S[i-1]+Sd[i-1]*T)*np.sin(PHI) + Cp[i-1]*np.cos(PHI))
        d_e = d_e + ((C[i-1]+Cd[i-1]*T)*np.cos(PHI) + Sp[i-1]*np.sin(PHI))

    e0 = 84381.406
    e = e0 - 46.836769*T - 0.0001831 *T**2 + 0.00200340 *T**3 - 0.000000576 *T**4 - 0.0000000434 *T**5
    e = e * s2r

    F = phi[12] * s2r
    D = phi[13] * s2r
    om = phi[14] * s2r
    return d_psi, d_e, e, F, D, om

    

def nutation(d_psi, d_e, e):
    #d_psi, d_e, e, _, _, _ = create_phi(T)
    S1 = np.sin(e)
    S2 = np.sin(-d_psi)
    S3 = np.sin(-e-d_e)
    C1 = np.cos(e)
    C2 = np.cos(-d_psi)
    C3 = np.cos(-e-d_e)
    row1 = [C2, S2*C1, S2*S1]
    row2 = [-S2*C3, C3*C2*C1-S1*S3, C3*C2*S1+C1*S3]
    row3 = [S2*S3, -S3*C2*C1-S1*C3, -S3*C2*S1+C3*C1]
    return np.array([row1, row2, row3])


def equation_of_equinoxes(T, d_psi, e, F, D, om):
    # Equation of equinoxes (arcsec)
    eq_eq = d_psi * np.cos(e) \
        + 0.00264096 * np.sin(om) \
        + 0.00006352 * np.sin(2*om) \
        + 0.00001175 * np.sin(2*F - 2*D + 3*om) \
        + 0.00001121 * np.sin(2*F - 2*D + om) \
        - 0.00000455 * np.sin(2*F - 2*D + 2*om) \
        + 0.00000202 * np.sin(2*F + 3*om) \
        + 0.00000198 * np.sin(2*F + om) \
        - 0.00000172 * np.sin(3*om) \
        - 0.00000087 * T * np.sin(om)
    return eq_eq


def equation_of_origins(T, eq_eq):
    par = 0.014506 + 4612.156534*T + 1.3915817* T**2 \
          - 0.00000044 *T**3 - 0.000029956 *T**4 - 0.0000000368 *T**5
    return -(par + eq_eq)

# B : frame bias matrix
mas2rad = 4.84813681109536e-09

da0 = -14.6 * mas2rad
e0 = -16.6170 * mas2rad
n0 = -6.8192 * mas2rad

B = np.array([[ 1-0.5*(da0**2 + e0**2), da0, -e0],
              [-da0-n0*e0, 1-0.5*(da0**2 + n0**2), -n0],
              [e0-n0*da0, n0+e0*da0, 1-0.5*(n0**2+e0**2)]])

def tete_rotmat(t):
    T = get_T(t)
    P = precession(T)
    d_psi, d_e, e, _, _, _ = create_phi(T)
    N = nutation(d_psi, d_e, e)
    rot_mat = np.matmul(np.matmul(N, P), B)
    return rot_mat

def gcrs2tete(pos, t):
    rot_mat = tete_rotmat(t)
    return np.matmul(rot_mat, pos)


def get_C(T):
    """
    Get C matrix

    T can be calculated with this formula:
    T = (datetime_to_jd(utc2tdb(t)) - 2451545) / 36525

    Using C matrix, one can tranfsorm from GCRS to CIRS:
    r_CIRS = C * r_GCRS

    Arguments:
    ----------
        T: number of centuries of TDB from J2000

    Return:
    -------
        C matrix    
    """
    P = precession(T)
    d_psi, d_e, e, F, D, om = create_phi(T)
    N = nutation(d_psi, d_e, e)
    
    tmp = np.matmul(np.matmul(B.T, P.T), N.T)
    n_GCRS = np.matmul(tmp, np.array([0,0,1]))
    eq_GCRS = np.matmul(tmp, np.array([1,0,0]))

    eq_eq = equation_of_equinoxes(T, d_psi, e, F, D, om)
    eq_or = equation_of_origins(T, eq_eq)
    eq_o_rad = eq_or * s2r

    # sig_GCRS: position of CIO in GCRS
    sig_GCRS = eq_GCRS * np.cos(eq_o_rad) - np.cross(n_GCRS, eq_GCRS) * np.sin(eq_o_rad)
    y_GCRS = np.cross(n_GCRS, sig_GCRS) # page 65

    C = np.array([sig_GCRS, y_GCRS, n_GCRS])
    # page 66 and 48 (MOHEM: combined method)
    return C


def loc2itrs(lon, lat, h):
    LON, LAT = lon*d2r, lat*d2r
    a = 6378137
    f = 1 / 298.257223563
    C = 1 / np.sqrt(np.cos(LAT)**2 + (1-f)**2 * np.sin(LAT)**2)
    S = (1 - f)**2 * C
    r = np.array([(a*C+h)*np.cos(LAT)*np.cos(LON),
                  (a*C+h)*np.cos(LAT)*np.sin(LON),
                  (a*S+h)*np.sin(LAT)])
    return r


def interpolate(t, array):
    mjd = datetime_to_jd(t) - 2400000.5
    m1 = mjd>=array[:,0]
    m2 = mjd<array[:,0]
    t1 = array[m1][-1] # previous
    t2 = array[m2][0] # next
    # interpolate
    cf = np.polyfit([t1[0],t2[0]], [t1[1],t2[1]], 1)
    f = np.poly1d(cf)
    return f(mjd)


def get_finals2000a():
    should_download = True
    if os.path.exists(FINALS2000A_FNAME):
        if is_recent(FINALS2000A_FNAME):
            should_download = False
    if should_download:
        print('Downloading '+FINALS2000A_FNAME+'...')
        urlretrieve(FINALS2000A_URL, FINALS2000A_FNAME)

    with open(FINALS2000A_FNAME, 'r') as f:
        data = f.read().split('\n')

    data = [i for i in data if len(i.strip())>68]
    mjd = [int(i[7:12]) for i in data]
    dut1 = [float(i[58:68]) for i in data]
    pm_x = [float(i[18:27]) for i in data]
    pm_y = [float(i[37:46]) for i in data]

    dut1_array = np.array(list(zip(mjd,dut1)))
    pm_x = np.array(list(zip(mjd, pm_x)))
    pm_y = np.array(list(zip(mjd, pm_y)))
    return dut1_array, pm_x, pm_y


def ut1_utc(t, dut1_array):
    mjd = datetime_to_jd(t) - 2400000.5
    m1 = mjd>=dut1_array[:,0]
    m2 = mjd<dut1_array[:,0]
    t1 = dut1_array[m1][-1] # previous
    t2 = dut1_array[m2][0] # next
    # interpolate
    cf = np.polyfit([t1[0],t2[0]], [t1[1],t2[1]], 1)
    f = np.poly1d(cf)
    return f(mjd)


def get_W(t, pm_x, pm_y):
    x = interpolate(t, pm_x) * s2r
    y = interpolate(t, pm_y) * s2r
    #T = get_T(t)
    #s_p = -47 * microascsec * T

    row1 = [np.cos(x), np.sin(x)*np.sin(y), -np.sin(x)*np.cos(y)]
    row2 = [0, np.cos(y), np.sin(y)]
    row3 = [np.sin(x), -np.sin(y), np.cos(x)*np.cos(y)]
    W = np.array([row1, row2, row3])
    # or: W = np.array([[1, 0, -x], [0, 1, y], [x, -y, 1]])
    return W

def R1(x):
    r1 = [1, 0, 0]
    r2 = [0, np.cos(x), np.sin(x)]
    r3 = [0, -np.sin(x), np.cos(x)]
    return np.array([r1, r2, r3])

def R2(x):
    r1 = [np.cos(x), 0, -np.sin(x)]
    r2 = [0, 1, 0]
    r3 = [np.sin(x), 0, np.cos(x)]
    return np.array([r1, r2, r3])

def R3(x):
    r1 = [np.cos(x), np.sin(x), 0]
    r2 = [-np.sin(x), np.cos(x), 0]
    r3 = [0, 0, 1]
    return np.array([r1, r2, r3])

def get_ERA(utc, dut1):
    #dut1_arr, _, _ =get_finals2000a()
    #dut1 = ut1_utc(utc, dut1_arr)
    ut1 = utc + timedelta(seconds=dut1)
    Du = datetime_to_jd(ut1) - 2451545
    ERA = 2*np.pi * (0.7790572732640 + 1.00273781191135448 * Du)
    era_deg = ERA * (180/np.pi)
    return era_deg % 360


def get_era_eqor(ra, dec, obs_loc, t):
    """
    Calculate hour angle of an object

    Arguments
    ---------
        obs_loc (tuple): (longtitude, latitude) of observer
        ra (float): RA of the object (equinox of date)
        dec (float): DEC of the object (equinox of date)
        t (datetime): time of observation in UTC        

    Returns
    -------
        hour angle (degrees)
    """
    LON, LAT = obs_loc
    T = get_T(t)
    
    # Calculate ERA
    # =============
    dut1_array, pm_x, pm_y = get_finals2000a()
    dut1 = ut1_utc(t, dut1_array)
    ERA = get_ERA(t, dut1)

    # Find equation of origins
    d_psi, d_e, e, F, D, om = create_phi(T)
    eq_eq = equation_of_equinoxes(T, d_psi, e, F, D, om)
    eq_or = equation_of_origins(T, eq_eq)

    # Real lon and lat
    # ================
    #lon , lat = LON, LAT
    xp = interpolate(t, pm_x)
    yp = interpolate(t, pm_y)
    lon = LON + (xp*np.sin(LON*d2r) + yp*np.cos(LON*d2r)) * np.tan(LAT*d2r) / 3600
    lat = LAT + (xp*np.cos(LON*d2r) - yp*np.sin(LON*d2r)) / 3600

    # Calculate hour angle of object
    # ==============================
    # ra     : RA of object wrt equinox of date
    # ra_sig : RA of object wrt CIO
    eq_or = eq_or / 3600 # convert to degrees
    #ra_sig = (ra + eq_or) % 360
    #ha = (ERA - ra_sig + lon) % 360
    return ERA, eq_or, (lon, lat)


def get_ha(ra, dec, obs_loc, t):
    """
    Calculate hour angle of an object

    Arguments
    ---------
        obs_loc (tuple): (longtitude, latitude) of observer
        ra (float): RA of the object (equinox of date)
        dec (float): DEC of the object (equinox of date)
        t (datetime): time of observation in UTC        

    Returns
    -------
        hz (deg), obs_loc (lon, lat)
    """
    
    era, eqor, obs_loc = get_era_eqor(ra, dec, obs_loc, t)
    lon = obs_loc[0]
    ra_sig = (ra + eqor) % 360
    ha = (era - ra_sig + lon) % 360
    return ha, obs_loc


"""
T = get_T(t)
GMST = 86400 * ERA + (0.014506 + 4612.156534*T + 1.3915817* T**2 \
         - 0.00000044 *T**3 - 0.000029956 *T**4 - 0.0000000368 *T**5)/15
GAST = GMST + equation_of_equinoxes/15
LAST = GAST + (3600/15)*lon # lon should be calculated with eq page 17
"""
