"""
Module iau2000A
===============
Transform coordinates from GCRS to True Equator True Equinox
Ref: https://arxiv.org/abs/astro-ph/0602086
"""

import csv
import os.path
from urllib.request import urlretrieve
import numpy as np
from .time import utc2tdb, datetime_to_jd

url = 'https://raw.githubusercontent.com/behrouzz/astrodatascience/main/data/iau2000A_nut.csv'
s2r = (1/3600) * (np.pi/180)


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



def nutation(T):
    fname = 'iau2000A_nut.csv'
    if not os.path.exists(fname):
        print('Downloading '+fname+'...')
        urlretrieve(url, fname)

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

    #ep = e + d_e

    S1 = np.sin(e)
    S2 = np.sin(-d_psi)
    S3 = np.sin(-e-d_e)

    C1 = np.cos(e)
    C2 = np.cos(-d_psi)
    C3 = np.cos(-e-d_e)

    row1 = [C2, S2*C1, S2*S1]
    row2 = [-S2*C3, C3*C2*C1-S1*S3, C3*C2*S1+C1*S3]
    row3 = [S2*S3, -S3*C2*C1-S1*C3, -S3*C2*S1+C3*C1]

    N = np.array([row1, row2, row3])
    return N

# B : frame bias matrix
mas2rad = 4.84813681109536e-09

da0 = -14.6 * mas2rad
e0 = -16.6170 * mas2rad
n0 = -6.8192 * mas2rad

B = np.array([[ 1-0.5*(da0**2 + e0**2), da0, -e0],
              [-da0-n0*e0, 1-0.5*(da0**2 + n0**2), -n0],
              [e0-n0*da0, n0+e0*da0, 1-0.5*(n0**2+e0**2)]])

def tete_rotmat(t):
    t = utc2tdb(t)
    T = (datetime_to_jd(t) - 2451545) / 36525
    P = precession(T)
    N = nutation(T)
    rot_mat = np.matmul(np.matmul(N, P), B)
    return rot_mat

def gcrs2tete(pos, t):
    rot_mat = tete_rotmat(t)
    return np.matmul(rot_mat, pos)

