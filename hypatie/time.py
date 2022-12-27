from datetime import datetime, timedelta
import math
import pickle, os
import numpy as np
import pandas as pd
from scipy import interpolate
from urllib.request import urlretrieve


d2r = math.pi/180

leaps = {
2272060800: 10, # 1 Jan 1972
2287785600: 11, # 1 Jul 1972
2303683200: 12, # 1 Jan 1973
2335219200: 13, # 1 Jan 1974
2366755200: 14, # 1 Jan 1975
2398291200: 15, # 1 Jan 1976
2429913600: 16, # 1 Jan 1977
2461449600: 17, # 1 Jan 1978
2492985600: 18, # 1 Jan 1979
2524521600: 19, # 1 Jan 1980
2571782400: 20, # 1 Jul 1981
2603318400: 21, # 1 Jul 1982
2634854400: 22, # 1 Jul 1983
2698012800: 23, # 1 Jul 1985
2776982400: 24, # 1 Jan 1988
2840140800: 25, # 1 Jan 1990
2871676800: 26, # 1 Jan 1991
2918937600: 27, # 1 Jul 1992
2950473600: 28, # 1 Jul 1993
2982009600: 29, # 1 Jul 1994
3029443200: 30, # 1 Jan 1996
3076704000: 31, # 1 Jul 1997
3124137600: 32, # 1 Jan 1999
3345062400: 33, # 1 Jan 2006
3439756800: 34, # 1 Jan 2009
3550089600: 35, # 1 Jul 2012
3644697600: 36, # 1 Jul 2015
3692217600: 37 # 1 Jan 2017
}


def sec_to_jd(seconds):
    """seconds since J2000 to jd"""
    return 2451545.0 + seconds / 86400.0


def jd_to_sec(jd):
    """jd to seconds since J2000"""
    return (jd - 2451545.0) * 86400.0

  
def datetime_to_jd(t):
    year = t.year
    month = t.month
    t_d = t.day
    t_H = t.hour
    t_M = t.minute
    t_S = t.second
    t_MS = t.microsecond
    day = t_d + t_H/24 + t_M/(24*60) + t_S/(24*60*60) + t_MS/(24*60*60*1000000)
    
    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month
    
    if ((year < 1582) or
        (year == 1582 and month < 10) or
        (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)
        
    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)
        
    D = math.trunc(30.6001 * (monthp + 1))
    jd = B + C + D + day + 1720994.5
    
    return jd


def jd_to_datetime(jd):

    jd = jd + 0.5
    F, I = math.modf(jd)
    I = int(I)
    A = math.trunc((I - 1867216.25)/36524.25)
    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    D = math.trunc((C - 122.1) / 365.25)
    E = math.trunc(365.25 * D)
    G = math.trunc((C - E) / 30.6001)
    day = C - E + F - math.trunc(30.6001 * G)
    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715

    d_ = day
    d = int(d_)
    h_ = (d_-d)*24
    h = int(h_)
    m_ = (h_-h)*60
    m = int(m_)
    s_ = (m_-m)*60
    s = int(s_)
    ms_ = (s_-s)*1000000
    ms = int(ms_)
    dt = datetime(year, month, d, h, m, s, ms)
        
    return dt


def get_gmst(t):
    """Greenwich mean sidereal time (in degrees)"""
    JD = datetime_to_jd(t)
    T = (JD - 2451545) / 36525
    GMST = 280.46061837 + 360.98564736629 * (JD - 2451545) + (0.000387933 * T**2) - (T**3 / 38710000)
    return GMST % 360


def get_lst(t, lon):
    GMST = get_gmst(t)
    LST = GMST + lon
    return LST % 360


def get_ha(t, lon, ra):
    LST = get_lst(t, lon)
    ha = LST - ra
    return ha % 360


def get_lp(t):
    s = (t-datetime(1900,1,1)).total_seconds()
    if t < datetime(1972,1,1):
        lp = 0
    else:
        before = [i for i in list(leaps.keys()) if i<=s]
        lp = leaps[sorted(before)[-1]]
    return lp


def utc2tt(t):
    lp = get_lp(t)
    TAI = t + timedelta(seconds=lp)
    TT = TAI + timedelta(seconds=32.184)
    return TT

def tt2utc(t):
    lp = get_lp(t)
    TAI = t - timedelta(seconds=32.184)
    utc = TAI - timedelta(seconds=lp)
    return utc

def tt2tdb(t):
    T = (datetime_to_jd(t) - 2451545) / 36525 # centuries from J2000
    g = (2*math.pi/360) * (357.528 + 35999.050*T)
    dt = 0.001658 * math.sin(g + 0.0167*math.sin(g))
    TDB = t + timedelta(seconds=dt)
    return TDB

def tdb2tt(tdb):
    g = 357.53 + 0.9856003*(datetime_to_jd(tdb)-2451545)
    dt = (0.001658*math.sin(g*d2r) + 0.000014*math.sin(2*g*d2r))
    tt = tdb - timedelta(seconds=dt)
    return tt

def utc2tdb(utc):
    tt = utc2tt(utc)
    tdb = tt2tdb(tt)
    return tdb

def tdb2utc(tdb):
    tt = tdb2tt(tdb)
    utc = tt2utc(tt)
    return utc

def tcb2utc(t):
    mjd = datetime_to_jd(t) - 2400000.5
    dt = (1.55051976772e-8 * (mjd-43144) * 86400 + 6.55e-5)
    tdb = t - timedelta(seconds=dt)
    return tdb2utc(tdb)


def tim2ord(t):
    d = (t - datetime(t.year, 1, 1)).total_seconds()/86400
    return d


def download_eot_file(path=None):
    if path is None:
        path = ''
    tbl_file = 'eot_2020_2050.csv'
    url = 'https://raw.githubusercontent.com/behrouzz/eot/main/'+tbl_file
    print(f'Downloading "{tbl_file}"...')
    urlretrieve(url, path+tbl_file)
    print(f'"{tbl_file}" downloaded.')
    print()
    return path+tbl_file


def eot_func(year=None, tbl_file=None, interp_kind='cubic'):
    if year is None:
        year = datetime.utcnow().year
    if tbl_file is None:
        tbl_file = 'eot_2020_2050.csv'
        if not os.path.isfile(tbl_file):
            tbl_file = download_eot_file()
            
    df = pd.read_csv(tbl_file)
    y = df[str(year)].dropna().values
    x = np.linspace(-0.5, len(y)-1.5, len(y))
    f = interpolate.interp1d(x, y, kind=interp_kind)
    return f


def get_eot(t, tbl_file=None, interp_kind='cubic'):
    """
    Equation of time for a given moment
    
    Arguments
    ---------
        t (datetime)     : time (UTC)
        tbl_file (str)   : path to csv file (table of daily values of EoT)
                           (Default is None, i.e. file will be downloaded)
        interp_kind (str): interpolation kind (linear, quadratic, cubic, etc.)

    Returns
    -------
        eot (float): equation of time in minutes

    Note: the csv file containing daily values of EoT from 2020 to 2050:
    https://raw.githubusercontent.com/behrouzz/eot/main/eot_2020_2050.csv
    """
    f = eot_func(year=t.year, tbl_file=tbl_file, interp_kind=interp_kind)
    equ_of_time = f(tim2ord(t)).flatten()[0]/60
    return equ_of_time


def eot_time_window(time_window, tbl_file=None, interp_kind='cubic'):
    """
    Equation of time for a given time window (interval)
    
    Arguments
    ---------
        time_window (array of datetimes) : list or numpy array of times (UTC)
        tbl_file (str)   : path to csv file (table of daily values of EoT)
                           (Default is None, i.e. file will be downloaded)
        interp_kind (str): interpolation kind (linear, quadratic, cubic, etc.)

    Returns
    -------
        eot_arr (array of floats): equation of time in minutes for all moments

    Note: the csv file containing daily values of EoT from 2020 to 2050:
    https://raw.githubusercontent.com/behrouzz/eot/main/eot_2020_2050.csv
    """
    time_window = pd.Series(time_window)
    years = time_window.dt.year.unique()
    eot_arr = []
    for yr in years:
        tw = time_window[time_window.dt.year==yr]
        ords = np.array([tim2ord(i) for i in tw])
        f = eot_func(yr, tbl_file, interp_kind)
        eot_arr = eot_arr + list(f(ords)/60)
    eot_arr = np.array(eot_arr)
    return eot_arr


def solar_time(t, lon, tbl_file=None):
    """
    Mean and True solar times
    
    Arguments
    ---------
        t (datetime)   : time in UTC
        lon (float)    : longtitude of observer
        tbl_file (str) : path to csv file (table of daily values of EoT)
                         (Default is None, i.e. file will be downloaded)
                             
    Returns
    -------
        mean_solar_time (datetime)
        true_solar_time (datetime)
    """
    mean_solar_time = t + timedelta(hours=(lon/15))
    equ_of_time = get_eot(t=t, tbl_file=tbl_file)
    equ_of_time = timedelta(minutes=equ_of_time)
    true_solar_time = mean_solar_time - equ_of_time
    return mean_solar_time, true_solar_time


def get_noon(t, lon, tbl_file=None):
    """
    Calculate the noon time for a given longtitude in UTC

    Arguments
    ---------
        t (str/datetime) : time (UTC)
        lon (float)      : longtitude
        tbl_file (str)   : path to csv file (table of daily values of EoT)
                           (Default is None, i.e. file will be downloaded)

    Returns
    -------
        noon (datetime): noon time in UTC
    """
    if isinstance(t, str):
        t = datetime.strptime(t[:10], '%Y-%m-%d')
    # mean solar time at noon (local in UTC):
    mean_st = datetime(t.year, t.month, t.day, 12) - \
              timedelta(hours=(lon/15))
    equ_of_time = get_eot(t=mean_st, tbl_file=tbl_file)
    equ_of_time = timedelta(minutes=equ_of_time)
    # true solar time at noon (local in UTC):
    true_st = mean_st - equ_of_time
    return true_st



