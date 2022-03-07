from datetime import datetime, timedelta
import math

d2r = math.pi/180
LP = 37


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
  

def utc2tt(t, lp=LP):
    TAI = t + timedelta(seconds=lp)
    TT = TAI + timedelta(seconds=32.184)
    return TT

def tt2utc(t, lp=LP):
    TAI = t - timedelta(seconds=32.184)
    utc = TAI - timedelta(seconds=lp)
    return utc

def tt2tdb(t, lp=LP):
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

def utc2tdb(utc, lp=LP):
    tt = utc2tt(utc, lp)
    tdb = tt2tdb(tt, lp)
    return tdb

def tdb2utc(tdb, lp=LP):
    tt = tdb2tt(tdb)
    utc = tt2utc(tt, lp=lp)
    return utc
