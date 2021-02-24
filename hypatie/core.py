import re
from datetime import datetime, timedelta
from urllib.request import urlopen
import numpy as np


#https://ssd.jpl.nasa.gov/horizons.cgi?show=1#results
#center : geocentric='500@399', ssb='500@0'

def xyz_url(body, t1, t2, step, center='500@0'):
    """Returns cartesian vector of a body"""
    base = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&'
    params = f"""COMMAND='{body}'
    CENTER='{center}'
    MAKE_EPHEM='YES'
    TABLE_TYPE='VECTORS'
    START_TIME='{t1}'
    STOP_TIME='{t2}'
    STEP_SIZE='{step}'
    CAL_FORMAT='CAL'
    TIME_DIGITS='FRACSEC'
    OUT_UNITS='KM-S'
    REF_PLANE='FRAME'
    REF_SYSTEM='J2000'
    VECT_CORR='NONE'
    VEC_LABELS='NO'
    VEC_DELTA_T='NO'
    CSV_FORMAT='YES'
    OBJ_DATA='NO'
    VEC_TABLE='1'"""
    params = params.replace('\n', '&').replace(' ', '')
    url = base + params
    return url

def apparent_radec_url(body, t1, t2, step, center='500@399'):
    """Returns apparent RA, DEC of body"""
    base = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&'
    params = f"""COMMAND='{body}'
    CENTER='{center}'
    MAKE_EPHEM='YES'
    TABLE_TYPE='OBSERVER'
    START_TIME='{t1}'
    STOP_TIME='{t2}'
    STEP_SIZE='{step}'
    CAL_FORMAT='CAL'
    TIME_DIGITS='FRACSEC'
    ANG_FORMAT='DEG'
    OUT_UNITS='KM-S'
    RANGE_UNITS='KM'
    APPARENT='AIRLESS'
    SUPPRESS_RANGE_RATE='YES'
    SKIP_DAYLT='YES'
    EXTRA_PREC='YES'
    R_T_S_ONLY='NO'
    REF_SYSTEM='J2000'
    CSV_FORMAT='YES'
    OBJ_DATA='NO'
    QUANTITIES='2'"""
    params = params.replace('\n', '&').replace(' ', '')
    url = base + params
    return url



def get_position(base, body, t1, t2, step, center=None):
    """Returns positions of a body between t1 and t2.

    Keyword arguments:
    base   -- possible values: 'obs' for observer, 'vec' for vector
    body   -- name or number of a major body
    t1     -- starting time; datetime or str in format: '%Y-%m-%d %H:%M:%S'
    t2     -- starting time; datetime or str in format: '%Y-%m-%d %H:%M:%S'
    step   -- number of equal intervals between t1 and t2
    center -- center of frame, defalult '500@0' for 'vec', '500@399' for 'obs'

    Returns: [dates, positions]
    dates     -- list of dates for each step
    positions -- a numpy array with two (ra,dec) or three (x,y,z) columns
    """
    
    if isinstance(t1, datetime):
        t1 = t1.isoformat()
        t2 = t2.isoformat()
    elif isinstance(t1, str) and bool(re.match("\d{4}-\d\d-\d\d \d\d:\d\d:\d\d", t1)):
        t1 = t1.replace(' ', 'T')
        t2 = t2.replace(' ', 'T')
    else:
        raise Exception("t1 & t2 should be datetime or a string in format: '%Y-%m-%d %H:%M:%S'")
        
    if base=='vec':
        if center is None:
            center = '500@0'
        url = xyz_url(body, t1, t2, step, center)
    elif base=='obs':
        if center is None:
            center='500@399'
        url = apparent_radec_url(body, t1, t2, step, center)
    else:
        raise Exception("possible values for base: 'vec' for vector, 'obs' for observer")
    
    with urlopen(url) as r:
        text = r.read().decode('utf-8')

    mark1 = text.find('$$SOE')
    text = text[mark1+6:]
    mark2 = text.find('$$EOE')
    text = text[:mark2]
    states = text.split('\n')[:-1]
    
    if base=='vec':
        dates = [i.split(',')[1].strip()[5:] for i in states]
        pS = [i.split(',')[2:5] for i in states]
        ls = []
        for p in pS:
            tmp_ls = [float(i.strip())*1000 for i in p]
            ls.append(tmp_ls)
    elif base=='obs':
        dates = [i.split(',')[0].strip() for i in states]
        pS = [i.split(',')[-3:-1] for i in states]
        ls = []
        for p in pS:
            tmp_ls = [float(i.strip()) for i in p]
            ls.append(tmp_ls)
    else:
        raise Exception("base invalid")
    
    dates = [datetime.strptime(i, '%Y-%b-%d %H:%M:%S.%f') for i in dates]
    return [dates, np.array(ls)]
