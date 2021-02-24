from datetime import datetime, timedelta
from urllib.request import urlopen
import numpy as np

# https://ssd.jpl.nasa.gov/horizons.cgi?show=1#results

def get_body_url(body, t1, t2, step):
    """Returns apparent RA and DEC of Sun in geocentric frame"""
    base = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&'
    params = f"""COMMAND='{body}'
    CENTER='500@399'
    MAKE_EPHEM='YES'
    TABLE_TYPE='VECTORS'
    START_TIME='{t1}'
    STOP_TIME='{t2}'
    STEP_SIZE='{step}'
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

def get_traject(body, t1, t2, step):
    url = get_body_url(body, t1, t2, step)
    with urlopen(url) as r:
        text = r.read().decode('utf-8')
    mark1 = text.find('$$SOE')
    text = text[mark1+6:]
    mark2 = text.find('$$EOE')
    text = text[:mark2]
    states = text.split('\n')[:-1]
    dates = [i.split(',')[1].strip()[5:] for i in states]
    pS = [i.split(',')[2:5] for i in states]
    ls = []
    for p in pS:
        tmp_ls = [float(i.strip())*1000 for i in p]
        ls.append(tmp_ls)
    return np.array(ls)
