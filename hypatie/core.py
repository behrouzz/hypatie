import re
from datetime import datetime, timedelta
from urllib.request import urlopen
import numpy as np

BASE_URL = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&'

#https://ssd.jpl.nasa.gov/horizons.cgi?show=1#results
#center : geocentric='500@399', ssb='500@0'

def _time_format(t):
    correct = True
    if isinstance(t, datetime):
        t = t.isoformat()
    elif isinstance(t, str) and bool(re.match("\d{4}-\d\d-\d\d \d\d:\d\d:\d\d", t)):
        t = t.replace(' ', 'T')
    else:
        correct = False
    return [correct, t]
    

class Vector:
    """
    Vector type coordinates (x, y, z)
    
    Parameters
    ----------
    target : str
        targer body
    t1 : datetime or str in format '%Y-%m-%d %H:%M:%S'
        time from which to request data
    t2 : datetime or str in format '%Y-%m-%d %H:%M:%S'
        time to which to request data
    steps: int
        number of time intervals
    center: str
        origin of coordinates. format: site@body (default '500@0')
    reg_plane: str
        refrence plane. can be 'FRAME', 'ECLIPTIC' or 'BODY EQUATOR' (default 'FRAME')
    vec_table: int:
        vector table type: enter 1 for position; 2 for position and velocity
    
    Attributes
    ----------
    time : datetime
        time at which the coordinates is presented
    pos : numpy array with shape (n, 3); n: number of steps
        coordinates array with three columns (x, y, z)
    x : numpy array with shape (n,); n: number of steps
        x of cartesian coordinates 
    y : numpy array with shape (n,); n: number of steps
        y of cartesian coordinates
    z : numpy array with shape (n,); n: number of steps
        z of cartesian coordinates
    url : str
        url of NASA's JPL Horizons that procuded the results 
    """

    def __init__(self, target, t1, t2=None, step=None, center='500@0', ref_plane='FRAME', vec_table=1):
        self.target = target
        self.center = center
        self.ref_plane = ref_plane
        self.vec_table = vec_table
        
        correct, t1 = _time_format(t1)
        if correct:
            self.t1 = t1
        else:
            raise Exception("t1 must be datetime or str: '%Y-%m-%d %H:%M:%S'")
        self.init = False
        if t2 is None:
            self.init = True # user just wants the initial state
            step = 1
            t2 = datetime.strptime(self.t1, '%Y-%m-%dT%H:%M:%S') + timedelta(seconds=1)
        
        correct, t2 = _time_format(t2)
        if correct:
            self.t2 = t2
        else:
            raise Exception("t2 must be datetime or str: '%Y-%m-%d %H:%M:%S'")
            
        self.step = step

        error_msg, time, pos = self.get_request()
        
        if len(error_msg)==0:
            if self.init:
                self.time = time[0]
                self.pos = pos[0]
                self.x, self.y, self.z = pos[0][0], pos[0][1], pos[0][2]
            else:
                self.time = time
                self.pos = pos
                self.x, self.y, self.z = pos[:,0], pos[:,1], pos[:,2]
        else:
            raise Exception('\n'+error_msg[:-2])

    def vector_url(self, target, t1, t2, step, center, ref_plane, vec_table):
        """Returns cartesian vector of a target body"""
        params = f"""COMMAND='{target}'
        CENTER='{center}'
        MAKE_EPHEM='YES'
        TABLE_TYPE='VECTORS'
        START_TIME='{t1}'
        STOP_TIME='{t2}'
        STEP_SIZE='{step}'
        CAL_FORMAT='CAL'
        TIME_DIGITS='FRACSEC'
        OUT_UNITS='KM-S'
        REF_PLANE='{ref_plane}'
        REF_SYSTEM='J2000'
        VECT_CORR='NONE'
        VEC_LABELS='NO'
        VEC_DELTA_T='NO'
        CSV_FORMAT='YES'
        OBJ_DATA='NO'
        VEC_TABLE='{vec_table}'"""
        params = params.replace('\n', '&').replace(' ', '')
        url = BASE_URL + params
        self.url = url
        return url

    def get_request(self):
        error_msg = ''
        url = self.vector_url(self.target, self.t1, self.t2, self.step,
                              self.center, self.ref_plane, self.vec_table)
        with urlopen(url) as r:
            text = r.read().decode('utf-8')
        if ('$$SOE' not in text) or ('$$EOE' not in text):
            error_msg = text[:text.find('$$SOF')]
            return [error_msg, np.array([]), np.array([])]
        mark1 = text.find('$$SOE')
        text = text[mark1+6:]
        mark2 = text.find('$$EOE')
        text = text[:mark2]
        rows = text.split('\n')[:-1]
        times = [i.split(',')[1].strip()[5:] for i in rows]
        times = [datetime.strptime(i, '%Y-%b-%d %H:%M:%S.%f') for i in times]
        pS = [i.split(',')[2:5] for i in rows]
        pos_list = []
        for p in pS:
            tmp_ls = [float(i.strip())*1000 for i in p]
            pos_list.append(tmp_ls)
        return [error_msg, np.array(times), np.array(pos_list)]


class Observer:
    """
    Observer type coordinates (RA and DEC)
    
    Parameters
    ----------
    target : str
        targer body
    t1 : datetime or str in format '%Y-%m-%d %H:%M:%S'
        time from which to request data
    t2 : datetime or str in format '%Y-%m-%d %H:%M:%S'
        time to which to request data
    steps: int
        number of time intervals
    center: str
        origin of coordinates. format: site@body (default '500@399')
    quantities: int:
        enter 1 for astrometric RA & DEC and 2 for apparent RA & DEC (default 2)
    
    Attributes
    ----------
    time : datetime
        time at which the coordinates is presented
    pos : numpy array with shape (n, 2); n: number of steps
        coordinates array with two columns (RA & DEC)
    ra : numpy array with shape (n,); n: number of steps
        right ascension
    dec : numpy array with shape (n,); n: number of steps
        declination
    url : str
        url of NASA's JPL Horizons that procuded the results 
    """
    def __init__(self, target, t1, t2=None, step=None, center='500@399', quantities=2):
        self.target = target
        self.center = center
        self.quantities = quantities
        
        correct, t1 = _time_format(t1)
        if correct:
            self.t1 = t1
        else:
            raise Exception("t1 must be datetime or str: '%Y-%m-%d %H:%M:%S'")
        self.init = False
        if t2 is None:
            self.init = True # user just wants the initial state
            step = 1
            t2 = datetime.strptime(self.t1, '%Y-%m-%dT%H:%M:%S') + timedelta(seconds=1)
        
        correct, t2 = _time_format(t2)
        if correct:
            self.t2 = t2
        else:
            raise Exception("t2 must be datetime or str: '%Y-%m-%d %H:%M:%S'")
            
        self.step = step

        error_msg, time, pos = self.get_request()
        
        if len(error_msg)==0:
            if self.init:
                self.time = time[0]
                self.pos = pos[0]
                self.ra, self.dec = pos[0][0], pos[0][1]
            else:
                self.time = time
                self.pos = pos
                self.ra, self.dec = pos[:,0], pos[:,1]
        else:
            raise Exception('\n'+error_msg[:-2])

    def observer_url(self, target, t1, t2, step, center, quantities):
        """Returns RA, DEC of target body"""
        params = f"""COMMAND='{target}'
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
        QUANTITIES='{quantities}'"""
        params = params.replace('\n', '&').replace(' ', '')
        url = BASE_URL + params
        self.url = url
        return url

    def get_request(self):
        error_msg = ''
        url = self.observer_url(self.target, self.t1, self.t2, self.step,
                                self.center, self.quantities)
        with urlopen(url) as r:
            text = r.read().decode('utf-8')
        if ('$$SOE' not in text) or ('$$EOE' not in text):
            error_msg = text[:text.find('$$SOF')]
            return [error_msg, np.array([]), np.array([])]
        mark1 = text.find('$$SOE')
        text = text[mark1+6:]
        mark2 = text.find('$$EOE')
        text = text[:mark2]
        rows = text.split('\n')[:-1]
        times = [i.split(',')[0].strip() for i in rows]
        times = [datetime.strptime(i, '%Y-%b-%d %H:%M:%S.%f') for i in times]
        pS = [i.split(',')[-3:-1] for i in rows]
        pos_list = []
        for p in pS:
            tmp_ls = [float(i.strip()) for i in p]
            pos_list.append(tmp_ls)
        return [error_msg, np.array(times), np.array(pos_list)]
