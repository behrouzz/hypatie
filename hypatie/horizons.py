"""
Module horizons
===============
This module has two classes (Vector and Observer) to return positions
of the solar sytem objects by querying NASA's JPL Horizons API.

The Vector class provides positions with three dimensions: x,y,z.

The Observer class provides positions with two dimensions: ra,dec or
alt,az; depending on the request.
"""
import re
from datetime import datetime, timedelta
from urllib.request import urlopen
import numpy as np
from .plots import plot_altaz, plot_xyz

BASE_URL = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&'

def _time_format(t): #check time format
    correct = True
    if isinstance(t, datetime):
        t = t.isoformat()[:19]
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
        target (str): targer body
        t1 (datetime or str): time from which to request data
        t2 (datetime or str): time to which to request data
        step (int): number of time intervals
        center (str): origin of coordinates. format: site@body
                      default '500@0'
        ref_plane (str): refrence plane: 'FRAME', 'ECLIPTIC' or 'BODY EQUATOR'
                         default 'FRAME' (equator plane)
        vec_table (int): vector table type: 1 for position; 2 for position and velocity
    
    Attributes
    ----------
        time (list): list of datetime objects; each element is the time
                     at which the coordinates is presented
        pos (np.array): coordinates array with three columns (x, y, z)
                        shape of array is (n, 3); n: number of steps
        x (np.array): x component of position vector
        y (np.array): y component of position vector
        z (np.array): z component of position vector
        vel (np.array): velocity array with three columns (vx, vy, vz)
                        shape of array is (n, 3); n: number of steps
        vx (np.array): x component of velocity vector
        vy (np.array): y component of velocity vector
        vz (np.array): z component of velocity vector
        url (str): url of NASA's JPL Horizons that produced the results
    
    Methods
    -------
        plot: plot position(s) (x,y,z) of the target body
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

        error_msg, time, pos, vel = self.get_request()
        
        if len(error_msg)==0:
            if self.init:
                self.time = time[0]
                self.pos = pos[0]
                self.x, self.y, self.z = pos[0][0], pos[0][1], pos[0][2]
                if vel.shape!=(0,):
                    self.vel = vel[0]
                    self.vx, self.vy, self.vz = vel[0,0], vel[0,1], vel[0,2]
                else:
                    self.vel = vel
                
            else:
                self.time = time
                self.pos = pos
                self.x, self.y, self.z = pos[:,0], pos[:,1], pos[:,2]
                self.vel = vel
                if vel.shape!=(0,):
                    self.vx, self.vy, self.vz = vel[:,0], vel[:,1], vel[:,2]
                else:
                    self.vx, self.vy, self.vz = vel, vel, vel
        else:
            raise Exception('\n'+error_msg[:-2])

    def vector_url(self, target, t1, t2, step, center, ref_plane, vec_table):
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
        n_columns = len(rows[0].split(',')[:-1])
        
        times = [i.split(',')[1].strip()[5:] for i in rows]
        times = [datetime.strptime(i, '%Y-%b-%d %H:%M:%S.%f') for i in times]

        # positions
        pS = [i.split(',')[2:5] for i in rows]
        pos_list = []
        for p in pS:
            tmp_ls = [float(i.strip())*1000 for i in p]
            pos_list.append(tmp_ls)

        # velocities
        if n_columns==8:
            vS = [i.split(',')[5:8] for i in rows]
            vel_list = []
            for v in vS:
                tmp_ls = [float(i.strip())*1000 for i in v]
                vel_list.append(tmp_ls)
        else:
            vel_list = []
        return [error_msg, times, np.array(pos_list), np.array(vel_list)]

    def plot(self, color='b', size=10):
        """
        plot position(s) (ra/dec or az/alt) of the target body
        
        Arguments
        ----------
            color (str): color of body
            size (int): size of body
        
        Returns
        -------
            matplotlib.axes.Axes object
        """
        if isinstance(self.time, list):
            return plot_xyz(self.x, self.y, self.z, color, size)
        else:
            return plot_xyz([self.x], [self.y], [self.z], color, size)


class Observer:
    """
    Observer type coordinates (ra/dec or az/alt)
    
    Parameters
    ----------
        ttarget (str): targer body
        t1 (datetime or str): time from which to request data
        t2 (datetime or str): time to which to request data
        step (int): number of time intervals
        center (str): origin of coordinates. format: site@body
                      default '500@0'
        quantities (int): 1 for ra/dec, 2 for az/alt; (default 2)
    
    Attributes
    ----------
        time (list): list of datetime objects; each element is the time
                     at which the coordinates is presented
        pos (np.array): coordinates array with two columns (ra/dec or az/alt)
                        shape of array is (n, 2); n: number of steps
        ra (np.array): right ascension OR azimuth (depending on quantities chosen)
        dec (np.array): declination OR altitude (depending on quantities chosen)
        url (str): url of NASA's JPL Horizons that produced the results
    
    Methods
    -------
        plot: plot position(s) (ra/dec or az/alt) of the target body
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

    def coord(self, center):        
        valid_center = True
        str_center = ''

        if len(center.split('@'))==2:
            body = center.split('@')[-1]
        elif len(center.split('@'))==1:
            body = '399'
        else:
            valid_center = False

        crd = center.split('@')[0]
        if len(crd.split(','))==3:
            crd = ','.join(crd.split(',')[:-1]) + ',' + str(int(crd.split(',')[-1])/1000)
        elif len(crd.split(','))==2:
            crd = crd + ',0'
        else:
            valid_center = False

        if valid_center:
            str_center = f"coord@{body}'&COORD_TYPE='GEODETIC'&SITE_COORD='{crd}"

        return [valid_center, str_center]
    
    def observer_url(self, target, t1, t2, step, center, quantities, skip_daylight=True):
        if len(center.split(',')) > 1: #coordinates
            valid_center, center = self.coord(center)
            skip_daylight = False
            quantities = 4 # apparent Ra & DEC
            if not valid_center:
                raise Exception('Center is not valid!')

        skip_daylight='YES' if skip_daylight==True else 'NO'
        
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
        SKIP_DAYLT='{skip_daylight}'
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
        return [error_msg, times, np.array(pos_list)]

    def plot(self, color='b', size=10):
        """
        plot position(s) (ra/dec or az/alt) of the target body
        
        Arguments
        ----------
            color (str): color of body
            size (int): size of body
        
        Returns
        -------
            matplotlib.axes.Axes object
        """
        plot_altaz(self.ra, self.dec, color, size)
