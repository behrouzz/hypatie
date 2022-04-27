import os
from datetime import datetime
import numpy as np


def is_recent(fname):
    lm = os.path.getmtime(fname)
    lm = datetime.utcfromtimestamp(lm)
    now = datetime.utcnow()
    dt = (now - lm).total_seconds()
    recent = True if dt < 86400 else False
    return recent


def _time(t):
    if isinstance(t, datetime):
        return t
    elif isinstance(t, str) and bool(re.match("\d{4}-\d\d-\d\d \d\d:\d\d:\d\d", t)):
        return datetime.strptime(t, '%Y-%m-%d %H:%M:%S')
    else:
        raise Exception("only datetime or str: '%Y-%m-%d %H:%M:%S'")


def mag(x):
    """Returns magnitude of a vector"""
    return np.linalg.norm(np.array(x))

def unit(x):
    """Returns unit vector of a vector"""
    return x / mag(x)

def rev(x):
    return x % 360


def _in_vec(inp):
    if len(inp.shape)>1:
        if inp.shape[-1]==3:
            vec_out = inp[:,0], inp[:,1], inp[:,2]
        elif inp.shape[-1]==2:
            vec_out = inp[:,0], inp[:,1]
        else:
            raise Exception('Input imension not sopported!')
    elif len(inp.shape)==1:
        vec_out = inp
    else:
        raise Exception('Dimensionality problem! 2')
    return vec_out


def _out_vec(vec):
    if len(vec.shape)>1:
        vec = vec.T
        if vec.shape[-1]==3:
            x, y, z = vec[:,0], vec[:,1], vec[:,2]
            vec_out = np.vstack((x, y, z)).T
        elif vec.shape[-1]==2:
            x, y = vec[:,0], vec[:,1]
            vec_out = np.vstack((x, y)).T
        else:
            raise Exception('Dimensionality problem! 3')
    elif len(vec.shape)==1:
        vec_out = vec
    else:
        raise Exception('Dimensionality problem! 4')
    return vec_out
