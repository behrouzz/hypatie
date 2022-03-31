# https://maia.usno.navy.mil/ser7/readme.finals2000A

import os.path
from urllib.request import urlretrieve
from hypatie.time import datetime_to_jd
import numpy as np

FINALS2000A_URL = 'https://maia.usno.navy.mil/ser7/finals2000A.all'
FINALS2000A_FNAME = 'finals2000A.all'

if not os.path.exists(FINALS2000A_FNAME):
    print('Downloading '+FINALS2000A_FNAME+'...')
    urlretrieve(FINALS2000A_URL, FINALS2000A_FNAME)


with open(FINALS2000A_FNAME, 'r') as f:
    data = f.read().split('\n')

data = [i for i in data if len(i.strip())>68]
mjd = [int(i[7:12]) for i in data]
dut1 = [float(i[58:68]) for i in data]

dut1_array = np.array(list(zip(mjd,dut1)))

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

