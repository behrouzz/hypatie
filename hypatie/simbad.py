"""
Module simbad
=============
This module deals with deep-sky objects by querying SIMBAD astronomical
database from Centre de données astronomiques de Strasbourg (CDS).

Note about otype
----------------
Several functions in this module have a parameter 'otype'. This is the
type of deep-sky object defined by SIMBAD in this url:
http://simbad.u-strasbg.fr/simbad/sim-display?data=otypes

This a hierarchical classification. For exact type of an object, you
shoud enter its numeric code as a string, such as '14.09.05.03' for
Cepheid variable Star. You can also enter '14.09.05' for all Pulsating
variable Stars (including cepheids and other Pulsating variable Stars
subclasses). If you enter '14.09.00.00' or '14.09' (they are same), it
means all types of variable stars.

For simplicity, we have also provided 7 string values ('radiation',
'gravitation', 'candidate', 'multiple', 'interstellar', 'star' and
'galaxy'). Theses are general categories which encompass many subcategories. 
"""
import numpy as np
import pandas as pd
from urllib.request import urlopen
import requests
import json, csv, re
from io import StringIO
from datetime import datetime
from .simbad_otypes import text
from .transform import altaz_to_radec, radec_to_altaz


data = StringIO(text)

code, short, long = [], [], []

for row in csv.reader(data):
    code.append(row[0])
    short.append(row[1])
    long.append(row[2])

_dc = {k:v for (k,v) in zip(code, short)}
otype_dict = {k:v for (k,v) in zip(code, long)}

def object_type(num):
    
    if num[0].isalpha():
        if num.lower()=='radiation':
            num = '999'
        elif num.lower()=='gravitation':
            num = '09'
        elif num.lower()=='candidate':
            num = '10'
        elif num.lower()=='multiple':
            num = '12'
        elif num.lower()=='interstellar':
            num = '13'
        elif num.lower()=='star':
            num = '14'
        elif num.lower()=='galaxy':
            num = '15'
            
    if bool(re.match("\d\d.\d\d.\d\d.\d\dZ", num+'Z')):
        otypes = [_dc[num]]
    elif bool(re.match("\d\d.\d\d.\d\dZ", num+'Z')):
        otypes = [_dc[i] for i in _dc.keys() if i[:8]==num]
    elif bool(re.match("\d\d.\d\dZ", num+'Z')):
        otypes = [_dc[i] for i in _dc.keys() if i[:5]==num]
    elif bool(re.match("\d\dZ", num+'Z')):
        otypes = [_dc[i] for i in _dc.keys() if i[:2]==num]
    elif num=='999':
        otypes = ['Rad','mR','cm','mm','smm','HI','rB','Mas','IR','FIR',
                  'MIR','NIR','red','ERO','blu','UV','X','UX?','ULX',
                  'gam', 'gB']
    else:
        raise Exception('Object type invalid!')
            
    return otypes


def bright_objects(otype=None):
    """
    Return Deep-sky objects with apparent magnitudes less than 4.5

    Arguments
    ---------
        otype (str): type of object. Default None of all objects

    Returns
    -------
        Two lists : data field names, rows of data
    """
    
    from .simbad_bright_objects import mag45
    f_mag45 = StringIO(mag45)
    rows_mag45 = []
    for row in csv.reader(f_mag45):
        rows_mag45.append(row)
    rows = [[i[0], i[1], float(i[2]), float(i[3]), float(i[4])] \
          for i in rows_mag45]
    if otype is not None:
        otype = object_type(otype)
        rows = [i for i in rows if i[1] in otype]
    field_names = ['main_id', 'otype_txt', 'V', 'ra', 'dec']
    return field_names, rows
    


BASE_SIMBAD = 'http://simbad.u-strasbg.fr/simbad/sim-tap/sync?\
request=doQuery&lang=adql&format=json&query='


def get_objects(otype=None, n_max=1000):
    """
    Return deep-sky objects

    Arguments
    ---------
        otype (str): type of object. Default None of all objects
        n_max (int): maximum number of rows to return; default 1000

    Returns
    -------
        Two lists : data field names, rows of data
    """
    
    if otype is not None:
        otype = str(tuple(object_type(otype)))
        add_otype = f" WHERE otype_txt in {otype}"
    else:
        add_otype = ""
    
    sql = f"""SELECT TOP {n_max} main_id, otypedef.otype_longname, allfluxes.V, ra, dec
    FROM basic LEFT JOIN allfluxes ON basic.oid=allfluxes.oidref
    LEFT JOIN otypedef ON basic.otype=otypedef.otype"""+add_otype+" ORDER BY V"
    sql = ' '.join([i.strip() for i in sql.split('\n')])
    url = (BASE_SIMBAD+sql).replace(' ', '%20')
    r = json.loads(urlopen(url).read().decode('utf-8'))
    field_names = [dc['name'] for dc in r['metadata']]
    rows = r['data']
    return field_names, rows


def explore_region(ra, dec, r, n_max=1000):
    """
    Explore a circular region in the sky with ICRS coordinates

    Arguments
    ---------
        ra (float): right ascension of the cirular region
        dec (float): declination of the cirular region
        r (float): radius of the cirular region in degrees
        n_max (int): maximum number of rows to return; default 1000

    Returns
    -------
        Two lists : data field names, rows of data
    """
    return search_region(ra=ra, dec=dec, r=r, otype=None, n_max=n_max)

def search_region(ra, dec, r, otype=None, n_max=1000):
    """
    Explore a circular region in the sky with ICRS coordinates

    Arguments
    ---------
        ra (float): right ascension of the cirular region
        dec (float): declination of the cirular region
        r (float): radius of the cirular region in degrees
        otype (str): type of object. Default None of all objects
        n_max (int): maximum number of rows to return; default 1000

    Returns
    -------
        Two lists : data field names, rows of data
    """
    if otype is not None:
        otype = str(tuple(object_type(otype)))
        add_otype = f" AND otype_txt in {otype}"
    else:
        add_otype = ""
    sql = f"""SELECT TOP {n_max} main_id, otypedef.otype_longname, allfluxes.V, ra, dec
    FROM basic LEFT JOIN allfluxes ON basic.oid=allfluxes.oidref
    LEFT JOIN otypedef ON basic.otype=otypedef.otype
    WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {r})) = 1
    AND ra IS NOT NULL AND dec IS NOT NULL"""+add_otype+" ORDER BY V"
    sql = ' '.join([i.strip() for i in sql.split('\n')])
    url = (BASE_SIMBAD+sql).replace(' ', '%20')
    r = json.loads(urlopen(url).read().decode('utf-8'))
    field_names = [dc['name'] for dc in r['metadata']]
    rows = r['data']
    return field_names, rows

def explore_sky(lon, lat, t=None, az=0, alt=90, r=1, n_max=1000):
    """
    Explore a circular region in the local sky with AltAz coordinates

    Arguments
    ---------
        lon (float): longtitude of observer location
        lat (float): latitude of observer location
        t (datetime or str): time of observation in UTC. default now.
        az (float): azimuth of center of the circle
        alt (float): altitude of center of the circle
        r (float): radius of the cirular region in degrees
        n_max (int): maximum number of rows to return; default 1000

    Returns
    -------
        Two lists : data field names, rows of data
    """
    return search_sky(lon, lat, t, az, alt, r, otype=None, n_max=1000)

def search_sky(lon, lat, t=None, az=0, alt=90, r=1, otype=None, n_max=1000):
    """
    search in a circular region in the local sky with AltAz coordinates

    Arguments
    ---------
        lon (float): longtitude of observer location
        lat (float): latitude of observer location
        t (datetime or str): time of observation in UTC. default now.
        az (float): azimouth of center of the circle
        alt (float): altitude of center of the circle
        r (float): radius of the cirular region in degrees
        otype (str): type of object. Default None of all objects
        n_max (int): maximum number of rows to return; default 1000

    Returns
    -------
        Two lists : data field names, rows of data
    """

    if t is None:
        t = datetime.utcnow()
    elif isinstance(t, datetime):
        pass
    elif isinstance(t, str) and bool(re.match("\d{4}-\d\d-\d\d \d\d:\d\d:\d\d", t)):
        t = datetime.strptime(t, '%Y-%m-%d %H:%M:%S')
    else:
        raise Exception("t should be a datetime or str: '%Y-%m-%d %H:%M:%S'")
    
    ra, dec = altaz_to_radec(lon, lat, az, alt, t)
    return search_region(ra=ra, dec=dec, r=r, otype=otype, n_max=n_max)

def sql2df(script):
    BASE = 'http://simbad.u-strasbg.fr/simbad/sim-tap/sync?request=doQuery&lang=adql&format=csv&query='
    script = ' '.join(script.strip().split('\n'))
    url = BASE+script.replace(' ', '%20') + '&format=csv'
    req = requests.request('GET', url)
    r = req.content.decode('utf-8')
    lines = r.splitlines()
    col = lines[0].split(',')
    data_lines = [i.split(',') for i in lines[1:]]
    return pd.DataFrame(data_lines, columns=col)
