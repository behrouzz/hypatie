"""
Module hipparcos
================
This module retrieve data from Hipparcos Catalogue via Vizier database.
Hipparcos is an acronym for HIgh Precision PARallax COllecting Satellite.
The satellite launched in 1989 and operated until 1993.

The catalogue has 118218 data rows, and this module return all of them,
but just for the following rows:

field names
-----------
hip      : Identifier (HIP number) (H1)
app_mag  : Apparent magnitude in Johnson V (H5)
abs_mag  : Absolute magnitude computed by hypatie
distance : Distance in parsec computed by hypatie
parallax : Trigonometric parallax (H11)
B_V      : Johnson B-V colour (H37)
period   : Variability period (days) (H51)
var_type : [CDMPRU] variability type (H52)
sp_type  : Spectral type (H76)
ra       : Right Ascension (ICRS, epoch 2000) computed by CDS
dec      : Declination (ICRS, epoch 2000) computed by CDS

Note on var_type:
-----------------
var_type (H52) is the Hipparcos-defined type of variability.
It can be one of the following values:

' ' : could not be classified as variable or constant
'C' : no variability detected ("constant")
'D' : duplicity-induced variability
'M' : possibly micro-variable (amplitude < 0.03mag)
'P' : periodic variable
'R' : V-I colour index was revised due to variability analysis
'U' : unsolved variable which does not fall in the other categories
"""

import numpy as np
from urllib.request import urlopen
import json

def _hipparcos(filter_str=None, order_by=None, n_max=None):
    BASE_HIP = """http://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync?\
                request=doQuery&lang=adql&format=json&query=""".replace(' ','')

    filter_str = "" if filter_str is None else " WHERE "+filter_str
    order_by = "" if order_by is None else " ORDER BY "+order_by
    n_max = "" if n_max is None else "TOP "+str(n_max)+" "

    sql = 'SELECT ' + n_max
    sql = sql + 'HIP,Vmag,Plx,"B-V",Period,HvarType,SpType,"_RA.icrs","_DE.icrs"'
    sql = sql + ' FROM "I/239/hip_main"' + filter_str + order_by

    url = (BASE_HIP+sql).replace(' ', '%20')
    r = json.loads(urlopen(url).read().decode('utf-8'))
    field_names = [dc['name'] for dc in r['metadata']]
    rows = r['data']
    return field_names, rows


def hipparcos(filter_str=None, order_by=None, n_max=None):
    """
    Return stars from Hipparcos Catalogue

    Arguments
    ---------
        filter_str (str): filter the results. Default None.
        order_by (str): name of the filed to order the results by;
                        Default None (order by hip field).
        n_max (int): maximum number of rows to return;
                     default None (all results)

    Returns
    -------
        Two lists : data field names, rows of data

    ex:
    >>> fn, rows = hipparcos(filter_str='app_mag < 1', order_by='app_mag')
    >>> print(fn)
    >>> ['hip', 'app_mag', 'abs_mag', 'distance', 'parallax', 'B_V', 'period', 'var_type', 'sp_type', 'ra', 'dec']
    >>> print(rows[0])
    >>> [32349, -1.44, 16.4544, 0.00264, 379.21, 0.009, None, 'U', 'A0m...', 101.28715539, -16.71611582]
    """
    # translate field names from hypatie to vizier
    dc = {'hip':'HIP', 'app_mag':'Vmag', 'parallax':'Plx', 'B_V':'"B-V"',
          'period':'Period', 'var_type':'HvarType', 'sp_type':'SpType',
          'ra':'"_RA.icrs"', 'dec':'"_DE.icrs"'}
    if filter_str is not None:
        for k in dc.keys():
            if k in filter_str:
                filter_str = filter_str.replace(k, dc[k])
    if order_by is not None:
        if order_by in dc.keys():
            order_by = dc[order_by]
        else:
            raise Exception('field name is not valid!')

    field_names, rows = _hipparcos(filter_str, order_by, n_max)
    app_mag = np.array([i[1] for i in rows])
    
    plx = np.array([i[2] for i in rows])
    dist = []
    for i in plx:
        if i is None:
            dist.append(644997784)
        elif i==0:
            dist.append(644997784)
        else:
            dist.append(abs(1/i))
    dist = np.array(dist)
    
    abs_mag = np.where(dist!=644997784, app_mag + 5 - 5*np.log10(dist), np.nan)
    dist = np.where(dist!=644997784, dist, np.nan)
    
    rows = [v[:2] + [float("{:.5f}".format(abs_mag[i]))] + \
            [float("{:.5f}".format(dist[i]))] + v[2:] \
            for i,v in enumerate(rows)]
    
    field_names = ['hip','app_mag','abs_mag','distance','parallax','B_V',
                   'period','var_type','sp_type','ra','dec']
    return field_names, rows
