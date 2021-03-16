"""
Module catalogues
=================
This module retrieves data from different catalogues via VizieR TAP Service.
You can find the current available catalogues with the dict:

>>> from hypatie import available_catalogues
>>> print(available_catalogues)
{'lqac5': 'Large Quasar Astrometric Catalogue (5th release)',
 'ucac4': 'U.S. Naval Observatory CCD Astrograph Catalog (4th release)',
 'hipparcos': 'HIgh Precision PARallax COllecting Satellite (Hipparcos)',
 'sdss12': 'SDSS Photometric Catalogue, (Release 12)',
 'glimpse': 'GLIMPSE Source Catalog (I + II + 3D) (IPAC 2008)',
 'gaia2': 'Gaia data release 2',
 'gaiae3': 'Gaia early data release 3'}
"""
from urllib.request import urlopen
import pandas as pd
import json

BASE = """http://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync?\
request=doQuery&lang=adql&format=json&query=""".replace(' ','')

# Catalogue names
lqac5_nm = 'Large Quasar Astrometric Catalogue (5th release)'
ucac4_nm = 'U.S. Naval Observatory CCD Astrograph Catalog (4th release)'
hipparcos_nm = 'HIgh Precision PARallax COllecting Satellite (Hipparcos)'
sdss12_nm = 'SDSS Photometric Catalogue, (Release 12)'
glimpse_nm = 'GLIMPSE Source Catalog (I + II + 3D) (IPAC 2008)' #GLIMPSE sources in Galactic Center (104240613 rows)
gaia2_nm = 'Gaia data release 2'
gaiae3_nm = 'Gaia early data release 3'
dfgrs2_nm = '2dF Galaxy Redshift Survey'
glade2_nm = 'GLADE v2.3 catalog (Galaxy List for the Advanced Detector Era)'

# Table names
lqac5_tb = 'J/A%2bA/624/A145/lqac5'
ucac4_tb = 'I/322A/out'
hipparcos_tb = 'I/239/hip_main'
sdss12_tb = 'V/147/sdss12'
glimpse_tb = 'II/293/glimpse'
gaia2_tb = 'I/345/gaia2'
gaiae3_tb = 'I/350/gaiaedr3'
dfgrs2_tb = 'VII/250/2dfgrs'
glade2_tb = 'VII/281/glade2'

# Field names
lqac5_fn = ['LQAC','RAJ2000','DEJ2000','Plx','pmRA*','pmDE','Ref','GAIA',
            'z','Gmag2','BPmag','RPmag','umag','Bmag','Vmag','gmag',
            'rmag','imag','zmag','Jmag','Kmag','W1mag','W2mag','W3mag',
            'W4mag','S1.4','S2.3','S5.0','S8.4','S23','BMAG','IMAG','uMAG',
            'gMAG','rMAG','iMAG','zMAG']
ucac4_fn = ['UCAC4','RAJ2000','DEJ2000','pmRA','pmDE','2Mkey','f.mag',
            'a.mag','Jmag','Hmag','Kmag','Bmag','Vmag','gmag','rmag',
            'imag','H']
hipparcos_fn = ['HIP','_RA.icrs','_DE.icrs','Vmag','Plx','e_Plx',
                'pmRA','pmDE','B-V','Period','HvarType','SpType']
sdss12_fn = ['RA_ICRS','DE_ICRS','class','objID','SDSS12','ObsDate','Q',
             'zsp','e_zsp','zph','e_zph','<zph>','spType','spCl']
glimpse_fn = ['S','GLIMPSE','2MASS','Glon','Glat','RAJ2000','DEJ2000',
              'Jmag','Hmag','Kmag','3.6mag','4.5mag','5.8mag','8.0mag']
gaia2_fn = ['designation','source_id','random_index','ra_epoch2000',
            'dec_epoch2000','l','b','parallax','lum_val','radius_val',
            'pmra','pmdec','radial_velocity','teff_val','phot_rp_mean_mag',
            'phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_flux',
            'phot_g_mean_flux','phot_bp_mean_flux','a_g_val',
            'phot_bp_rp_excess_factor']
gaiae3_fn = ['EDR3Name','RAJ2000','DEJ2000','Epoch','Source','RandomI',
             'Plx','e_Plx','pmRA','pmDE','NAL','NAC','MatchObs','FG','e_FG',
             'FGCorr','Gmag','e_Gmag','GmagCorr','FBP','e_FBP','BPmag',
             'e_BPmag','FRP','e_FRP','RPmag','e_RPmag','RVDR2','e_RVDR2',
             'Tefftemp','GLON','GLAT']
dfgrs2_fn = ['SeqNum','Name','RAJ2000','DEJ2000','Bjmag','Bjsel','z','q_z',
             'n_z','z.em','o_z.em','SNR']
glade2_fn = ['PGC','GWGC','HyperLEDA','2MASS','SDSS-DR12','Flag1','RAJ2000',
             'DEJ2000','Dist','e_Dist','z','Bmag','e_Bmag','BMAG','Jmag',
             'e_Jmag','Hmag','e_Hmag','Kmag','e_Kmag','Flag2','Flag3']

dc = {'lqac5'     : [lqac5_tb, lqac5_fn],
      'ucac4'     : [ucac4_tb, ucac4_fn],
      'hipparcos' : [hipparcos_tb, hipparcos_fn],
      'sdss12'    : [sdss12_tb, sdss12_fn],
      'glimpse'   : [glimpse_tb, glimpse_fn],
      'gaia2'     : [gaia2_tb, gaia2_fn],
      'gaiae3'    : [gaiae3_tb, gaiae3_fn],
      'dfgrs2'    : [dfgrs2_tb, dfgrs2_fn],
      'glade2'    : [glade2_tb, glade2_fn]}

long_names = [lqac5_nm, ucac4_nm, hipparcos_nm, sdss12_nm,
              glimpse_nm, gaia2_nm, gaiae3_nm, dfgrs2_nm, glade2_nm]

# dictionary containing names of available catalogues
available_catalogues = {k:v for (k,v) in [*zip(dc.keys(), long_names)]}

class Catalogue:
    """
    Catalogue class

    Parameters
    ----------
        name (str)     : name of the catalogue
        columns (list) : columns of data to be requested
        where (str)    : condition (WHERE clause in SQL)
        order_by (str) : column to order (ORDER BY clause in SQL)
        n_max (int)    : maximum number of rows to return; dfault 100.
        random_sampling (bool): return a rondom sample of data.
            Note that huge datasets can create Gateway timeout error
            using this parameter; try to reduce the size of data using
            'where' in the case of Gateway timeout error.

    Attributes
    ----------
        tbl_name        : name of the table
        columns         : columns of the table
        random_sampling : whether or not 'random_sampling' is True
        sql             : ADQL script to produce data
        url             : URL of VizieR to produce data

    Methods
    -------
        available_columns : return all available possible columns in the
                            table of the catalogue
        download : download data from the catalogue
    """
    def __init__(self,
                 name,
                 columns=None,
                 where=None,
                 order_by=None,
                 n_max=100,
                 random_sampling=False):
        
        self.tbl_name = dc[name][0]
        self.random_sampling = random_sampling
        
        if columns is None:
            self.columns = dc[name][1]
        elif isinstance(columns, list):
            self.columns = columns
        elif columns=='all':
            self.columns = list(self.available_columns().index)

        if where is not None:
            where_words = where.split(' ')
            where_cols = [i for i in where_words if i in \
                          list(self.available_columns().index)]
            if len(where_cols)==0:
                raise Exception("column(s) defined is 'where' is/are not valid!")
            _where = []
            for word in where_words:
                if word in where_cols:
                    _where.append('"'+word+'"')
                else:
                    _where.append(word)
            where = ' WHERE ' + ' '.join(_where)
        else:
            where = ''

        if order_by is None:
            order_by = ''
        else:
            if ' desc' in order_by.lower():
                order_by = ' ORDER BY "'+order_by[:-5]+'"' + " DESC"
            else:
                order_by = ' ORDER BY "'+order_by+'"'
        
        if self.random_sampling and isinstance(n_max, int):
            rnd = 'RAND(1) as rnd_sample,'
            order_by = ' ORDER BY rnd_sample'
        else:
            rnd = ''
        if isinstance(n_max, int):
            n_max = f" TOP {n_max} "
        elif n_max=='all':
            n_max = ''
        else:
            raise Exception("'n_max' should be integer or 'all'")

        sql_columns = [f'"{i}"' for i in self.columns]
        self.sql = 'SELECT' + n_max + rnd + ','.join(sql_columns) + \
                   ' FROM "' + self.tbl_name+'"' + where + order_by
        self.url = (BASE+self.sql).replace(' ', '%20')

    def available_columns(self):
        """
        Get all available columns in the main table of the catalogue

        Returns
        -------
            pandas DataFrame of available columns and descriptions

        Ex:
        >>> cat = Catalogue('gaia2')
        >>> print(cat.available_columns())
        """
        
        sql = 'SELECT column_name, description, unit FROM tap_schema.columns'
        sql = sql + f" WHERE table_name='{self.tbl_name}'"
        url = (BASE+sql).replace(' ', '%20')
        r = json.loads(urlopen(url).read().decode('utf-8'))
        cols = [i[0] for i in r['data']]
        desc = [i[1] for i in r['data']]
        df = pd.DataFrame(r['data'],
                columns=['column','description','unit'])
        return df.set_index('column')

    def download(self):
        """
        Download data from the main table of the catalogue

        Returns
        -------
            data : pandas DataFrame of the main data table
            meta : pandas DataFrame of the meta-data table

        Ex:
        >>> cat = Catalogue('gaia2')
        >>> data, meta = cat.download()
        """
        r = json.loads(urlopen(self.url).read().decode('utf-8'))
        cols = [i['name'] for i in r['metadata']]
        data = pd.DataFrame(r['data'], columns=cols)
        meta = pd.DataFrame(r['metadata'])
        if self.random_sampling:
            data = data.drop('rnd_sample', axis=1)
            meta = meta.drop([0])
        meta = meta[['name', 'description', 'unit']].set_index('name')
        return data, meta


