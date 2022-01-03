"""
Module plot
================
This module supplies several functions for plotting reasons.
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections.abc import Iterable
from urllib.request import urlopen
from .simbad import bright_objects, object_type, sql2df
from .transform import radec_to_altaz, altaz_to_radec, angular_separation

def _equalize_scale(X,Y,Z, ax):
    """
    Equalize scales of x, y and z axis

    Arguments
    ---------
        X (np.array): x values
        Y (np.array): y values
        Z (np.array): z values
        ax (axes): axes object

    Returns
    -------
        matplotlib axes object
    """
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
    for xb, yb, zb in zip(Xb, Yb, Zb):
       ax.plot([xb], [yb], [zb], 'w')
    return ax

def plot_xyz(x, y, z, color, size, au=False):
    """
    Draw 3d cartesian plot of a body

    Arguments
    ---------
        x (np.array): x values
        y (np.array): y values
        z (np.array): z values
        color (str): color
        size (int): size
        au (bool): whether or not in AU unit; default False.

    Returns
    -------
        matplotlib axes object
    """
    fig = plt.figure(figsize=plt.figaspect(0.5)*1.2)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='z', labelsize=8)
    au = 149597870700.0 if au else 1
    ax.scatter(x/au, y/au, z/au, c=color, s=size)
    ax = _equalize_scale(x/au, y/au, z/au, ax)
    return ax

def plot_altaz(az, alt, mag=None, size=None, color='k', alpha=1, marker='o', ax=None):
    """
    Plot positions of bodies based on altitude/azimuth coordinates

    Arguments
    ---------
        az (iter of float): azimuth values
        alt (iter of float): altitude values
        mag (iter of float): apparent magnitudes; default None.
        size (int): size; default None.
        color (str): color; default 'k'.
        alpha (float): alpha value (transparency), between 0 and 1; default 1.
        marker (str): marker shape; default 'o'.
        ax (axes): axes object; default None.

    Returns
    -------
        matplotlib axes object
    """

    if isinstance(az, Iterable):
        az = np.array(az)
        alt = np.array(alt)
        mag = np.array(mag) if mag is not None else None
    else:
        az = np.array([az])
        alt = np.array([alt])
        mag = None
        
    az  = az*(np.pi/180)

    if size is None:
        if mag is None:
            size = [5] * len(az)
        else:
            size = (0.1 + max(mag)-mag)**1.8
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.8,0.8], polar=True)
        ax.set_theta_zero_location('N')
        ax.set_rlim(90, 0, 1)
        ax.set_yticks(np.arange(0,91,30))
        ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi, 5*np.pi/4, 6*np.pi/4, 7*np.pi/4])
        ax.set_xticklabels(['N','NE','E','SE','S','SW','W','NW'])
        ax.tick_params(axis=u'both', which=u'both',length=0)
    if matplotlib.__version__ < '3.0.0':
        alt = [90-i for i in alt]
    ax.scatter(az, alt, c=color, s=size, alpha=alpha, marker=marker)
    ax.grid(True, alpha=0.7)
    return ax

def plot_radec(ra, dec, mag=None, size=None):
    """
    Plot positions of bodies based on RA/Dec coordinates

    Arguments
    ---------
        ra (iter of float): Right Ascension values
        dec (iter of float): Declination values
        mag (iter of float): apparent magnitudes; default None.
        size (int): size; default None.

    Returns
    -------
        matplotlib axes object
    """
        
    if isinstance(ra, Iterable):
        ra = np.array(ra)
        dec = np.array(dec)
        mag = np.array(mag) if mag is not None else None
    else:
        ra = np.array([ra])
        dec = np.array([dec])
        mag = None

    ra = (np.mod(ra-180, 360)-180) * (np.pi/180)
    dec = dec * (np.pi/180)

    if size is None:
        if mag is None:
            size = [5] * len(ra)
        else:
            size = (0.1 + max(mag)-mag)**1.8

    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111, projection="aitoff")
    ax.grid(True)
    ax.scatter(ra, dec, s=size, c='k')
    plt.subplots_adjust(top=0.95,bottom=0.0)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    return ax

def star_chart(lon, lat, t=None, otype=None):
    """
    Plot the star chart for a location and time on earth

    Arguments
    ---------
        lon (float): longitude of observer
        lat (float): latitude of observer
        t (datetime or str): time of observation; default None.
        otype (str): object type to be presented in the chart; default None.

    Returns
    -------
        matplotlib axes object
    """
    _, rows = bright_objects(otype)
    mag = [i[2] for i in rows]
    ra  = [i[3] for i in rows]
    dec = [i[4] for i in rows]
    alt, az = radec_to_altaz(lon, lat, ra, dec, t)
    ax = plot_altaz(az, alt, mag=mag)
    return ax

class Telescope:
    """
    Plot image of a location in the sky by querying NASA's Virtual
    Observatory (SkyView)

    Arguments
    ---------
        target_loc (tuple): location of a target point in the sky;
                            (ra, dec) or (az, alt)
        obs_loc (tuple): location of observer on Earth; format: (lon, lat)
        fov (float) : Field of View in degrees; default 1
        t (datetime or str): time of observation; default now
        survey (str): name of survey; default 'digitized sky survey'

    Attributes
    ----------
        data (numpy array): image data
        url (str) : image url

    Methods
    -------
        plot : plots the image and returns fig and ax
        show : shows the image
    """
    def __init__(self,
                 target_loc,
                 obs_loc=None,
                 fov=1,
                 t=None,
                 survey='digitized sky survey'):
        
        survey = survey.replace(' ', '+')

        if obs_loc is not None:
            lon, lat = obs_loc
            az , alt = target_loc
            ra, dec = altaz_to_radec(lon, lat, az, alt, t)
            ra, dec = float(ra), float(dec)
            position = str((ra, dec)).replace(' ', '')[1:-1]
        else:
            position = str(target_loc).replace(' ', '')[1:-1]

        BASE = 'https://skyview.gsfc.nasa.gov/cgi-bin/images?'
        params = f'Survey={survey}&position={position}&Size={fov}&Return=JPEG'
        self.url = BASE + params
        f = urlopen(self.url)
        
        self.data = plt.imread(f, format='jpeg')

    def plot(self):
        """
        Plots the image

        Returns
        -------
            fig, ax
        """
        fig, ax = plt.subplots()
        ax.imshow(self.data, cmap='gray')
        ax.set_xticks([])
        ax.set_yticks([])
        return fig, ax

    def show(self):
        """Shows the image"""
        fig, ax = self.plot()
        plt.show()

def plot_pm(ra, dec, pmra, pmdec, alpha=0.5, color=None):
    """
    Plot the proper motion of objects

    Arguments
    ---------
        ra (iterable of floats): Array of RA
        dec (iterable of floats): Array of DEC
        pmra (iterable of floats): Array of pmRA
        pmdec (iterable of floats): Array of pmDEC
        color (iterable of trings): Array of colors

    Returns
    -------
        fig, ax
    """
    origin = np.array([ra, dec])
    fig, ax = plt.subplots()
    ax.quiver(*origin, pmra, pmdec, alpha=alpha)
    ax.scatter(ra, dec, color=color)
    ax.set_aspect('equal')
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    ax.ticklabel_format(useOffset=False)
    return fig, ax

def explore_pm(ra, dec, r, otype=None, pm_valid=True, alpha=0.5, mag_max=None, n_max=1000):
    """
    Explore a circular region showing objects and their proper motion

    Arguments
    ---------
        ra (float): RA of center of cirular region
        dec (float): DEC of center of cirular region
        r (float): radius of the cirular region in degrees
        otype (str): type of object ('star','galaxy',etc.). Default None of all objects
        pm_valid (bool): only return objects with valid pmRA and pmDEC
        alpha (float): alpha (between 0 and 1) of the flesh marks
        mag_max (float): maximum apparent value. Default None of all objects
        n_max (int): maximum number of rows to return; default 1000

    Returns
    -------
        df, fig, ax

    Ex:
    ---
        ra, dec = 266.41681662499997, -29.00782497222222 # Sgr A*
        df, fig, ax = explore_pm(ra, dec, r=0.001, otype='star')
    """
    if otype is not None:
        otype = str(tuple(object_type(otype)))
        add_otype = f" AND otype_txt in {otype}"
    else:
        add_otype = ""
    if pm_valid:
        add_pm = " AND pmra IS NOT NULL"
    else:
        add_pm = ""
    if mag_max is not None:
        add_mag = f" AND V <= {mag_max}"
    else:
        add_mag = ""
    script = f"""SELECT TOP {n_max} main_id, otypedef.otype_longname, allfluxes.V, ra, dec, pmra, pmdec
    FROM basic LEFT JOIN allfluxes ON basic.oid=allfluxes.oidref
    LEFT JOIN otypedef ON basic.otype=otypedef.otype
    WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {r})) = 1""" + add_otype + add_pm + add_mag
    df = sql2df(script)
    cols = ['V', 'ra', 'dec', 'pmra', 'pmdec']
    for i in cols:
        df.loc[df[i]=='', i] = np.nan
    df[cols] = df[cols].astype(float)
    df['ang_sep'] = [angular_separation(ra, dec, r, d) for r, d in zip(df['ra'], df['dec'])]
    df = df.sort_values('ang_sep')
    if len(df)>0:
        fig, ax = plot_pm(df['ra'], df['dec'], df['pmra'], df['pmdec'], alpha=alpha, color=None)
        ax.scatter([ra], [dec], marker='+', c='r')
    else:
        fig, ax = None, None
    return df, fig, ax
