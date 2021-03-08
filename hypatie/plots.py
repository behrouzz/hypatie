import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections.abc import Iterable
from .simbad import bright_objects
from .transform import radec_to_altaz

def _equalize_scale(X,Y,Z, ax):
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
    for xb, yb, zb in zip(Xb, Yb, Zb):
       ax.plot([xb], [yb], [zb], 'w')
    return ax

def plot_xyz(x, y, z, color, size, au=False):
    """cartesian plot"""
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
    """plot altitude/azimuth coordinates"""

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
        ax.set_yticklabels([])
        ax.tick_params(axis=u'both', which=u'both',length=0)
        #ax.set_yticklabels(ax.get_yticks()[::-1])
    if matplotlib.__version__ < '3.0.0':
        alt = [90-i for i in alt]
    ax.scatter(az, alt, c=color, s=size, alpha=alpha, marker=marker)
    ax.grid(True, alpha=0.7)
    return ax

def plot_radec(ra, dec, mag=None, size=None):
    """plot ra dec coordinates"""
        
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
            size = 0.1 + max(mag)-mag

    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111, projection="aitoff")
    ax.grid(True)
    for i in range(len(ra)):
        ax.plot(ra[i], dec[i], 'o', markersize=size[i], c='k')
    plt.subplots_adjust(top=0.95,bottom=0.0)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    return ax

def star_chart(lon, lat, t=None, otype=None):
    """plot the star chart for a location and time on earth"""
    _, rows = bright_objects(otype)
    mag = [i[2] for i in rows]
    ra  = [i[3] for i in rows]
    dec = [i[4] for i in rows]
    alt, az = radec_to_altaz(lon, lat, ra, dec, t)
    ax = plot_altaz(az, alt, mag=mag)
    return ax
