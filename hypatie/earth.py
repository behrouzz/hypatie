import numpy as np

d2r = np.pi/180
r2d = 180/np.pi


re = 6378.1366 # equatorial radius
rp = 6356.7519 # polar radius
f = (re-rp)/re # flattening factor

e = (2*f - f**2) ** 0.5 # eccentricity
om = 7.292115e-5 # angular velocity (rad/s)
ge = 9.780327 # surface gravity at equator (m/s2)


def surface_gravity(lat):
    """Surface gravity of Earth at a given geodetic latitude at sea level"""
    ge = 9.780327 # surface gravity at equator (m/s2)
    g = ge * (1 + 0.0053024*np.sin(lat*d2r)**2 - 0.0000058*np.sin(2*lat*d2r)**2)
    return g


def radius(lat):
    """Radius of Earth at a given latitude"""
    re = 6378.1366 # equatorial radius
    rp = 6356.7519 # polar radius
    f = (re-rp)/re # flattening factor
    return re * (1 - f*(np.sin(lat*d2r)**2))


def geocentric_latitude(lat):
    """
    Calculate geocentric latitude for a given geodetic latitude

    Arguments
    ---------
        lat (deg): geodetic latitude
        
    Returns
    -------
        geocentric latitude

    Another formula: np.arctan((1-e**2)*np.tan(lat*d2r)) * r2d
    """
    tmp = (692.74 * np.sin(2*lat*d2r) - 1.16 * np.sin(4*lat*d2r))
    return lat - tmp/3600


def h_min_max(lat, dec):
    h_min = -90 + abs(lat+dec)
    h_max = 90 - abs(lat-dec)
    return np.array([h_min, h_max])


def geodetic_to_geocentric(obs_loc):
    """
    Convert geodetic coordinates (lon, lat, h) to geocentric coordinates (x, y, z)

    Arguments
    ---------
        obs_loc : Geodetic (geographic) coordinates as (lon (deg), lat (deg), alt (m))

    Returns
    -------
        xyz : geocentric coordinates as (x, y, z) in km

    Ref: Astrophysical formulae, 1999, vol.2, p.5
    """
    if len(obs_loc)==2:
        lon, lat = obs_loc
        h = 0
    else:
        lon, lat, h = obs_loc
        h = h/1000 # convert m to km
        
    re = 6378.1366 # equatorial radius
    rp = 6356.7519 # polar radius
    f = (re-rp)/re # flattening factor
    C = (np.cos(lat*d2r)**2 + (1-f)**2 * np.sin(lat*d2r)**2) ** (-0.5)
    S = (1-f)**2 * C
    x = (re*C + h) * np.cos(lat*d2r) * np.cos(lon*d2r)
    y = (re*C + h) * np.cos(lat*d2r) * np.sin(lon*d2r)
    z = (re*S + h) * np.sin(lat*d2r)
    return np.array([x,y,z])


def geocentric_to_geodetic(xyz):
    """
    Convert geocentric coordinates (x, y, z) to geodetic coordinates (lon, lat, h)

    Arguments
    ---------
        xyz : geocentric coordinates as (x, y, z) in km

    Returns
    -------
        obs_loc : geodetic (geographic) coordinates as (lon (deg), lat (deg), alt (m))

    Ref: Astrophysical formulae, 1999, vol.2, p.6-7
    """
    x, y, z = xyz
    re = 6378.1366 # equatorial radius
    rp = 6356.7519 # polar radius
    rp = np.where(z<0, -rp, rp).ravel()[0]
    r = (x**2 + y**2) ** 0.5
    E = (rp*z - (re**2 - rp**2)) / (re*r)
    F = (rp*z + (re**2 - rp**2)) / (re*r)
    P = (4/3)*(E*F+1)
    Q = 2*(E**2 - F**2)
    D = P**3 + Q**2
    if D < 0:
        v = 2*np.sqrt(-P)*np.cos((1/3)*np.arccos((Q/P)*(-P)**(-0.5)))
    else:
        v = ((D**0.5 - Q)**(1/3)) - ((D**0.5 + Q)**(1/3))
    if abs(z)<1: # or z==0:
        v = -(v**3 + 2*Q)/(3*P)
    G = 0.5 * ((E**2 + v)**0.5 + E)
    t = ((G**2 + (F-v*G)/(2*G-E)) ** 0.5) - G
    lat = np.arctan(re * (1-t**2)/(2*rp*t))
    h = (r-re*t)*np.cos(lat) + (z-rp)*np.sin(lat)
    lon = np.arctan(y/x)
    obs_loc = lon*r2d, lat*r2d, h*1000
    return obs_loc
