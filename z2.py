from hypatie.transform import to_tete
import numpy as np
from datetime import datetime

t = datetime(2022, 3, 18)

# GCRS coordinates
pos = np.array([0.73859258, 0.13935437, 0.65959182])

# True Equator and True equinox of t
pos_tete = to_tete(pos, t)

print(pos_tete)
#[0.73649269 0.14295327 0.66116782]




from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy import units as u
from hypatie.transform import sph2car, car2sph

ra, dec, r = car2sph(pos)

tt = Time(t)
c = SkyCoord(ra=ra, dec=dec, unit='deg')
cc = c.tete
pos_tete2 = sph2car(np.array([cc.ra.value, cc.dec.value, r]))
print(pos_tete2)
