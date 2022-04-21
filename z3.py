from hypatie.iau import get_ha
from datetime import datetime
from hypatie.transform import radec_to_tete

t = datetime(2022, 4, 21, 21, 49, 55)
lon, lat = 7, 48

RA, DEC = 15, 35

ra, dec = radec_to_tete(RA, DEC, t)

ha, _ = get_ha(ra, dec, (lon, lat), t)


from astropy.coordinates import SkyCoord, EarthLocation, AltAz, HADec
from astropy.time import Time
from astropy import units as u

tt = Time(t)
c = SkyCoord(ra=ra, dec=dec, unit='deg')
#c = c.tete
me = EarthLocation(lon=lon, lat=lat)
hadec = HADec(location=me, obstime=tt)
cc = c.transform_to(hadec)

print('Behrouz Hour angle:', ha)
print('Astropy Hour angle:', cc.ha.deg)

from hypatie.transform import radec_to_altaz, hadec_to_altaz

az, alt = radec_to_altaz(ra, dec, (lon, lat), t)
print(az, alt)


aa = AltAz(location=me, obstime=tt)
ccc= c.transform_to(aa)
print(ccc.az.value, ccc.alt.value)
