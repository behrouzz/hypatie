from hypatie.iau import get_ha
from datetime import datetime
from hypatie.transform import radec_to_tete

t = datetime.utcnow()
lon, lat = 7, 48
RA, DEC = 15, 35

ra, dec = radec_to_tete(RA, DEC, t)

ha = get_ha(ra, dec, (lon, lat), t)


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
