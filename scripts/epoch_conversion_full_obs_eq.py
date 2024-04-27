
#from ephem import compat as ephem
from astropy.coordinates import FK5
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import numpy as np
import sys
from novas import compat as novas
from astropy.wcs import WCS
from astropy.io import fits
import glob


filename = sys.argv[1]
data_raw, header_raw = fits.getdata(filename,0,header=True)

date_obs = header_raw['DATE-OBS']


ast_ra,ast_dec,x,y=[],[],[],[]
with open("../extras/astrometric_xy.txt") as ip:
	crd=ip.read().splitlines(True)
for line in crd[2:]:
    row = line.split()
    ast_ra.append(float(row[3]))
    ast_dec.append(float(row[4]))
    x.append(float(row[1]))
    y.append(float(row[2]))


def precess_new(ra,dec,epoch_old,epoch_new):
	""" Precesses one epoch to another
	epoch input format : ( ra_in_degrees, dec_in_degrees, (date,month, year), (date,month, year))
	"""
	old_year = epoch_old[2];old_month = epoch_old[1]; old_date = epoch_old[0]
	a = np.floor((14-old_month)/12.0)
	y = old_year + 4800 - a
	m = old_month + 12*a - 3
	JD_old = old_date + np.floor((153*m + 2)/5) + 365*y + np.floor(y/4.0) - np.floor(y/100.0) + np.floor(y/400.0) - 32045

	new_year = epoch_new[2];new_month = epoch_new[1]; new_date = epoch_new[0]
	a = np.floor((14-new_month)/12.0)
	y = new_year + 4800 - a
	m = new_month + 12*a - 3
	JD_new = new_date + np.floor((153*m + 2)/5) + 365*y + np.floor(y/4.0) - np.floor(y/100.0) + np.floor(y/400.0) - 32045


	ra_n = [];dec_n = []

	for i in range(len(ra)):
		x = np.cos(np.deg2rad(ra[i]))*np.cos(np.deg2rad(dec[i]))
		y = np.sin(np.deg2rad(ra[i]))*np.cos(np.deg2rad(dec[i]))
		z = np.sin(np.deg2rad(dec[i]))
		x_n,y_n,z_n = novas.precession(JD_old,(x,y,z),JD_new)
		dec_n.append(np.rad2deg(np.arcsin(z_n)))
		r = np.rad2deg(np.arctan(y_n/x_n))
		if x_n < 0:
			ra_n.append(r + 180)
		elif y_n < 0 and x_n > 0:
			ra_n.append(r + 360)
		else:
			ra_n.append(r)
	return np.asarray(ra_n), np.asarray(dec_n)

#corrected_radec = precess_new(ast_ra,ast_dec,(24,10,2022),(1,1,2000))
corrected_radec = precess_new(ast_ra,ast_dec,(int(date_obs[8:10]),int(date_obs[5:7]),int(date_obs[0:4])),(1,1,2000))
#print(corrected_radec)


with open("../extras/astro_epoch_corrected.txt", "w+") as d:
	d.write('id\tx\ty\t\tAst_ra\t\tAst_dec\t\tCorrected_ra\tCorrected_dec')
	for x1 in range(len(x)):
		d.write("\n%.2d\t "%(x1+1))
		d.write("%4.2f\t "%(x[x1]))
		d.write("%4.2f\t "%(y[x1]))
		d.write("%.8f\t"%(ast_ra[x1]))
		d.write("%.8f\t "%(ast_dec[x1]))
		d.write("%.8f\t "%(corrected_radec[0][x1]))
		d.write("%.8f\t "%(corrected_radec[1][x1]))

#print('\n')


