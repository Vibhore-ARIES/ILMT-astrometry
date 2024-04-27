
import numpy as np
#from lmfit import Model
from astropy.io.fits import getdata
import sep
#from astropy.wcs import WCS
#from astropy import wcs
#from astropy.coordinates import SkyCoord

import sys
#from astropy import units as u
from astroquery.gaia import Gaia
import warnings
warnings.filterwarnings('ignore')#,category=DeprecationWarning)
#**************************************************************************

x,y,ast_ra,ast_dec,cor_ra,cor_dec=[],[],[],[],[],[]
z='../extras/astro_epoch_corrected.txt'
with open(z) as ip:
	f = ip.read().splitlines(True)
	for line in f[1:]:  
		row = line.split()
		x.append(float(row[1]))
		y.append(float(row[2]))
		ast_ra.append(float(row[3]))
		ast_dec.append(float(row[4]))
		cor_ra.append(float(row[5]))
		cor_dec.append(float(row[6]))

xy = list(zip(x,y))
ra_dec_det = list(zip(cor_ra,cor_dec))
ra_min, ra_max = np.min(cor_ra),np.max(cor_ra)

from astropy.table import Table
gaia_dat = Table.read('/Users/rubiniucaa/Cat/GAIA_DR3_ILMT_full.fits', format='fits')
gaia_all_ra,gaia_all_dec,gaia_all_mag = gaia_dat['ra'].data, gaia_dat['dec'].data,gaia_dat['phot_g_mean_mag'].data
gaia_all_pmra, gaia_all_pmdec, gaia_all_parallax = gaia_dat['pmra'].data,gaia_dat['pmdec'].data,gaia_dat['parallax'].data
gaia_all_bp_rp = gaia_dat['bp_rp']

id_range = np.where( (gaia_all_ra > ra_min-0.1) & (gaia_all_ra < ra_max+0.1))[0]
gaia_all_ra,gaia_all_dec,gaia_all_mag = gaia_all_ra[id_range],gaia_all_dec[id_range],gaia_all_mag[id_range]
gaia_all_pmra, gaia_all_pmdec, gaia_all_parallax = gaia_all_pmra[id_range], gaia_all_pmdec[id_range], gaia_all_parallax[id_range]
gaia_all_bp_rp = gaia_all_bp_rp[id_range]

def sdss_query(ra,dec,rad = 1.0):
	tar_ra,tar_dec = ra,dec
	sep_ra = abs(gaia_all_ra - tar_ra ) * np.cos(np.radians((gaia_all_dec + tar_dec)/2))
	sep_dec  = abs(gaia_all_dec - tar_dec)
	sep = np.sqrt( (sep_ra)**2 + (sep_dec)**2)
	#print(np.min(sep))
	id_cntrprt = np.where(sep <= (rad/3600))
	gaia_ra,gaia_dec,gaia_mag = gaia_all_ra[id_cntrprt],gaia_all_dec[id_cntrprt],gaia_all_mag[id_cntrprt]
	gaia_pmra,gaia_pmdec,gaia_parallax = gaia_all_pmra[id_cntrprt], gaia_all_pmdec[id_cntrprt], gaia_all_parallax[id_cntrprt]
	gaia_bp_rp = gaia_all_bp_rp[id_cntrprt]
	#print(tar_ra,tar_dec,'\t',gaia_ra,gaia_dec, len(gaia_ra),len(gaia_dec))
	n_comp = len(gaia_ra)
	return [gaia_ra[0], gaia_dec[0],gaia_mag[0],gaia_pmra[0],gaia_pmdec[0],gaia_parallax[0],n_comp,gaia_bp_rp[0]]

ra_dec_cat,ra_dec_ref,xy_ref =[],[],[]

#for i in range(20):
for i in range(len(xy)):
    try:
        ra_dec_cat.append(sdss_query(ra_dec_det[i][0],ra_dec_det[i][1]))
        ra_dec_ref.append(ra_dec_det[i])
        xy_ref.append(xy[i])
    except:
        pass


print(len(ra_dec_cat),' sources matched with GAIA ',' out of ',len(ra_dec_det))


with open("../extras/full_gaia_query_pm_1arcsec_no_companion.txt", "w+") as d:
	d.write('id\tx\ty\tGaia_ra\t\tGAIA_dec\t\tAstro_ra\tAstro_dec\tSDSS_gmag\tPM_RA\tPM_DEC\tParallax\tN_comp\tBP_RP\n')
	for x in range(len(ra_dec_cat)):
		d.write("\n%.2d\t "%(x+1))
		d.write("%4.2f\t "%(xy_ref[x][0]))
		d.write("%4.2f\t "%(xy_ref[x][1]))
		d.write("%.9f\t"%(ra_dec_cat[x][0]))
		d.write("%.9f\t "%(ra_dec_cat[x][1]))
		d.write("%.9f\t "%(ra_dec_ref[x][0]))
		d.write("%.9f\t "%(ra_dec_ref[x][1]))
		d.write("%.9f\t "%(ra_dec_cat[x][2]))
		d.write("%.9f\t "%(ra_dec_cat[x][3]))
		d.write("%.9f\t "%(ra_dec_cat[x][4]))
		d.write("%.9f\t "%(ra_dec_cat[x][5]))
		d.write("%d\t "%(ra_dec_cat[x][6]))
		d.write("%.5f\t "%(ra_dec_cat[x][7]))




#*****************************************************************************************
