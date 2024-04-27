
import numpy as np
from lmfit import Model
from astropy.io.fits import getdata
import sep
from astropy.wcs import WCS
from astropy import wcs
from astropy.coordinates import SkyCoord

import sys
from astropy import units as u
from astroquery.gaia import Gaia
import warnings
warnings.filterwarnings('ignore')#,category=DeprecationWarning)




data, header = getdata('../extras/chunk_image1.fits',0,header=True)

data = np.asarray(data,dtype=float)
#data = data.byteswap().newbyteorder()
x_matrix,y_matrix=[],[]
bkg=sep.Background(data)#,bw=1000,bh=1000,fw=3,fh=3)
bkg_image=bkg.back()
bkg_rms = bkg.rms()
back_sub_data = data - bkg_image
objects=sep.extract(back_sub_data,3,err=bkg.globalrms,maskthresh=20000,minarea=5)
x_matrix.append(objects['x'])
y_matrix.append(objects['y'])
xy = list(zip(x_matrix[0],y_matrix[0]))

w = wcs.WCS('../extras/new_image_chunk1.fits',naxis=0)
pixcrd1 = np.array(xy,dtype=np.float64)
ra_dec_det = w.all_pix2world(pixcrd1,0)


from astropy.table import Table
gaia_dat = Table.read('/Users/rubiniucaa/Cat/GAIA_DR3_ILMT_full.fits', format='fits')
gaia_all_ra,gaia_all_dec,gaia_all_mag = gaia_dat['ra'].data, gaia_dat['dec'].data,gaia_dat['phot_g_mean_mag'].data
gaia_all_pmra, gaia_all_pmdec, gaia_all_parallax = gaia_dat['pmra'].data,gaia_dat['pmdec'].data,gaia_dat['parallax'].data
gaia_all_bp_rp = gaia_dat['bp_rp']


ra_min,ra_max = np.min(ra_dec_det[:,0]),np.max(ra_dec_det[:,0])
id_range = np.where( (gaia_all_ra > ra_min-0.1) & (gaia_all_ra < ra_max+0.1))[0]
gaia_all_ra,gaia_all_dec,gaia_all_mag = gaia_all_ra[id_range],gaia_all_dec[id_range],gaia_all_mag[id_range]
gaia_all_pmra, gaia_all_pmdec, gaia_all_parallax = gaia_all_pmra[id_range], gaia_all_pmdec[id_range], gaia_all_parallax[id_range]
gaia_all_bp_rp = gaia_all_bp_rp[id_range]

def sdss_query(ra,dec,rad = 2.0):
    tar_ra,tar_dec = ra,dec
    sep_ra = abs(gaia_all_ra - tar_ra ) * np.cos(np.radians((gaia_all_dec + tar_dec)/2))
    sep_dec  = abs(gaia_all_dec - tar_dec)
    sep = np.sqrt( (sep_ra)**2 + (sep_dec)**2)
    #print(np.min(sep))
    id_cntrprt = np.where(sep <= (rad/3600))
    gaia_ra,gaia_dec = gaia_all_ra[id_cntrprt],gaia_all_dec[id_cntrprt]
    #print(tar_ra,tar_dec,'\t',gaia_ra,gaia_dec)
    return [gaia_ra[0], gaia_dec[0]]

ra_dec_cat,ra_dec_ref,xy_ref =[],[],[]

#for i in range(100):
for i in range(len(xy)):
    try:
        ra_dec_cat.append(sdss_query(ra_dec_det[i][0],ra_dec_det[i][1]))
        ra_dec_ref.append(ra_dec_det[i])
        xy_ref.append(xy[i])
    except:
        pass


print(len(ra_dec_cat),' sources matched with SDSS ',' out of ',len(ra_dec_det))

with open("../extras/chunk_gaia_query.txt", "w+") as d:
	d.write('id\tx\ty\t\tSDSS_ra\t\tSDSS_dec\t\tAstro_ra\tAstro_dec\n')
	for x in range(len(ra_dec_cat)):
		d.write("\n%.2d\t "%(x+1))
		d.write("%4.2f\t "%(xy_ref[x][0]+1024))
		d.write("%4.2f\t "%(xy_ref[x][1]))
		d.write("%.9f\t"%(ra_dec_cat[x][0]))
		d.write("%.9f\t "%(ra_dec_cat[x][1]))
		d.write("%.9f\t "%(ra_dec_ref[x][0]))
		d.write("%.9f\t "%(ra_dec_ref[x][1]))
