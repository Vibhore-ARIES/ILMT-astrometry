
import numpy as np
from lmfit import Model
from astropy.io.fits import getdata
import sep
from astropy.wcs import WCS
from astropy import wcs
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from scipy.optimize import curve_fit
from scipy.stats import norm
from astropy.stats import sigma_clip
import matplotlib.patches as mpatches
import glob
from astropy import units as u
#from astroquery.gaia import Gaia
import warnings
warnings.filterwarnings('ignore')#,category=DeprecationWarning)
import glob
# import matplotlib.pyplot as plt
from astropy.io import ascii, fits
import sys
from astropy.wcs import WCS
from pathlib import Path
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

import time

start_time = time.time()
print('start_time: ',start_time)

#*******Read GAIA catalog in ILMT strip ######################
from astropy.table import Table
gaia_dat = Table.read('../multiprocessing_test/GAIA_DR3_ILMT_full.fits', format='fits')
#gaia_dat = Table.read('../GAIA_DR3_ILMT_full.fits', format='fits')

gaia_cat_ra,gaia_cat_dec,gaia_cat_mag = gaia_dat['ra'].data, gaia_dat['dec'].data,gaia_dat['phot_g_mean_mag'].data
gaia_cat_pmra, gaia_cat_pmdec, gaia_cat_parallax = gaia_dat['pmra'].data,gaia_dat['pmdec'].data,gaia_dat['parallax'].data
gaia_cat_bp_rp = gaia_dat['bp_rp']
#print('Length of Gaia stars in ILMT strip: ',len(gaia_cat_ra))
########################################################################################

def create_chunks(filename):
    data, header = getdata(filename,0,header=True)
    data = np.asarray(data,dtype=float)
    #******************************************************
    chunk1 = data[0:2048,1024:3072]
    fits.writeto('chunk1.fits',chunk1,overwrite=True)
    chunk2 = data[34816:,1024:3072]
    fits.writeto('chunk2.fits',chunk2,overwrite=True)
    return None
#********************************************************************************
######### Chunks generated, now perform astrometry on the two chunks ############
#********************************************************************************

from astroquery.astrometry_net import AstrometryNet

def do_astrometry(infile,outfile):
    ast = AstrometryNet()
    ast.api_key = 'ezduozpdlmrxxisu'
    wcs_header = ast.solve_from_image(infile, force_image_upload=True)
    #print(type(wcs_header))
    data = fits.getdata(infile)
    fits.writeto(outfile,data,header=wcs_header,overwrite=True)
    outfile = wcs_header#outfile
    return outfile

#************************************************************************************************
######### Chunks astrometry done, now gaia/sdss query of all stars in the two chunks ############
#************************************************************************************************


def gaia_query(ra,dec,rad = 1.0):

    gaia_all_ra,gaia_all_dec,gaia_all_mag = gaia_cat_ra,gaia_cat_dec,gaia_cat_mag
    gaia_all_pmra, gaia_all_pmdec, gaia_all_parallax = gaia_cat_pmra, gaia_cat_pmdec, gaia_cat_parallax
    gaia_all_bp_rp = gaia_cat_bp_rp

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


def get_gaia_coordinates(img):


    gaia_all_ra,gaia_all_dec,gaia_all_mag = gaia_cat_ra,gaia_cat_dec,gaia_cat_mag
    gaia_all_pmra, gaia_all_pmdec, gaia_all_parallax = gaia_cat_pmra, gaia_cat_pmdec, gaia_cat_parallax
    gaia_all_bp_rp = gaia_cat_bp_rp
    data, header = getdata(img,0,header=True)
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

    w = wcs.WCS('new_image_%s'%(img),naxis=0)
    pixcrd1 = np.array(xy,dtype=np.float64)
    ra_dec_det = w.all_pix2world(pixcrd1,0)

    ra_min,ra_max = np.min(ra_dec_det[:,0]),np.max(ra_dec_det[:,0])
    #print('Minimum and Max RA of chunk: ',ra_min,ra_max )
    id_range = np.where( (gaia_all_ra > ra_min-0.1) & (gaia_all_ra < ra_max+0.1))[0]
    gaia_all_ra,gaia_all_dec,gaia_all_mag = gaia_all_ra[id_range],gaia_all_dec[id_range],gaia_all_mag[id_range]
    gaia_all_pmra, gaia_all_pmdec, gaia_all_parallax = gaia_all_pmra[id_range], gaia_all_pmdec[id_range], gaia_all_parallax[id_range]
    gaia_all_bp_rp = gaia_all_bp_rp[id_range]
    print('Length of Gaia stars in the chunk : ',len(gaia_all_ra))
    ra_dec_cat,ra_dec_ref,xy_ref =[],[],[]

    for i in range(len(xy)):
        try:
            ra_dec_cat.append(gaia_query(ra_dec_det[i][0],ra_dec_det[i][1]))
            ra_dec_ref.append(ra_dec_det[i])
            xy_ref.append(xy[i])
        except:
            pass

    print(len(ra_dec_cat),' sources matched with GAIA',' out of ',len(ra_dec_det))

    if(img=='chunk1.fits'):
        with open("chunk_gaia_query.txt", "w+") as d:
            d.write('id\tx\ty\t\tSDSS_ra\t\tSDSS_dec\t\tAstro_ra\tAstro_dec\n')
            for x in range(len(ra_dec_cat)):
                d.write("\n%.2d\t "%(x+1))
                d.write("%4.2f\t "%(xy_ref[x][0]+1024))
                d.write("%4.2f\t "%(xy_ref[x][1]))
                d.write("%.9f\t"%(ra_dec_cat[x][0]))
                d.write("%.9f\t "%(ra_dec_cat[x][1]))
                d.write("%.9f\t "%(ra_dec_ref[x][0]))
                d.write("%.9f\t "%(ra_dec_ref[x][1]))

    elif(img=='chunk2.fits'):
        with open("chunk_gaia_query.txt", "a") as d:
            for x in range(len(ra_dec_cat)):
                d.write("\n%.2d\t "%(x+1))
                d.write("%4.2f\t "%(xy_ref[x][0]+1024))
                d.write("%4.2f\t "%(xy_ref[x][1]+34816))
                d.write("%.9f\t"%(ra_dec_cat[x][0]))
                d.write("%.9f\t "%(ra_dec_cat[x][1]))
                d.write("%.9f\t "%(ra_dec_ref[x][0]))
                d.write("%.9f\t "%(ra_dec_ref[x][1]))

    return None


#********************************************************************************
######### Chunks gaia/sdss query done, now epoch conversion ############
#********************************************************************************

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

from novas import compat as novas
from astropy.io import fits
import glob

def get_chunk_epoch_corrected():
    data_raw, header_raw = fits.getdata(filename,0,header=True)
    date_obs = header_raw['DATE-OBS']

    gaia_ra,gaia_dec,x,y=[],[],[],[]
    with open("chunk_gaia_query.txt") as ip:
        crd=ip.read().splitlines(True)
    for line in crd[2:]:
        row = line.split()
        gaia_ra.append(float(row[3]))
        gaia_dec.append(float(row[4]))
        x.append(float(row[1]))
        y.append(float(row[2]))

    corrected_radec = precess_new(gaia_ra,gaia_dec,(1,1,2000),(int(date_obs[8:10]),int(date_obs[5:7]),int(date_obs[0:4])))

    with open("list_selected_epoch_corrected.txt", "w+") as d:
        d.write('id\tx\ty\t\tGaia_ra\t\tGaia_dec\t\tCorrected_ra\tCorrected_dec')
        for x1 in range(len(x)):
            d.write("\n%.2d\t "%(x1+1))
            d.write("%4.2f\t "%(x[x1]))
            d.write("%4.2f\t "%(y[x1]))
            d.write("%.8f\t"%(gaia_ra[x1]))
            d.write("%.8f\t "%(gaia_dec[x1]))
            d.write("%.8f\t "%(corrected_radec[0][x1]))
            d.write("%.8f\t "%(corrected_radec[1][x1]))


#********************************************************************************
######### Chunks epoch conversion done, now chunk transfomrations ############
#********************************************************************************

def chunk_transformations():
    x,y,gaia_ra,gaia_dec,cor_ra,cor_dec=[],[],[],[],[],[]
    z='list_selected_epoch_corrected.txt'
    with open(z) as ip:
        f = ip.read().splitlines(True)
        for line in f[1:]:  
            row = line.split()
            x.append(float(row[1]))
            y.append(float(row[2]))
            gaia_ra.append(float(row[3]))
            gaia_dec.append(float(row[4]))
            cor_ra.append(float(row[5]))
            cor_dec.append(float(row[6]))

    #**********************************************************************************************************************************
    x1,y1 = np.array(x),np.array(y)
    A, D= cor_ra[0], cor_dec[0]
    x0,y0 = x[0], y[0]

    def func(X,g1,f1,g2,f2,g3,f3,g4,f4,g5,f5):
        x = X[0]
        y = X[1]
        ra  = f1 + (x * f2) + (y * f3)+ (x*x*f4) + (y*y*f5)
        dec = g1 + (y * g2) + (x * g3)+ (y*y*g4) + (x*x*g5)
        return np.concatenate((ra,dec))



    fmodel=Model(func)
    params= fmodel.make_params(g1=10,f1=10,g2=10,f2=10,g3=10,f3=10,g4=10,f4=10,g5=10,f5=10)
    result=fmodel.fit(np.concatenate((cor_ra,cor_dec)),params,X=(x1,y1))

    dict = result.best_values
    values = [float(x11) for x11 in list(dict.values())]
    popt=[]
    for value in values:
        popt.append(value)

    p5 = func((x1,y1), popt[0], popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7],popt[8],popt[9])
    #p5 = func((x_matrix[0],y_matrix[0]), popt[0], popt[1])
    hf = int(len(p5)/2)


    #*********************************************************************
    data, header = getdata(filename,0,header=True)

    data = np.asarray(data,dtype=float)
    #data = data.byteswap().newbyteorder()
    x_matrix,y_matrix=[],[]
    bkg=sep.Background(data)#,bw=1000,bh=1000,fw=3,fh=3)
    bkg_image=bkg.back()
    bkg_rms = bkg.rms()
    back_sub_data = data - bkg_image

    try:
        objects=sep.extract(back_sub_data,3,err=bkg.globalrms,maskthresh=20000,minarea=5)
    except Exception as e:
        buff_message = 'internal pixel buffer full'
        if str(e)[0:26] == buff_message:
            sep.set_extract_pixstack(600000)
        try:
            objects=sep.extract(back_sub_data,3,err=bkg.globalrms,maskthresh=20000,minarea=5)
        except Exception as e:
            if str(e)[0:26] == buff_message:
                sep.set_extract_pixstack(900000)
                objects=sep.extract(back_sub_data,3,err=bkg.globalrms,maskthresh=20000,minarea=5)

    x_matrix.append(objects['x'])
    y_matrix.append(objects['y'])
    x_all = x_matrix[0]
    y_all = y_matrix[0]
    xy = list(zip(x_all,y_all))
    xy = np.array(xy)

    radec_all = func((x_all,y_all), popt[0], popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7],popt[8],popt[9])
    hf = int(len(radec_all)/2)
    all_ra=radec_all[0:hf]
    all_dec=radec_all[hf:]


    with open("astrometric_xy.txt", "w+") as d:
        d.write('id\tx\t\y\tAstro_ra\tAstro_dec\n')
        for x in range(len(all_ra)):
            d.write("\n%.2d\t "%(x+1))
            d.write("%4.2f\t "%(xy[x][0]))
            d.write("%4.2f\t "%(xy[x][1]))
            d.write("%.9f\t"%(all_ra[x]))
            d.write("%.9f\t "%(all_dec[x]))

    return None

#********************************************************************************
######### Chunk transfomrations done and fitted, now convert back to J2000 epoch ############
#********************************************************************************

def get_full_J2000_coordinates():
    ast_ra,ast_dec,x,y=[],[],[],[]
    with open("astrometric_xy.txt") as ip:
        crd=ip.read().splitlines(True)
    for line in crd[2:]:
        row = line.split()
        ast_ra.append(float(row[3]))
        ast_dec.append(float(row[4]))
        x.append(float(row[1]))
        y.append(float(row[2]))

    data_raw, header_raw = fits.getdata(filename,0,header=True)
    date_obs = header_raw['DATE-OBS']
    #corrected_radec = precess_new(ast_ra,ast_dec,(24,10,2022),(1,1,2000))
    corrected_radec = precess_new(ast_ra,ast_dec,(int(date_obs[8:10]),int(date_obs[5:7]),int(date_obs[0:4])),(1,1,2000))
    #print(corrected_radec)


    with open("astro_epoch_corrected.txt", "w+") as d:
        d.write('id\tx\ty\t\tAst_ra\t\tAst_dec\t\tCorrected_ra\tCorrected_dec')
        for x1 in range(len(x)):
            d.write("\n%.2d\t "%(x1+1))
            d.write("%4.2f\t "%(x[x1]))
            d.write("%4.2f\t "%(y[x1]))
            d.write("%.8f\t"%(ast_ra[x1]))
            d.write("%.8f\t "%(ast_dec[x1]))
            d.write("%.8f\t "%(corrected_radec[0][x1]))
            d.write("%.8f\t "%(corrected_radec[1][x1]))


#********************************************************************************************************************
######### Full frame coordinates converted back to J2000 epoch, now gaia/sdss query of full frame ###############
#********************************************************************************************************************
def full_gaia_query():

    gaia_all_ra,gaia_all_dec,gaia_all_mag = gaia_cat_ra,gaia_cat_dec,gaia_cat_mag
    gaia_all_pmra, gaia_all_pmdec, gaia_all_parallax = gaia_cat_pmra, gaia_cat_pmdec, gaia_cat_parallax
    gaia_all_bp_rp = gaia_cat_bp_rp

    x,y,ast_ra,ast_dec,cor_ra,cor_dec=[],[],[],[],[],[]
    z='astro_epoch_corrected.txt'
    with open(z) as ip:
        f = ip.read().splitlines(True)
        for line in f[1:]:  
            row = line.split()
            x.append(float(row[1]))
            y.append(float(row[2]))
            #ast_ra.append(float(row[3]))
            #ast_dec.append(float(row[4]))
            cor_ra.append(float(row[5]))
            cor_dec.append(float(row[6]))

    xy = list(zip(x,y))
    ra_dec_det = list(zip(cor_ra,cor_dec))
    ra_min, ra_max = np.min(cor_ra),np.max(cor_ra)

    id_range = np.where( (gaia_all_ra > ra_min-0.1) & (gaia_all_ra < ra_max+0.1))[0]
    gaia_all_ra,gaia_all_dec,gaia_all_mag = gaia_all_ra[id_range],gaia_all_dec[id_range],gaia_all_mag[id_range]
    gaia_all_pmra, gaia_all_pmdec, gaia_all_parallax = gaia_all_pmra[id_range], gaia_all_pmdec[id_range], gaia_all_parallax[id_range]
    gaia_all_bp_rp = gaia_all_bp_rp[id_range]


    ra_dec_cat,ra_dec_ref,xy_ref =[],[],[]

    #for i in range(20):
    for i in range(len(xy)):
        try:
            ra_dec_cat.append(gaia_query(ra_dec_det[i][0],ra_dec_det[i][1]))
            ra_dec_ref.append(ra_dec_det[i])
            xy_ref.append(xy[i])
        except:
            pass


    print(len(ra_dec_cat),' sources matched with GAIA in full frame (trial 1)',' out of ',len(ra_dec_det))


    with open("full_gaia_query_pm_1arcsec_no_companion.txt", "w+") as d:
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


#********************************************************************************************************************
######### Full gaia/sdss query done, again epoch conversion for fitting ###############
#********************************************************************************************************************

def get_epoch_corrected_full():
    data_raw, header_raw = fits.getdata(filename,0,header=True)
    date_obs = header_raw['DATE-OBS']

    ast_ra,ast_dec,x,y,g_mag,pm_ra,pm_dec,parallax,N_comp=[],[],[],[],[],[],[],[],[]
    with open("full_gaia_query_pm_1arcsec_no_companion.txt") as ip:
        crd=ip.read().splitlines(True)
    for line in crd[2:]:
        row = line.split()
        ast_ra.append(float(row[3]))
        ast_dec.append(float(row[4]))
        x.append(float(row[1]))
        y.append(float(row[2]))
        pm_ra.append(float(row[8]))
        pm_dec.append(float(row[9]))
        parallax.append(float(row[10]))
        N_comp.append(int(row[11]))


    corrected_radec = precess_new(ast_ra,ast_dec,(1,1,2000),(int(date_obs[8:10]),int(date_obs[5:7]),int(date_obs[0:4])))

    with open("gaia_epoch_corrected_final_1arcsec_no_companion.txt", "w+") as d:
        d.write('id\tx\ty\t\tAst_ra\t\tAst_dec\t\tCorrected_ra\tCorrected_dec')
        for x1 in range(len(x)):
            d.write("\n%.2d\t "%(x1+1))
            d.write("%4.2f\t "%(x[x1]))
            d.write("%4.2f\t "%(y[x1]))
            d.write("%.8f\t"%(ast_ra[x1]))
            d.write("%.8f\t "%(ast_dec[x1]))
            d.write("%.8f\t "%(corrected_radec[0][x1]))
            d.write("%.8f\t "%(corrected_radec[1][x1]))
            d.write("%d\t "%(N_comp[x1]))
    return None




#********************************************************************************************************************
######### Full frame epoch conversion done, now fit transformations ###############
#********************************************************************************************************************

def plot_astrometry(xpix,ypix,RA_new,DEC_new,RA_g,DEC_g,save_file_name):
    fnt_sz,leg_sz=12,10
    #######
    plt.figure(1, figsize=(14,8))
    plt.subplots_adjust(wspace=0.35,hspace=0.45)
    #plt.text(6000, 2.8, 'q1 = (RA - RA0) - ( (y-y0)*f1 + f2)', ha='center',fontsize=12)
    #########
    DEC_diff_arcsec = (DEC_new - DEC_g)*3600.0
    RA_diff_arcsec = (RA_new - RA_g)*3600.0
    RA_diff_arcsec = sigma_clip(RA_diff_arcsec, sigma=3, maxiters=None,masked=False, copy=False)
    DEC_diff_arcsec = sigma_clip(DEC_diff_arcsec, sigma=3, maxiters=None, masked=False, copy=False)
    ####################
    ax1 = plt.subplot(2, 2, 1)
    cpdf_1 = (RA_new - RA_g)*3600.0
    #print(cpdf_1)
    plt.xlabel('Along RA',fontsize=fnt_sz)
    plt.ylabel( 'RA gaia - RA cal (arcsec)' ,fontsize=fnt_sz)
    plt.scatter(ypix, cpdf_1,c=xpix,marker='.',s=2)
    plt.colorbar(label="Dec pixel", orientation="vertical")
    ##############
    ax2 = plt.subplot(2, 2, 2)

    nbins = 75
    n, bins, patches = ax2.hist(RA_diff_arcsec,nbins, density=True, facecolor = 'grey', alpha = 0.5); 

    centers = (0.5*(bins[1:]+bins[:-1]))
    pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), centers, n, p0=[0,1])

    ax2.plot(centers, norm.pdf(centers,*pars), 'k--',linewidth = 2) 




    #sns.distplot(RA_diff_arcsec ,hist=True  , hist_kws={'color': 'r'}, fit=stats.norm, kde=False)
    #mu, sigma = stats.norm.fit(RA_diff_arcsec)    
    plt.xlabel('RA gaia - RA cal (arcsec)',fontsize=fnt_sz); plt.ylabel( 'PDF' ,fontsize=fnt_sz)
    handles, labels = plt.gca().get_legend_handles_labels()
    empty_patch = mpatches.Patch(color='none'); handles.append(empty_patch)  
    #emp_label1 = u'$\mu$='+ str('%.2E' % mu); emp_label2 = u'$\sigma$='+ str('%.2E' % sigma)
    emp_label1 = u'$\mu$='+ str('%.2E' % pars[0]); emp_label2 = u'$\sigma$='+ str('%.2E' % pars[1])
    labels.append(emp_label1+'\n'+emp_label2); plt.legend(handles, labels,frameon=False,loc='upper left',prop={'size': leg_sz})
    #ax2.set_title('$\mu={:.4f}\pm{:.4f}$, $\sigma={:.4f}\pm{:.4f}$'.format(pars[0],np.sqrt(cov[0,0]), pars[1], np.sqrt(cov[1,1 ])))
    ###################

    ax3 = plt.subplot(2, 2, 3)
    cpdf_1 = (DEC_new - DEC_g)*3600.0
    plt.xlabel('Along DEC',fontsize=fnt_sz)
    plt.ylabel( 'DEC gaia - DEC cal (arcsec)' ,fontsize=fnt_sz)
    #plt.scatter(xpix, cpdf_1,c=ypix,marker='.',s=2)
    plt.scatter(xpix, cpdf_1,c=ypix,marker='.',s=2)
    plt.colorbar(label="RA pixel", orientation="vertical")
    ##############
    ax4 = plt.subplot(2, 2, 4)

    n, bins, patches = ax4.hist(DEC_diff_arcsec,nbins, density=True, facecolor = 'grey', alpha = 0.5); 

    centers = (0.5*(bins[1:]+bins[:-1]))
    pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), centers, n, p0=[0,1])

    ax4.plot(centers, norm.pdf(centers,*pars), 'k--',linewidth = 2) 


    #sns.distplot(DEC_diff_arcsec ,hist=True  , hist_kws={'color': 'r'}, fit=stats.norm, kde=False)
    #mu, sigma = stats.norm.fit(DEC_diff_arcsec)    
    plt.xlabel('DEC gaia - DEC cal (arcsec)',fontsize=fnt_sz); plt.ylabel( 'PDF' ,fontsize=fnt_sz)
    handles, labels = plt.gca().get_legend_handles_labels()
    empty_patch = mpatches.Patch(color='none'); handles.append(empty_patch)  
    #emp_label1 = u'$\mu$='+ str('%.2E' % mu); emp_label2 = u'$\sigma$='+ str('%.2E' % sigma)
    emp_label1 = u'$\mu$='+ str('%.2E' % pars[0]); emp_label2 = u'$\sigma$='+ str('%.2E' % pars[1])
    labels.append(emp_label1+'\n'+emp_label2); plt.legend(handles, labels,frameon=False,loc='upper left',prop={'size': leg_sz})    
    ############# 
    plt.suptitle('Astrometric offsets in one of the TDI frames observed on 2022-10-24 ')
    plt.savefig(save_file_name)
    #plt.show()
    plt.clf()
    plt.cla()
    plt.close()

def final_transformations():
    x,y,gaia_ra,gaia_dec,cor_ra,cor_dec,g_mag=[],[],[],[],[],[],[]
    pm_ra,pm_dec,parallax,N_comp=[],[],[],[]
    z='gaia_epoch_corrected_final_1arcsec_no_companion.txt'
    with open(z) as ip:
        f = ip.read().splitlines(True)
        for line in f[1:]:  
            row = line.split()
            x.append(float(row[1]))
            y.append(float(row[2]))
            gaia_ra.append(float(row[3]))
            gaia_dec.append(float(row[4]))
            cor_ra.append(float(row[5]))
            cor_dec.append(float(row[6]))
            N_comp.append(int(row[7]))


    id_nans = np.argwhere(~np.isnan(cor_ra))
    x,y,gaia_ra,gaia_dec,cor_ra,cor_dec,N_comp = np.array(x)[id_nans],np.array(y)[id_nans],np.array(gaia_ra)[id_nans],np.array(gaia_dec)[id_nans],np.array(cor_ra)[id_nans],np.array(cor_dec)[id_nans],np.array(N_comp)[id_nans]

    id_isolated = np.where(np.array(N_comp)==1)[0]
    x,y,gaia_ra,gaia_dec,cor_ra,cor_dec,N_comp = np.array(x)[id_isolated],np.array(y)[id_isolated],np.array(gaia_ra)[id_isolated],np.array(gaia_dec)[id_isolated],np.array(cor_ra)[id_isolated],np.array(cor_dec)[id_isolated],np.array(N_comp)[id_isolated]


    #**********************************************************************************************************************************
    x1,y1 = np.array(x),np.array(y)

    y0,x0 = 18432.0, 2048.0

    def func(X,g1,f1,g2,f2,g3,f3,g5):
        x = X[0]
        y = X[1]

        ra  = f1 + ((x-x0) * f2) + ((y-y0) * f3)
        dec = g1 + ((y-y0) * g2) + ((x-x0) * g3) + ((x-x0)*(x-x0)*g5)

        return np.concatenate((ra,dec))

    fmodel=Model(func)
    params= fmodel.make_params(g1=10,f1=10,g2=10,f2=10,g3=10,f3=10,g5=10)
    result=fmodel.fit(np.concatenate((cor_ra,cor_dec)),params,X=(x1,y1))

    dict = result.best_values
    values = [float(x) for x in list(dict.values())]
    popt=[]
    for value in values:
        popt.append(value)


    print(popt)
    p5 = func((x1,y1), popt[0], popt[1],popt[2],popt[3],popt[4],popt[5],popt[6])
    #p5 = func((x_matrix[0],y_matrix[0]), popt[0], popt[1])
    hf = int(len(p5)/2)
    ras=p5[0:hf]
    decs=p5[hf:]

    RA_deg_all = ras
    DEC_deg_all = decs



    plot_astrometry(x,y,RA_deg_all,DEC_deg_all,cor_ra,cor_dec,'new_plot_4d_1arcsec_no_companion.pdf')

    fitsname = filename
    data, header = getdata(fitsname,0,header=True)
    header.append('f1')
    header['f1']= popt[1]
    header.append('f2')
    header['f2']= popt[3]
    header.append('f3')
    header['f3']= popt[5]
    header.append('g1')
    header['g1']= popt[0]
    header.append('g2')
    header['g2']= popt[2]
    header.append('g3')
    header['g3']= popt[4]
    header.append('g5')
    header['g5']= popt[6]
    header.append('relation')
    header['relation']= 'RA = f1 + ((x-x0) * f2) + ((y-y0) * f3)  and DEC = g1 + ((y-y0) * g2) + ((x-x0) * g3) + ((x-x0)*(x-x0)*g5)'
    fits.writeto('%s_ast.fits'%(fitsname[:-5]),data,header,overwrite=True)
    return None


#********************************************************************************************************************
######### Final transformations done, now header update###############
#********************************************************************************************************************

'''
Functions that are used later on
'''
def pixel_to_sky(x, y, f1, f2, f3, f4, f5, g1, g2, g3, g4, g5):
    x0 = 2048
    y0 = 36864/2
    X = x - x0
    Y = y - y0
    ra = f1 + f2*X + f3*Y + f4*X**2 + f5*Y**2
    dec = g1 + g2*X + g3*Y + g4*X**2 + g5*Y**2
    # print(ra, f1, dec, g1)
    return ra, dec

def removewcs(parent_header):
    '''
    Adapted from https://github.com/Astro-Sean/autophot 
    '''

    if 'UPWCS' in parent_header.keys():
        del parent_header['UPWCS']

    keywords = ['CD1_1','CD1_2', 'CD2_1','CD2_2', 'CRVAL1','CRVAL2', 'CRPIX1','CRPIX2',
                'CUNIT1','CUNIT2', 'CTYPE1','CTYPE2', 'WCSAXES','EQUINOX', 'LONPOLE','LATPOLE',
                'CDELT1','CDELT2', 'A_ORDER', 'A_0_0', 'A_0_1','A_0_2', 'A_1_0','A_1_1',
                'A_2_0', 'B_ORDER', 'B_0_0','B_0_1', 'B_0_2','B_1_0', 'B_1_1','B_2_0',
                'AP_ORDER', 'AP_0_0','AP_0_1', 'AP_0_2','AP_1_0', 'AP_1_1','AP_2_0',
                'BP_ORDER', 'BP_0_0','BP_0_1', 'BP_0_2','BP_1_0', 'BP_1_1','BP_2_0',
                'PROJP1','PROJP3', 'RADECSYS', 'PV1_1','PV1_2', 'PV2_1','PV2_2', 'LTV1','LTV2',
                'LTM1_1','LTM2_2', 'PC1_1','PC1_2', 'PC2_1','PC2_2', 'RADESYS']

    for i in keywords:
        if i in parent_header.keys():
            del parent_header[i]

    return parent_header
    
def generate_wcs_header(header,equinox):
    '''
    Updating header with WCS information
    '''
    
    f1, f2, f3, f4, f5 = header['F1'], header['F2'], header['F3'], 0, 0
    g1, g2, g3, g4, g5 = header['G1'], header['G3'], header['G2'], header['G5'], 0
    if f1 > 180:
        f1 -= 360

    x0 = 2048
    y0 = 36864/2

    M = np.array([
        [f2, f3],
        [g2, g3]
    ])
      
    old_header, header = header, removewcs(header)

    M_inv = np.linalg.inv(M)

    A1, B1 = ( M_inv @ np.array([[f4], [g4]]) ).reshape(2)
    A2, B2 = ( M_inv @ np.array([[f5], [g5]]) ).reshape(2)


    header.set('CTYPE1', 'RA---CAR-SIP')
    header.set('CTYPE2', 'DEC--CAR-SIP')

    header.set('CUNIT1', 'deg')
    header.set('CRVAL1', f1)
    header.set('CRPIX1', x0+1)

    header.set('CUNIT2', 'deg')
    header.set('CRVAL2', g1)
    header.set('CRPIX2', y0+1)

    header.set('CD1_1' , f2)
    header.set('CD1_2' , f3)
    header.set('CD2_1' , g2)
    header.set('CD2_2' , g3)
    header.set('EQUINOX', equinox)
    header.set('RADESYS', 'FK5')

    header.set('PV1_0', 1)
    header.set('PV1_1', f1)
    header.set('PV1_2', g1)

    header.set('A_ORDER', 2)
    header.set('A_2_0', A1)
    header.set('A_0_2', A2)
    header.set('B_ORDER', 2)
    header.set('B_2_0', B1)
    header.set('B_0_2', B2)
    
    header.set('LONPOLE', 0)
    header.set('LATPOLE', 90)
    
    return header


def update_wcs_header():
    foldername = Path('.')
    #foldername = Path('ast_new')
    flist = list(foldername.glob('*ast.fits'))
    nfiles = len(flist)

    savefolder_ast = foldername.parent / Path('.')
    savefolder_wcs = foldername.parent / Path('.')#Path('/home/ilmt/wcs_rename')


    orig_name_list = []
    save_name_list = []
    for i, f in enumerate(flist, 1):
        save_file = True
        print('\n')
        print(f'File {i} of {nfiles}')
        print(f'Filename: {f.name}') 
        
        img, header = fits.getdata(f, header=True)
            
            
        '''
        Reading Information from Header
        '''

        # print(header['RELATION'])

        from astropy.time import Time
        t = '{}T{}'.format(header['DATE-OBS'], header['UTSTART'])
        t = Time(t)
        equinox = t.jyear
        print('The utc of frame is: ',t.utc ,equinox)  
        
            
        header = generate_wcs_header(header,equinox)

        '''
        Testing the WCS information against the fitted relation
        '''
        f1, f2, f3, f4, f5 = header['F1'], header['F2'], header['F3'], 0, 0
        g1, g2, g3, g4, g5 = header['G1'], header['G3'], header['G2'], header['G5'], 0
        x0 = 2048
        y0 = 36864/2
        
        w = WCS(header)
        if not w.wcs.get_pv():
            print(f'Your astropy version still has a bug in WCS module.')
            w.wcs.set_pv([(1, 0, 1.0), (1, 1, f1), (1, 2, g1), (2, 3, 0), (2, 4, 90)])
            w.wcs.set()

        from astropy.coordinates import FK5

        test_pixels = [(2048, 18432), (0, 0), (4096, 0), (4096, 36800), (0, 36800)]
        ra1, dec1 = [], []

        for x, y in test_pixels:
            c = w.pixel_to_world(x, y)
            try:
                c.transform_to('icrs')
            except:
                print('nan values found.')
                save_file = False
            ra1.append(c.ra.deg)
            dec1.append(c.dec.deg)
            
        # ### Vibhore's relation
        ra2, dec2 = [], []
        for x, y in test_pixels:
            c = pixel_to_sky(x, y, f1, f2, f3, f4, f5, g1, g2, g3, g4, g5)
            ra2.append(c[0])
            dec2.append(c[1])


        print(f'Maximum offset in ra is : {max([abs(a - b)*3600 for a, b in zip(ra1, ra2)]):E} arcsec')
        print(f'Maximum offset in dec is : {max([abs(a - b)*3600 for a, b in zip(dec1, dec2)]):E} arcsec')
        
        
        '''
        Saving the new fits file
        '''
        if max( [abs(a - b)*3600 for a, b in zip(ra1, ra2)] + [abs(a - b)*3600 for a, b in zip(dec1, dec2)] ) > 1e-9:
            save_file = False
            print(f'\033[91m Deviation from the fitting relation is > 1e-9. Not saving the file. \033[00m')

        print(np.array(ra1))
        
        if len(sys.argv) <= 1:
            pass
        elif sys.argv[1] == 'nosave':
            save_file = False
        
        if save_file:    
            new_hdu = fits.PrimaryHDU(img, header)
            if not savefolder_ast.is_dir():
                savefolder_ast.mkdir()
            if not savefolder_wcs.is_dir():
                savefolder_wcs.mkdir()
            
            print(savefolder_ast / f.name)
            f.rename(savefolder_ast / f.name)
            
            flt = header['FILTER']
            start_coord = w.pixel_to_world(2048, 0)
            c = start_coord.transform_to('icrs')
            start_lst = c.ra.to_string(unit='hourangle', fields=2)
            save_name = savefolder_wcs / f'{f.stem[:-9]}_{flt}_{start_lst}{f.suffix}'
            print(save_name)
            orig_name_list.append(f.stem)
            save_name_list.append(save_name.stem)
            new_hdu.writeto(save_name, overwrite=True)

    print(orig_name_list)
    print(save_name_list)
    namelist = Table([orig_name_list, save_name_list])
    namelist.write('namelist.dat', format='ascii', overwrite=True)
    return None





def main(filename):
    create_chunks(filename)
    do_astrometry('chunk1.fits','new_image_chunk1.fits')
    do_astrometry('chunk2.fits','new_image_chunk2.fits')
    get_gaia_coordinates('chunk1.fits')
    get_gaia_coordinates('chunk2.fits')
    get_chunk_epoch_corrected()
    chunk_transformations()
    get_full_J2000_coordinates()
    full_gaia_query()
    get_epoch_corrected_full()
    final_transformations()
    update_wcs_header()


if __name__ == '__main__':
    files = glob.glob('*37.fits')
    for filename in files:
        main(filename)


end_time = time.time()
elapsed_time = end_time - start_time

print('Elapsed time ',elapsed_time)
