
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

from astropy import units as u
#from astroquery.gaia import Gaia
import warnings
warnings.filterwarnings('ignore')#,category=DeprecationWarning)
import glob


filename = sys.argv[1]
data, header = getdata(filename,0,header=True)

data = np.asarray(data,dtype=float)
#data = data.byteswap().newbyteorder()
x_matrix,y_matrix=[],[]
bkg=sep.Background(data)#,bw=1000,bh=1000,fw=3,fh=3)
bkg_image=bkg.back()
bkg_rms = bkg.rms()
back_sub_data = data - bkg_image

objects=sep.extract(back_sub_data,3,err=bkg.globalrms,maskthresh=20000,minarea=25)
'''try:
    objects=sep.extract(back_sub_data,15,err=bkg.globalrms,maskthresh=20000,minarea=25)
except:
    sep.set_extract_pixstack(120000)
    try:
        objects=sep.extract(back_sub_data,15,err=bkg.globalrms,maskthresh=20000,minarea=25)
    except:
        sep.set_extract_pixstack(180000)
        try:
            objects=sep.extract(back_sub_data,15,err=bkg.globalrms,maskthresh=20000,minarea=25)
        except:
            try:
                objects=sep.extract(back_sub_data,15,err=bkg.globalrms,maskthresh=20000,minarea=25)
            except:
                try:
                    objects=sep.extract(back_sub_data,15,err=bkg.globalrms,maskthresh=20000,minarea=25)
                except:
                    objects=sep.extract(back_sub_data,15,err=bkg.globalrms,maskthresh=20000,minarea=25)'''
            
x_matrix.append(objects['x'])
y_matrix.append(objects['y'])
xy = list(zip(x_matrix[0],y_matrix[0]))

#******************************************************


chunk1 = data[0:2048,1024:3072]
'''
UT_chunk_centre = start_UT_hr + offset + (exp_time_hr)*300/10000.0
LST_chunk_centre = UT2LST(obs_date,UT_chunk_centre)[0]
Dec_chunk_centre  = latitude
print('RA and DEC for chunk centre: ',LST_chunk_centre*15,Dec_chunk_centre)
'''
fits.writeto('../extras/chunk_image1.fits',chunk1,overwrite=True)

chunk2 = data[34816:,1024:3072]
fits.writeto('../extras/chunk_image2.fits',chunk2,overwrite=True)
