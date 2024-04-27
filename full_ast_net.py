
import numpy as np
from astropy.io import fits
import sep
import warnings
warnings.filterwarnings('ignore')#,category=DeprecationWarning)
from astroquery.astrometry_net import AstrometryNet
import time
start = time.time()

img = '20230522-0003.fits'

def do_astrometry(infile,outfile):
    ast = AstrometryNet()
    ast.api_key = 'ezduozpdlmrxxisu'
    data, header = fits.getdata(infile,0,header=True)
    data = np.asarray(data,dtype=float)
    x_matrix,y_matrix=[],[]
    bkg=sep.Background(data)#,bw=1000,bh=1000,fw=3,fh=3)
    bkg_image=bkg.back()
    bkg_rms = bkg.rms()
    back_sub_data = data - bkg_image
    objects=sep.extract(back_sub_data,3,err=bkg.globalrms,maskthresh=20000,minarea=5)
    objects = np.sort(objects, order='flux')
    objects = objects[::-1]
    x_matrix.append(objects['x'])
    y_matrix.append(objects['y'])
    xy = list(zip(x_matrix[0],y_matrix[0]))
    #print(x_matrix)
    #print(objects['x'])
    wcs_header = ast.solve_from_source_list(x_matrix[0], y_matrix[0],data.shape[1],data.shape[0],solve_timeout=300)
    #wcs_header = ast.solve_from_image(infile, force_image_upload=True)
    data = fits.getdata(infile)
    fits.writeto(outfile,data,header=wcs_header,overwrite=True)
    outfile = wcs_header#outfile
    return outfile
'''
def do_astrometry(infile,outfile):
    ast = AstrometryNet()
    ast.api_key = 'ezduozpdlmrxxisu'
    wcs_header = ast.solve_from_image(infile, force_image_upload=True)
    print(type(wcs_header))
    data = fits.getdata(infile)
    fits.writeto(outfile,data,header=wcs_header,overwrite=True)
    outfile = wcs_header#outfile
    return outfile
'''
do_astrometry(img,'%s_astnet.fits'%(img[0:-5]))

end = time.time()
print('Executed in ',end-start, 'seconds')
