20
# for checking 25 stars and plotting histogram

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from astropy.modeling import models, fitting
from lmfit import Model
from scipy.optimize import curve_fit
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
import math
from lmfit import Model
from itertools import combinations as comb
from scipy.stats import norm
from astropy.io.fits import getdata
import sep
import glob
import sys

filename = sys.argv[1]


x,y,gaia_ra,gaia_dec,cor_ra,cor_dec=[],[],[],[],[],[]
z='../extras/list_selected_epoch_corrected.txt'
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
x1 = np.array(x)
y1 = np.array(y)
A=cor_ra[0]
D=cor_dec[0]
x0 = x[0]
y0 = y[0]
#*******************************************************************************************
#**************************** Using LMFit **************************************************

#*******************************************************************************************
'''
def func(X,g1,f1,g2,f2):
    x = X[0]
    y = X[1]
    #print(x,y)
    ra = A +  (f1*y-y0)/f2 
    dec = D +  (g1*x-x0)/g2
    return np.concatenate((ra,dec))
'''

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
ras=p5[0:hf]
decs=p5[hf:]


#*******************************************************************************************

#*****************Plots******************
'''
fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,5))#,sharey=True)
fig.suptitle('Image: 19990312_1210_V.fits, relation used: RA = f1 + (x * f2) + (y * f3) and DEC = g1 + (y * g2) + (x * g3) ',color='r')

ax1.title.set_text('RA offset')
ax1.plot(y,(np.array(cor_ra) - np.array(ras))*3600,'bo',label='RA_gaia - RA_fitted')
ax1.set_xlabel('ccd y_axis')
ax1.set_ylabel('RA')
#ax1.set_ylim(-1.0,1.0)
ax1.legend()
ax2.title.set_text('DEC offset')
ax2.plot(x,(np.array(cor_dec) - np.array(decs))*3600,'bo',label='DEC_gaia - DEC_fitted')
ax2.set_xlabel('ccd x_axis')
ax2.set_ylabel('DEC')
#ax2.set_ylim(-1.0,1.0)
ax2.legend()
plt.savefig('chunk_offsets.pdf')
plt.close()




RA_deg_all = ras
DEC_deg_all = decs


RA_diff = np.array(cor_ra) - np.array(RA_deg_all)
DEC_diff = np.array(cor_dec) - np.array(DEC_deg_all)
RA_diff_arcsec = RA_diff*3600
DEC_diff_arcsec = DEC_diff*3600

(mu_1, sigma_1) = norm.fit(RA_diff_arcsec)
(mu_2, sigma_2) = norm.fit(DEC_diff_arcsec)
idx = np.where(DEC_diff_arcsec > -1.5)
fig, (ax1,ax2) = plt.subplots(1, 2,figsize=(10,5))#,sharey=True)
#fig.suptitle('Image: 19990312_1210_V.fits, relation used: (RA-RA0)*Cos(DEC-DEC0) = (y-y0)/f1 +f2 and DEC-DEC0 = (x-x0)/g1 +g2',color='r')
ax1.title.set_text('Histogram showing offset in RA')
#ax1.text(0.5, 25, r'mu= %f'%(mu_1), fontsize=10)
#ax1.text(0.5, 24, r' sigma= %f'%(sigma_1), fontsize=10)
ax1.hist(RA_diff_arcsec,label='RA_gaia - RA_fitted')
ax1.set_xlabel('Offset(arc-sec)')
ax1.set_ylabel('No. of sources')
ax1.set_xlim(-1.0,1.0)
ax1.grid()
ax1.minorticks_on()
ax1.grid(which='major', linestyle='-', linewidth='0.5', color='red')
ax1.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
ax1.legend()


ax2.title.set_text('Histogram showing offset in DEC')
#ax2.text(-1, 40, r'mu= %f'%(mu_2), fontsize=10)
#ax2.text(-1, 38, r' sigma= %f'%(sigma_2), fontsize=10)
ax2.hist(DEC_diff_arcsec[idx],label='DEC_gaia - DEC_fitted')
ax2.set_xlabel('Offset(arc-sec)')
ax2.set_ylabel('No. of sources')
ax2.grid()
ax2.set_xlim(-1.0,1.0)
ax2.minorticks_on()
ax2.grid(which='major', linestyle='-', linewidth='0.5', color='red')
ax2.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
ax2.legend()
plt.savefig('chunk_offsets_hist.pdf')
plt.close()
'''

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


with open("../extras/astrometric_xy.txt", "w+") as d:
	d.write('id\tx\t\y\tAstro_ra\tAstro_dec\n')
	for x in range(len(all_ra)):
		d.write("\n%.2d\t "%(x+1))
		d.write("%4.2f\t "%(xy[x][0]))
		d.write("%4.2f\t "%(xy[x][1]))
		d.write("%.9f\t"%(all_ra[x]))
		d.write("%.9f\t "%(all_dec[x]))



