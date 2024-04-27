
# for checking 25 stars and plotting histogram

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lmfit import Model
import math
from lmfit import Model
from itertools import combinations as comb
from scipy.stats import norm
from astropy.io.fits import getdata
import sep
from astropy.io import fits

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
x1 = np.array(x)
y1 = np.array(y)
#A=cor_ra[8]
#D=cor_dec[8]
y0 = 18432.0#x[8]
x0 = 2048.0#y[8]

#*******************************************************************************************
#**************************** Using LMFit **************************************************

#*******************************************************************************************

def func(X,g1,f1,g2,f2,f3):
    x = X[0]
    y = X[1]

    ra  = f1 + ((x-x0) * f2) + ((y-y0) * f3)
    dec = g1 + ((x-x0) * g2) #+ ((x-x0)*(x-x0)*g3)

    return np.concatenate((ra,dec))

fmodel=Model(func)
params= fmodel.make_params(g1=10,f1=10,g2=10,f2=10,f3=10)
result=fmodel.fit(np.concatenate((cor_ra,cor_dec)),params,X=(x1,y1))

dict = result.best_values
values = [float(x) for x in list(dict.values())]
popt=[]
for value in values:
    popt.append(value)


print(popt)
p5 = func((x1,y1), popt[0], popt[1],popt[2],popt[3],popt[4])
#p5 = func((x_matrix[0],y_matrix[0]), popt[0], popt[1])
hf = int(len(p5)/2)
ras=p5[0:hf]
decs=p5[hf:]


#*******************************************************************************************

#*****************Plots******************

RA_deg_all = ras
DEC_deg_all = decs


#RA_diff = np.array(cor_ra) - np.array(RA_deg_all)
#DEC_diff = np.array(cor_dec) - np.array(DEC_deg_all)
#RA_diff_arcsec = RA_diff*3600
#DEC_diff_arcsec = DEC_diff*3600

from astropy.stats import sigma_clip

import seaborn as sns
from scipy import stats
import matplotlib.patches as mpatches
import glob

files = glob.glob('202*.fits')
filename = files[0]


def plot_scatter(xpix,ypix,RA_new,DEC_new,RA_g,DEC_g,save_file_name):
    fnt_sz,leg_sz=12,10
    #######
    plt.figure(1, figsize=(10,4))
    plt.subplots_adjust(wspace=0.4,hspace=0.5)
    #plt.text(6000, 2.8, 'q1 = (RA - RA0) - ( (y-y0)*f1 + f2)', ha='center',fontsize=12)
    #########
    DEC_diff_arcsec = (DEC_new - DEC_g)*3600.0
    RA_diff_arcsec = (RA_new - RA_g)*3600.0
    RA_diff_arcsec = sigma_clip(RA_diff_arcsec, sigma=3, maxiters=None,masked=False, copy=False)
    DEC_diff_arcsec = sigma_clip(DEC_diff_arcsec, sigma=3, maxiters=None, masked=False, copy=False)
    ####################
    ax1 = plt.subplot(1, 2, 1)
    cpdf_1 = (RA_new - RA_g)*3600.0
    #print(cpdf_1)
    plt.xlabel('Pixels along RA',fontsize=fnt_sz)
    plt.ylabel( 'RA gaia - RA cal (arcsec)' ,fontsize=fnt_sz)
    plt.scatter(ypix, cpdf_1,c=xpix,marker='.',s=2)
    plt.colorbar(label="Dec pixel", orientation="vertical")
    ##############s

    ax3 = plt.subplot(1, 2, 2)
    cpdf_1 = (DEC_new - DEC_g)*3600.0
    plt.xlabel('Pixels along DEC',fontsize=fnt_sz)
    plt.ylabel( 'DEC gaia - DEC cal (arcsec)' ,fontsize=fnt_sz)
    #plt.scatter(xpix, cpdf_1,c=ypix,marker='.',s=2)
    plt.scatter(xpix, cpdf_1,c=ypix,marker='.',s=2)
    plt.colorbar(label="RA pixel", orientation="vertical")
    ##############
    ############# 
    #plt.suptitle('Astrometric offsets in one of the TDI frames observed on 2022-10-24 ')
    plt.savefig(save_file_name,bbox_inches='tight',pad_inches=0.2)
    #plt.show()
    plt.clf()
    plt.cla()
    plt.close()

def plot_histogram(xpix,ypix,RA_new,DEC_new,RA_g,DEC_g,save_file_name):
    fnt_sz,leg_sz=12,10
    #######
    plt.figure(1, figsize=(10,5))
    plt.subplots_adjust(wspace=0.25,hspace=0.25)
    #plt.text(6000, 2.8, 'q1 = (RA - RA0) - ( (y-y0)*f1 + f2)', ha='center',fontsize=12)
    #########
    DEC_diff_arcsec = (DEC_new - DEC_g)*3600.0
    RA_diff_arcsec = (RA_new - RA_g)*3600.0
    RA_diff_arcsec = sigma_clip(RA_diff_arcsec, sigma=3, maxiters=None,masked=False, copy=False)
    DEC_diff_arcsec = sigma_clip(DEC_diff_arcsec, sigma=3, maxiters=None, masked=False, copy=False)
    ####################
    ax2 = plt.subplot(1, 2, 1)

    nbins = 75
    n, bins, patches = ax2.hist(RA_diff_arcsec,nbins, density=True, facecolor = 'grey', alpha = 0.5); 

    centers = (0.5*(bins[1:]+bins[:-1]))
    pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), centers, n, p0=[0,1])

    ax2.plot(centers, norm.pdf(centers,*pars), 'k--',linewidth = 2) 

  
    plt.xlabel('RA gaia - RA cal (arcsec)',fontsize=fnt_sz); plt.ylabel( 'PDF' ,fontsize=fnt_sz)
    handles, labels = plt.gca().get_legend_handles_labels()
    empty_patch = mpatches.Patch(color='none'); handles.append(empty_patch)  
    #emp_label1 = u'$\mu$='+ str('%.2E' % mu); emp_label2 = u'$\sigma$='+ str('%.2E' % sigma)
    emp_label1 = u'$\mu$='+ str('%.2f' % abs(pars[0])); emp_label2 = u'$\sigma$='+ str('%.2f' % pars[1])
    labels.append(emp_label1+'\n'+emp_label2); plt.legend(handles, labels,frameon=False,loc='upper left',prop={'size': leg_sz})
    #ax2.set_title('$\mu={:.4f}\pm{:.4f}$, $\sigma={:.4f}\pm{:.4f}$'.format(pars[0],np.sqrt(cov[0,0]), pars[1], np.sqrt(cov[1,1 ])))
    ###################


    ax4 = plt.subplot(1, 2, 2)

    n, bins, patches = ax4.hist(DEC_diff_arcsec,nbins, density=True, facecolor = 'grey', alpha = 0.5); 

    centers = (0.5*(bins[1:]+bins[:-1]))
    pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), centers, n, p0=[0,1])

    ax4.plot(centers, norm.pdf(centers,*pars), 'k--',linewidth = 2) 

    plt.xlabel('DEC gaia - DEC cal (arcsec)',fontsize=fnt_sz); plt.ylabel( 'PDF' ,fontsize=fnt_sz)
    handles, labels = plt.gca().get_legend_handles_labels()
    empty_patch = mpatches.Patch(color='none'); handles.append(empty_patch)  
    #emp_label1 = u'$\mu$='+ str('%.2E' % mu); emp_label2 = u'$\sigma$='+ str('%.2E' % sigma)
    emp_label1 = u'$\mu$='+ str('%.2f' % pars[0]); emp_label2 = u'$\sigma$='+ str('%.2f' % pars[1])
    labels.append(emp_label1+'\n'+emp_label2); plt.legend(handles, labels,frameon=False,loc='upper left',prop={'size': leg_sz})    
    ############# 
    #plt.suptitle('Astrometric offsets in one of the TDI frames observed on 2022-10-24 ')
    plt.savefig(save_file_name,bbox_inches='tight',pad_inches=0.2)
    #plt.show()
    plt.clf()
    plt.cla()
    plt.close()

plot_scatter(x,y,RA_deg_all,DEC_deg_all,cor_ra,cor_dec,'plot_scatter_4f.pdf')
plot_histogram(x,y,RA_deg_all,DEC_deg_all,cor_ra,cor_dec,'plot_histogram_4f.pdf')

