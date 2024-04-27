#!/usr/bin/env python
# coding: utf-8

import numpy as np
# import matplotlib.pyplot as plt
from astropy.io import ascii, fits
import sys
from astropy.wcs import WCS
import os
from pathlib import Path

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

foldername = Path('.')
#foldername = Path('ast_new')
#flist = [PosixPath(sys.argv[1][0:-5]+'_ast.fits')]
flist = list(foldername.glob('../extras/'+sys.argv[1][0:-5]+'*ast.fits'))
nfiles = len(flist)

savefolder_ast = foldername.parent / Path('../extras')
savefolder_wcs = foldername.parent / Path('../ast_calibrated')

'''
Functions that are used later on
'''
def pixel_to_sky(x, y, f1, f2, f3, f4, f5, g1, g2, g3, g4, g5):
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
    
def generate_wcs_header(header):
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
    
        
    header = generate_wcs_header(header)

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
        
    # ### Vibhore relation
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
from astropy.table import Table
namelist = Table([orig_name_list, save_name_list])
namelist.write('namelist.dat', format='ascii', overwrite=True)
