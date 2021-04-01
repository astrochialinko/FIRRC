#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: ClassStacking.py
Name: Chia-Lin Ko
Create Date: March 29, 2021
------------------------
This program aims to do the stacking
"""
import numpy as np
import astropy.io.fits as pyfits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import gvar as gv

# my own packing
from ClassStatistic import MyStatistic

class Stacking:

    def __init__(self, filename):
        self.fn   = filename
        self.f    = pyfits.getdata(self.fn)
        self.f_hd = pyfits.getheader(self.fn)
    
    def make_stacking_1d(self, img_rms_fn, img_rms, coord_pix_ds9_arr, stacking_stat):
        
        # set the rms map
        if img_rms_fn is not None:
            img_rms_f   = pyfits.getdata(img_rms_fn)
        else:
            img_rms_f   = self.f*0 + img_rms

        # stacking in 1d 
        row_num, column_num = np.shape(coord_pix_ds9_arr)
        stacking_1d_arr     = np.zeros(row_num)
        rms_1d_arr          = np.zeros(row_num)
        # 1-based (as in the FITS convention, for example coordinates coming from ds9)
        ra_pix_ds9_arr      = coord_pix_ds9_arr.T[0]    # ra    (pixel, 1-based)
        dec_pix_ds9_arr     = coord_pix_ds9_arr.T[1]    # dec   (pixel, 1-based)

        for i, coord_pix in enumerate(coord_pix_ds9_arr):    
            # 0-based (as in Numpy arrays)
            ra_pix_np           = int(np.float(ra_pix_ds9_arr[i])-1)  # ra    (pixel, 0-based)
            dec_pix_np          = int(np.float(dec_pix_ds9_arr[i])-1)  # dec   (pixel, 0-based)
            stacking_1d_arr[i]  = self.f[0][0][dec_pix_np][ra_pix_np]
            rms_1d_arr[i]       = img_rms_f[0][0][dec_pix_np][ra_pix_np]

        # stacking for different statistic method
        if stacking_stat == 'mean':
            stacking_1d     = np.nanmean(stacking_1d_arr)
        elif stacking_stat == 'weighted_mean':
            stacking_1d     = MyStatistic.nanaverage(stacking_1d_arr, 1/rms_1d_arr, axis=0)
        elif stacking_stat == 'medium':
            stacking_1d     = np.nanmedian(stacking_1d_arr)
        else:
            print('Please select the stacking statistic from: (1) mean (2) weighted_mean (3) medium.')

        # error

        return stacking_1d

    def wcs2pix(self, coord_arr, coord_name_lst, is_savetxt=False, outname='coord_pix.txt'):
        '''
        Transfor the world coordinates to pixel from the input FITS file
        
        Input Parameter
            fits_fn         [str]       : the filename of a FITs file
            coord_arr       [2D-ndarray]: right ascension and declination   (deg in wcs)
        Return
            coord_pix_arr   [2D-ndarray]: right ascension and declination   (pixel)
        '''   
        
        # loading WCS information from the FITS file
        wcs             = WCS(self.f_hd)
        coord_pix_arr   = np.zeros((len(coord_arr), len(coord_name_lst)))
              
        for i in range(len(coord_arr)):

            # get the info of world coordinate (RA, Dec) (deg)
            crd_ra      = coord_arr.T[0][i]
            crd_dec     = coord_arr.T[1][i]
            
            # transforms world coordinates to pixel coordinates
            # 1-based (as in the FITS convention, for example coordinates coming from DS9)
            pixcrd                  = wcs.wcs_world2pix([[crd_ra, crd_dec, 0, 0]], 1)[0] 
            coord_pix_arr.T[0][i]   = int(np.around(pixcrd[0]))  # ra_pix,  get the closest pixel by round the number
            coord_pix_arr.T[1][i]   = int(np.around(pixcrd[1]))  # dec_pix, get the closest pixel by round the number

        if is_savetxt:
            # save the coordinate in pixel to the txt file
            columns     = " ".join(coord_name_lst)
            fmt_lst     = ['%d' for i in coord_name_lst]
            fmt_str     = ' '.join(fmt_lst)
            np.savetxt(outname, coord_pix_arr, fmt=fmt_str, header=columns)
        
        return coord_pix_arr