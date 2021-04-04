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
        self.fn     = filename
        self.f      = pyfits.getdata(self.fn)
        self.f_hd   = pyfits.getheader(self.fn)

        
        # convert the image to 2 dimension
        self.f = self.make_image_2d(self.f)

    def make_image_2d(self, fitsfiles):
        """
        convert the image to 2 dimension
        """
        # check the dimension of the image
        f_shape = np.shape(fitsfiles)    # the shape of the input data( array)
        f_dim   = len(f_shape)   # the dimension of the input data


        if f_dim == 4:
            fitsfiles = fitsfiles[0][0]
            print('Covert the 4D image to 2D')
        elif f_dim == 2:
            print('The input fitfile is a 2D image.')
        else:
            print('Please check the dimension of the input fitsfile')

        return fitsfiles


    def make_stacking(self, img_rms_fn, img_rms, coord_pix_ds9_arr, stacking_stat='mean', 
                      stacking_half_width='20', stackging_2d_fn='stacking.fits'):
        """
        stack the image based on the input coordinate info

        Input Parameter
            img_rms_fn          [str]           : filename of the rms image
            img_rms             [float]         : rms value of the original image 
            coord_pix_ds9_arr   [2d-ndarray]    : 2d array contains the ra and dec arr
            stacking_stat       [str]           : statistic method for stacking 
                                                  (e.g., mean, weighted mean, median)
            stacking_half_width [int]           : half side width of the stacking square (pixel)
            stackging_2d_fn     [str]           : filename for the output stacking 2d image

        Return
            stacking_1d         [float]         : value of the stacking at the postion
            stacking_2d         [2d-ndarray]    : 2d image of the stacking

        """
        # set the range for stacking (2d array)  
        stacking_width = 2*stacking_half_width + 1      # the width (side length) of stacking area [pixel]

        # set the rms map
        img_rms_f           = self.make_rms_map(img_rms_fn, img_rms)
        img_rms_f           = self.make_image_2d(img_rms_f)
        
        # stacking in 2d
        row_num, column_num = np.shape(coord_pix_ds9_arr)
        stacking_1d_arr     = np.zeros(row_num)
        stacking_3d_arr     = np.zeros((row_num, stacking_width, stacking_width))
        rms_1d_arr          = np.zeros(row_num)
        rms_3d_arr          = np.zeros((row_num, stacking_width, stacking_width))
        # 1-based (as in the FITS convention, for example coordinates coming from ds9)
        ra_pix_ds9_arr      = coord_pix_ds9_arr.T[0]    # ra    (pixel, 1-based)
        dec_pix_ds9_arr     = coord_pix_ds9_arr.T[1]    # dec   (pixel, 1-based)


        for i, coord_pix in enumerate(coord_pix_ds9_arr):    
            # 0-based (as in Numpy arrays)
            ra_pix_np           = int(np.float(ra_pix_ds9_arr[i]) -1)  # ra    (pixel, 0-based)
            dec_pix_np          = int(np.float(dec_pix_ds9_arr[i])-1)  # dec   (pixel, 0-based)
            ra_pix_np_start     = ra_pix_np  - stacking_half_width
            ra_pix_np_end       = ra_pix_np  + stacking_half_width
            dec_pix_np_start    = dec_pix_np - stacking_half_width
            dec_pix_np_end      = dec_pix_np + stacking_half_width

            # stacking in 1d 
            stacking_1d_arr[i]  = self.f[dec_pix_np][ra_pix_np]
            rms_1d_arr[i]       = img_rms_f[dec_pix_np][ra_pix_np]

            # stacking in 2d 
            # split the 2d area
            stacking_2d_arr     = self.f[dec_pix_np_start:dec_pix_np_end+1].T[ra_pix_np_start:ra_pix_np_end+1].T
            rms_2d_arr          = img_rms_f[dec_pix_np_start:dec_pix_np_end+1].T[ra_pix_np_start:ra_pix_np_end+1].T
            
            # save the stacking image to a 3d array
            stacking_3d_arr[i] = stacking_2d_arr
            rms_3d_arr[i]      = rms_2d_arr

        # stacking for different statistic method
        if stacking_stat == 'mean':
            stacking_1d     = np.nanmean(stacking_1d_arr)
            stacking_2d     = np.nanmean(stacking_3d_arr, axis=0)
        elif stacking_stat == 'weighted_mean':
            stacking_1d     = MyStatistic.nanaverage(stacking_1d_arr, 1/rms_1d_arr, axis=0)

            stacking_2d     = MyStatistic.nanaverage(stacking_3d_arr, 1/rms_3d_arr, axis=0)
        elif stacking_stat == 'mediam':
            stacking_1d     = np.nanmedian(stacking_1d_arr)
            stacking_2d     = np.nanmedian(stacking_3d_arr, axis=0)
        else:
            print('Please select the stacking statistic from: (1) mean (2) weighted_mean (3) mediam.')

        # save the 2d stacking to FITS file
        pyfits.writeto(stackging_2d_fn, stacking_2d, self.f_hd, overwrite=True)

        # calculate the rms



        return stacking_1d, stacking_2d


    def make_rms_map(self, img_rms_fn, img_rms):
        # set the rms map
        if img_rms_fn is not None:
            img_rms_f   = pyfits.getdata(img_rms_fn)
        else:
            img_rms_f   = self.f*0 + img_rms

        return img_rms_f


    def wcs2pix(self, coord_arr, coord_name_lst, is_savetxt=False, outname='coord_pix.txt'):
        """
        Transfor the world coordinates to pixel from the input FITS file
        
        Input Parameter
            fits_fn         [str]       : the filename of a FITs file
            coord_arr       [2D-ndarray]: right ascension and declination   (deg in wcs)
        Return
            coord_pix_arr   [2D-ndarray]: right ascension and declination   (pixel)
        """   
        
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