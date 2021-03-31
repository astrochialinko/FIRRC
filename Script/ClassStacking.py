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

class Stacking:

    def __init__(self, filename):
        self.fn   = filename
        self.f    = pyfits.getdata(self.fn)
        self.f_hd = pyfits.getheader(self.fn)

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
            # 1-based (as in the FITS convention, for example coordinates coming from DS9).
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