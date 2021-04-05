#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: ClassStacking.py
Name: Chia-Lin Ko
Create Date: March 29, 2021
------------------------
This program aims to do the stacking
"""
import numpy as np
import gvar as gv
import astropy.io.fits as pyfits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from photutils import aperture_photometry
from photutils import EllipticalAperture

# my own packing
from ClassStatistic import MyStatistic

class Stacking:

    def __init__(self, filename):
        self.fn     = filename
        self.f      = pyfits.getdata(self.fn)
        self.f_hd   = pyfits.getheader(self.fn)
        self.f_hd_d = {}
        self.set_beam_info()
       
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
        stacking_num, column_num = np.shape(coord_pix_ds9_arr)
        stacking_1d_arr     = np.zeros(stacking_num)
        stacking_3d_arr     = np.zeros((stacking_num, stacking_width, stacking_width))
        rms_1d_arr          = np.zeros(stacking_num)
        rms_3d_arr          = np.zeros((stacking_num, stacking_width, stacking_width))
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

        # calculate the error
        mask_half_width_pix = int(self.f_hd_d['bmaj_pix'])
        # theoretical error = orginal rms/ sqrt(N)
        stacking_err_theor  = np.divide(img_rms, np.sqrt(stacking_num))
        stacking_err_corr   = self.cal_correlated_error(stacking_2d, mask_half_width_pix, num_times=1000)
        # correlated error
        print('theoretical error =', stacking_err_theor)
        print('correlated error =', stacking_err_corr)

        return stacking_1d, stacking_2d

    def cal_correlated_error(self, img_stacking, mask_half_width, num_times=1000):
        """
        Calculate the correlated flux error by randomly measure the aperature 
        from the central-object-masked image and measure stddev

        Input Parameter
            img_stacking_arr    [2d-ndarray]: image in 2d array of 
            mask_half_width     [int]       : half length of the central mask area (pixel)
            num_times           [int]       : number of repeating times of calcualte the error

        Return
            stacking_err_corr   [float]     : correlated error
        """
        # convert the unit of image (uJy/beam to uJy/pixel)
        img_stacking_pix        = np.divide( img_stacking, self.f_hd_d['beam_area_pix'])
        img_stacking_mask_pix   = self.make_mask_image(img_stacking_pix, mask_half_width)
        pos_sample_pix          = self.make_random_position(img_stacking_pix, mask_half_width, num_times)

        # measure the correlated error
        pos_num     = len(pos_sample_pix)
        aper_arr    = np.ones(pos_num)
        for i, pos_pix in enumerate(pos_sample_pix):
            aperture    = EllipticalAperture(pos_pix,   a=self.f_hd_d['bmaj_pix']/2, b=self.f_hd_d['bmin_pix']/2, 
                                                        theta=self.f_hd_d['bpa_deg'])
            phot_table  = aperture_photometry(img_stacking_mask_pix, aperture)
            aper_arr[i] = phot_table[0][3]
        stacking_err_corr    = np.std(aper_arr)

        return stacking_err_corr

    def make_random_position(self, img_org_arr, mask_half_width, num_times):
        """
        mask the outer and central masked image for randomly calcualte the rms

        Input Parameter
            img_org_arr         [2d-ndarray]: original image for masking
            mask_half_width     [int]       : half length of the central mask area (pixel)
            num_times           [int]       : number of repeating times of calcualte the error
        Return
            pos_sample_pix_arr  [2d-dnarray]: postion sample
        """
        ny, nx  = np.shape(img_org_arr)
        nx_c    = nx//2
        ny_c    = ny//2
        bm_pix  = int(self.f_hd_d['bmaj_pix'])
        # outer layer
        nx_o1   = bm_pix
        nx_o2   = nx - bm_pix
        ny_o1   = bm_pix
        ny_o2   = ny - bm_pix
        # inner layer
        nx_i0   = int(nx_c-mask_half_width-bm_pix)
        nx_i1   = int(nx_c+mask_half_width+bm_pix)
        ny_i0   = int(ny_c-mask_half_width-bm_pix)
        ny_i1   = int(ny_c+mask_half_width+bm_pix)

        # create the mask array within certain range of x and y 
        mask_arr = np.full((ny, nx), True, dtype=bool)
        mask_arr[:ny_o1]    = False # mask the outer upper region
        mask_arr[ny_o2:]    = False # mask the outer bottom region
        mask_arr.T[:nx_o1]  = False # mask the outer left region
        mask_arr.T[nx_o2:]  = False # mask the outer right region
        mask_arr[ny_i0:ny_i1+1].T[nx_i0:nx_i1+1] = False # mask the inner central region

        # make the masked position (pixel)
        x       = np.arange(0, nx)
        y       = np.arange(0, ny)
        xx, yy  = np.meshgrid(x, y)
        xx_mask = xx[mask_arr]
        yy_mask = yy[mask_arr]

        # generate the random position (pixel) from the masked position
        pos_pix_arr         = np.stack((xx[mask_arr], yy[mask_arr]), axis=-1) # [2d-ndarray]
        sample              = np.random.randint(len(pos_pix_arr), size=num_times)
        pos_sample_pix_arr  = pos_pix_arr[sample]

        return pos_sample_pix_arr

    def make_mask_image(self, img_org_arr, mask_half_width):
        """
        make the central masked image

        Input Parameter
            img_org_arr     [2d-ndarray]: original image for masking
            maks_length     [flaot]     : length of the mask area
        Return
            img_mask_arr    [2d-dnarray]: masked image
        """
        ny, nx  = np.shape(img_org_arr)
        nx_c    = nx//2
        ny_c    = ny//2
        x0      = int(nx_c-mask_half_width) # the start value in the x-axis for masking (pixel)
        x1      = int(nx_c+mask_half_width) # the end value in the x-axis for masking (pixel)
        y0      = int(ny_c-mask_half_width) # the start value in the y-axis for masking (pixel)
        y1      = int(ny_c-mask_half_width) # the end value in the y-axis for masking (pixel)

        # create the mask array within certain range of x and y 
        mask_arr = np.full((ny, nx), False, dtype=bool)
        mask_arr[y0:y1+1].T[x0:x1+1] = True
        # mask the image based on the mask array (Ture value)
        image_mask_arr = np.where(mask_arr, np.nan, img_org_arr)

        return image_mask_arr

    def make_rms_map(self, img_rms_fn, img_rms):
        """
        Set the rms map

        Input Parameter
            img_rms_fn  [str]       : filename of rms image. None if no rmg image
            img_rms     [float]     : rms value of the image
        Return
            img_rms_f   [2d-ndarray]: rms image
        """
        if img_rms_fn is not None:
            img_rms_f   = pyfits.getdata(img_rms_fn)
        else:
            img_rms_f   = self.f*0 + img_rms

        return img_rms_f

    def set_beam_info(self):
        """set the beam information from the header to the f_hd_d dictionaty"""

        self.f_hd_d['pix_deg']          = self.f_hd['CDELT2']                   # 1 pixel length           [deg]
        self.f_hd_d['pix_arcs']         = self.deg2arcs(self.f_hd_d['pix_deg']) # 1 pixel length           [arcs]
        self.f_hd_d['bmaj_deg']         = self.f_hd['BMAJ']                     # major axis of beam size  [deg]
        self.f_hd_d['bmin_deg']         = self.f_hd['BMIN']                     # minor axis of beam size  [deg]
        self.f_hd_d['bpa_deg']          = self.f_hd['BPA']                      # angle of the beam        [deg]
        # major axis of beam size [pixel]
        self.f_hd_d['bmaj_pix']         = self.deg2pixel(self.f_hd_d['bmaj_deg'], self.f_hd_d['pix_deg'] )
        # minor axis of beam size [pixel]  
        self.f_hd_d['bmin_pix']         = self.deg2pixel(self.f_hd_d['bmin_deg'], self.f_hd_d['pix_deg'] )
        # beam area [pixel]
        self.f_hd_d['beam_area_pix']    = np.divide(np.pi * self.f_hd_d['bmaj_pix'] * self.f_hd_d['bmin_pix'], 
                                                    4 * np.log(2))

    @staticmethod
    def deg2arcs(input_deg):
        return input_deg * 36e2

    @staticmethod
    def deg2pixel(input_deg, pix_deg):
        return np.divide(input_deg, pix_deg)

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