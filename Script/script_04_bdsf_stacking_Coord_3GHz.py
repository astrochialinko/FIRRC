#! /anaconda3/bin/python
"""
File: script_bdsf_stacking_binning_irac.py
Name: Chia-Lin Ko
Create Date: Jun 11, 2021
Last Modified Date: Jun 15, 2021
------------------------
This program aims to calculate the flux density
"""
import os
import bdsf
import numpy as np
import gvar as gv
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

##### Global paramters and keywords ######################
# path
PATH_DATA_IN        = '../Data/COSMOS/Image/Stacking/XS_3GHzDet/'
# PATH_DATA_IN        = '../Data/COSMOS/Image/Stacking/XS_IRACDet/'
# PATH_DATA_IN        = '../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/'


# list
IMG_LABEL_LIST      = ['3GHz', '1d4GHzXS']
BEAM_FACTOR = [10]
THLD_ISL = [5]
THLD_PIX = [5]
##########################################################

def main():
    stacking_stat   = 'mediam'
    for beam_factor in BEAM_FACTOR:
        beam_factor_str = '%.1f'%(beam_factor)
        for thld_isl in THLD_ISL:
            for thld_pix in THLD_PIX:

                folder_lst = listdir_nohidden(PATH_DATA_IN)
                for folder in folder_lst:
                    # create a directory for saving Pybdsf results
                    path_data_out      = '%s%s/Pybdsf_bm%s_thIsl%s_thPix%s/'%(PATH_DATA_IN, folder, beam_factor_str, thld_isl, thld_pix)
                    # remove the previous results
                    command = 'rm -r %s'%(path_data_out) 
                    os.system(command)
                    # create a new directory
                    command = 'mkdir %s'%(path_data_out) 
                    os.system(command)

                    img_fn_lst = listdir_nohidden(PATH_DATA_IN+folder)
                    for img_fn in img_fn_lst:
                        img_fn = img_fn.split('.fits')[0]
                        input_fits_fn = '%s%s/%s.fits'%(PATH_DATA_IN, folder, img_fn)
                        try:
                            x_min, x_max, y_min, y_max = get_trim_box(input_fits_fn, beam_factor=beam_factor)
                            
                            img = bdsf.process_image(input_fits_fn, 
                                thresh='hard', thresh_isl=thld_isl, thresh_pix=thld_pix,
                                advanced_opts = True, ini_gausfit= 'default', trim_box=(x_min, x_max, y_min, y_max),
                                output_opts=True, plot_allgaus=True, plot_islands=True
                                )
                                 
                            # Write the source list catalog. File is named automatically.
                            img.write_catalog(outfile='%s%s.pybdsm.gaul.csv'%(path_data_out, img_fn), 
                                format='csv', catalog_type='gaul', clobber=True)
                            
                            # Write the residual image. File is named automatically.
                            img.export_image(outfile='%s%s.pybdsm_gaus_resid.fits'%(path_data_out, img_fn),
                                img_type='gaus_resid', img_format='fits', clobber=True)
                            
                            # Write the model image. Filename is specified explicitly.
                            img.export_image(outfile='%s%s.pybdsm_gaus_model.fits'%(path_data_out, img_fn),
                                img_type='gaus_model', img_format='fits', clobber=True)

                            # move the results into folder
                            command = 'mv %s%s/%s_pybdsm %s'%(PATH_DATA_IN, folder, img_fn, path_data_out) 
                            os.system(command)
                            command = 'mv %s%s/%s.fits.pybdsf.log %s'%(PATH_DATA_IN, folder, img_fn, path_data_out) 
                            os.system(command)

                        except:
                            print('File %s not found'%(input_fits_fn))


def listdir_nohidden(path):
    return [f for f in os.listdir(path) if not f.startswith('.')]

def get_trim_box(fits_fn, beam_factor = 3):
    img_hd = pyfits.getheader(fits_fn)
    img_hd_dict = {}
    img_hd_dict['pix_deg']          = img_hd['CDELT2']
    img_hd_dict['pix_arcs']         = img_hd_dict['pix_deg']*36e2
    img_hd_dict['bmaj_deg']         = img_hd['BMAJ']
    img_hd_dict['bmaj_pix']         = img_hd_dict['bmaj_deg']/img_hd_dict['pix_deg']
    img_hd_dict['img_size_pix']     = img_hd['NAXIS1']
    img_hd_dict['img_size_arcs']    = img_hd_dict['img_size_pix']*img_hd_dict['pix_arcs']
    if img_hd_dict['img_size_pix']%2 == 0:
        print('Warning! The length of image is an even number.')
        print('Please check the trim_box value carefully.')
    center_pix = int(img_hd_dict['img_size_pix']/2)
    trim_radius_pix = round(beam_factor*img_hd_dict['bmaj_pix']/2)
    x_min = y_min = center_pix - trim_radius_pix
    x_max = y_max = center_pix + trim_radius_pix+1

    return x_min, x_max, y_min, y_max

if __name__ == '__main__':
    main()
