#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: Script_measure_flux_density.py
Name: Chia-Lin Ko
Create Date: March 7, 2021
------------------------
This program aims to measure the
(1) peak flux (2) integrated flux density and (3) the aperture flux
"""
import numpy as np
import astropy.io.fits as pyfits
from photutils import aperture_photometry
from photutils import EllipticalAperture
from flux_density import FluxDensity

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
from matplotlib import rc, rcParams, style
from matplotlib.colors import ListedColormap
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle

# Constant

# path
PATH_DATA      = '../Data/COSMOS/Image/Stacking/'
PATH_DATA_VLA  = '../Data/COSMOS/Image/VLA/'
PATH_TABLE     = '../Data/Tables/'
PATH_FIGURE    = '../Figures/'

RANDOM_SEED    = 6
INPUT_FILENAME_LIST = [ 'COSMOS_3GHz_stacking_1d4GHz_dp_rand%d_median'%(RANDOM_SEED),
                        'COSMOS_1d4GHz_lp_stacking_1d4GHz_dp_rand%d_median'%(RANDOM_SEED),
                        'COSMOS_1d4GHz_stacking_1d4GHz_dp_rand%d_median'%(RANDOM_SEED)
                      ]

POS_FILENAME = [
                'Coord_pix_3GHz_radio_COSMOS_RadioCat_nAGN_1d4GHz_dp_rand%d'%(RANDOM_SEED),
                'Coord_pix_1d4GHz_lp_radio_COSMOS_RadioCat_nAGN_1d4GHz_dp_rand%d'%(RANDOM_SEED),
                'Coord_pix_1d4GHz_radio_COSMOS_RadioCat_nAGN_1d4GHz_dp_rand%d'%(RANDOM_SEED)
                ]

IMG_FILENAME = [
                '%scosmos_vla_3GHz_2017image_uJy.fits'%(PATH_DATA_VLA),
                '%scosmos_vla_1d4GHz_lp_2010image_uJy.fits'%(PATH_DATA_VLA),
                '%scosmos_vla_1d4GHz_dp_2010image_uJy.fits'%(PATH_DATA_VLA)
                ]

ERROR_LIST = [0.33, 1.88, 2.59]
LABEL_LIST = ['3GHz', '1d4GHz_lp', '1d4GHz_dp']

# switch
IS_FIT_VLA = False
IS_PLOT    = True
IS_BOX     = True

# Global variable



def main():

    for i, input_filename in enumerate(INPUT_FILENAME_LIST):


        if IS_FIT_VLA:
            pos_filename = '%s%s.txt'%(PATH_TABLE, POS_FILENAME[i])
            img_filename = IMG_FILENAME[i]
            #make_fits_image(fd, img_filename, pos_filename)
            sub_num = 2
            subimagename = (img_filename.split('.fits')[0] + '_sub_%d.fits'%(sub_num))
            imagename = subimagename
        else:
            imagename = '%s%s.fits' % (PATH_DATA, input_filename)
        

        fd = FluxDensity(imagename)
        fitsfile    = fd.f
        fitsfile_hd = fd.f_hd

        # 2D Gaussian fitting
        initial_guess   = fd.get_2d_gaussian_ini()
        #initestname = '%s%s_est.txt' % (PATH_TABLE, input_filename)


        # make the box for the gaussian fitting        
        if IS_BOX:
            box_bm_list    = [0.5, 1.0, 1.5, 2.0, 2.5]   # the length of the recentangular region to select for fitting [beam]
            dfit_shape_x = np.ones(len(box_bm_list))
            dfit_shape_y = np.ones(len(box_bm_list))
            box_data_fitted = []

            for j, box_bm in enumerate(box_bm_list):
                box_txtfn = '%s%s_gaufit_box%s.txt'%(PATH_TABLE, input_filename, box_bm)
                box_list  = fd.make_box_list(beam_multiply=box_bm, ismake_txt=True, txtfilename=box_txtfn)
                dfit_shape_x[j] = box_list[2]-box_list[0]+1
                dfit_shape_y[j] = box_list[3]-box_list[1]+1
                gaussian_result, curvefit_result = fd.fit_2d_gaussian(isbox_clip=True, box=box_list)
                box_data_fitted.append(fd.data_fitted)
        else:
            dfit_shape_y = fd.f_shape_y
            dfit_shape_x = fd.f_shape_x
            gaussian_result, curvefit_result = fd.fit_2d_gaussian()



        fd.print_fitting_result(gaussian_result, initial_guess)
        #print("")

        parm_plot_dict = {}
        parm_plot_dict['half_length_arcsec'] =  5
        xy = fd.make_2d_xy()
        x, y = xy
        parm_plot_dict['width_arcs']       = (x.max() - x.min()) * fd.f_pix_arcs
        parm_plot_dict['height_arcs']      = (y.max() - y.min()) * fd.f_pix_arcs
        parm_plot_dict['dfit_shape_x']     = dfit_shape_x
        parm_plot_dict['dfit_shape_y']     = dfit_shape_y
        parm_plot_dict['dfit_x_center']    = dfit_shape_x/2
        parm_plot_dict['dfit_y_center']    = dfit_shape_y/2
        parm_plot_dict['dfit_width_arcs']  = dfit_shape_x * fd.f_pix_arcs
        parm_plot_dict['dfit_height_arcs'] = dfit_shape_y * fd.f_pix_arcs
        if IS_BOX:
            parm_plot_dict['box_list']        = box_list
            parm_plot_dict['box_bm_list']     = box_bm_list
            parm_plot_dict['box_data_fitted'] = box_data_fitted
        
        if IS_PLOT:
            # plot the image
            #plot_image(fd, parm_plot_dict, input_filename, gaussian_result, curvefit_result)
            # plot the fitting result
            #plot_2d_gaussian_fit(fd, parm_plot_dict, input_filename, gaussian_result, curvefit_result)
            # plot the fitting reuslt (cruve of growth)
            plot_curve_of_growth(i, fd, parm_plot_dict, input_filename, gaussian_result, curvefit_result)


        
def make_fits_image(fd, imgfilename, posfilename):

    fitsimage    = pyfits.getdata(imgfilename)[0][0]
    fitsimage_hd = pyfits.getheader(imgfilename)
    pos_txtfile  = np.loadtxt(posfilename, dtype = str)
    lines_num, item_num = np.shape(pos_txtfile)
    half_length_pix     = 50
    half_length_beam    = 5
    half_length_arcsec  = 10

    for i in range(lines_num):
        # 0:RA, 2:Dec, 4: z_m
        # 0-based (as in Numpy arrays)
        pix_ra  = int(np.float(pos_txtfile[i][0])-1)  # ra  (deg)
        pix_dec = int(np.float(pos_txtfile[i][2])-1)  # dec (deg)

        total_length_pix = 2*half_length_pix + 1   # [pixel]
        pix_ra_start  = pix_ra  - half_length_pix
        pix_ra_end    = pix_ra  + half_length_pix
        pix_dec_start = pix_dec - half_length_pix
        pix_dec_end   = pix_dec + half_length_pix

        #----------------------------------------------------------
        # subimg
        subimage     = fitsimage[pix_dec_start:pix_dec_end+1].T[pix_ra_start:pix_ra_end+1].T
        subimagename = (imgfilename.split('.fits')[0] + '_sub_%d.fits'%(i))
        # temp!
        #print(subimage)
        pyfits.writeto(subimagename, subimage , fitsimage_hd, overwrite=True)

        #print(txtfile)



def plot_image(fd, parm_plot_dict, input_filename, gaussian_result, curvefit_result):
    # ---------------------------------
    # plot image
    #data_fitted = fd.func_2d_gaussian(xy, *curvefit_result[0])

    # parameter
    half_length_arcsec = parm_plot_dict['half_length_arcsec']
    width_arcs    = parm_plot_dict['width_arcs']
    height_arcs   = parm_plot_dict['height_arcs']
    dfit_shape_x  = parm_plot_dict['dfit_shape_x']
    dfit_shape_y  = parm_plot_dict['dfit_shape_y']
    dfit_width_arcs  = parm_plot_dict['dfit_width_arcs']
    dfit_height_arcs = parm_plot_dict['dfit_height_arcs']


    # print(-1*round(width_arcs/2), -1*round(height_arcs/2))

    fig, ax = plt.subplots(1, 1)
    im = ax.imshow(fd.f, cmap='viridis', origin='lower',
                   extent=(-1 * (width_arcs / 2), width_arcs / 2,
                           -1 * (height_arcs / 2), height_arcs / 2)
                   )


    ax.contour(fd.data_fitted.reshape(fd.f_shape_x, fd.f_shape_y), 5, colors='w',
               extent=(-1 * (width_arcs / 2), width_arcs / 2,
                       -1 * (height_arcs / 2), height_arcs / 2)
               )
    # beam
    ellipse = Ellipse(np.array([half_length_arcsec - (fd.f_bmaj_deg * 36e2),
                                -half_length_arcsec + (fd.f_bmaj_deg * 36e2)]),
                      width=fd.f_bmaj_deg * 36e2, height=fd.f_bmin_deg * 36e2, angle=fd.f_bpa_deg,
                      color='white', alpha=0.6, \
                      )
    ax.add_patch(ellipse)

    ax.set_title(input_filename)
    ax.set_xlim(-1 * half_length_arcsec, half_length_arcsec)
    ax.set_ylim(-1 * half_length_arcsec, half_length_arcsec)
    ax.set_xlabel('[arcsec]', fontsize=12)
    ax.set_ylabel('[arcsec]', fontsize=12)

    fig.colorbar(im, pad=0.03)
    plt.savefig('%s%s_gaufit_image.pdf' % (PATH_FIGURE, input_filename))
    plt.show()
    plt.close()


def plot_2d_gaussian_fit(fd, parm_plot_dict, input_filename, gaussian_result, curvefit_result):

    half_length_arcsec = parm_plot_dict['half_length_arcsec']
    width_arcs         = parm_plot_dict['width_arcs']
    height_arcs        = parm_plot_dict['height_arcs']
    dfit_shape_x       = parm_plot_dict['dfit_shape_x']
    dfit_shape_y       = parm_plot_dict['dfit_shape_y']
    dfit_x_center      = parm_plot_dict['dfit_x_center']
    dfit_y_center      = parm_plot_dict['dfit_y_center']
    dfit_width_arcs    = parm_plot_dict['dfit_width_arcs']
    dfit_height_arcs   = parm_plot_dict['dfit_height_arcs']
    if IS_BOX:
        box_list           = parm_plot_dict['box_list']
    coeff, coeff_std, chi2, dof, chi2_r = gaussian_result


    # PLOT
    fig1 = plt.figure(1)
    # Plot Data-model
    frame1 = fig1.add_axes((.1, .3, .8, .6))

    # xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
    ind      = np.arange(fd.f_shape_x) * fd.f_pix_arcs - width_arcs / 2
    dfit_ind = (np.arange(dfit_shape_x)+fd.f_pix_arcs) * fd.f_pix_arcs - dfit_width_arcs / 2

    plt.plot(ind, fd.f[int(fd.x_center)], label='stacked (h)', color='tab:blue')
    plt.plot(ind, fd.f.T[int(fd.x_center)], label='stacked (v)', color='tab:purple')
    if IS_BOX:
        print(box_list)
        plt.plot(dfit_ind, fd.f[int(fd.x_center)][box_list[0]:box_list[2]+1], label='stacked box (h)', color='blue')
        plt.plot(dfit_ind, fd.f.T[int(fd.x_center)][box_list[1]:box_list[3]+1], label='stacked box (v)', color='purple')
    
    plt.plot(ind, fd.data_fitted.reshape(fd.f_shape_x, fd.f_shape_y)[int(fd.x_center)], \
             label='gaussian fitting (h)', color='tab:red')  # Best fit model
    plt.plot(ind, fd.data_fitted.reshape(fd.f_shape_x, fd.f_shape_y).T[int(fd.x_center)], \
             label='gaussian fitting (v)', color='tab:orange')  # Best fit model
    frame1.set_xticklabels([])  # Remove x-tic labels for the first frame

    # FWHM
    y_list = np.linspace(fd.f.min(), fd.f.max(), 20)
    plt.plot(((y_list * 0 + coeff[1] + fd.fwhm_x.mean) * fd.f_pix_arcs - width_arcs / 2) / 2, y_list,
             '--', color='lightblue', alpha=0.8, label='FWHM')
    plt.plot(((y_list * 0 + coeff[1] - fd.fwhm_x.mean) * fd.f_pix_arcs - width_arcs / 2) / 2, y_list,
             '--', color='lightblue', alpha=0.8)

    plt.plot(y_list * 0 + fd.f_bmaj_deg * 36e2 / 2, y_list,
             '--', color='gray', alpha=0.5, label='Beamsize')
    plt.plot(y_list * 0 - fd.f_bmaj_deg * 36e2 / 2, y_list,
             '--', color='gray', alpha=0.5)

    plt.title(fd.fn)
    plt.ylabel('Flux density [uJy/beam]', fontsize=12)
    plt.xlim(-1 * half_length_arcsec, half_length_arcsec)
    plt.grid(alpha=0.5)
    plt.legend()

    # Residual plot
    difference = fd.f - fd.data_fitted.reshape(fd.f_shape_x, fd.f_shape_y)

    frame2 = fig1.add_axes((.1, .1, .8, .2))
    plt.plot(ind, difference[int(fd.x_center)]  , label='residual (h)', color='tab:blue')
    plt.plot(ind, difference.T[int(fd.y_center)], label='residual (v)', color='tab:purple')
    
    plt.xlabel('Offset [arcsec]', fontsize=12)
    plt.xlim(-1 * half_length_arcsec, half_length_arcsec)

    plt.legend(loc='upper right')
    plt.grid(alpha=0.5)
    plt.savefig('%s%s_gaufit_compared.pdf' % (PATH_FIGURE, input_filename))
    plt.show()
    plt.close()
    


def plot_curve_of_growth(i, fd, parm_plot_dict, input_filename, gaussian_result, curvefit_result):
    dfit_shape_x       = parm_plot_dict['dfit_shape_x']
    dfit_shape_y       = parm_plot_dict['dfit_shape_y']
    dfit_x_center      = parm_plot_dict['dfit_x_center']
    dfit_y_center      = parm_plot_dict['dfit_y_center']
    dfit_width_arcs    = parm_plot_dict['dfit_width_arcs']
    dfit_height_arcs   = parm_plot_dict['dfit_height_arcs']
    ape_radius = np.linspace(0.2 * (fd.f_bmaj_pix / 2), 4 * (fd.f_bmaj_pix / 2), 19 + 1)

    if IS_BOX:
        box_list       = parm_plot_dict['box_list']
        box_bm_list    = parm_plot_dict['box_bm_list']
        box_data_fitted = parm_plot_dict['box_data_fitted']
        gau_ape_arr = np.ones((len(ape_radius), len(box_bm_list)))
    else:
        gau_ape_arr = np.ones(len(ape_radius))

    
    pos_target = np.array([50, 50])

    ape_arr = np.ones(len(ape_radius))
    ape_err_arr = np.ones(len(ape_radius))
    # gau_ape_err_arr = np.ones(len(ape_radius))

    for j, radius_pix in enumerate(ape_radius):
        
        # stacked image
        aperture = EllipticalAperture(pos_target,
                                      a=radius_pix, b=radius_pix, theta=fd.f_bpa_deg)
        phot_table = aperture_photometry(fd.f / fd.f_beam_area_pix, aperture, error=fd.f * 0 + ERROR_LIST[i])
        ape_arr[j]    = phot_table[0][3]
        ape_err_arr[j] = phot_table[0][4]

        # gaufit
        if IS_BOX:
            for k, box_bm in enumerate(box_bm_list):
                
                gau_phot_table = aperture_photometry(np.divide(box_data_fitted[k].reshape(fd.f_shape_x, fd.f_shape_y),
                                                     fd.f_beam_area_pix),
                                                     aperture)
                gau_ape_arr[j][k] = gau_phot_table[0][3]
        else:
            gau_phot_table = aperture_photometry(np.divide(fd.data_fitted.reshape(fd.f_shape_x, fd.f_shape_y),
                                                 fd.f_beam_area_pix),
                                                 aperture)
            gau_ape_arr[j] = gau_phot_table[0][3]

    plt.errorbar(ape_radius / (fd.f_bmaj_pix/2), ape_arr, yerr=ape_err_arr, label=LABEL_LIST[i])
    # if is_box:
    #     ape_arr = np.ones(len(ape_radius))
    #     plt.errorbar(ape_radius / (fd.f_bmaj_pix / 2), ape_arr, yerr=ape_err_arr, label=LABEL_LIST[i]+'box')
    
    if IS_BOX:
        for k, box_bm in enumerate(box_bm_list):
            plt.errorbar(ape_radius / (fd.f_bmaj_pix/2), gau_ape_arr.T[k].T, label='gaufit'+'(d=%.1f beam)'%(box_bm_list[k]))
    else:
        plt.errorbar(ape_radius / (fd.f_bmaj_pix/2), gau_ape_arr, label='gaufit')
    # plt.plot(ape_radius/(pix_deg*36e2), ape_arr, label=label_list[i])
    # plt.plot(y_list*0+bmaj_deg*36e2/2, y_list)
    # plt.xlim(0, 6)
    # plt.ylim(0,250)
    plt.xlabel('Aperture diameter [beam]', fontsize=12)
    plt.ylabel('Total flux [uJy]', fontsize=12)
    plt.legend()
    plt.grid(alpha=0.5)
    plt.savefig('%s%s_gaufit_cog.pdf' % (PATH_FIGURE, input_filename))
    plt.show()
    plt.close()


if __name__ == '__main__':
    main()