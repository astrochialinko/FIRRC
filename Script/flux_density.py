#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: flux_density.py
Name: Chia-Lin Ko
Create Date: March 7, 2021
------------------------
This program aims to measure the
(1) peak flux (2) integrated flux density and (3) the aperture flux
"""
import numpy as np
from scipy import optimize
import gvar as gv
import astropy.io.fits as pyfits
from photutils import aperture_photometry
from photutils import EllipticalAperture


# Constant

# Global variable


class FluxDensity:

    def __init__(self, filename):
        self.fn   = filename
        self.f    = pyfits.getdata(self.fn)
        self.f_hd = pyfits.getheader(self.fn)

        # header
        self.f_pix_deg  = self.f_hd['CDELT2']            # 1 pixel length [deg]
        self.f_pix_arcs = self.deg2arcs(self.f_pix_deg)  # 1 pixel length [arcs]
        self.f_bmaj_deg = self.f_hd['BMAJ']              # major axis of beam size [deg]
        self.f_bmin_deg = self.f_hd['BMIN']              # minor axis of beam size [deg]
        self.f_bpa_deg  = self.f_hd['BPA']               # angle of the beam [deg]
        self.f_bmaj_pix = self.deg2pixel(self.f_bmaj_deg, self.f_pix_deg)  # major axis of beam size [pixel]
        self.f_bmin_pix = self.deg2pixel(self.f_bmin_deg, self.f_pix_deg)  # minor axis of beam size [pixel]
        self.f_beam_area_pix = np.divide(np.pi * self.f_bmaj_pix * self.f_bmin_pix, 4 * np.log(2))

        # check the dimension of the image
        self.f_shape = np.shape(self.f)    # the shape of the input data(array)
        self.f_dim   = len(self.f_shape)   # the dimension of the input data
        # convert the image to 2 dimension
        if self.f_dim == 4:
            self.f = self.f[0][0]
            print('Covert the 4D image to 2D')
        elif self.f_dim == 2:
            print('The input "%s" a 2D image.' % self.fn)
        else:
            print('Please check the dimension of the input fitsfile')
        self.f_shape_y, self.f_shape_x = self.f_shape

        #
        self.x_center = self.f_shape_x/2  # x axis of the center position [pixel]
        self.y_center = self.f_shape_x/2  # y axis of the center position [pixel]

        self.fwhm_x = 0
        self.fwhm_y = 0


    def get_2d_gaussian_ini(self):
        """
        This function estimate the initial guess for the 2D Gaussian fitting from the image
        :return  initial_guess: tuple, the initial guess of the 2D Gaussian fitting
        """

        # initial guess
        amp = self.f[self.f_shape_y // 2][self.f_shape_x // 2]  # amplitude of the center position [pixel]
        #self.x_center = self.f_shape_x/2                       # x axis of the center position [pixel]
        #self.y_center = self.f_shape_x/2                       # y axis of the center position [pixel]
        sigma_x = self.f_bmaj_pix                               # std in x axis [pixel]
        sigma_y = self.f_bmin_pix                               # std in y axis [pixel]
        theta = self.f_bpa_deg                                  # angle [pixel]
        initial_guess = (amp, self.x_center, self.y_center, sigma_x, sigma_y, theta)

        return initial_guess

    @staticmethod
    def func_2d_gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta):
        x, y = xy
        xo = float(xo)
        yo = float(yo)
        a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
        b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
        c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
        g = amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo)
                                  + c * ((y - yo) ** 2)))
        return g.ravel()

    def make_2d_xy(self):
        # Create x and y indices
        x = np.linspace(0, self.f_shape_x - 1, self.f_shape_x)
        y = np.linspace(0, self.f_shape_y - 1, self.f_shape_y)
        x, y = np.meshgrid(x, y)
        xy_meshgrid = (x, y)

        return xy_meshgrid

    def make_box_list(self, beam_multiply=1, ismake_txt=False, txtfilename='box.txt'):

        x_start = int( self.x_center - self.f_bmaj_pix*beam_multiply/2 )
        x_end   = int( self.x_center + self.f_bmaj_pix*beam_multiply/2 )
        y_start = int( self.y_center - self.f_bmaj_pix*beam_multiply/2 )
        y_end   = int( self.y_center + self.f_bmaj_pix*beam_multiply/2 )

        box_list = [x_start, y_start, x_end, y_end]

        if ismake_txt:
            np.savetxt(txtfilename, box_list, fmt='%d')

        return box_list

    def fit_2d_gaussian(self, isbox_clip = False, box = None):

        # Parameters: xpos, ypos, sigmaX, sigmaY, amp

        # Create x and y indices
        xy = self.make_2d_xy()
        x, y = xy

        if isbox_clip:
            ###
            pix_ra_start  = box[0] # x_start
            pix_ra_end    = box[2] # x_end
            pix_dec_start = box[1] # y_start
            pix_dec_end   = box[3] # y_end

            data_x  = x[pix_dec_start:pix_dec_end + 1].T[pix_ra_start:pix_ra_end + 1].T
            data_y  = y[pix_dec_start:pix_dec_end + 1].T[pix_ra_start:pix_ra_end + 1].T
            data_xy = data_x, data_y
            data_f  = self.f[pix_dec_start:pix_dec_end + 1].T[pix_ra_start:pix_ra_end + 1].T
        else:
            data_x  = x 
            data_y  = y
            data_xy = xy
            data_f  = self.f

        ###

        # initial guess
        initial_guess = self.get_2d_gaussian_ini()

        # 2D Gaussian fitting
        curvefit_result = optimize.curve_fit(self.func_2d_gaussian, 
                                                (data_x, data_y), data_f.flatten(),
                                                p0=initial_guess,
                                                sigma=data_f.flatten()*0 + 0.33,
                                                absolute_sigma=True, full_output=True)

        # coeff, coeff_std, chi2, dof, chi2_r
        gaussian_result = self.cal_chi2_from_curvefit(curvefit_result)
        self.data_fitted = self.func_2d_gaussian(xy, *curvefit_result[0])

        return gaussian_result, curvefit_result

    @staticmethod
    def cal_chi2_from_curvefit(curvefit_result):

        popt, pcov, infodict, errmsg, ier = curvefit_result
        coeff  = popt                               # optimal values for the coefficients
        perr   = np.sqrt(np.diag(pcov))             # standard error of coefficients
        chi2   = (infodict['fvec'] ** 2).sum()      # chi square
        dof    = len(infodict['fvec']) - len(popt)  # degree of freedom
        chi2_r = np.divide(chi2, dof)               # reduced chi square

        return coeff, perr, chi2, dof, chi2_r

    def print_fitting_result(self, gaussian_result, initial_guess):
        parm_list = ['amplitude (uJy/beam)', 'x_center  (pixel)', 'y_center  (pixel)',
                     'sigma_x', 'sigma_y', 'theta     (degree) '
                    ]

        coeff, coeff_std, chi2, dof, chi2_r = gaussian_result

        for j, coef in enumerate(coeff):
            if j in [3, 4]:
                print(parm_list[j] + '   (pixel) \t = %2g +/- %2g \t[%2g]' % (coeff[j], coeff_std[j], initial_guess[j]))
                print(parm_list[j] + '   (arcsec)\t = %2g +/- %2g \t[%2g]' % (coeff[j] * self.deg2arcs(self.f_pix_deg),
                                                                              coeff_std[j] * self.deg2arcs(self.f_pix_deg),
                                                                              initial_guess[j] * self.deg2arcs(self.f_pix_deg)))
            else:
                print('%s \t = %2g +/- %2g \t[%2g]' % (parm_list[j], coeff[j], coeff_std[j], initial_guess[j]))

        self.fwhm_x = self.get_fwhm_gaussian(gv.gvar(coeff[3], coeff_std[3]))
        self.fwhm_y = self.get_fwhm_gaussian(gv.gvar(coeff[4], coeff_std[4]))

        print('FWHM_x    (pixel) \t = %2g +/- %2g' % (self.fwhm_x.mean, self.fwhm_x.sdev))
        print('FWHM_y    (pixel) \t = %2g +/- %2g' % (self.fwhm_y.mean, self.fwhm_y.sdev))
        print('FWHM_x    (arcsec) \t = %2g +/- %2g' % (self.fwhm_x.mean * self.deg2arcs(self.f_pix_deg),
                                                       self.fwhm_x.sdev * self.deg2arcs(self.f_pix_deg)))
        print('FWHM_y    (arcsec) \t = %2g +/- %2g' % (self.fwhm_y.mean * self.deg2arcs(self.f_pix_deg),
                                                       self.fwhm_y.sdev * self.deg2arcs(self.f_pix_deg)))

        # print('vulumn = %2g'%(2*np.pi*coeff[0]*coeff[3]*pix_deg*36e2*coeff[4]*pix_deg*36e2) )
        print('chi2 \t\t\t = %2g' % (chi2))
        print('dof \t\t\t = %2g' % (dof))
        print('chi2_r \t\t\t = %2g' % (chi2_r))

    @staticmethod
    def get_fwhm_gaussian(sigma):
        return 2 * sigma * np.sqrt(2 * np.log(2))

    @staticmethod
    def deg2arcs(input_deg):
        return input_deg * 36e2

    @staticmethod
    def deg2pixel(input_deg, pix_deg):
        return np.divide(input_deg, pix_deg)


