#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: ClassStatistic.py
Name: Chia-Lin Ko
Create Date: March 29, 2021
------------------------
This program aims to measure the
(1) peak flux (2) integrated flux density and (3) the aperture flux
"""
import numpy as np
from scipy.optimize import curve_fit
import gvar as gv
from astropy.stats import bootstrap


class MyStatistic:
    def __init__(self):
        self.n = ''

    def mybootstrap(self, data, num_samples, statistic, alpha):
        """Returns bootstrap estimate of 100.0*(1-alpha) CI (Confidence Interval) for statistic."""
        n = len(data)
        np.random.seed()
        idx = np.random.randint(0, n, (num_samples, n))
        samples = data[idx]
        stat = np.sort(statistic(samples, 1))
        return (stat[int((alpha/2.0)*num_samples)], stat[int((1.-alpha/2.0)*num_samples)])

    def bootstrap_error(self, arr, number=1000):
        return np.nanmedian(arr), self.mybootstrap( arr, number, np.nanmedian, 1.-0.683) - np.nanmedian(arr)

    @staticmethod
    def gv2arr(gv_arr):   
        if gv_arr.ndim == 1:
            gv_value_arr = np.array([gv_arr[i].mean for i in range(len(gv_arr))])
            gv_err_arr   = np.array([gv_arr[i].sdev for i in range(len(gv_arr))])
        elif gv_arr.ndim == 2:
            y_ax, x_ax = np.shape(gv_arr)
            gv_value_arr = np.zeros((y_ax, x_ax))
            gv_err_arr = np.zeros((y_ax, x_ax))
            for i in range(y_ax):
                for j in range(x_ax):
                    gv_value_arr[i][j] = gv_arr[i][j].mean
                    gv_err_arr[i][j] = gv_arr[i][j].sdev
        else:
            print('Please check the dimension of your input array')
        return gv_value_arr, gv_err_arr

    @staticmethod
    def nanaverage(arr, weights, axis):
        return np.divide( np.nansum(arr*weights, axis=axis), ((~np.isnan(arr))*weights).sum(axis=axis) )

    @staticmethod
    def binning(x_arr, y_arr, bin_num, stat='median'):
        #
        x_add_arr       = np.empty(bin_num - len(x_arr)%bin_num)
        x_add_arr[:]    = np.nan
        x_arr           = np.concatenate((x_arr, x_add_arr), axis=None)
        
        y_add_arr       = np.empty(bin_num - len(y_arr)%bin_num)
        y_add_arr[:]    = np.nan
        y_arr           = np.concatenate((y_arr, y_add_arr), axis=None)

        
        xy_2d_arr       = np.vstack((x_arr, y_arr)).T
        x_ord_y_2d_arr  = xy_2d_arr[xy_2d_arr[:,0].argsort()]
        x_ord_x         = x_ord_y_2d_arr.T[0]
        y_ord_x         = x_ord_y_2d_arr.T[1] 
        x_ord_2darr     = x_ord_x[:(len(x_ord_x)//bin_num)*bin_num].reshape(-1,bin_num)
        y_ord_2darr     = y_ord_x[:(len(y_ord_x)//bin_num)*bin_num].reshape(-1,bin_num)

        # mean binning
        if stat=='mean':
            xbin_m   = np.nanmean(x_ord_2darr, axis=1)
            xbin_std = np.nanstd(x_ord_2darr, axis=1)
            ybin_m   = np.nanmean(y_ord_2darr, axis=1)
            ybin_std = np.nanstd(y_ord_2darr, axis=1)

            return xbin_m, xbin_std, ybin_m, ybin_std

        elif stat=='median':
            q1_x = np.nanpercentile(x_ord_2darr, 25, interpolation='midpoint', axis=1)
            q2_x = np.nanpercentile(x_ord_2darr, 50, interpolation='midpoint', axis=1)
            q3_x = np.nanpercentile(x_ord_2darr, 75, interpolation='midpoint', axis=1)

            q1_y = np.nanpercentile(y_ord_2darr, 25, interpolation='midpoint', axis=1)
            q2_y = np.nanpercentile(y_ord_2darr, 50, interpolation='midpoint', axis=1)
            q3_y = np.nanpercentile(y_ord_2darr, 75, interpolation='midpoint', axis=1)

            xbin_mid = q2_x
            xbin_iqr = [q2_x-q1_x, q3_x-q2_x]
            ybin_mid = q2_y
            ybin_iqr = [q2_y-q1_y, q3_y-q2_y]
            
            return xbin_mid, xbin_iqr, ybin_mid, ybin_iqr

        else:
            print("Please select the statistic method from mean or medain")





