#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: ClassSpecIndex.py
Name: Chia-Lin Ko
Create Date: April 7, 2021
------------------------
This program aims to calculate the spectral index
"""
import numpy as np
import astropy.io.fits as pyfits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import gvar as gv
import matplotlib.pyplot as plt

# my own packing
from ClassStatistic import MyStatistic

class SpecIndex:
    """
    Class for calculate the spectral index
    """
    def __init__(self, cat_name_lst):
        """
        Input Parameter:
            cat_name_lst    [list of str]: a list for the catalog filename in FITS format
        """
        self.cat_name_lst   = cat_name_lst
        self.cat_dict       = {}

    def read_fits_catalog(self):
        """
        Read the FITS format catalog and return the column and data info

        Return:
            self.cat_dict   [dict]: a dict for the column and data info 
                {catalog_name, key: [{column_name, key: column_unit}, data]}
        """
        for i, cat_name in enumerate(self.cat_name_lst):
            with pyfits.open(cat_name) as hdul:
                hdul.verify('fix')
                tb_data     = hdul[1].data
                tb_columns  = hdul[1].columns
                tb_hd = {}
                for i, name in enumerate(tb_columns.names):
                    tb_hd[name] = tb_columns.units[i]
                self.cat_dict[cat_name] = [tb_hd, tb_data]
        return self.cat_dict

    def get_slt_data_dict(self, cat_data, slt_dict, column_name_lst):
        """
        get the selcted data based on the input column_namd, catalog, and selected criteria

        Input Parameter:
            cat_data            [fitsrec]   : catalog data
            slt_dict            [dict]      : selected dictionary
            column_name_lst     [list]      : column name
        Return:
            slt_data_dict       [dict]      : selected data dictionary
                {slt_name, key: {column_name, key: selected data} }
        """
        slt_data_dict = {}
        for slt_name, slt_arr in slt_dict.items():
            slt_data_dict[slt_name] = {}
            for column_name in column_name_lst:
                slt_data_dict[slt_name][column_name]= self.get_catalog_array(cat_data, column_name, slt_crit=slt_arr)
        return slt_data_dict

    def get_catalog_array(self, catalog, column_name, slt_crit=None):
        """
        Return the column data (ndarray) based on the input selection criteria

        Input Parameter:
            catalog     [fitsrec]           : catalog
            column_name [str]               : column name
            slt_crit    [ndarray of bool]   : selection criteria

        Return:
            column_arr  [ndarray]    : the selected column data
        """
        if slt_crit is None:
            column_arr = catalog[column_name]
        else:
            column_arr = self.get_slt_array(catalog[column_name], slt_crit)

        return column_arr

    @staticmethod
    def get_slt_array(ini_arr, slt_crit):
        return ini_arr[np.where(slt_crit)]

    @staticmethod
    def cal_spcetral_index(S1, freq1, S2, freq2):
        return np.divide(np.log(S1)-np.log(S2), np.log(freq1)-np.log(freq2))