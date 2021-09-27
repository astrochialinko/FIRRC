#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: ClassSpecIndex.py
Name: Chia-Lin Ko
Create Date: April 7, 2021
------------------------
This program aims to calculate the spectral index
"""
import numpy as np
import pandas as pd 
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
    def get_coord_wcs(slt_arr_dict, column_name_lst, is_savetxt=False, outname='coord_wcs.txt'):
        """
        Get the wvs coordinate from the input dict

        Input Parameter:
            slt_data_dict       [dict]          : 
                slt_data_dict: {slt_name, key: {column_name, key: selected data} }
            column_name_lst     [list of str]   : list of the column name
        Return:
            column_arr          [ndarray]       : the data 
        """
        row_num         = len(slt_arr_dict[column_name_lst[0]])
        column_num      = len(column_name_lst)
        column_arr      = np.zeros((row_num, column_num+1))
        
        for i, column_name in enumerate(column_name_lst):
            column_arr.T[i] = slt_arr_dict[column_name]
        column_arr.T[column_num] = np.array([int(i) for i in range(row_num)])

            
        if is_savetxt:
            # save the coordinate to the txt file
            column_name_lst.append('index')
            columns     = " ".join(column_name_lst)
            fmt_lst     = ['%f' for i in column_name_lst]
            fmt_str     = ' '.join(fmt_lst)
            np.savetxt(outname, column_arr, fmt=fmt_str, header=columns)


            # save csv
            csv_fn = outname.split('.txt')[0] + '.csv'
            column_arr
            df = pd.DataFrame(column_arr, columns= column_name_lst)
            df.to_csv (csv_fn, header=True)

            # individual_csv_fn = '%sFluxDensity_%s_bm_%s.csv'%(PATH_TABLE, slt_name, BEAM_FACTOR)
            # txt_arr = np.c_[range(SUB_NUM), 
            #     St_arr.T[0][1], St_arr.T[1][1], Sp_arr.T[0][1], Sp_arr.T[1][1],
            #     S_1d4GHz_res_unres_arr, S_1d4GHz_res_unres_err_arr, is_resolve_1d4GHz, rms_1d4GHz_arr, 
            #     St_Sp_1d4GHz_arr, St_Sp_1d4GHz_err_arr, SN_ratio_1d4GHz_arr, SN_ratio_1d4GHz_err_arr,
            #     St_arr.T[0][0], St_arr.T[1][0], Sp_arr.T[0][0], Sp_arr.T[1][0],
            #     S_3GHz_res_unres_arr, S_3GHz_res_unres_err_arr, is_resolve_3GHz, rms_3GHz_arr,
            #     St_Sp_3GHz_arr, St_Sp_3GHz_err_arr, SN_ratio_3GHz_arr, SN_ratio_3GHz_err_arr,
            #     spec_ind_3_1d4_arr, spec_ind_3_1d4_err_arr
            #     ]
            # columns     = ['ID', 'St_1d4GHz', 'St_err_1d4GHz', 'Sp_1d4GHz', 'Sp_err_1d4GHz', 
            #                 'S_1d4GHz', 'S_err_1d4GHz', 'res_1d4GHz', 'rms_1d4GHz',
            #                 'St_Sp_1d4GHz', 'St_Sp_1d4GHz_err', 'SN_ratio_1d4GHz', 'SN_ratio_1d4GHz_err',
            #                 'St_3GHz', 'St_err_3GHz', 'Sp_3GHz', 'Sp_err_3GHz', 
            #                 'S_3GHz', 'S_err_3GHz', 'res_3GHz', 'rms_3GHz',
            #                 'St_Sp_3GHz', 'St_Sp_3GHz_err', 'SN_ratio_3GHz', 'SN_ratio_3GHz_err',
            #                 'spec_index', 'spec_index_err']
            # df = pd.DataFrame(txt_arr, columns= columns)
            # df.to_csv (individual_csv_fn, header=True)



        return column_arr

    @staticmethod
    def make_region(input_txt_name, output_reg_name='output.reg', radius_arcs=1, reg_color='green', reg_width=1):
        """
        Generate a ds9 region file

        Input Parameter:
            input_txt_name      [str]       : Input txt filename that contain the information of RA and Dec
            output_reg_name     [str]       : Output region filename.       Default: 'output.reg'
            radius_arcs         [float]     : Radius of circle (arcsec).    Defualt: 1
            reg_color           [str]       : Color of region.              Defualt: 'green'
            reg_width           [int]       : Line width of region.         Defualt: 1
        """
        # load the input txt file (ra dec info)
        radec_list = np.loadtxt(input_txt_name, dtype = str, delimiter=' ')
        lines_num, item_num = np.shape(radec_list)

        # Open the output regoin file and write the header
        reg_file = open(output_reg_name, "w")
        reg_file.write(
            '# Region file format: DS9 version 4.1\n'
            'global color={0} dashlist=8 3 width={1} font="helvetica 10 normal roman" '
            'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 '
            'source=1\nfk5\n'.format(reg_color, reg_width))

        for i in range(lines_num):
            ra      = np.float(radec_list[i][0])    # ra        (deg)
            dec     = np.float(radec_list[i][1])    # dec       (deg)
            radius  = '%s"'%(radius_arcs)           # radius    (arcsec)
            # write the region
            reg_file.write("circle({0},{1},{2})\n".format(
                            ra, dec, radius))
        reg_file.close()



    @staticmethod
    def get_slt_array(ini_arr, slt_crit):
        return ini_arr[np.where(slt_crit)]

    @staticmethod
    def cal_spcetral_index(S1, freq1, S2, freq2):
        return np.divide(np.log(S1)-np.log(S2), np.log(freq1)-np.log(freq2))