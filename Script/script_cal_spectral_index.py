#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: script_cal_spectral_index.py
Name: Chia-Lin Ko
Create Date: March 7, 2021
------------------------
This program aims to calculate the spectral index
"""
import numpy as np
import astropy.io.fits as pyfits
import gvar as gv
import matplotlib.pyplot as plt

# my own packing
from ClassStacking import Stacking
from ClassStatistic import MyStatistic

# Constant
# path
PATH_DATA               = '../Data/COSMOS/Image/Stacking/'
PATH_DATA_VLA           = '../Data/COSMOS/Image/VLA/'
PATH_DATA_STACKING      = '../Data/COSMOS/Image/Stacking/'
PATH_CATALOG_CROSSMATCH = '../Data/COSMOS/Catalog/CrossMatch/'
PATH_TABLE              = '../Data/Tables/'
PATH_IMFIT              = '../Data/Log_imfit/'
PATH_FIGURE             = '../Figures/'

CATALOG_FILENAME_LIST   =   [ 
    '%scosmos_jcmt_450um_2021cat_2020cat_3GHz_1d4GHz_irac_agnYYC_agn3GHz_deblended.fits'%(PATH_CATALOG_CROSSMATCH) 
                            ]
IMG_FILENAME_LIST       =   [
    '%scosmos_vla_3GHz_2017image_uJy.fits'%(PATH_DATA_VLA),
    '%scosmos_vla_1d4GHz_lp_2010image_uJy.fits'%(PATH_DATA_VLA),
    '%scosmos_vla_1d4GHz_dp_2010image_uJy.fits'%(PATH_DATA_VLA)
                            ]

S_3GHZ_DET_LIMIT        = 5*2.3e-3      # 5 sigma detection limit (mJy/beam), from Smolčić+17
S_1d4GHZ_DET_LIMIT      = 4*12e-3       # 4 sigma detection limit (mJy/beam), from Schinnerer+10

class ClassSpecIndex:
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


        
def main():

    
    csi = ClassSpecIndex(CATALOG_FILENAME_LIST)
    catalog_dict = csi.read_fits_catalog()

    #####################################################################################
    # individual source
    #####################################################################################   
    for cat_name, cat_info in catalog_dict.items():
        cat_hd, cat_data = cat_info

        column_name_lst = ['flux_3GHz', 'flux_err_3GHz', 'flux_1d4GHz', 'ferr_1d4GHz', 'ra_3GHz', 'dec_3GHz'] 

        # slt_dict:         {slt_name, key: ndarray(bool)}
        slt_dict        = set_slt_dict(cat_data)    
        # slt_data_dict:    {slt_name, key: {column_name, key: selected data} }
        slt_data_dict   = csi.get_slt_data_dict(cat_data, slt_dict, column_name_lst)


        for slt_name, slt_arr in slt_dict.items():

            # total flux
            S_3GHz          = slt_data_dict[slt_name]['flux_3GHz']/1e3      # uJy -> mJy
            S_3GHz_err      = slt_data_dict[slt_name]['flux_err_3GHz']/1e3  # uJy -> mJy
            S_1d4GHz        = slt_data_dict[slt_name]['flux_1d4GHz']        # mJy
            S_1d4GHz_err    = slt_data_dict[slt_name]['ferr_1d4GHz']        # mJy
            S_3GHz_gv       = gv.gvar(S_3GHz,     S_3GHz_err  )
            S_1d4GHz_gv     = gv.gvar(S_1d4GHz,   S_1d4GHz_err)

            # spectral index
            spec_ind_3_1d4_gv  = csi.cal_spcetral_index(S_3GHz_gv, 3, S_1d4GHz_gv, 1.4)
            spec_ind_3_1d4, spec_ind_3_1d4_err = MyStatistic.gv2arr(spec_ind_3_1d4_gv)

            # plot the figure (individual source only)
            x_arr       = np.array([S_3GHz], dtype=float)
            y_arr       = np.array([spec_ind_3_1d4], dtype=float)
            x_err_arr   = np.array([S_3GHz_err], dtype=float)
            y_err_arr   = np.array([spec_ind_3_1d4_err], dtype=float)
            xlabel      = r'S$_{\rm 3 GHz}$ [mJy]'
            c1_list     = ['navy']
            c2_list     = ['royalblue']
            marker_list = ['o']
            line_list   = ['Gao(>4SNR)']
            xlim_low, xlim_up = 0, 0.15
            ylim_low, ylim_up = -3.5, 0.5
            slt_thld    = S_3GHZ_DET_LIMIT*3
            slt_arr     = x_arr
            # plot_spec_ind(  csi, x_arr, y_arr, x_err_arr, y_err_arr, 
            #                 xlabel, c1_list, c2_list, marker_list, line_list,
            #                 xlim_low, xlim_up, ylim_low, ylim_up,
            #                 is_thld=True, slt_thld=slt_thld, slt_arr=slt_arr)
    

    #####################################################################################
    # stacking
    #####################################################################################
    
    # get ra and dec
    coord_name_lst  = ['ra_3GHz', 'dec_3GHz']
    coord_arr       = get_coord_wcs(slt_data_dict, coord_name_lst, is_savetxt=True)
 
    # create the region file (optional)


    # transfer the wcs to pixel coordinate
    for i, img_filename in enumerate(IMG_FILENAME_LIST):
        csk = Stacking(img_filename)
    
    # start stacking



    # plot the figure


#####################################################################################
# function
#####################################################################################
def get_coord_wcs(slt_data_dict, column_name_lst, is_savetxt=False):
    """
    Input Parameter:
        slt_data_dict       [dict]          : 
            slt_data_dict: {slt_name, key: {column_name, key: selected data} }
        column_name_lst     [list of str]   : list of the column name
    Return:
        column_arr          [ndarray]       : the data 
    """
    for slt_name, slt_arr in slt_data_dict.items():
        
        row_num         = len(slt_arr[column_name_lst[0]])
        column_num      = len(column_name_lst)
        column_arr      = np.zeros((row_num,column_num))
        
        for i, column_name in enumerate(column_name_lst):
            column_arr.T[i] = slt_arr[column_name]
            
        if is_savetxt:
            # save the coordinate to the txt file
            filename    = '%sCoord_wcs_%s.txt'%(PATH_TABLE, slt_name) 
            columns     = ", ".join(column_name_lst)
            np.savetxt(filename, column_arr, fmt='%f %f', header=columns)

        return column_arr

def plot_spec_ind(csi, x_arr, y_arr, x_err_arr, y_err_arr, xlabel, c1_list, c2_list, marker_list, line_list,
                  xlim_low, xlim_up, ylim_low, ylim_up, is_thld=False, slt_thld=None, slt_arr=None):

    fig = plt.figure()
    plt.rcParams["font.family"] = "serif"
    ax  = fig.add_axes([0.12,0.12,0.8,0.75]) # left, bottom, width, height (range 0 to 1)

    for i in range(len(x_arr)):
        if is_thld:
            slt_crit    = (slt_arr[i]>slt_thld) & (~np.isnan(slt_arr[i]))
            x_arr       = csi.get_slt_array(x_arr[i], slt_crit)
            y_arr       = csi.get_slt_array(y_arr[i], slt_crit)
            x_err_arr   = csi.get_slt_array(x_err_arr[i], slt_crit)
            y_err_arr   = csi.get_slt_array(y_err_arr[i], slt_crit)
        else:
            x_arr       = x_arr[i]
            y_arr       = y_arr[i]
            x_err_arr   = x_err_arr[i]
            y_err_arr   = y_err_arr[i]

        ax.errorbar(x_arr, y_arr, xerr=x_err_arr, yerr=y_err_arr,
                    fmt=marker_list[i], color=c2_list[i], ecolor =c2_list[i], alpha=0.6,
                    markersize=8, markeredgecolor=c1_list[i], markeredgewidth=1,
                    label='%s'%(line_list[i])
                    )

        # plot the medium
        x_line = np.linspace(xlim_low, xlim_up, 20)
        median, median_bs = bootstrap_error(arr=y_arr, number=1000)
        bs_low, bs_up     = median_bs[0], median_bs[1]
        ax.plot(x_line, x_line*0 + median, 
                linestyle='--', linewidth=2, color = c1_list[i], alpha = 0.8, zorder=10,\
                label=r'$\alpha_{\rm median}$=$%.2f^{+%.2f}_{%.2f}$'%(median, bs_up, bs_low),\
               )
        ax.fill_between(x_line, median + bs_up, median + bs_low, color=c2_list[i], alpha=0.2)

    # plot the 1.4 GHz detection limit
    S_3GHz_arr  = np.linspace(xlim_low, xlim_up, 100)
    alpha_lim   = csi.cal_spcetral_index(S_3GHz_arr, 3, S_1d4GHZ_DET_LIMIT, 1.4)
    ax.plot(S_3GHz_arr, alpha_lim, color='tab:red', 
            label=r'1.4 GHz detection limit ($S_{\rm 1.4GHz}=%d$uJy)'%(S_1d4GHZ_DET_LIMIT*1e3))

    # plot the 3 GHz threshold
    y_plot_arr = np.linspace(ylim_low, ylim_up, 50)   
    ax.plot(y_plot_arr*0+S_3GHZ_DET_LIMIT, y_plot_arr, color='silver', linestyle='-.', 
            label=r'3 GHz detection limit ($S_{\rm 3GHz}=%d$uJy)'%(S_3GHZ_DET_LIMIT*1e3))

    # plot the 3 GHz cutoff
    y_plot_arr = np.linspace(ylim_low, ylim_up, 50)   
    ax.plot(y_plot_arr*0+slt_thld, y_plot_arr, color='gray', linestyle='-.', 
            label=r'3 GHz cutoff ($S_{\rm 3GHz}=%d$uJy)'%(slt_thld*1e3))


    ax.set_xlim(xlim_low, xlim_up)
    ax.set_ylim(ylim_low, ylim_up)
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(r'spectral index $\alpha$', fontsize=14)
    ax.tick_params(which='major', direction='in', length=6, width=1, labelsize=12, top=False, right=True)
    ax.tick_params(which='minor', direction='in', length=4, width=1, labelsize=12, top=False, right=True)
    plt.legend(loc='lower right', fontsize=9)
    plt.show()

### temp
def mybootstrap(data, num_samples, statistic, alpha):
    """Returns bootstrap estimate of 100.0*(1-alpha) CI (Confidence Interval) for statistic."""
    n = len(data)
    np.random.seed()
    idx = np.random.randint(0, n, (num_samples, n))
    samples = data[idx]
    stat = np.sort(statistic(samples, 1))
    return (stat[int((alpha/2.0)*num_samples)], stat[int((1.-alpha/2.0)*num_samples)])

def bootstrap_error(arr, number=1000):
    return np.nanmedian(arr), mybootstrap( arr, number, np.nanmedian, 1.-0.683) - np.nanmedian(arr)
### temp

def set_slt_dict(cat_data):
    """
    Set the selction criteria

    Input Parameter:
        cat_data    [fitsrec]   : catalog data
    Return:
        slt_dict    [dict]      : selection criteria dictionary
            {slt_name, key: ndarray(bool)}
    """
    slt_dict = {}

    # selection criteria
    slt_4SNR    = (cat_data['SNR_450']>4)
    slt_5SNR    = (cat_data['SNR_450']>5)
    slt_nAGN    = (cat_data['XrayAGN']!='T') & (cat_data['MIRAGN']!='T')  & (cat_data['SEDAGN']!='T') & \
                  (cat_data['RExcess']!='T') & (~cat_data['agn_xxx'])     & (~cat_data['agn_c17b'])
    # flux
    slt_3GHz    = (~np.isnan(cat_data['flux_3GHz']))
    slt_1d4GHz  = (~np.isnan(cat_data['flux_1d4GHz']))
    
    # 4 SNR
    #slt_dict['slt_4SNR_nAGN_1d4_or_3GHz']   = slt_4SNR & slt_nAGN & (slt_3GHz   | slt_1d4GHz)
    slt_dict['slt_4SNR_nAGN_both_1d4_3GHz'] = slt_4SNR & slt_nAGN & slt_3GHz    & slt_1d4GHz
    #slt_dict['slt_4SNR_nAGN_3GHzOnly']      = slt_4SNR & slt_nAGN & slt_3GHz    & (~slt_1d4GHz)
    #slt_dict['slt_4SNR_nAGN_nRadio']        = slt_4SNR & slt_nAGN & (~slt_3GHz) & (~slt_1d4GHz)       
    
    # 5 SNR
    # slt_dict['slt_5SNR_nAGN_both_1d4_3GHz'] = slt_5SNR & slt_nAGN & slt_3GHz    & slt_1d4GHz
    # slt_dict['slt_5SNR_nAGN_1d4_or_3GHz']   = slt_5SNR & slt_nAGN & (slt_3GHz   | slt_1d4GHz)
    # slt_dict['slt_5SNR_nAGN_3GHzOnly']      = slt_5SNR & slt_nAGN & slt_3GHz    & (~slt_1d4GHz)
    #slt_dict['slt_5SNR_nAGN_nRadio']        = slt_5SNR & slt_nAGN & (~slt_3GHz) & (~slt_1d4GHz)

    return slt_dict


if __name__ == '__main__':
    main()