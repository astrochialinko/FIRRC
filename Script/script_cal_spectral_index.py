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
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import gvar as gv
import pandas as pd 
import matplotlib.pyplot as plt

# my own packing
from ClassSpecIndex import SpecIndex
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
PATH_REGION             = '../Data/Regions/'

BINNING_NUM = 6
BEAM_FACTOR = 10
THLD_ISL = 5
THLD_PIX = 5
PATH_DATA_LIST          = [ '../Data/COSMOS/Image/Stacking/NoMulti_bin%s/Pybdsf_bm%.1f_thIsl%s_thPix%s/'%(BINNING_NUM, BEAM_FACTOR, THLD_ISL, THLD_PIX),
                            '../Data/COSMOS/Image/Stacking/NoMulti_binAll/Pybdsf_bm%.1f_thIsl%s_thPix%s/'%(BEAM_FACTOR, THLD_ISL, THLD_PIX),                        
                        ]
CATALOG_FILENAME_LIST   =   [ 
    '%scosmos_jcmt_450um_2021cat_2020cat_3GHz_1d4GHz_irac_agnYYC_agn3GHz_deblended.fits'%(PATH_CATALOG_CROSSMATCH) 
                            ]
IMG_FILENAME_LIST       =   [
    '%scosmos_vla_3GHz_2017image_uJy.fits'%(PATH_DATA_VLA),
    #'%scosmos_vla_1d4GHz_lp_2010image_uJy.fits'%(PATH_DATA_VLA),
    '%scosmos_vla_1d4GHz_dp_2010image_uJy.fits'%(PATH_DATA_VLA)
                            ]
IMG_RMS_FILENAME_LIST   =   [
    '%scosmos_vla_3GHz_2017rms_uJy.fits'%(PATH_DATA_VLA),
    #None,
    None
                            ]
SLT_NAME_LIST       = [ 
    # 'slt_nAGN_both_1d4_3GHz_nomulti',
    # 'slt_nAGN_both_1d4_3GHz_15SNR3GHz_nomulti',
    # 'slt_nAGN_3GHzOnly_nomulti',
    'slt_4SNR_nAGN_3GHz_nomulti',
    # 'slt_4SNR_nAGN_3GHz_nomulti',
    # 'slt_4SNR_nAGN_both_1d4_3GHz_nomulti',
    # 'slt_4SNR_nAGN_both_1d4_3GHz_lt15SNR3GHz_nomulti',
    # 'slt_4SNR_nAGN_both_1d4_3GHz_15SNR3GHz_nomulti',
    # 'slt_4SNR_nAGN_3GHzOnly_nomulti',
                       ]
SLT_LABEL_LIST = [  '3 GHz detected (binning)', '3 GHz detected',
                    'both 1.4 and 3 GHz', 'both 1.4 and 3 GHz (<15 SNR)', 'both 1.4 and 3 GHz (>15 SNR)', '3 GHz only', 
                    # 'both 1.4 and 3 GHz', 'both 1.4 and 3 GHz (<15 SNR)', 'both 1.4 and 3 GHz (>15 SNR)', '3 GHz only'
                    ]
IMG_LABEL_LIST          = ['3GHz', '1d4GHz_dp']
IMG_RMS_LIST            = [ 2.3, 12]    # rms (uJy/beam)
# IMG_LABEL_LIST          = ['3GHz', '1d4GHz_lp', '1d4GHz_dp']
#IMG_RMS_LIST            = [ 2.3, 15, 12]    # rms (uJy/beam)

S_3GHZ_DET_LIMIT        = 5*2.3e-3      # 5 sigma detection limit (mJy/beam), from Smolčić+17
S_1d4GHZ_DET_LIMIT      = 4*12e-3       # 4 sigma detection limit (mJy/beam), from Schinnerer+10


BWS_FACTOR  = 0.9
        
def main():

    
    csi = SpecIndex(CATALOG_FILENAME_LIST)
    catalog_dict = csi.read_fits_catalog()

    is_indiv_spec_ind = True
    is_stacking = False
    
    
    for cat_name, cat_info in catalog_dict.items():
        cat_hd, cat_data = cat_info

        #####################################################################################
        # individual source
        #####################################################################################   

        column_name_lst = [ 'flux_3GHz', 'flux_err_3GHz', 'res_3GHz',
                            'flux_1d4GHz', 'ferr_1d4GHz', 'res_1d4GHz',
                            'ra_3GHz', 'dec_3GHz',
                            ] 

        # slt_dict:         {slt_name, key: ndarray(bool)}
        slt_dict        = set_slt_dict(cat_data)    
        # slt_data_dict:    {slt_name, key: {column_name, key: selected data} }
        slt_data_dict   = csi.get_slt_data_dict(cat_data, slt_dict, column_name_lst)

        if is_indiv_spec_ind:
            #for slt_name, slt_arr in slt_dict.items():

            slt_name = 'slt_4SNR_nAGN_both_1d4_3GHz_nomulti'
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

            res_3GHz        = slt_data_dict[slt_name]['res_3GHz']
            res_1d4GHz      = slt_data_dict[slt_name]['res_1d4GHz']
            isres_3GHz      = res_3GHz == 1
            unres_3GHz      = res_3GHz == 0
            isres_1d4GHz    = res_1d4GHz > 0
            unres_1d4GHz    = res_1d4GHz <=0

            is_res_all = isres_3GHz & isres_1d4GHz
            unres_all  = unres_3GHz & unres_1d4GHz

            # print(is_res_all.sum())
            # print(unres_all.sum())
            # print(is_res_all)
            is_plot_all = True

            # plot the figure (individual source only)
            if is_plot_all:
                x_arr           = np.array([S_3GHz], dtype=float)
                y_arr           = np.array([spec_ind_3_1d4], dtype=float)
                x_err_arr       = np.array([S_3GHz_err], dtype=float)
                y_err_arr       = np.array([spec_ind_3_1d4_err], dtype=float)
                xlabel          = r'S$_{\rm 3 GHz}$ [mJy]'
                c1_list         = ['navy']
                c2_list         = ['royalblue']
                marker_list     = ['o', 'o']
                line_list       = ['Gao(>4SNR)']
            else:
                slt_1 = is_res_all
                slt_2 = unres_all
                # [spec_ind_3_1d4[slt] for slt in slt_lst]


                x_arr           = [S_3GHz[slt_1], S_3GHz[slt_2]]
                y_arr           = [spec_ind_3_1d4[slt_1], spec_ind_3_1d4[slt_2]]
                x_err_arr       = [S_3GHz_err[slt_1], S_3GHz_err[slt_2]]
                y_err_arr       = [spec_ind_3_1d4_err[slt_1], spec_ind_3_1d4_err[slt_2]]
                xlabel          = r'S$_{\rm 3 GHz}$ [mJy]'
                c1_list         = ['navy']*2
                c2_list         = ['royalblue']*2
                marker_list     = ['o', 'o']*2
                line_list       = ['Gao(>4SNR)']*2
            xlim_low, xlim_up = 0   , 0.15
            ylim_low, ylim_up = -3.5, 0.5
            #xlim_low, xlim_up = 0   , 0.3
            #ylim_low, ylim_up = 0   , 0.4
            slt_thld        = S_3GHZ_DET_LIMIT*3
            slt_arr         = x_arr
            plot_spec_ind(  csi, x_arr, y_arr, x_err_arr, y_err_arr, 
                            xlabel, c1_list, c2_list, marker_list, line_list,
                            xlim_low, xlim_up, ylim_low, ylim_up,
                            is_thld=True, slt_thld=slt_thld, slt_arr=slt_arr)

    # plot the figure


#####################################################################################
# function
#####################################################################################

def plot_spec_ind(csi, x_arr, y_arr, x_err_arr, y_err_arr, xlabel, c1_list, c2_list, marker_list, line_list,
                  xlim_low, xlim_up, ylim_low, ylim_up, is_thld=False, slt_thld=None, slt_arr=None):

    fig = plt.figure()
    plt.rcParams["font.family"] = "serif"
    ax  = fig.add_axes([0.12,0.12,0.8,0.75]) # left, bottom, width, height (range 0 to 1)

    for i in range(len(x_arr)):
        if is_thld:

            x_all_arr       = x_arr[i]
            y_all_arr       = y_arr[i]
            x_all_err_arr   = x_err_arr[i]
            y_all_err_arr   = y_err_arr[i]

            ax.errorbar(x_all_arr, y_all_arr, xerr=x_all_err_arr, yerr=y_all_err_arr,
                    fmt=marker_list[i], color=c2_list[i], ecolor =c2_list[i], alpha=0.2,
                    markersize=8, markeredgecolor=c1_list[i], markeredgewidth=1,
                    #label='%s (all, #=%s)'%(line_list[i], len(x_all_arr))
                    )

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
                    label='%s (clipped, #=%s)'%(line_list[i], len(x_arr))
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

    #####################
    color_lst = ['gray', 'm', 'orange', 'green','cyan', 'red']
    for j, slt_name in enumerate(SLT_NAME_LIST):
        
        # set up
        stacking_stat = 'mediam' 
        is_binning = True if j == 0 else False
        path_data       = PATH_DATA_LIST[0] if j == 0 else PATH_DATA_LIST[1]
        if is_binning:
            arr_num = BINNING_NUM
        else:
            arr_num = 1
        img_num = len(IMG_LABEL_LIST)
        
        St_arr      = np.zeros((arr_num, img_num, 2))
        Sp_arr      = np.zeros((arr_num, img_num, 2))
        rms_arr     = np.zeros((arr_num, img_num))

        for l, img_label in enumerate(IMG_LABEL_LIST):

            for k in range(arr_num):
                if is_binning:
                    csv_fn   = '%sStacking_%s_%s_%s_bin%s_%s.pybdsm.gaul.csv'%(path_data, img_label, 
                                slt_name, stacking_stat, BINNING_NUM, k)
                else:
                    csv_fn   = '%sStacking_%s_%s_%s.pybdsm.gaul.csv'%(path_data, img_label, 
                                slt_name, stacking_stat)
                try:
                    df = pd.read_csv(csv_fn, skipinitialspace=True, skiprows=4)
                    ####
                    St_s = df['Total_flux']
                    slt_crit = St_s.argmax()

                    St = df['Total_flux'][slt_crit]
                    Sp = df['Peak_flux'][slt_crit]
                    St_err = df['E_Total_flux'][slt_crit]
                    Sp_err = df['E_Peak_flux'][slt_crit]
                    rms = df['Isl_rms'][slt_crit]

                    ####

                    St_arr[k][l] = np.array([St, St_err])
                    #St_arr[k][l] = np.array([float(df['Total_flux'].sum()), float(df['E_Total_flux'][0])])
                    if img_label == '3GHz':
                        Sp_arr[k][l] = np.array([Sp, Sp_err])
                    elif img_label == '1d4GHz_dp': # BWS correction
                        Sp_arr[k][l]    = np.array([Sp/BWS_FACTOR , Sp_err/BWS_FACTOR])
                    rms_arr[k][l]   = rms
                except FileNotFoundError:
                    St_arr[k][l]    = np.array([np.nan, np.nan])
                    Sp_arr[k][l]    = np.array([np.nan, np.nan])
                    rms_arr[k][l]   = np.nan
                    print('Not found: ', csv_fn)

        St_gv_arr       = gv.gvar(St_arr.T[0].T,  St_arr.T[1].T)
        Sp_gv_arr       = gv.gvar(Sp_arr.T[0].T,  Sp_arr.T[1].T)
        St_Sp_gv_arr    = np.divide(St_gv_arr,  Sp_gv_arr)
        SN_ratio_gv_arr = np.divide(Sp_gv_arr,  rms_arr)
        St_Sp_arr   , St_Sp_err_arr     = MyStatistic.gv2arr(St_Sp_gv_arr)
        SN_ratio_arr, SN_ratio_err_arr  = MyStatistic.gv2arr(SN_ratio_gv_arr)
        

        def resolve_curve_1d4GHz(SNratio):
            return np.power(0.35, -11/(SNratio)**1.45) 
        def resolve_curve_3GHz(SNratio):
            return 1+6*np.power(SNratio, -1.44) 

        is_resolve_3GHz     = (St_Sp_arr-resolve_curve_3GHz(SN_ratio_arr)).T[0].T > 0
        is_resolve_1d4GHz   = (St_Sp_arr-resolve_curve_1d4GHz(SN_ratio_arr)).T[1].T > 0

        S_3GHz_gv_arr   = St_gv_arr.T[0]
        S_1d4GHz_gv_arr = St_gv_arr.T[1]
        S_3GHz_arr, S_3GHz_err_arr      = MyStatistic.gv2arr(S_3GHz_gv_arr)
        S_1d4GHz_arr, S_1d4GHz_err_arr  = MyStatistic.gv2arr(S_1d4GHz_gv_arr)

        # 1.4 GHz
        S_1d4GHz_res_unres_arr      = np.zeros(arr_num)
        S_1d4GHz_res_unres_err_arr  = np.zeros(arr_num)
        for i in range(arr_num):
            if is_resolve_1d4GHz[i]:
                S_1d4GHz_res_unres_arr[i]       = St_arr.T[0][1][i] # S_1.4 GHz, St
                S_1d4GHz_res_unres_err_arr[i]   = St_arr.T[1][1][i] # Serr_1.4 GHz, St
            else:
                S_1d4GHz_res_unres_arr[i]       = Sp_arr.T[0][1][i] # S_1d4GHz, Sp
                S_1d4GHz_res_unres_err_arr[i]   = Sp_arr.T[1][1][i] # Serr_1d4GHz, Sp

        # 3 GHz
        S_3GHz_res_unres_arr      = np.zeros(arr_num)
        S_3GHz_res_unres_err_arr  = np.zeros(arr_num)
        for i in range(arr_num):
            if is_resolve_3GHz[i]:
                S_3GHz_res_unres_arr[i]       = St_arr.T[0][0][i]
                S_3GHz_res_unres_err_arr[i]   = St_arr.T[1][0][i]
            else:
                S_3GHz_res_unres_arr[i]       = Sp_arr.T[0][0][i]
                S_3GHz_res_unres_err_arr[i]   = Sp_arr.T[1][0][i]

        # S_1d4GHz_res_unres_arr       = Sp_arr.T[0][1]
        # S_1d4GHz_res_unres_err_arr   = Sp_arr.T[1][1]

        # print(S_3GHz_res_unres_arr)
        # print(is_resolve_3GHz)
        # print(S_3GHz_gv_arr)
        # print(St_gv_arr)
        # print(Sp_gv_arr)

        S_1d4GHz_res_unres_gv_arr = gv.gvar(S_1d4GHz_res_unres_arr, S_1d4GHz_res_unres_err_arr)
        S_3GHz_res_unres_gv_arr = gv.gvar(S_3GHz_res_unres_arr, S_3GHz_res_unres_err_arr)
        

        #spec_ind_3_1d4_gv_arr  = cal_spcetral_index(S_3GHz_gv_arr, 3, S_1d4GHz_gv_arr, 1.4)
        spec_ind_3_1d4_gv_arr  = cal_spcetral_index(S_3GHz_res_unres_gv_arr, 3, S_1d4GHz_res_unres_gv_arr, 1.4)
        spec_ind_3_1d4_arr, spec_ind_3_1d4_err_arr     = MyStatistic.gv2arr(spec_ind_3_1d4_gv_arr)


        print('-------------------')
        print(img_label, slt_name)
        # print('S_t = ', St_arr.T[0].T)
        # print('S_p = ', Sp_arr.T[0].T)
        # print('S = ', Sp_arr.T[0].T)
        # print('St_Sp_arr=', St_Sp_arr)
        # print('SN_ratio_arr=', SN_ratio_arr)
        # print('S_1d4GHz_arr = ', S_1d4GHz_arr)
        print('S_1d4GHz_res_unres_arr =', S_1d4GHz_res_unres_arr)
        # print('S_3GHz_arr = ', S_3GHz_arr)
        print('S_3GHz_res_unres_arr =', S_3GHz_res_unres_arr)
        print('alpha =', spec_ind_3_1d4_arr)

        is_plot_res_unres = True
        # slt_freq = is_resolve_3GHz
        slt_freq = is_resolve_1d4GHz

        if not is_plot_res_unres:
            xx_arr     = np.array([S_3GHz_res_unres_arr], dtype=float)/1e3
            xx_err_arr = np.array([S_3GHz_res_unres_err_arr], dtype=float)/1e3
            yy_arr     = np.array([spec_ind_3_1d4_arr], dtype=float)
            yy_err_arr = np.array([spec_ind_3_1d4_err_arr], dtype=float)
            color = [color_lst[j]]
            ecolor = ['k']
            mcolor = ['k']
            label_2 = ['binning']
        else:
            # slt_freq = is_resolve_3GHz
            # label_2 = ['res(3 GHz)', 'unres(3 GHz)']

            slt_freq = is_resolve_1d4GHz
            label_2 = ['res(1.4 GHz)', 'unres(1.4 GHz)']

            xx_arr     = np.array([S_3GHz_res_unres_arr[slt_freq], S_3GHz_res_unres_arr[~slt_freq]])/1e3
            xx_err_arr = np.array([S_3GHz_res_unres_err_arr[slt_freq], S_3GHz_res_unres_err_arr[~slt_freq]])/1e3
            yy_arr     = np.array([spec_ind_3_1d4_arr[slt_freq], spec_ind_3_1d4_arr[~slt_freq]])
            yy_err_arr = np.array([spec_ind_3_1d4_err_arr[slt_freq], spec_ind_3_1d4_err_arr[~slt_freq]])
            mcolor = ['#4F4F4F', 'darkviolet']
            ecolor = ['silver', 'violet']
            color  = ['darkgray', 'violet']

        for k in range(len(xx_arr)):
            if j == 0:
                # plot the medium
                x_line = np.linspace(xlim_low, xlim_up, 20)
                median, median_bs = bootstrap_error(arr=yy_arr[k], number=1000)
                bs_low, bs_up     = median_bs[0], median_bs[1]
                ax.plot(x_line, x_line*0 + median, 
                        linestyle='--', linewidth=2, color = color[k], alpha = 0.8, zorder=10,\
                        label=r'$\alpha_{\rm median}$(%s)=$%.2f^{+%.2f}_{%.2f}$'%(label_2[k], median, bs_up, bs_low),\
                       )
                ax.fill_between(x_line, median + bs_up, median + bs_low, color=color_lst[j], alpha=0.2)
                label_alpha = r'%s ($\alpha$=$%.2f^{+%.2f}_{%.2f}$)'%(SLT_LABEL_LIST[j], median, bs_up, bs_low)
            else:
                label_alpha = r'%s ($\alpha$=%.2f$\pm$%.2f)'%(SLT_LABEL_LIST[j], 
                            np.nanmedian(yy_arr[k]), np.nanmedian(yy_err_arr[k]))

            # plot the binned alpha
            ax.errorbar(xx_arr[k], yy_arr[k], xerr=xx_err_arr[k], yerr=yy_err_arr[k], 
                        fmt='s', color=color[k], ecolor =ecolor[k], alpha=0.8, zorder = 20,
                        markersize=9, markeredgecolor=mcolor[k], markeredgewidth=2,
                        label=label_alpha
                        )
            print(SLT_LABEL_LIST[j], np.nanmedian(yy_arr[k]), np.nanmedian(yy_err_arr[k]))

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
    #ax.set_ylabel(r'S$_{\rm 1.4 GHz}$', fontsize=14)
    ax.set_ylabel(r'spectral index $\alpha$', fontsize=14)
    ax.tick_params(which='major', direction='in', length=6, width=1, labelsize=12, top=False, right=True)
    ax.tick_params(which='minor', direction='in', length=4, width=1, labelsize=12, top=False, right=True)
    plt.legend(loc='lower right', fontsize=7)
    if is_binning:
        fig_fn = '%sCOSMOS_alpha_Gao4SNR_nAGN_bin%s_bm%.1f_thIsl%s_thPix%s.pdf'%(PATH_FIGURE, BINNING_NUM, BEAM_FACTOR, THLD_ISL, THLD_PIX)
        ax.set_title(fig_fn[len(PATH_FIGURE):-4], fontsize = 14)
    else:
        fig_fn = '%sCOSMOS_alpha_Gao4SNR_nAGN_bin%s_all_bm%.1f_thIsl%s_thPix%s.pdf'%(PATH_FIGURE, BINNING_NUM, BEAM_FACTOR, THLD_ISL, THLD_PIX)
        ax.set_title(fig_fn[len(PATH_FIGURE):-4], fontsize = 12)
    plt.savefig(fig_fn)
    # plt.show()

### temp
def temp_spec_index(S_3GHz, S_3GHz_err, S_1d4GHz, S_1d4GHz_err):
    S_3GHz_gv       = gv.gvar(S_3GHz,     S_3GHz_err  )
    S_1d4GHz_gv     = gv.gvar(S_1d4GHz,   S_1d4GHz_err)
    spec_ind_3_1d4_gv  = cal_spcetral_index(S_3GHz_gv, 3, S_1d4GHz_gv, 1.4)
    spec_ind_3_1d4      = spec_ind_3_1d4_gv.mean
    spec_ind_3_1d4_err  = spec_ind_3_1d4_gv.sdev
    #spec_ind_3_1d4, spec_ind_3_1d4_err = MyStatistic.gv2arr(spec_ind_3_1d4_gv)

    return spec_ind_3_1d4, spec_ind_3_1d4_err

def cal_spcetral_index(S1, freq1, S2, freq2):
        return np.divide(np.log(S1)-np.log(S2), np.log(freq1)-np.log(freq2))

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
                  (cat_data['RExcess']!='T') 
                   #& (~cat_data['agn_xxx'])     & (~cat_data['agn_c17b'])
    slt_nomulti = ((cat_data['multi_3GHz']==0) | (cat_data['mult_1d4GHz']==0))
    
    # flux
    slt_3GHz            = (~np.isnan(cat_data['flux_3GHz']))
    slt_3GHz_15SNR3GHz  = (cat_data['flux_3GHz']>S_3GHZ_DET_LIMIT*3*1e3)
    slt_1d4GHz          = (~np.isnan(cat_data['flux_1d4GHz']))
    
    # 4 SNR
    #slt_dict['slt_4SNR_nAGN_1d4_or_3GHz']   = slt_4SNR & slt_nAGN & (slt_3GHz   | slt_1d4GHz)
    #slt_dict['slt_4SNR_nAGN_both_1d4_3GHz'] = slt_4SNR & slt_nAGN & slt_3GHz    & slt_1d4GHz & slt_nomulti
    slt_dict['slt_4SNR_nAGN_both_1d4_3GHz_nomulti'] = slt_4SNR & slt_nAGN & slt_3GHz    & slt_1d4GHz & slt_nomulti
    #slt_dict['slt_4SNR_nAGN_both_1d4_3GHz_15SNR3GHz'] = slt_4SNR & slt_nAGN & slt_3GHz    & slt_1d4GHz & slt_3GHz_15SNR3GHz
    slt_dict['slt_4SNR_nAGN_3GHzOnly_nomulti']      = slt_4SNR & slt_nAGN & slt_3GHz    & (~slt_1d4GHz) & slt_nomulti
    #slt_dict['slt_4SNR_nAGN_nRadio']        = slt_4SNR & slt_nAGN & (~slt_3GHz) & (~slt_1d4GHz)       
    
    # 5 SNR
    #slt_dict['slt_5SNR_nAGN_both_1d4_3GHz'] = slt_5SNR & slt_nAGN & slt_3GHz    & slt_1d4GHz
    #slt_dict['slt_5SNR_nAGN_1d4_or_3GHz']   = slt_5SNR & slt_nAGN & (slt_3GHz   | slt_1d4GHz)
    #slt_dict['slt_5SNR_nAGN_3GHzOnly']      = slt_5SNR & slt_nAGN & slt_3GHz    & (~slt_1d4GHz)
    #slt_dict['slt_5SNR_nAGN_nRadio']        = slt_5SNR & slt_nAGN & (~slt_3GHz) & (~slt_1d4GHz)

    return slt_dict


if __name__ == '__main__':
    main()