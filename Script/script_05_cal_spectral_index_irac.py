#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: script_cal_spectral_index.py
Name: Chia-Lin Ko
Create Date: Jun 11, 2021
Last Modified Date: Jun 17, 2021
------------------------
This program aims to calculate the spectral index
"""
import os
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

# binning_num = int(folder.split('bin')[-1])
BEAM_FACTOR = 10
THLD_ISL = 5
THLD_PIX = 5

# PATH_DATA_LIST          = [ '../Data/COSMOS/Image/Stacking/IRAC/irac_bin%s/Pybdsf_bm%.1f_thIsl%s_thPix%s/'%(binning_num, BEAM_FACTOR, THLD_ISL, THLD_PIX),
#                             '../Data/COSMOS/Image/Stacking/IRAC/Pybdsf_bm%.1f_thIsl%s_thPix%s/'%(BEAM_FACTOR, THLD_ISL, THLD_PIX),                        
#                         ]
# PATH_DATA_LIST          = [ '../Data/COSMOS/Image/Stacking/Coord_3GHz/'#%s/Pybdsf_bm%.1f_thIsl%s_thPix%s/'%(folder, BEAM_FACTOR, THLD_ISL, THLD_PIX),
#                             # '../Data/COSMOS/Image/Stacking/3GHz/Pybdsf_bm%.1f_thIsl%s_thPix%s/'%(BEAM_FACTOR, THLD_ISL, THLD_PIX),                        
#                         ]
CATALOG_FILENAME_LIST   =   [ 
    '%scosmos_match_Ugne_comb.fits'%(PATH_CATALOG_CROSSMATCH)
                            ]
# CATALOG_FILENAME_CSV_LIST   =   [ 
#     '%scosmos_match_Ugne_radec.csv'%(PATH_CATALOG_CROSSMATCH)
#                                 ]  
IMG_FILENAME_LIST       =   [
    '%scosmos_vla_3GHz_2017image_uJy.fits'%(PATH_DATA_VLA),
    '%scosmos_vla_1d4GHz_XS_2021image_uJy.fits'%(PATH_DATA_VLA)
                            ]
IMG_RMS_FILENAME_LIST   =   [
    '%scosmos_vla_3GHz_2017rms_uJy.fits'%(PATH_DATA_VLA),
    '%scosmos_vla_1d4GHz_XS_2021rms_uJy.fits'%(PATH_DATA_VLA),
                            ]
SLT_NAME_LIST       = [ 
    '3GHz',
                       ]
SLT_LABEL_LIST = [  '3GHz (binning)', 
                    ]
IMG_LABEL_LIST          = ['3GHz', '1d4GHzXS']
IMG_RMS_LIST            = [ 2.3, 1.8]    # rms (uJy/beam)

S_3GHZ_DET_LIMIT        = 5*2.3e-3      # 5 sigma detection limit (mJy/beam), from Smolčić+17
# S_1d4GHZ_DET_LIMIT      = 4*12e-3       # 4 sigma detection limit (mJy/beam), from Schinnerer+10
S_1d4GHZ_DET_LIMIT      = 5*1.8e-3      # 4 sigma detection limit (mJy/beam), from Hiddo+21


BWS_FACTOR  = 1.
        
def main():
    
    # spcetral_index(path_data_in='../Data/COSMOS/Image/Stacking/Coord_3GHz/') 
    # spcetral_index(path_data_in='../Data/COSMOS/Image/Stacking/Coord_IRAC/')
    # spcetral_index(path_data_in='../Data/COSMOS/Image/Stacking/Det_3GHzIRAC/') 
    # spcetral_index(path_data_in='../Data/COSMOS/Image/Stacking/XS_3GHzDet/') 
    # spcetral_index(path_data_in='../Data/COSMOS/Image/Stacking/XS_IRACDet/') 
    spcetral_index(path_data_in='../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/') 
    # spcetral_index(path_data_in='../Data/COSMOS/Image/Stacking/Det_3GHz_nRadioAGN/')
    # spcetral_index(path_data_in='../Data/COSMOS/Image/Stacking/Det_IRAC_nRadioAGN/')


def spcetral_index(path_data_in='./'):
    
    folder_lst = listdir_nohidden(path_data_in)
    for folder in folder_lst:

        cat_bin_fn  = folder.split('mediam')[0]+'bin'
        binning_num = int(folder.split('gp')[-1])
        binning_gp  = int(folder.split('bin')[-1].split('gp')[0])
        csi = SpecIndex(CATALOG_FILENAME_LIST)
        catalog_dict = csi.read_fits_catalog()
        df_cat_bin  = pd.read_csv(path_data_in+folder+'/'+cat_bin_fn+'.csv')

        is_indiv_spec_ind = False
        # is_stacking = True
        
        
        for cat_name, cat_info in catalog_dict.items():
            cat_hd, cat_data = cat_info

            #####################################################################################
            # individual source
            #####################################################################################   

            color_lst = ['gray', 'm', 'orange', 'green','cyan', 'red']
            if is_indiv_spec_ind:
                #for slt_name, slt_arr in slt_dict.items():

                # slt_name = 'slt_4SNR_nAGN_both_1d4_3GHz_nomulti'
                # total flux
                S_3GHz          = slt_data_dict[slt_name]['flux_3GHz']/1e3      # uJy -> mJy
                S_3GHz_err      = slt_data_dict[slt_name]['flux_err_3GHz']/1e3  # uJy -> mJy
                S_1d4GHz        = slt_data_dict[slt_name]['Flux_corr_1d4GHz']*1e3         # mJy
                S_1d4GHz_err    = slt_data_dict[slt_name]['Eflux_corr_1d4GHz']*1e3       # mJy
                # S_1d4GHz        = slt_data_dict[slt_name]['flux_1d4GHz']        # mJy
                # S_1d4GHz_err    = slt_data_dict[slt_name]['ferr_1d4GHz']        # mJy
                            
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
            else:

                fig = plt.figure()
                plt.rcParams["font.family"] = "serif"
                ax  = fig.add_axes([0.12,0.12,0.8,0.75]) # left, bottom, width, height (range 0 to 1)
                xlabel          = r'S$_{\rm 3 GHz}$ [mJy]'
                xlim_low, xlim_up = 0   , 0.15
                ylim_low, ylim_up = -3.5, 0.5

                for j, slt_name in enumerate(SLT_NAME_LIST):   
                    # set up
                    stacking_stat = 'mediam' 
                    is_binning = True if j == 0 else False
                    # path_data       = PATH_DATA_LIST[0] if j == 0 else PATH_DATA_LIST[1]
                    path_data       = path_data_in+folder+'/' if j == 0 else path_data
                    path_data += 'Pybdsf_bm%.1f_thIsl%s_thPix%s/'%(BEAM_FACTOR, THLD_ISL, THLD_PIX)
                    if is_binning:
                        arr_num = binning_num
                    else:
                        arr_num = 1
                    img_num = len(IMG_LABEL_LIST)
                    
                    St_arr      = np.zeros((arr_num, img_num, 2))
                    Sp_arr      = np.zeros((arr_num, img_num, 2))
                    rms_arr     = np.zeros((arr_num, img_num))
                    pd_value    = df_cat_bin['PB_1d4XS'].to_numpy()

                    for l, img_label in enumerate(IMG_LABEL_LIST):
                        for k in range(arr_num):
                            if is_binning:
                                csv_fn   = '%sStacking_img%s_%s_%s.pybdsm.gaul.csv'%(path_data, img_label, folder, k)
                            else:
                                csv_fn   = '%sStacking_img%s_%s.pybdsm.gaul.csv'%(path_data, img_label, folder)
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

                                
                                #St_arr[k][l] = np.array([float(df['Total_flux'].sum()), float(df['E_Total_flux'][0])])
                                if img_label == '3GHz':
                                    Sp_arr[k][l]    = np.array([Sp, Sp_err])
                                    St_arr[k][l]    = np.array([St, St_err])
                                elif img_label == '1d4GHzXS': # BWS correction
                                    Sp_arr[k][l]    = np.array([Sp/BWS_FACTOR , Sp_err/BWS_FACTOR])#/pd_value[k]
                                    St_arr[k][l]    = np.array([St, St_err])#/pd_value[k]
                                rms_arr[k][l]   = rms
                            except FileNotFoundError:
                                St_arr[k][l]    = np.array([np.nan, np.nan])
                                Sp_arr[k][l]    = np.array([np.nan, np.nan])
                                rms_arr[k][l]   = np.nan
                                print('Not found: ', csv_fn)

                    St_gv_arr       = gv.gvar(St_arr.T[0].T,  St_arr.T[1].T)
                    Sp_gv_arr       = gv.gvar(Sp_arr.T[0].T,  Sp_arr.T[1].T)
                    print(Sp_gv_arr)
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

                    St_1d4GHz_arr     = St_arr.T[0][1]
                    St_1d4GHz_err_arr = St_arr.T[1][1]
                    Sp_1d4GHz_arr     = Sp_arr.T[0][1]
                    Sp_1d4GHz_err_arr = Sp_arr.T[1][1]
                    rms_S1d4GHz_arr   = rms_arr.T[1]
                    print('rms=',rms_S1d4GHz_arr)

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

                    St_3GHz_arr     = St_arr.T[0][0]
                    St_3GHz_err_arr = St_arr.T[1][0]
                    Sp_3GHz_arr     = Sp_arr.T[0][0]
                    Sp_3GHz_err_arr = Sp_arr.T[1][0]
                    rms_S3GHz_arr   = rms_arr.T[0]

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
                    print('St_1d4GHz_arr = ', St_1d4GHz_arr)
                    print('Sp_1d4GHz_arr = ', Sp_1d4GHz_arr)
                    print('is_resolve_1d4GHz =', is_resolve_1d4GHz)
                    print('S_1d4GHz_res_unres_arr =', S_1d4GHz_res_unres_arr)
                    print('St_3GHz_arr = ', St_3GHz_arr)
                    print('Sp_3GHz_arr = ', Sp_3GHz_arr)
                    print('is_resolve_3GHz =', is_resolve_3GHz)
                    print('S_3GHz_res_unres_arr =', S_3GHz_res_unres_arr)
                    print('alpha =', spec_ind_3_1d4_arr)

                    # save csv file
                    df_bin = pd.DataFrame({
                        "S_1d4GHz_res_unres_arr"    : S_1d4GHz_res_unres_arr,
                        "S_1d4GHz_res_unres_err_arr": S_1d4GHz_res_unres_err_arr,
                        "St_1d4GHz_arr"             : St_1d4GHz_arr,
                        "St_1d4GHz_err_arr"         : St_1d4GHz_err_arr,
                        "Sp_1d4GHz_arr"             : Sp_1d4GHz_arr,
                        "Sp_1d4GHz_err_arr"         : Sp_1d4GHz_err_arr, 
                        "rms_S1d4GHz_arr"           : rms_S1d4GHz_arr,         
                        "S_3GHz_res_unres_arr"      : S_3GHz_res_unres_arr,
                        "S_3GHz_res_unres_err_arr"  : S_3GHz_res_unres_err_arr,
                        "St_3GHz_arr"               : St_3GHz_arr,
                        "St_3GHz_err_arr"           : St_3GHz_err_arr,
                        "Sp_3GHz_arr"               : Sp_3GHz_arr,
                        "Sp_3GHz_err_arr"           : Sp_3GHz_err_arr,
                        "rms_S3GHz_arr"             : rms_S3GHz_arr,    
                        "is_resolve_1d4GHz"         : is_resolve_1d4GHz, 
                        "is_resolve_3GHz"           : is_resolve_3GHz,
                        "spec_ind_3_1d4_arr"        : spec_ind_3_1d4_arr,
                        "spec_ind_3_1d4_err_arr"    : spec_ind_3_1d4_err_arr
                        })
                    df_bin.to_csv("%sS_%s_bm%.1f_thIsl%s_thPix%s.csv"%(PATH_TABLE, folder, BEAM_FACTOR, THLD_ISL, THLD_PIX), index=True)
                    
                    
                    is_plot_res_unres = False
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
            # y_plot_arr = np.linspace(ylim_low, ylim_up, 50)   
            # ax.plot(y_plot_arr*0+slt_thld, y_plot_arr, color='gray', linestyle='-.', 
            #         label=r'3 GHz cutoff ($S_{\rm 3GHz}=%d$uJy)'%(slt_thld*1e3))
            
            ax.set_xlim(xlim_low, xlim_up)
            ax.set_ylim(ylim_low, ylim_up)
            ax.set_xlabel(xlabel, fontsize=14)
            #ax.set_ylabel(r'S$_{\rm 1.4 GHz}$', fontsize=14)
            ax.set_ylabel(r'spectral index $\alpha$', fontsize=14)
            ax.tick_params(which='major', direction='in', length=6, width=1, labelsize=12, top=False, right=True)
            ax.tick_params(which='minor', direction='in', length=4, width=1, labelsize=12, top=False, right=True)
            plt.legend(loc='lower right', fontsize=7)
            if is_binning:
                # fig_fn = '%sCOSMOS_alpha_Gao4SNR_nAGN_bin%s_bm%.1f_thIsl%s_thPix%s.pdf'%(PATH_FIGURE, binning_num, BEAM_FACTOR, THLD_ISL, THLD_PIX)
                fig_fn = '%s%s_bm%.1f_thIsl%s_thPix%s.pdf'%(PATH_FIGURE, folder, BEAM_FACTOR, THLD_ISL, THLD_PIX)
                ax.set_title(fig_fn[len(PATH_FIGURE):-4], fontsize = 11)
            else:
                fig_fn = '%sCOSMOS_alpha_Gao4SNR_nAGN_bin%sgp%s_all_bm%.1f_thIsl%s_thPix%s.pdf'%(PATH_FIGURE, binning_num, binning_gp, BEAM_FACTOR, THLD_ISL, THLD_PIX)
                ax.set_title(fig_fn[len(PATH_FIGURE):-4], fontsize = 12)
            plt.savefig(fig_fn)

        # plot the figure


#####################################################################################
# function
#####################################################################################
def listdir_nohidden(path):
    return [f for f in os.listdir(path) if not f.startswith('.')]

def listdir_csv(path):
    return [f for f in os.listdir(path) if f.endswith('.csv')]

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
        path_data       = path_data_in+folder+'/' if j == 0 else path_data_in
        path_data += 'Pybdsf_bm%.1f_thIsl%s_thPix%s/'%(BEAM_FACTOR, THLD_ISL, THLD_PIX)
        # path_data       = PATH_DATA_LIST[0] if j == 0 else PATH_DATA_LIST[1]
        if is_binning:
            arr_num = binning_num
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
                                slt_name, stacking_stat, binning_num, k)
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
        # fig_fn = '%sCOSMOS_alpha_Gao4SNR_nAGN_bin%s_bm%.1f_thIsl%s_thPix%s.pdf'%(PATH_FIGURE, binning_num, BEAM_FACTOR, THLD_ISL, THLD_PIX)
        fig_fn = '%s%s_bm%.1f_thIsl%s_thPix%s.pdf'%(PATH_FIGURE, folder, BEAM_FACTOR, THLD_ISL, THLD_PIX)
        ax.set_title(fig_fn[len(PATH_FIGURE):-4], fontsize = 14)
    else:
        fig_fn = '%sCOSMOS_alpha_Gao4SNR_nAGN_bin%sgp%s_all_bm%.1f_thIsl%s_thPix%s.pdf'%(PATH_FIGURE, binning_num, binning_gp, BEAM_FACTOR, THLD_ISL, THLD_PIX)
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

if __name__ == '__main__':
    main()