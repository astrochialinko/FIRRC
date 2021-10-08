#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: script_plot_stacking_result.py
Name: Chia-Lin Ko
Create Date: Jun 15, 2021
Last Modified Date: Jun 15, 2021
------------------------
This program aims to plot the results of stacking
 
"""
import os
import numpy as np
import pandas as pd 
import gvar as gv
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import FlatLambdaCDM

# my own packing
from ClassStatistic import MyStatistic

# Constant
# path
PATH_DATA               = '../Data/'
PATH_DATA_VLA           = '../Data/COSMOS/Image/VLA/'
# PATH_DATA_STACKING      = '../Data/COSMOS/Image/Stacking/Coord_IRAC/'
PATH_CATALOG_CROSSMATCH = '../Data/COSMOS/Catalog/CrossMatch/'
PATH_TABLE              = '../Data/Tables/'
PATH_FIGURE             = '../Figures/'

L_SUN = 3.9e33 # (cgs)
BEAM_FACTOR = 10
THLD_ISL = 5
THLD_PIX = 5
S_3GHZ_DET_LIMIT        = 5*2.3e-3      # 5 sigma detection limit (mJy/beam), from Smolčić+17
S_1d4GHZ_DET_LIMIT      = 5*1.8e-3       # 4 sigma detection limit (mJy/beam), from Schinnerer+10


def main():

    # plot_stacking_result(path_data_in='../Data/COSMOS/Image/Stacking/Coord_3GHz/')
    # plot_stacking_result(path_data_in='../Data/COSMOS/Image/Stacking/Coord_IRAC/')
    # plot_stacking_result(path_data_in='../Data/COSMOS/Image/Stacking/Det_3GHzIRAC/')
    # plot_stacking_result(path_data_in='../Data/COSMOS/Image/Stacking/XS_3GHzDet/')
    # plot_stacking_result(path_data_in='../Data/COSMOS/Image/Stacking/XS_IRACDet/')
    plot_stacking_result(path_data_in='../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/')
    # plot_stacking_result(path_data_in='../Data/COSMOS/Image/Stacking/Det_3GHz_nRadioAGN/') 
    # plot_stacking_result(path_data_in='../Data/COSMOS/Image/Stacking/Det_IRAC_nRadioAGN/') 

def plot_stacking_result(path_data_in):
    
    slt_label_lst = listdir_nohidden(path_data_in)

    for slt_label in slt_label_lst:

        binning_num = int(slt_label.split('gp')[-1])
        binning_gp  = int(slt_label.split('bin')[-1].split('gp')[0])

        fn_flux = "%sS_%s_bm%.1f_thIsl%s_thPix%s.csv"%(PATH_TABLE, slt_label, BEAM_FACTOR, THLD_ISL, THLD_PIX)
        df_flux = pd.read_csv(fn_flux) 

        fn_bin = '%s%s/%s_bin.csv'%(path_data_in, slt_label, '_'.join(slt_label.split('_')[0:-2])) 
        df_bin = pd.read_csv(fn_bin) 
        fn_ind = '%s%s/%s.csv'%(path_data_in, slt_label, '_'.join(slt_label.split('_')[0:-2]))
        df_ind = pd.read_csv(fn_ind)

        df_m = df_flux.merge(df_bin, left_on='Unnamed: 0', right_on='Unnamed: 0', 
            how='inner', suffixes=['_1', '_2'])
        # if 'sortReshift' in slt_label:
        # if 'sortRedshift' in slt_label:
        if 'sortredshift' in slt_label:
            fig_fn = '%sCosmos_alpha_%s_bm%.1f_thIsl%s_thPix%s.pdf'%(PATH_FIGURE, slt_label, BEAM_FACTOR, THLD_ISL, THLD_PIX)
            # plot_spectral_index(df_m, df_ind, fig_fn)
            fig_fn = '%sCosmos_res_unres_%s_bm%.1f_thIsl%s_thPix%s'%(PATH_FIGURE, slt_label, BEAM_FACTOR, THLD_ISL, THLD_PIX)
            # plot_res_unres(df_m, df_ind, fig_fn)
            fig_fn = '%sCosmos_qIR_%s_bm%.1f_thIsl%s_thPix%s.pdf'%(PATH_FIGURE, slt_label, BEAM_FACTOR, THLD_ISL, THLD_PIX)
            plot_qIR(df_m, df_ind, fig_fn, fn_flux)


def plot_res_unres(df, df_ind, fig_fn):
    img_lst = ['3GHz', '1d4GHz']

    #-----
    # binning
    # is_res_arr          = df['is_resolve_1d4GHz']
    # S_3GHz_bin          = df['S_3GHz_res_unres_arr']/1e3
    # S_3GHz_bin_err      = df['S_3GHz_res_unres_err_arr']/1e3
    # rms_3GHz_bin_err    = df['rms_3GHz_arr']/1e3
    

    # S_1d4GHz_bin        = df['S_1d4GHz_res_unres_arr']
    # S_1d4GHz_bin_err    = df['S_1d4GHz_res_unres_err_arr']
    # rms_1d4GHz_bin_err  = df['rms_1d4GHz_arr']

    fmt = ['s'] *2
    color = ['gray']*2
    ecolor = ['k'] *2 
    mcolor = ['k'] *2
    alpha = [0.8]*2

    for i, img_str in enumerate(img_lst):
        fig = plt.figure()
        plt.rcParams["font.family"] = "serif"
        ax  = fig.add_axes([0.12,0.12,0.8,0.75]) # left, bottom, width, height (range 0 to 1)

        if img_str == '3GHz':
            St_gv_arr       = gv.gvar(df['St_3GHz_arr'].to_numpy(),  df['St_3GHz_err_arr'].to_numpy())/1e3
            Sp_gv_arr       = gv.gvar(df['Sp_3GHz_arr'].to_numpy(),  df['Sp_3GHz_err_arr'].to_numpy())/1e3
            rms_arr         = df['rms_S3GHz_arr'].to_numpy()/1e3

        elif img_str == '1d4GHz':
            St_gv_arr       = gv.gvar(df['St_1d4GHz_arr'].to_numpy(),  df['St_1d4GHz_err_arr'].to_numpy())
            Sp_gv_arr       = gv.gvar(df['Sp_1d4GHz_arr'].to_numpy(),  df['Sp_1d4GHz_err_arr'].to_numpy())
            rms_arr          = df['rms_S1d4GHz_arr'].to_numpy()


        St_Sp_gv_arr    = np.divide(St_gv_arr,  Sp_gv_arr)
        SN_ratio_gv_arr = np.divide(Sp_gv_arr,  rms_arr)     

        St_Sp_arr   , St_Sp_err_arr     = MyStatistic.gv2arr(St_Sp_gv_arr)
        SN_ratio_arr, SN_ratio_err_arr  = MyStatistic.gv2arr(SN_ratio_gv_arr)

        ax.errorbar(SN_ratio_arr, St_Sp_arr, xerr=SN_ratio_err_arr, yerr=St_Sp_err_arr, 
                     fmt=fmt[i], color=color[i], ecolor =ecolor[i], alpha=alpha[i], markersize=9, markeredgecolor=mcolor[i], 
                     markeredgewidth=2, label='%s'%(img_str)
                     )

        # line
        SNratio = np.logspace(0.3,3)
        St_Sp_1d4GHz     = np.power(0.35, -11/(SNratio)**1.45) 
        St_Sp_1d4GHz_mir = 2-St_Sp_1d4GHz
        St_Sp_3GHz       = 1+6*np.power(SNratio, -1.44) 
        St_Sp_3GHz_mir   = 2-St_Sp_3GHz

        if img_str == '3GHz': # 3 GHz
            ax.plot(SNratio, St_Sp_3GHz, c='C1', label = '3 GHz')
            ax.plot(SNratio, St_Sp_3GHz_mir, c='C1')
            ax.set_ylabel(r'$S_{\rm total}$/ $S_{\rm peak}$', fontsize = 14)
        elif img_str == '1d4GHz': # 1.4 GHz
            ax.plot(SNratio, St_Sp_1d4GHz, c='C0', label = '1.4 GHz')
            ax.plot(SNratio, St_Sp_1d4GHz_mir, c='C0') 
            ax.set_ylabel(r'$S_{\rm total}$/ $S_{\rm peak, corr}$', fontsize = 14)

        ax.set_xlim(2, 1000)
        ax.set_ylim(0.5, 10)
        ax.set_xlabel('S/N', fontsize = 14)
        
        plt.loglog()
        ax.tick_params(which='major', direction='in', length=6, width=1, labelsize=12, top=False, right=True)
        ax.tick_params(which='minor', direction='in', length=4, width=1, labelsize=12, top=False, right=True)
        plt.legend(loc='upper right', fontsize=7)
        ax.set_title(fig_fn[len(PATH_FIGURE):-4], fontsize = 11)
        plt.savefig(fig_fn+'_'+img_str+'.pdf')
        plt.close()


def plot_spectral_index(df, df_ind, fig_fn):

    fig = plt.figure()
    plt.rcParams["font.family"] = "serif"
    ax  = fig.add_axes([0.12,0.12,0.8,0.75]) # left, bottom, width, height (range 0 to 1)

    #-----
    # individual source
    S_3GHz        = df_ind['flux_3GHz'].to_numpy()/1e3      # uJy -> mJy
    S_3GHz_err    = df_ind['flux_err_3GHz'].to_numpy()/1e3  # uJy -> mJy
    S_1d4GHz      = df_ind['flux_1d4GHz'].to_numpy()      # mJy  # flux_1d4GHz
    S_1d4GHz_err  = df_ind['ferr_1d4GHz'].to_numpy()   # mJy  # ferr_1d4GHz
    # S_1d4GHz      = df_ind['VLA_20cm_flux'].to_numpy()      # mJy  # flux_1d4GHz
    # S_1d4GHz_err  = df_ind['VLA_20cm_fluxerr'].to_numpy()   # mJy  # ferr_1d4GHz

    S_3GHz_gv       = gv.gvar(S_3GHz,     S_3GHz_err  )
    S_1d4GHz_gv     = gv.gvar(S_1d4GHz,   S_1d4GHz_err)
    # spectral index
    spec_ind_3_1d4_gv  = cal_spcetral_index(S_3GHz_gv, 3, S_1d4GHz_gv, 1.4)
    spec_ind_3_1d4, spec_ind_3_1d4_err = MyStatistic.gv2arr(spec_ind_3_1d4_gv)

    #-----
    # binning
    is_res_arr = df['is_resolve_1d4GHz']
    S_3GHz_bin          = df['S_3GHz_res_unres_arr']/1e3
    S_3GHz_bin_err      = df['S_3GHz_res_unres_err_arr']/1e3
    S_1d4GHz_bin        = df['S_1d4GHz_res_unres_arr']
    S_1d4GHz_bin_err    = df['S_1d4GHz_res_unres_err_arr']

    if 'sort3GHz' in fig_fn:
        x_arr           = S_3GHz_bin
        x_err_arr       = S_3GHz_bin_err
        x_ind_arr       = S_3GHz
        x_ind_err_arr   = S_3GHz_err
        xlim_low, xlim_up = 0   , 0.15
        ax.set_xlabel(r'S$_{\rm 3 GHz}$ [mJy]', fontsize=14)
    elif 'sortIRAC' in fig_fn:
        x_arr           = S_3GHz_bin
        x_err_arr       = S_3GHz_bin_err
        x_ind_arr       = S_3GHz
        x_ind_err_arr   = S_3GHz_err
        xlim_low, xlim_up = 0   , 0.15
        ax.set_xlabel(r'S$_{\rm 3 GHz}$ [mJy]', fontsize=14)
    elif ('sortRedshift' in fig_fn ) or ('sortReshift' in fig_fn):
        if 'Lim' in fig_fn:
            x_arr           = df['redshift']
            x_err_arr       = df['redshift_err']
            x_ind_arr       = df_ind['redshift']
            x_ind_err_arr   = df_ind['redshift_err']
        elif 'Ugne' in fig_fn:
            x_arr           = df['z_m_Ugne']
            x_err_arr       = 0.5*(df['z_16_Ugne']+df['z_84_Ugne']- 2*x_arr)
            x_ind_arr       = df_ind['z_m_Ugne']
            x_ind_err_arr   = 0.5*(df_ind['z_16_Ugne']+df_ind['z_84_Ugne']- 2*x_ind_arr)
        xlim_low, xlim_up = 0   , 5
        ax.set_xlabel('redshift', fontsize=14)
    elif 'sortLIR' in fig_fn:
        if 'Lim' in fig_fn:
            x_arr           = df['logLIR']
            x_err_arr       = df['logLIR_err']
            x_ind_arr       = df_ind['logLIR']
            x_ind_err_arr   = df_ind['logLIR_err']
        elif 'Ugne' in fig_fn:
            x_arr           = df['Ldust_m_Ugne']
            x_err_arr       = 0.5*(df['Ldust_16_Ugne']+df['Ldust_84_Ugne']- 2*x_arr)
            x_ind_arr       = df_ind['Ldust_m_Ugne']
            x_ind_err_arr   = 0.5*(df_ind['Ldust_16_Ugne']+df_ind['Ldust_84_Ugne']- 2*x_ind_arr)
        xlim_low, xlim_up = 10   , 14

    y_arr       = df['spec_ind_3_1d4_arr']
    y_err_arr   = df['spec_ind_3_1d4_err_arr']

    mcolor = ['#4F4F4F', 'darkviolet', 'k']
    ecolor = ['silver', 'violet', 'k']
    color  = ['darkgray', 'violet', 'k']
    label_alpha = ['res(1.4 GHz)', 'unres(1.4 GHz)', 'all']
    
    ylim_low, ylim_up = -3.5, 0.5


    # plot
    y_ind_arr     = spec_ind_3_1d4
    y_ind_err_arr = spec_ind_3_1d4_err

    x_line = np.linspace(xlim_low, xlim_up, 20)
    med_arr = y_ind_arr
    median, median_bs = bootstrap_error(arr=med_arr, number=1000)
    bs_low, bs_up     = median_bs[0], median_bs[1]
    ax.plot(x_line, x_line*0 + median, 
            linestyle='--', linewidth=2, color = 'royalblue', alpha = 0.8, zorder=10,\
            label=r'$\alpha_{\rm median}$=$%.2f^{+%.2f}_{%.2f}$'%(median, bs_up, bs_low),\
           )
    ax.fill_between(x_line, median + bs_up, median + bs_low, color='royalblue', alpha=0.1)

    ax.errorbar(x_ind_arr, y_ind_arr, 
                xerr=x_ind_err_arr, 
                yerr=y_ind_err_arr,
                fmt='o', color='royalblue', ecolor ='lightblue', alpha=0.3, zorder = 20,
                markersize=9, markeredgecolor='navy', markeredgewidth=2,
                label='individual source')

    slt_crit_lst = [is_res_arr, ~is_res_arr]
    for i, slt_crit in enumerate(slt_crit_lst):

        # plot the medium
        x_line = np.linspace(xlim_low, xlim_up, 20)
        med_arr = y_arr[slt_crit].to_numpy()
        median, median_bs = bootstrap_error(arr=med_arr, number=1000)
        bs_low, bs_up     = median_bs[0], median_bs[1]
        ax.plot(x_line, x_line*0 + median, 
                linestyle='--', linewidth=2, color = color[i], alpha = 0.8, zorder=10,\
                label=r'$\alpha_{\rm median}$(%s)=$%.2f^{+%.2f}_{%.2f}$'%(label_alpha[i], median, bs_up, bs_low),\
               )
        ax.fill_between(x_line, median + bs_up, median + bs_low, color=color[i], alpha=0.1)

    
        ax.errorbar(x_arr[slt_crit], y_arr[slt_crit], 
                    xerr=x_err_arr[slt_crit], 
                    yerr=y_err_arr[slt_crit],
                    fmt='s', color=color[i], ecolor =ecolor[i], alpha=0.8, zorder = 20,
                    markersize=9, markeredgecolor=mcolor[i], markeredgewidth=2,
                    label=label_alpha[i])

    y_arr_all = np.append(y_arr[is_res_arr].to_numpy(), y_arr[~is_res_arr].to_numpy())

    # plot the medium
    i+=1
    med_arr = y_arr_all
    median, median_bs = bootstrap_error(arr=med_arr, number=1000)
    bs_low, bs_up     = median_bs[0], median_bs[1]
    ax.plot(x_line, x_line*0 + median, 
            linestyle='--', linewidth=2, color = color[i], alpha = 0.8, zorder=10,\
            label=r'$\alpha_{\rm median}$(%s)=$%.2f^{+%.2f}_{%.2f}$'%(label_alpha[i], median, bs_up, bs_low),\
           )
    ax.fill_between(x_line, median + bs_up, median + bs_low, color=color[i], alpha=0.1)

    if 'sort3GHz' in fig_fn or 'sortIRAC' in fig_fn:
        # plot the 1.4 GHz detection limit
        S_3GHz_arr  = np.linspace(xlim_low, xlim_up, 100)
        alpha_lim   = cal_spcetral_index(S_3GHz_arr, 3, S_1d4GHZ_DET_LIMIT, 1.4)
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
    # alpha_dict = {
    #     'x_ind_arr': x_ind_arr,
    #     'x_ind_err_arr': x_ind_err_arr,
    #     'alpha_ind': spec_ind_3_1d4,
    #     'alpha_err_ind': spec_ind_3_1d4_err,
    #     'x_arr':x_arr,
    #     'x_err_arr': x_err_arr,
    #     'slt_crit_lst':slt_crit_lst
    # }


    ax.set_xlim(xlim_low, xlim_up)
    ax.set_ylim(ylim_low, ylim_up)
    ax.set_ylabel(r'spectral index $\alpha$', fontsize=14)
    ax.tick_params(which='major', direction='in', length=6, width=1, labelsize=12, top=False, right=True)
    ax.tick_params(which='minor', direction='in', length=4, width=1, labelsize=12, top=False, right=True)
    plt.legend(loc='lower right', fontsize=7)
    ax.set_title(fig_fn[len(PATH_FIGURE):-4], fontsize = 11)
    plt.savefig(fig_fn)

def plot_qIR(df, df_ind, fig_fn, fn_flux):

    # def lum_limit_sample(df_in):
    #     df = df_in.copy()
    #     df = df.drop(df[df.z_m_Ugne> 4.0].index)
    #     df = df.drop(df[df.spec_ind_3_1d4_arr< 1.55].index)
    #     # df = df.drop(df[df.z_m > 4.0].index)
    #     return df

    # df_ind = lum_limit_sample(df_ind)


    # individual source
    S_3GHz_ind        = df_ind['flux_3GHz'].to_numpy()/1e6      # uJy -> Jy
    S_3GHz_ind_err    = df_ind['flux_err_3GHz'].to_numpy()/1e6  # uJy -> Jy
    S_1d4GHz_ind      = df_ind['Flux_corr_1d4GHz'].to_numpy()    # Jy
    S_1d4GHz_ind_err  = df_ind['Eflux_corr_1d4GHz'].to_numpy()   # Jy
    # S_1d4GHz_ind      = df_ind['flux_1d4GHz'].to_numpy()/1e3    # mJy -> Jy
    # S_1d4GHz_ind_err  = df_ind['ferr_1d4GHz'].to_numpy()/1e3    # mJy -> Jy
    res_1d4GHz_arr    = df_ind['res_1d4GHz'].to_numpy()   # 1 is resolved, 0 is unresolved

    S_3GHz_ind_gv       = gv.gvar(S_3GHz_ind,     S_3GHz_ind_err  )
    S_1d4GHz_ind_gv     = gv.gvar(S_1d4GHz_ind,   S_1d4GHz_ind_err)
    # spectral index
    alpha_ind_gv  = cal_spcetral_index(S_3GHz_ind_gv, 3, S_1d4GHz_ind_gv, 1.4)
    alpha_ind, alpha_ind_err = MyStatistic.gv2arr(alpha_ind_gv)

    # bin
    plot_ind_res_unres =  False
    plot_bin_res_unres = False
    is_res_arr = df['is_resolve_1d4GHz']
    # S_3GHz_bin          = df['S_3GHz_res_unres_arr'].to_numpy()/1e6
    # S_3GHz_bin_err      = df['S_3GHz_res_unres_err_arr'].to_numpy()/1e6
    # S_1d4GHz_bin        = df['S_1d4GHz_res_unres_arr'].to_numpy()/1e6
    # S_1d4GHz_bin_err    = df['S_1d4GHz_res_unres_err_arr'].to_numpy()/1e6
    
    S_3GHz_bin          = df['St_3GHz_arr'].to_numpy()/1e6
    S_3GHz_bin_err      = df['St_3GHz_err_arr'].to_numpy()/1e6
    # S_1d4GHz_bin        = df['Sp_1d4GHz_arr'].to_numpy()/1e6
    # S_1d4GHz_bin_err    = df['Sp_1d4GHz_err_arr'].to_numpy()/1e6
    # i = 1
    S_1d4GHz_bin        = df['St_1d4GHz_arr'].to_numpy()/1e6
    S_1d4GHz_bin_err    = df['St_1d4GHz_err_arr'].to_numpy()/1e6
    i = 0

    alpha_bin           = df['spec_ind_3_1d4_arr'].to_numpy()
    alpha_bin_err       = df['spec_ind_3_1d4_err_arr'].to_numpy()

    fig = plt.figure()
    plt.rcParams["font.family"] = "serif"
    ax  = fig.add_axes([0.12,0.12,0.8,0.75]) # left, bottom, width, height (range 0 to 1)

    # Local value
    x = np.linspace(0, 5, 20)
    y1_Bell_03 = x*0 + 2.64 + 0.26
    y2_Bell_03 = x*0 + 2.64 - 0.26
    # ax.plot(x, y1_Bell_03, '-', color='navy', linewidth=1)
    # ax.plot(x, y2_Bell_03, '-', color='navy', linewidth=1)
    # ax.fill_between(x, y1_Bell_03, y2_Bell_03, color='tab:blue', alpha=0.1, \
    #                 label = r'q$_{IR}$(local, Bell 2003)',\
    #                )

    qIR_Del = 2.88*np.power(1+x, -0.19)
    ax.plot(x, qIR_Del, '-', color='tab:red', linewidth=2,\
           label = r'q$_{IR}$(Delhaize+ 2017)')


    ax = plot_850um_qIR_z(ax)
    # if 'sort3GHz' in fig_fn:
    #     x_arr           = S_3GHz_bin
    #     x_err_arr       = S_3GHz_bin_err
    #     x_ind_arr       = S_3GHz_ind
    #     x_ind_err_arr   = S_3GHz_ind_err
    #     xlim_low, xlim_up = 0   , 0.15
    #     ax.set_xlabel(r'S$_{\rm 3 GHz}$ [mJy]', fontsize=14)
    # elif 'sortReshift' in fig_fn:

    is_use_Ugne_SED = True
    if ('Ugne' in fig_fn) and is_use_Ugne_SED:
        print('Ugne!')
        x_arr           = df['z_m_Ugne']
        x_err_arr       = 0.5*(df['z_16_Ugne']+df['z_84_Ugne']- 2*x_arr)
        x_ind_arr       = df_ind['z_m_Ugne']
        x_ind_err_arr   = 0.5*(df_ind['z_16_Ugne']+df_ind['z_84_Ugne']- 2*x_ind_arr)
    else: # if 'Lim' in fig_fn:
        x_arr           = df['redshift']
        x_err_arr       = df['redshift_err']
        x_ind_arr       = df_ind['redshift']
        x_ind_err_arr   = df_ind['redshift_err']

    xlim_low, xlim_up = 0   , 5
    ax.set_xlabel('redshift', fontsize=14)      

    # x_ind_arr           = alpha_dict['x_ind_arr']
    # x_ind_err_arr       = alpha_dict['x_ind_err_arr']
    # spec_ind_3_1d4      = alpha_dict['alpha_ind']
    # spec_ind_3_1d4_err  = alpha_dict['alpha_err_ind']

    # x_bin_arr           = alpha_dict['x_arr']
    # x_bin_err_arr       = alpha_dict['x_err_arr']
    # slt_crit_lst        = alpha_dict['slt_crit_lst']
    
    
    if ('Ugne' in fig_fn) and is_use_Ugne_SED:
        print('Ugne!')
        LIR_ind_arr     = df_ind['Ldust_m_Ugne'].to_numpy()
        LIR_ind_err_arr = 0.5*(df_ind['Ldust_16_Ugne'].to_numpy()+df_ind['Ldust_84_Ugne'].to_numpy()- 2*LIR_ind_arr)
        z_ind           =  df_ind['z_m_Ugne'].to_numpy()
        z_err_ind       =  0.5*(df_ind['z_16_Ugne'].to_numpy()+df_ind['z_84_Ugne'].to_numpy()- 2*z_ind)
        
        LIR_bin_arr     = df['Ldust_m_Ugne'].to_numpy()
        LIR_bin_err_arr = 0.5*(df['Ldust_16_Ugne'].to_numpy()+df['Ldust_84_Ugne'].to_numpy()- 2*LIR_bin_arr)
        z_bin           = df['z_m_Ugne'].to_numpy()
        z_err_bin       =  0.5*(df['z_16_Ugne'].to_numpy()+df['z_84_Ugne'].to_numpy()- 2*z_bin)
    else: # if 'Lim' in fig_fn:
        LIR_ind_arr     = df_ind['logLIR'].to_numpy()
        LIR_ind_err_arr = df_ind['logLIR_err'].to_numpy()
        z_ind           =  df_ind['redshift'].to_numpy()
        z_err_ind       =  df_ind['redshift_err'].to_numpy()
        LIR_bin_arr     = df['logLIR'].to_numpy()
        LIR_bin_err_arr = df['logLIR_err'].to_numpy()
        z_bin           =  df['redshift'].to_numpy()
        z_err_bin       =  df['redshift_err'].to_numpy()

    # calculate qIR error
    z_ind_gv        = gv.gvar(z_ind,     z_err_ind)
    S_3GHz_ind_gv   = gv.gvar(S_3GHz_ind,     S_3GHz_ind_err  )
    S_1d4GHz_ind_gv = gv.gvar(S_1d4GHz_ind,   S_1d4GHz_ind_err)
    alpha_ind_gv    = cal_spcetral_index(S_3GHz_ind_gv, 3, S_1d4GHz_ind_gv, 1.4)
    LIR_ind_arr_gv  = gv.gvar(LIR_ind_arr,     LIR_ind_err_arr)


    q_IR_ind_gv     = func_cal_FIRRC(z_ind_gv, S_1d4GHz_ind_gv, LIR_ind_arr_gv, alpha=alpha_ind_gv)
    q_IR_ind, q_IR_ind_err = MyStatistic.gv2arr(q_IR_ind_gv)

    z_bin_gv        = gv.gvar(z_bin,     z_err_bin)
    S_3GHz_bin_gv   = gv.gvar(S_3GHz_bin,     S_3GHz_bin_err  )
    S_1d4GHz_bin_gv = gv.gvar(S_1d4GHz_bin,   S_1d4GHz_bin_err)
    alpha_bin_gv    = cal_spcetral_index(S_3GHz_bin_gv, 3, S_1d4GHz_bin_gv, 1.4)
    LIR_bin_arr_gv  = gv.gvar(LIR_bin_arr,     LIR_bin_err_arr)

    q_IR_bin_gv     = func_cal_FIRRC(z_bin_gv, S_1d4GHz_bin_gv, LIR_bin_arr_gv, alpha=alpha_bin_gv)
    q_IR_bin, q_IR_bin_err = MyStatistic.gv2arr(q_IR_bin_gv)

    # plot
    
    y_ind_arr     = q_IR_ind
    y_ind_err_arr = q_IR_ind_err 

    if plot_ind_res_unres:
        slt_crit_lst = [np.ma.make_mask(res_1d4GHz_arr), ~np.ma.make_mask(res_1d4GHz_arr)]
        mcolor = ['#4F4F4F', 'darkviolet', 'k']
        ecolor = ['silver', 'violet', 'k']
        color  = ['darkgray', 'violet', 'k']
        label_alpha = ['res(1.4 GHz)', 'unres(1.4 GHz)', 'all']

        for i, slt_crit in enumerate(slt_crit_lst):

            # plot the medium
            x_line = np.linspace(xlim_low, xlim_up, 20)
            med_arr = y_ind_arr[slt_crit]
            median, median_bs = bootstrap_error(arr=med_arr, number=1000)
            bs_low, bs_up     = median_bs[0], median_bs[1]
            ax.plot(x_line, x_line*0 + median, 
                    linestyle='--', linewidth=2, color = color[i], alpha = 0.8, zorder=10,\
                    label=r'$\alpha_{\rm median}$(%s)=$%.2f^{+%.2f}_{%.2f}$'%(label_alpha[i], median, bs_up, bs_low),\
                   )
            # ax.fill_between(x_line, median + bs_up, median + bs_low, color=color[i], alpha=0.2)
        
            ax.errorbar(x_ind_arr[slt_crit], y_ind_arr[slt_crit], 
                        xerr=x_ind_err_arr[slt_crit], 
                        yerr=y_ind_err_arr[slt_crit],
                        fmt='o', color=color[i], ecolor =ecolor[i], alpha=0.3, zorder = 20,
                        markersize=9, markeredgecolor=mcolor[i], markeredgewidth=2,
                        label=label_alpha[i])

    else:
        y_ind_arr     = q_IR_ind
        y_ind_err_arr = q_IR_ind_err 

        x_line = np.linspace(xlim_low, xlim_up, 20)
        med_arr = y_ind_arr
        median, median_bs = bootstrap_error(arr=med_arr, number=1000)
        bs_low, bs_up     = median_bs[0], median_bs[1]
        ax.plot(x_line, x_line*0 + median, 
                linestyle='--', linewidth=2, color = 'royalblue', alpha = 0.7, zorder=10,\
                label=r'$q_{\rm IR}$=$%.2f^{+%.2f}_{%.2f}$'%(median, bs_up, bs_low),\
               )


        ax.fill_between(x_line, median + bs_up, median + bs_low, color='royalblue', alpha=0.1)

        ax.errorbar(x_ind_arr, y_ind_arr, 
                    xerr=x_ind_err_arr, 
                    yerr=y_ind_err_arr,
                    fmt='o', color='royalblue', ecolor ='lightblue', alpha=0.1, zorder = 20,
                    markersize=5, markeredgecolor='navy', markeredgewidth=1,
                    # label='individual source'
                    )

        xbin_m, xbin_std, ybin_m, ybin_std = func_binning(x_ind_arr, y_ind_arr, bin_num=18)
        xbin_std = np.nan_to_num(xbin_std)
        ybin_std = np.nan_to_num(ybin_std)

        ax.errorbar(xbin_m, ybin_m,\
                    xerr=np.abs(xbin_std).T,\
                    yerr=np.abs(ybin_std).T,\
                    fmt='o', color='royalblue', ecolor = 'royalblue', alpha=0.8,\
                    markersize=8, markeredgecolor='navy', markeredgewidth=1, capsize=3,\
                    label=r'450$\mu$m selection', 
                    zorder=20
                   )

    # plot bin
    x_arr = z_bin
    y_arr = q_IR_bin
    x_err_arr = z_err_bin
    y_err_arr = q_IR_bin_err
    slt_crit_lst = [is_res_arr, ~is_res_arr]
    mcolor = ['k', 'darkviolet', 'k']
    ecolor = ['#4F4F4F', 'violet', 'k']
    color  = ['#4F4F4F', 'violet', 'k']
    label_alpha = [r'450$\mu$m (stacks)', 'unres(1.4 GHz)', 'all']
    # mcolor = ['#4F4F4F', 'darkviolet', 'k']
    # ecolor = ['silver', 'violet', 'k']
    # color  = ['darkgray', 'violet', 'k']
    # label_alpha = ['res(1.4 GHz)', 'unres(1.4 GHz)', 'all']
    
    
    ylim_low, ylim_up = -3.5, 0.5

    # save the qIR and redshit to csv table
    df_flux = pd.read_csv(fn_flux)
    if 'Unnamed: 0' in df_flux.columns:
        df_flux.drop(df_flux.columns[df_flux.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)

    def append_nan(x):
        x = np.append(x, np.nan)
        return x
    # print(fn_flux)
    # print(df_flux)
    try:
        df_flux['qIR']          = q_IR_bin
    except ValueError:
        q_IR_bin = append_nan(q_IR_bin)
        q_IR_bin_err = append_nan(q_IR_bin_err)
        z_bin = append_nan(z_bin)
        z_err_bin = append_nan(z_err_bin)
        df_flux['qIR']          = q_IR_bin
    df_flux['qIR_err']      = q_IR_bin_err
    df_flux['z_bin']        = z_bin
    df_flux['z_bin_err']    = z_err_bin
    # save csv
    df_flux.to_csv(fn_flux, index=True, header=True)
    print('Save matched catalog %s'%(fn_flux))

    if plot_bin_res_unres:
        for i, slt_crit in enumerate(slt_crit_lst):

            # plot the medium
            x_line = np.linspace(xlim_low, xlim_up, 20)
            med_arr = y_arr[slt_crit]
            median, median_bs = bootstrap_error(arr=med_arr, number=1000)
            bs_low, bs_up     = median_bs[0], median_bs[1]
            ax.plot(x_line, x_line*0 + median, 
                    linestyle='--', linewidth=2, color = color[i], alpha = 0.8, zorder=10,\
                    # label=r'$q_{\rm IR}_{\rm median}$(%s)=$%.2f^{+%.2f}_{%.2f}$'%(label_alpha[i], median, bs_up, bs_low),\
                    label=r'$q_{\rm IR}$=$%.2f^{+%.2f}_{%.2f}$'%(median, bs_up, bs_low),\
                    )
            # ax.fill_between(x_line, median + bs_up, median + bs_low, color=color[i], alpha=0.2)

        
            ax.errorbar(x_arr[slt_crit], y_arr[slt_crit], 
                        xerr=x_err_arr[slt_crit], 
                        yerr=y_err_arr[slt_crit],
                        fmt='s', color=color[i], ecolor =ecolor[i], alpha=0.8, zorder = 20,
                        markersize=9, markeredgecolor=mcolor[i], markeredgewidth=2, capsize=3,
                        label=label_alpha[i])

        i+=1
    else:
        ax.errorbar(x_arr, y_arr, 
                        xerr=x_err_arr, 
                        yerr=y_err_arr,
                        fmt='s', color=color[i], ecolor =ecolor[i], alpha=0.8, zorder = 20,
                        markersize=9, markeredgecolor=mcolor[i], markeredgewidth=2, capsize=3,
                        label=label_alpha[i])

    y_arr_all = np.append(y_arr[is_res_arr], y_arr[~is_res_arr])

    # plot the medium
    
    # med_arr = y_arr_all
    med_arr = q_IR_bin
    median, median_bs = bootstrap_error(arr=med_arr, number=1000)
    bs_low, bs_up     = median_bs[0], median_bs[1]
    ax.plot(x_line, x_line*0 + median, 
            linestyle='--', linewidth=2, color = color[i], alpha = 0.8, zorder=10,\
            label=r'$q_{\rm IR}$=$%.2f^{+%.2f}_{%.2f}$'%(median, bs_up, bs_low),
            # label=r'$q{\rm IR}_{\rm median}$(%s)=$%.2f^{+%.2f}_{%.2f}$'%(label_alpha[i], median, bs_up, bs_low),\
           )
    ax.fill_between(x_line, median + bs_up, median + bs_low, color=color[i], alpha=0.1)


    ax.set_ylabel(r'q$_{IR}$', fontsize=14)
    ax.set_xlim(xlim_low, xlim_up)
    # ax.set_ylim(0.3, 3.2)
    ax.set_ylim(1.4, 3.0)
    ax.tick_params(which='major', direction='in', length=6, width=1, labelsize=12, top=False, right=True)
    ax.tick_params(which='minor', direction='in', length=4, width=1, labelsize=12, top=False, right=True)
    plt.legend(loc='upper right', fontsize=7)
    ax.set_title(fig_fn[len(PATH_FIGURE):-4], fontsize = 11)
    plt.savefig(fig_fn)


def plot_850um_qIR_z(ax):

    df_850um_Ugne   = pd.read_csv('%sAS2UDS_match_r1d6arcs_detected_1d4GHz.csv'%(PATH_DATA))

    def lum_limit_sample(df_in):
        df = df_in.copy()
        df = df.drop(df[df.z_m < 1.5].index)
        df = df.drop(df[df.z_m > 4.0].index)
        return df

    df_850um_Ugne = lum_limit_sample(df_850um_Ugne)

    S_1d4GHz_850     = df_850um_Ugne['Flux'].to_numpy()/1e3 # mJy
    S_1d4GHz_err_850 = df_850um_Ugne['Flux_err']/1e3 # mJy
    S_1d4GHz_850_gv  = gv.gvar(S_1d4GHz_850,   S_1d4GHz_err_850)
    z_850            = df_850um_Ugne['z_m'].to_numpy()
    z_err_850        = 0.5*(df_850um_Ugne['z_16'].to_numpy() + 
                            df_850um_Ugne['z_84'].to_numpy() - 
                            2*z_850)  
    z_850_gv        = gv.gvar(z_850,     z_err_850)
    L_IR_850         = df_850um_Ugne['Ldust_m'].to_numpy()      # [Lsun]
    L_IR_err_850     = 0.5*(df_850um_Ugne['Ldust_16'].to_numpy() + 
                            df_850um_Ugne['Ldust_84'].to_numpy() - 
                            2*L_IR_850)
    L_IR_850_gv  = gv.gvar(L_IR_850,     L_IR_err_850)
    q_IR_850_gv, L_1d4GHz_log10_850_gv  = func_cal_FIRRC_850um(z_850_gv, S_1d4GHz_850_gv, L_IR_850_gv, alpha=None)
    q_IR_850, q_IR_err_850 = MyStatistic.gv2arr(q_IR_850_gv)
    L_1d4GHz_log10_850, L_1d4GHz_log10_err_850 = MyStatistic.gv2arr(L_1d4GHz_log10_850_gv)

    q_IR_850[q_IR_850<1.55] = np.nan


    x_850_arr = z_850
    xerr_850_arr = z_err_850
    y_850_arr = q_IR_850
    yerr_850_arr = q_IR_err_850

    ax.errorbar(x_850_arr, y_850_arr, 
                xerr=xerr_850_arr, 
                yerr=yerr_850_arr,
                fmt='o', color='palegreen', ecolor ='lightgreen', alpha=0.1, #zorder = 20,
                markersize=5, markeredgecolor='limegreen', markeredgewidth=1,
#                         label=r'850$\mu$m selection (#=%d)'%(np.shape(y_850_arr)[0])
#                         label='Ugne (850 um) (#=%d)'%(np.shape(y_850_arr)[0])
                )

    x_line = np.linspace(0, 5, 100)
    try:
        med_arr = y_850_arr.to_numpy()
    except AttributeError:
        med_arr = y_850_arr
    median, median_bs = bootstrap_error(arr=med_arr, number=1000)
    bs_low, bs_up     = median_bs[0], median_bs[1]
    ax.plot(x_line, x_line*0 + median, 
            linestyle='--', linewidth=2, color = 'limegreen', alpha = 0.8, zorder=10,\
            # label=r'$qIR_{\rm median}$=$%.2f^{+%.2f}_{%.2f}$'%(median, bs_up, bs_low),\
            label=r'$q_{\rm IR}$=$%.2f^{+%.2f}_{%.2f}$'%(median, bs_up, bs_low)
           )
    
    xbin_m, xbin_std, ybin_m, ybin_std = func_binning(x_850_arr, y_850_arr, bin_num=40)
    xbin_std = np.nan_to_num(xbin_std)
    ybin_std = np.nan_to_num(ybin_std)

    ax.errorbar(xbin_m, ybin_m,\
                xerr=np.abs(xbin_std).T,\
                yerr=np.abs(ybin_std).T,\
                fmt='o', color='lime', ecolor = 'limegreen', alpha=0.8,\
                markersize=8, markeredgecolor='k', markeredgewidth=1, capsize=3,\
                label=r'850$\mu$m selection', 
                zorder=20
               )

    return ax

def func_binning(x_arr, y_arr, bin_num):

    # x_arr = z_arr
    # y_arr = q_IR_IR
    # bins_axis = x
    
    #
    x_add_arr    = np.empty(bin_num - len(x_arr)%bin_num)
    x_add_arr[:] = np.nan
    x_arr        = np.concatenate((x_arr, x_add_arr), axis=None)
    
    y_add_arr    = np.empty(bin_num - len(y_arr)%bin_num)
    y_add_arr[:] = np.nan
    y_arr        = np.concatenate((y_arr, y_add_arr), axis=None)

    # mean binning
    xy_2d_arr      = np.vstack((x_arr, y_arr)).T
    x_ord_y_2d_arr = xy_2d_arr[xy_2d_arr[:,0].argsort()]
    x_ord_x        = x_ord_y_2d_arr.T[0]
    y_ord_x        = x_ord_y_2d_arr.T[1] 

    #bin_num = 10
#     xbin_m   = np.nanmean(x_ord_x[:(len(x_ord_x)//bin_num)*bin_num].reshape(-1,bin_num), axis=1)
#     xbin_std = np.nanstd(x_ord_x[:(len(x_ord_x)//bin_num)*bin_num].reshape(-1,bin_num), axis=1)
#     ybin_m   = np.nanmean(y_ord_x[:(len(y_ord_x)//bin_num)*bin_num].reshape(-1,bin_num), axis=1)
#     ybin_std = np.nanstd(y_ord_x[:(len(y_ord_x)//bin_num)*bin_num].reshape(-1,bin_num), axis=1)
    xbin_m   = np.nanmedian(x_ord_x[:(len(x_ord_x)//bin_num)*bin_num].reshape(-1,bin_num), axis=1)
    xbin_std = np.nanstd(x_ord_x[:(len(x_ord_x)//bin_num)*bin_num].reshape(-1,bin_num), axis=1)
    ybin_m   = np.nanmedian(y_ord_x[:(len(y_ord_x)//bin_num)*bin_num].reshape(-1,bin_num), axis=1)
    ybin_std = np.nanstd(y_ord_x[:(len(y_ord_x)//bin_num)*bin_num].reshape(-1,bin_num), axis=1)

#     xbin_m, xbin_m_bs = bootstrap_error(arr=med_arr, number=1000)
#     bs_low, bs_up     = median_bs[0], median_bs[1]



    xbin_m, ybin_m, xbin_std, ybin_std = \
    remove_xynan(xbin_m, ybin_m, xerrdata=xbin_std, yerrdata=ybin_std)
    
    return xbin_m, xbin_std, ybin_m, ybin_std

def remove_xynan(xdata, ydata, xerrdata=None, yerrdata=None, return_mask=False):
    mask = ~np.isnan(xdata) & ~np.isnan(ydata)
    xdata_nonan = remove_mask(xdata, mask)
    ydata_nonan = remove_mask(ydata, mask)
    xerrdata_nonan = np.nan
    yerrdata_nonan = np.nan
    if xerrdata is not None:
        xerrdata_nonan = remove_mask(xerrdata, mask)
    if yerrdata is not None:
        yerrdata_nonan = remove_mask(yerrdata, mask)
    return xdata_nonan, ydata_nonan, xerrdata_nonan, yerrdata_nonan

def remove_mask(data, mask):
    data_nonan = data.copy()
    data_nonan = data[mask]
    return data_nonan

def listdir_nohidden(path):
    return [f for f in os.listdir(path) if not f.startswith('.')]

def cal_L_1d4GHz(S_1d4GHz_Jy, z, alpha = -0.8):    
    '''
    Purpose:   
      Calculate luminosity at a rest frequency of 1.4 GHz from the observer-frame frequency of 1.4 GH.
    
    Input:
    - S_1d4GHz_uJy [float]: flux density at the observer-frame frequency of 1.4 GHz 
                            (uJy = 1e-6 * 1e-23 * erg s-1 cm-2 Hz-1)
    - z            [float]: redshift (unitless)
    - alpha        [float]: spectral index (unitless)    
    - D_L          [float]: luminosity distance (cm)
        
    Return:    
    - L_1d4GHz [float]: luminosity at a rest frequency of 1.4 GHz (erg s-1 Hz-1)   
    '''
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    try:
        z_m, z_err = MyStatistic.gv2arr(z)
    except AttributeError:
        z_m = z
    u_Mpc2cm     = 3.08567758e24 # 1 Mpc in cm unit
    D_L_cm       = cosmo.luminosity_distance(z_m).value * u_Mpc2cm # luminosity distance (cm)
    S_1d4GHz_cgs = S_1d4GHz_Jy * 1e-23 # flux density at 1.4 GHz (cgs unit)
    L_1d4GHz     = np.divide( 4*np.pi * D_L_cm**2 * S_1d4GHz_cgs, np.power(1+z, 1+alpha) )  

    return L_1d4GHz

def convert_flux(S_org, freq_org, freq_fin, alpha):
    S_fin = S_org * np.power(np.divide(freq_fin, freq_org), alpha)
    return S_fin

def func_cal_FIRRC(z, S_1d4GHz_Jy, L_IR_log10, alpha=-0.8):
    '''
    Prupose:  
      Calculate the far-infrared radio correlation (q_IR) from the catalog    

    Input:
    - z           [float]: redshift
    - L_IR_log10  [float]: IR-luminosity (erg s-1 Hz-1)
    - S_radio_Jy  [float]: luminosity at a rest-frame frequency of 1.4 GHz (1e-23 * 1e-3 * erg s-1 Hz-1)
    
    Return:   
    - q_IR        [float]: far-infrared radio correlation (unitless)
    '''
    L_SUN = 3.9e33 # (cgs)
    L_1d4GHz_cgs = cal_L_1d4GHz(S_1d4GHz_Jy, z, alpha)
    L_IR_cgs     = np.divide(10**(L_IR_log10), 3.75e12)*L_SUN
    
    q_IR         = np.log(L_IR_cgs)/np.log(10) - np.log(L_1d4GHz_cgs)/np.log(10)     
    return q_IR

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

def cal_L_1d4GHz_850um(S_1d4GHz, z, alpha = -0.8):
    
    '''
    Purpose: 
    
      Calculate luminosity at a rest frequency of 1.4 GHz from the observer-frame frequency of 1.4 GH.
    
    
    Input:
    
    - S_1d4GHz [float]: flux density at the observer-frame frequency of 1.4 GHz 
                        (uJy = 1e-6 * 1e-23 * erg s-1 cm-2 Hz-1)
    - z        [float]: redshift (unitless)
    - alpha    [float]: spectral index (unitless)
    
    - D_L      [float]: luminosity distance (cm)
    
    
    Return:
    
    - L_1d4GHz [float]: luminosity at a rest frequency of 1.4 GHz (erg s-1 Hz-1)
    
    '''
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    try:
        z_m, z_err = MyStatistic.gv2arr(z)
    except AttributeError:
        z_m = z
    u_Mpc2cm     = 3.08567758e24                                   # 1 Mpc in cm unit
    D_L_cm       = cosmo.luminosity_distance(z_m).value * u_Mpc2cm   # luminosity distance (cm)
    S_1d4GHz_uJy = S_1d4GHz                                        # flux density at 1.4 GHz (uJy)
    S_1d4GHz_cgs = S_1d4GHz_uJy * 1e-6 * 1e-23                     # flux density at 1.4 GHz (cgs unit)
    
    # luminosity (cgs unit)
    L_1d4GHz     = np.divide( 4*np.pi * D_L_cm**2 * S_1d4GHz_cgs, np.power(1+z, 1+alpha) )  
    #L_1d4GHz     = 4*np.pi * D_L_cm**2 * S_1d4GHz_cgs * np.power(1+z, alpha -1)
    
    return L_1d4GHz
    
def cal_q_IR(L_FIR, L_1d4GHz):
    '''
    Prupose:
    
      Calculate the far-infrared radio correlation (q_IR)
    
    
    Input:
    
    - L_FIR    [float]: FIR-luminosity (erg s-1 Hz-1)
    - L_1d4GHz [float]: luminosity at a rest-frame frequency of 1.4 GHz (erg s-1 Hz-1)
    
    
    Return:
    
    - q_IR     [float]: far-infrared radio correlation (unitless)
    '''
    #return np.log10(np.divide( L_FIR, 3.75e12)) - np.log10(L_1d4GHz) 
    return np.log10(L_FIR) - np.log10(L_1d4GHz) 

def convert_L_log(L_org, freq_org, freq_fin, alpha):
    L_fin = L_org * np.power(np.divide(freq_fin, freq_org), alpha)
    #L_log_fin = L_log_org + alpha * (np.log10(freq_fin) - np.log10(freq_org))
    return L_fin


def func_cal_FIRRC_850um(z, S_radio_mJy, L_IR, alpha=None):
    '''
    Prupose:
    
      Calculate the far-infrared radio correlation (q_IR) from the catalog
    
    
    Input:
    
    - z           [float]: redshift
    - L_IR        [float]: IR-luminosity (erg s-1 Hz-1)
    - S_radio_mJy [float]: luminosity at a rest-frame frequency of 1.4 GHz (1e-23 * 1e-3 * erg s-1 Hz-1)
    
    
    Return:
    
    - q_IR        [float]: far-infrared radio correlation (unitless)
    '''
    
    if alpha == None:
        # 1.4 GHz detected     
        S_1d4GHz_uJy = S_radio_mJy*1e3        # convert mJy to uJy

    else:
        # 1.4 GHz undetected
        # calculate the 1.4 GHz flux from 3 GHz flux
        S_3GHz_uJy   = S_radio_mJy*1e3        # convert mJy to uJy
        S_1d4GHz_uJy = convert_L_log(S_3GHz_uJy, 3, 1.4, alpha)
        
  
    # calculate the luminosity at a rest frequency of 1.4 GHz
    L_1d4GHz_cgs = cal_L_1d4GHz_850um(S_1d4GHz_uJy, z)
    L_1d4GHz_log10 = np.log(L_1d4GHz_cgs/L_SUN)/np.log(10)
    # calculate the q_IR
    L_IR_log10  = L_IR
    L_IR_cgs    = np.divide(10**(L_IR_log10), 3.75e12)*L_SUN
    q_IR         = np.log(L_IR_cgs)/np.log(10) - np.log(L_1d4GHz_cgs)/np.log(10) 
#     q_IR         = cal_q_IR(L_IR_cgs, L_1d4GHz_cgs)
    
    return q_IR, L_1d4GHz_log10

if __name__ == '__main__':
    main()