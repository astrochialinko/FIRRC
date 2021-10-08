#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: script_make_radec.py
Name: Chia-Lin Ko
Create Date: Jun 14, 2021
Last Modified Date: Jun 14, 2021
------------------------
This program aims to make the stacking image
 
"""
import numpy as np
import gvar as gv
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import FlatLambdaCDM

# my own packing
import path
from ClassStatistic import MyStatistic

# Constant
# path
PATH_CATALOG_CROSSMATCH         = path.catalog_crossmatch
PATH_DATA_VLA                   = path.image_vla
IMG_FILENAME_1d4GHzXS           = '%scosmos_vla_1d4GHz_XS_2021image_uJy.fits'%(PATH_DATA_VLA)
IMG_FILENAME_PB                 = '%scosmos_vla_1d4GHz_XS_2021pb.fits'%(PATH_DATA_VLA)
IMG_FILENAME_1d4GHzXS_PBCORR    = '%scosmos_vla_1d4GHz_XS_2021image_uJy_pbcorr.fits'%(PATH_DATA_VLA)
ALPHA = -0.8 # spectral index
QIR_AGN_THLD = 1.55

def main():


    make_1d4GHz_pbcorr()

    # filename
    fn_match_Ugne_Lim_3GHzlp_irac       = '%scosmos_match_450um_Ugne_Lim_1d4GHzXS_3GHzlp_irac'%(PATH_CATALOG_CROSSMATCH)
    fn_match_Ugne_Lim_3GHzlp_mips_irac  = '%scosmos_match_450um_Ugne_Lim_1d4GHzXS_3GHzlp_mips_irac'%(PATH_CATALOG_CROSSMATCH)

    fn_match_Lim_3GHzlp_irac            = '%scosmos_match_450um_Lim_1d4GHzXS_3GHzlp_irac'%(PATH_CATALOG_CROSSMATCH)
    fn_match_Lim_3GHzlp_mips_irac       = '%scosmos_match_450um_Lim_1d4GHzXS_3GHzlp_mips_irac'%(PATH_CATALOG_CROSSMATCH)

    fn_match_Gao_Lim_3GHzlp_irac        = '%scosmos_match_450um_Gao_Lim_1d4GHzXS_3GHzlp_irac'%(PATH_CATALOG_CROSSMATCH)
    fn_match_Gao_Lim_3GHzlp_mips_irac   = '%scosmos_match_450um_Gao_Lim_1d4GHzXS_3GHzlp_mips_irac'%(PATH_CATALOG_CROSSMATCH)

    # combine catalog
    fn_Ugne = '%scosmos_match_Ugne_comb'%(PATH_CATALOG_CROSSMATCH)
    fn_Lim  = '%scosmos_match_Lim_comb'%(PATH_CATALOG_CROSSMATCH)
    fn_Gao  = '%scosmos_match_Gao_comb'%(PATH_CATALOG_CROSSMATCH)
    combine_catalog(fn_match_Ugne_Lim_3GHzlp_irac, 
                    fn_match_Ugne_Lim_3GHzlp_mips_irac, 
                    'ra_radio',  fn_Ugne)
    combine_catalog(fn_match_Lim_3GHzlp_irac, 
                    fn_match_Lim_3GHzlp_mips_irac, 
                    'ra_radio',  fn_Lim)
    combine_catalog(fn_match_Gao_Lim_3GHzlp_irac, 
                    fn_match_Gao_Lim_3GHzlp_mips_irac, 
                    'ra_radio',  fn_Gao)

    # Gao SNR>4
    clip_SNR(fn_Gao,            fn_Gao+'_4SNR',     snr_thld=4)
    clip_SNR(fn_Gao,            fn_Gao+'_5SNR',     snr_thld=5)

    # remove AGN
    remove_AGN(fn_Ugne,         fn_Ugne+'_nAGN') 
    remove_AGN(fn_Lim,          fn_Lim+'_nAGN')
    remove_AGN(fn_Gao,          fn_Gao+'_nAGN')
    remove_AGN(fn_Gao+'_4SNR',  fn_Gao+'_4SNR_nAGN')
    remove_AGN(fn_Gao+'_5SNR',  fn_Gao+'_5SNR_nAGN')

def make_1d4GHz_pbcorr():
    img_1d4GHzXS        = pyfits.getdata(IMG_FILENAME_1d4GHzXS)
    img_pd              = pyfits.getdata(IMG_FILENAME_PB)
    hd_1d4GHzXS         = pyfits.getheader(IMG_FILENAME_1d4GHzXS)

    img_1d4GHzXS_pdcorr = np.divide(img_1d4GHzXS, img_pd)
    pyfits.writeto(IMG_FILENAME_1d4GHzXS_PBCORR, img_1d4GHzXS_pdcorr, hd_1d4GHzXS, overwrite=True)


def clip_SNR(df_fn, df_SNR_fn, snr_thld):
    df = pd.read_csv(df_fn+'.csv')
    slt_SNR = (df['SNR_450'] > snr_thld)    

    df_SNR = df[slt_SNR].reset_index(drop=True)
    df_SNR = drop_unname(df_SNR)
    df_SNR.to_csv(df_SNR_fn+'.csv', index=True, header=True)

def remove_AGN(df_fn, df_nAGN_fn):
    df = pd.read_csv(df_fn+'.csv')
    # slt_AGN = (df['XrayAGN_3GHz_multiAGN'] != True) & (df['MIRAGN_3GHz_multiAGN'] != True) &\
    #       (df['SEDAGN_3GHz_multiAGN']  != True) & (df['RExcess_3GHz_multiAGN']!= True)
    # slt_nomulti = (df['multi_3GHz']!=1)
    # df_nAGN = df[slt_nomulti & slt_AGN].reset_index(drop=True)
    for col in df.columns:
        if 'is_radio_AGNs' in col:
            slt_nRadioAGN = df[col] != True
            df_nAGN = df[slt_nRadioAGN].reset_index(drop=True)   
            df_nAGN = drop_unname(df_nAGN)
            df_nAGN.to_csv(df_nAGN_fn+col.split('is_radio_AGNs')[-1]+'.csv', index=True, header=True)

def combine_catalog(df1_fn, df2_fn, det_band, df_comb_fn):
    df1 = pd.read_csv(df1_fn+'.csv')
    df2 = pd.read_csv(df2_fn+'.csv')

    df1_det   = df1[~df1[det_band].isna()]
    df2_undet = df2[df2[det_band].isna()]

    df3 = df1_det.append(df2_undet, ignore_index=True)
    df3 = cal_catalog(df3) # cal alpha, qIR
    df3 = drop_unname(df3)
    df3.to_csv(df_comb_fn+'.csv', index=True, header=True)

def drop_unname(df):
    if 'Unnamed: 0' in df.columns:
        df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
    return df

def cal_catalog(df_ind):

    def cal_alpha_qIR(  S_3GHz_ind,     S_3GHz_ind_err, S_1d4GHz_ind,   S_1d4GHz_ind_err,\
                        LIR_ind_arr,     LIR_ind_err_arr, z_ind,     z_ind_err):
        S_3GHz_ind_gv       = gv.gvar(S_3GHz_ind,     S_3GHz_ind_err  )
        S_1d4GHz_ind_gv     = gv.gvar(S_1d4GHz_ind,   S_1d4GHz_ind_err)
        LIR_ind_arr_gv  = gv.gvar(LIR_ind_arr,     LIR_ind_err_arr)
        z_ind_gv        = gv.gvar(z_ind,     z_ind_err)

        # spectral index
        alpha_ind_gv  = cal_spcetral_index(S_3GHz_ind_gv, 3, S_1d4GHz_ind_gv, 1.4)
        alpha_ind, alpha_ind_err = MyStatistic.gv2arr(alpha_ind_gv)
        # replace nan with -0.8
        alpha_ind[np.isnan(alpha_ind)] = ALPHA
        alpha_ind_err[np.isnan(alpha_ind_err)] = 0
        alpha_ind_gv = gv.gvar(alpha_ind,     alpha_ind_err)

        # S_1d4GHz_gv     = cal_S1d4GHz_from_S3GHz(S_3GHz_ind_gv, alpha_ind_gv)
        S_1d4GHz_from_3GHz_gv = cal_S1d4GHz_from_S3GHz(S_3GHz_ind_gv, alpha_ind_gv)
        S_1d4GHz_from_3GHz, S_1d4GHz_from_3GHz_err = MyStatistic.gv2arr(S_1d4GHz_from_3GHz_gv)
        S_1d4GHz        = np.where(~np.isnan(S_1d4GHz_ind), S_1d4GHz_ind, S_1d4GHz_from_3GHz)
        S_1d4GHz_err    = np.where(~np.isnan(S_1d4GHz_ind_err), S_1d4GHz_ind_err, S_1d4GHz_from_3GHz_err)
        S_1d4GHz_gv     = gv.gvar(S_1d4GHz,     S_1d4GHz_err)
        q_IR_ind_gv, L_1d4GHz_log10_ind_gv = func_cal_FIRRC(z_ind_gv, S_1d4GHz_gv, LIR_ind_arr_gv, alpha=alpha_ind_gv)
        q_IR_ind, q_IR_ind_err = MyStatistic.gv2arr(q_IR_ind_gv)
        L_1d4GHz_log10_ind, L_1d4GHz_log10_ind_err = MyStatistic.gv2arr(L_1d4GHz_log10_ind_gv)

        is_radio_AGNs = q_IR_ind<QIR_AGN_THLD

        return alpha_ind, alpha_ind_err, L_1d4GHz_log10_ind, L_1d4GHz_log10_ind_err, \
                q_IR_ind, q_IR_ind_err, is_radio_AGNs

    def wcs2pix(f_hd, ra_deg_arr, dec_deg_arr):
        wcs     = WCS(f_hd)
        len_arr = len(ra_deg_arr)
        ra_pix_arr  = np.zeros(len_arr)
        dec_pix_arr = np.zeros(len_arr)

        for i in range(len_arr):
            try:
                pixcrd          = wcs.wcs_world2pix([[ra_deg_arr[i], dec_deg_arr[i], 0, 0]], 1)[0] 

                ra_pix_arr[i]   = int(np.around(pixcrd[0]))  # ra_pix,  get the closest pixel by round the number
                dec_pix_arr[i]  = int(np.around(pixcrd[1]))
            except ValueError:
                ra_pix_arr[i]   = np.nan
                dec_pix_arr[i]  = np.nan

        return ra_pix_arr, dec_pix_arr

    def make_image_2d(fitsfiles):
        """
        convert the image to 2 dimension
        """
        # check the dimension of the image
        f_shape = np.shape(fitsfiles)    # the shape of the input data( array)
        f_dim   = len(f_shape)   # the dimension of the input data


        if f_dim == 4:
            fitsfiles = fitsfiles[0][0]
            #print('Covert the 4D image to 2D')
        elif f_dim == 2:
            pass
            #print('The input fitfile is a 2D image.')
        else:
            print('Please check the dimension of the input fitsfile')

        return fitsfiles

    # individual source
    S_3GHz_lp           = df_ind['flux_3GHz'].to_numpy()/1e6      # uJy -> Jy
    S_3GHz_lp_err       = df_ind['flux_err_3GHz'].to_numpy()/1e6  # uJy -> Jy
    S_1d4GHz_XS         = df_ind['Flux_corr_1d4GHz'].to_numpy()   # Jy
    S_1d4GHz_XS_err     = df_ind['Eflux_corr_1d4GHz'].to_numpy()  # Jy
    S_1d4GHz_dp         = df_ind['flux_1d4GHz'].to_numpy()/1e3    # mJy -> Jy
    S_1d4GHz_dp_err     = df_ind['ferr_1d4GHz'].to_numpy()/1e3    # mJy -> Jy

    # primary beam for 1.4 GHz XS
    ra_1d4GHz_XS_deg    = df_ind['RA_1d4GHz'].to_numpy()
    dec_1d4GHz_XS_deg   = df_ind['DEC_1d4GHz'].to_numpy()
    pb_img  = pyfits.getdata(IMG_FILENAME_PB)
    pb_img  = make_image_2d(pb_img)
    pb_hd   = pyfits.getheader(IMG_FILENAME_PB)
    ra_1d4GHz_XS_pix, dec_1d4GHz_XS_pix = wcs2pix(pb_hd, ra_1d4GHz_XS_deg, dec_1d4GHz_XS_deg)
    df_ind['RA_1d4XS_pix']  = ra_1d4GHz_XS_pix  #(1-based)
    df_ind['DEC_1d4XS_pix'] = dec_1d4GHz_XS_pix #(1-based)

    pb_value = np.zeros(len(ra_1d4GHz_XS_pix))
    for i in range(len(pb_value)):
        try:
            pb_value[i] = pb_img[int(dec_1d4GHz_XS_pix[i]-1)][int(ra_1d4GHz_XS_pix[i]-1)] #(0-based)
        except ValueError:
            pb_value[i] = np.nan
    df_ind['PB_1d4XS']  = pb_value

    try:
        LIR_Ugne     = df_ind['Ldust_m_Ugne'].to_numpy()
        LIR_Ugne_err = 0.5*(df_ind['Ldust_16_Ugne'].to_numpy()+df_ind['Ldust_84_Ugne'].to_numpy()-2*LIR_Ugne)      
        z_Ugne       = df_ind['z_m_Ugne'].to_numpy()
        z_Ugne_err   = 0.5*(df_ind['z_16_Ugne'].to_numpy()+df_ind['z_84_Ugne'].to_numpy()- 2*z_Ugne)
        isUgne_exist = True
    except KeyError: # 
        isUgne_exist = False
    LIR_Lim     = df_ind['logLIR'].to_numpy()
    LIR_Lim_err = df_ind['logLIR_err'].to_numpy()
    z_Lim       = df_ind['redshift'].to_numpy()
    z_Lim_err   = df_ind['redshift_err'].to_numpy()

    # 1: Lim, 1.4 GHz XS
    alpha_3lp_1d4XS, alpha_3lp_1d4XS_err, L_1d4GHz_log10_3lp_1d4XS, L_1d4GHz_log10_3lp_1d4XS_err,\
    qIR_3lp1d4XS_Lim, qIR_3lp1d4XS_Lim_err, is_radio_AGNs_3lp1d4XS_Lim = \
    cal_alpha_qIR(  S_3GHz_lp,     S_3GHz_lp_err,   S_1d4GHz_XS,    S_1d4GHz_XS_err,\
                    LIR_Lim,  LIR_Lim_err,    z_Lim,         z_Lim_err)
    df_ind['alpha_3lp_1d4XS']              = alpha_3lp_1d4XS
    df_ind['alpha_3lp_1d4XS_err']          = alpha_3lp_1d4XS_err
    df_ind['L_1d4GHz_log10_3lp_1d4XS']     = L_1d4GHz_log10_3lp_1d4XS
    df_ind['L_1d4GHz_log10_3lp_1d4XS_err'] = L_1d4GHz_log10_3lp_1d4XS_err    
    df_ind['qIR_3lp1d4XS_Lim']             = qIR_3lp1d4XS_Lim
    df_ind['qIR_3lp1d4XS_Lim_err']         = qIR_3lp1d4XS_Lim_err
    df_ind['is_radio_AGNs_3lp1d4XS_Lim']   = is_radio_AGNs_3lp1d4XS_Lim

    # 2: Lim, 1.4 GHz dp
    alpha_3lp_1d4dp, alpha_3lp_1d4dp_err, L_1d4GHz_log10_3lp_1d4dp, L_1d4GHz_log10_3lp_1d4dp_err,\
    qIR_3lp1d4dp_Lim, qIR_3lp1d4dp_Lim_err, is_radio_AGNs_3lp1d4dp_Lim = \
    cal_alpha_qIR(  S_3GHz_lp,     S_3GHz_lp_err,   S_1d4GHz_dp,    S_1d4GHz_dp_err,\
                    LIR_Lim,  LIR_Lim_err,    z_Lim,         z_Lim_err)
    df_ind['alpha_3lp_1d4dp']              = alpha_3lp_1d4dp
    df_ind['alpha_3lp_1d4dp_err']          = alpha_3lp_1d4dp_err
    df_ind['L_1d4GHz_log10_3lp_1d4dp']      =  L_1d4GHz_log10_3lp_1d4dp
    df_ind['L_1d4GHz_log10_3lp_1d4dp_err']  = L_1d4GHz_log10_3lp_1d4dp_err    
    df_ind['qIR_3lp1d4dp_Lim']             = qIR_3lp1d4dp_Lim
    df_ind['qIR_3lp1d4dp_Lim_err']         = qIR_3lp1d4dp_Lim_err
    df_ind['is_radio_AGNs_3lp1d4dp_Lim']   = is_radio_AGNs_3lp1d4dp_Lim

    if isUgne_exist:
        # 3: Ugne, 1.4 GHz XS
        alpha_3lp_1d4XS, alpha_3lp_1d4XS_err, L_1d4GHz_log10_3lp_1d4XS, L_1d4GHz_log10_3lp_1d4XS_err,\
        qIR_3lp1d4XS_Ugne, qIR_3lp1d4XS_Ugne_err, is_radio_AGNs_3lp1d4XS_Ugne = \
        cal_alpha_qIR(  S_3GHz_lp,     S_3GHz_lp_err,   S_1d4GHz_XS,    S_1d4GHz_XS_err,\
                        LIR_Ugne,  LIR_Ugne_err,    z_Ugne,         z_Ugne_err)
        # df_ind['alpha_3lp_1d4XS']               = alpha_3lp_1d4XS
        # df_ind['alpha_3lp_1d4XS_err']           = alpha_3lp_1d4XS_err
        df_ind['qIR_3lp1d4XS_Ugne']             = qIR_3lp1d4XS_Ugne
        df_ind['qIR_3lp1d4XS_Ugne_err']         = qIR_3lp1d4XS_Ugne_err
        df_ind['is_radio_AGNs_3lp1d4XS_Ugne']   = is_radio_AGNs_3lp1d4XS_Ugne

        # 4: Ugne, 1.4 GHz dp
        alpha_3lp_1d4dp, alpha_3lp_1d4dp_err, L_1d4GHz_log10_3lp_1d4dp, L_1d4GHz_log10_3lp_1d4dp_err,\
        qIR_3lp1d4dp_Ugne, qIR_3lp1d4dp_Ugne_err, is_radio_AGNs_3lp1d4dp_Ugne = \
        cal_alpha_qIR(  S_3GHz_lp,     S_3GHz_lp_err,   S_1d4GHz_dp,    S_1d4GHz_dp_err,\
                        LIR_Ugne,  LIR_Ugne_err,    z_Ugne,         z_Ugne_err)
        # df_ind['alpha_3lp_1d4dp']               = alpha_3lp_1d4dp
        # df_ind['alpha_3lp_1d4dp_err']           = alpha_3lp_1d4dp_err
        df_ind['qIR_3lp1d4dp_Ugne']             = qIR_3lp1d4dp_Ugne
        df_ind['qIR_3lp1d4dp_Ugne_err']         = qIR_3lp1d4dp_Ugne_err
        df_ind['is_radio_AGNs_3lp1d4dp_Ugne']   = is_radio_AGNs_3lp1d4dp_Ugne

    return df_ind

def cal_spcetral_index(S1, freq1, S2, freq2):
    return np.divide(np.log(S1)-np.log(S2), np.log(freq1)-np.log(freq2))

def cal_S1d4GHz_from_S3GHz(S_3GHz, alpha = -0.8):
    return np.power(1.4/3, alpha) * S_3GHz

def cal_L_1d4GHz (S_1d4GHz_Jy, z, alpha = -0.8):    
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
    L_sun = 3.9e33 # (cgs)
    L_1d4GHz_cgs = cal_L_1d4GHz(S_1d4GHz_Jy, z, alpha)
    L_1d4GHz_log10 = np.log(L_1d4GHz_cgs/L_sun)/np.log(10)
    L_IR_cgs     = np.divide(10**(L_IR_log10), 3.75e12)*L_sun
    q_IR         = np.log(L_IR_cgs)/np.log(10) - np.log(L_1d4GHz_cgs)/np.log(10)     
    return q_IR, L_1d4GHz_log10

if __name__ == '__main__':
    main()