#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: script_make_stacking_studies.py
Name: Chia-Lin Ko
Create Date: Jun 11, 2021
Last Modified Date: Jul 22, 2021
------------------------
This program aims to make the stacking image
 
"""
import os
import numpy as np
import pandas as pd
import astropy.io.fits as pyfits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from photutils import aperture_photometry
from photutils import EllipticalAperture

# my own packing
import path
# from ClassStacking import Stacking

# Constant
# path
PATH_CATALOG_CROSSMATCH = path.catalog_crossmatch
PATH_DATA_VLA           = path.image_vla
PATH_TABLE              = path.table
PATH_FIGURE             = path.figure

# PATH_DATA_VLA           = '../Data/COSMOS/Image/VLA/'
# PATH_CATALOG_CROSSMATCH = '../Data/COSMOS/Catalog/CrossMatch/'
# PATH_TABLE              = '../Data/Tables/'
# PATH_FIGURE             = '../Figures/'

# filename for catalog
FN_UGNE         = '%scosmos_match_Ugne_comb'%(PATH_CATALOG_CROSSMATCH)
FN_LIM          = '%scosmos_match_Lim_comb'%(PATH_CATALOG_CROSSMATCH)
FN_GAO4SNR      = '%scosmos_match_Gao_comb_4SNR'%(PATH_CATALOG_CROSSMATCH)
FN_GAO5SNR      = '%scosmos_match_Gao_comb_5SNR'%(PATH_CATALOG_CROSSMATCH)
FN_UGNE_NAGN    = FN_UGNE + '_nAGN_3lp1d4XS_Ugne'
FN_LIM_NAGN     = FN_LIM + '_nAGN_3lp1d4XS_Lim'
FN_GAO4SNR_NAGN = FN_GAO4SNR + '_nAGN_3lp1d4XS_Lim'
FN_GAO5SNR_NAGN = FN_GAO5SNR + '_nAGN_3lp1d4XS_Lim'

FN_GAO          = FN_GAO5SNR
FN_GAO_NAGN     = FN_GAO5SNR_NAGN

# list for images                      
IMG_FILENAME_LIST       =   [
    '%scosmos_vla_3GHz_2017image_uJy.fits'%(PATH_DATA_VLA),
    '%scosmos_vla_1d4GHz_XS_2021image_uJy_pbcorr.fits'%(PATH_DATA_VLA)
    # '%scosmos_vla_1d4GHz_XS_2021image_uJy.fits'%(PATH_DATA_VLA)
                            ]
IMG_RMS_FILENAME_LIST   =   [
    '%scosmos_vla_3GHz_2017rms_uJy.fits'%(PATH_DATA_VLA),
    '%scosmos_vla_1d4GHz_XS_2021rms_uJy.fits'%(PATH_DATA_VLA),
                            ]


IMG_LABEL_LIST          = ['3GHz', '1d4GHzXS']
IMG_RMS_LIST            = [ 2.3, 1.8]    # rms (uJy/beam)

# value
S_3GHZ_DET_LIMIT        = 5*2.3e-3      # 5 sigma detection limit (mJy/beam), from Smolčić+17
# S_1d4GHZ_DET_LIMIT      = 4*12e-3       # 4 sigma detection limit (mJy/beam), from Schinnerer+10
S_1d4GHZ_DET_LIMIT      = 5*1.8e-3      # 4 sigma detection limit (mJy/beam), from Hiddo+21

IS_BINNING = True
# BINNING_ELE_NUM = 20
BINNING_ELE_NUM = 18

def main():


    # stack_det3GHz(is_sort_3GHz=True,  is_sort_redshift=True, is_sort_LIR=False)
    # stack_detIRAC(is_sort_3GHz=True,  is_sort_1d4GHz=True, is_sort_IRAC=False, is_sort_redshift=True, is_sort_LIR=False)
    stack_det3GHzIRAC(  is_Coord3GHz_sort3GHz=True,     is_CoordIRAC_sort3GHz=True,
                        is_Coord3GHz_sort1d4GHz=True,   is_CoordIRAC_sort1d4GHz=True,
                        is_Coord3GHz_sortIRAC=True,     is_CoordIRAC_sortIRAC=True, 
                        is_Coord3GHz_sortredshift=True, is_CoordIRAC_sortredshift=True,
                        is_Coord3GHz_sortLIR=False,      is_CoordIRAC_sortLIR=False,
                        )


#####################################################################################
# functions
#####################################################################################

def stack_det3GHzIRAC(  is_Coord3GHz_sort3GHz=False,        is_CoordIRAC_sort3GHz=False,
                        is_Coord3GHz_sort1d4GHz=False,      is_CoordIRAC_sort1d4GHz=False,
                        is_Coord3GHz_sortIRAC=False,        is_CoordIRAC_sortIRAC=False, 
                        is_Coord3GHz_sortredshift=False,    is_CoordIRAC_sortredshift=False,
                        is_Coord3GHz_sortLIR=False,         is_CoordIRAC_sortLIR=False,
                        ):
    subset = 'det3GHzIRAC'    
    
    # Coordinate: VLA 3 GHz
    # sorted with 3 GHz
    if is_Coord3GHz_sort3GHz:  
        coord = '3GHz' 
        sort  = '3GHz'     
        fn_lst          = [FN_UGNE, FN_LIM, FN_GAO4SNR, FN_UGNE_NAGN, FN_LIM_NAGN, FN_GAO4SNR_NAGN]
        col_sort_lst    = ['flux_3GHz']*6
        out_label_lst   = [ 'Ugne_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s'%(subset, coord, sort),
                            'Ugne_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s_nAGN'%(subset, coord, sort)
                            ]
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   &
                        (df['ra_irac'].notnull())   ]   # only 3 GHz & IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ '3GHz', 'ra_3GHz', 'dec_3GHz', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/'
                        )

    # sorted with 1.4 GHz
    if is_Coord3GHz_sort1d4GHz:  
        coord = '3GHz' 
        sort  = '1d4GHz'     
        fn_lst          = [FN_UGNE, FN_LIM, FN_GAO4SNR, FN_UGNE_NAGN, FN_LIM_NAGN, FN_GAO4SNR_NAGN]
        col_sort_lst    = ['Flux_corr_1d4GHz']*6
        out_label_lst   = [ 'Ugne_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s'%(subset, coord, sort),
                            'Ugne_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s_nAGN'%(subset, coord, sort)
                            ]
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   &
                        (df['ra_irac'].notnull())   ]   # only 3 GHz & IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ '3GHz', 'ra_3GHz', 'dec_3GHz', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/'
                        )

    
    # sorted with IRAC
    if is_Coord3GHz_sortIRAC:
        coord = '3GHz' 
        sort  = 'IRAC'     
        fn_lst          = [FN_UGNE, FN_LIM, FN_GAO4SNR, FN_UGNE_NAGN, FN_LIM_NAGN, FN_GAO4SNR_NAGN]
        col_sort_lst    = ['flux_c1_1']*6
        out_label_lst   = [ 'Ugne_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s'%(subset, coord, sort),
                            'Ugne_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s_nAGN'%(subset, coord, sort)
                            ]
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   &
                        (df['ra_irac'].notnull())   ]   # only 3 GHz & IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ '3GHz', 'ra_3GHz', 'dec_3GHz', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/'
                        )
    
    # sorted with redshift
    if is_Coord3GHz_sortredshift:
        coord = '3GHz' 
        sort  = 'redshift'     
        fn_lst          = [FN_UGNE, FN_LIM, FN_UGNE_NAGN, FN_LIM_NAGN]
        col_sort_lst    = ['z_m_Ugne', 'redshift']*2
        out_label_lst   = [ 'Ugne_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Ugne_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            ]
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   &
                        (df['ra_irac'].notnull())   ]   # only 3 GHz & IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ '3GHz', 'ra_3GHz', 'dec_3GHz', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/'
                        )

    # sorted with redshift
    if is_Coord3GHz_sortLIR:
        coord = '3GHz' 
        sort  = 'LIR'     
        fn_lst          = [FN_UGNE, FN_LIM, FN_UGNE_NAGN, FN_LIM_NAGN]
        col_sort_lst    = ['Ldust_m_Ugne', 'logLIR']*2
        out_label_lst   = [ 'Ugne_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Ugne_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            ]
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   &
                        (df['ra_irac'].notnull())   ]   # only 3 GHz & IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ '3GHz', 'ra_3GHz', 'dec_3GHz', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/'
                        )

    #-------------------------
    # Coordinate: Spitzer IRAC
    
    # sorted with 3 GHz
    if is_CoordIRAC_sort3GHz:
        coord = 'IRAC' 
        sort  = '3GHz'     
        fn_lst          = [FN_UGNE, FN_LIM, FN_GAO4SNR, FN_UGNE_NAGN, FN_LIM_NAGN, FN_GAO4SNR_NAGN]
        col_sort_lst    = ['flux_3GHz']*6
        out_label_lst   = [ 'Ugne_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s'%(subset, coord, sort),
                            'Ugne_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s_nAGN'%(subset, coord, sort)
                            ]
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   &
                        (df['ra_irac'].notnull())   ]   # only 3 GHz & IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ 'IRAC', 'ra_irac', 'dec_irac', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/'
                        )

    # sorted with 1.4 GHz
    if is_CoordIRAC_sort1d4GHz:
        coord = 'IRAC' 
        sort  = '1d4GHz'     
        fn_lst          = [FN_UGNE, FN_LIM, FN_GAO4SNR, FN_UGNE_NAGN, FN_LIM_NAGN, FN_GAO4SNR_NAGN]
        col_sort_lst    = ['Flux_corr_1d4GHz']*6
        out_label_lst   = [ 'Ugne_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s'%(subset, coord, sort),
                            'Ugne_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s_nAGN'%(subset, coord, sort)
                            ]
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   &
                        (df['ra_irac'].notnull())   ]   # only 3 GHz & IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ 'IRAC', 'ra_irac', 'dec_irac', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/'
                        )

    # sorted with IRAC
    if is_CoordIRAC_sortIRAC:
        coord = 'IRAC' 
        sort  = 'IRAC'     
        fn_lst          = [FN_UGNE, FN_LIM, FN_GAO4SNR, FN_UGNE_NAGN, FN_LIM_NAGN, FN_GAO4SNR_NAGN]
        col_sort_lst    = ['flux_c1_1']*6
        out_label_lst   = [ 'Ugne_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s'%(subset, coord, sort),
                            'Ugne_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Gao4SNR_%s_coord%s_sort%s_nAGN'%(subset, coord, sort)
                            ]
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   &
                        (df['ra_irac'].notnull())   ]   # only 3 GHz & IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ 'IRAC', 'ra_irac', 'dec_irac', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/'
                        )

    # sorted with redshift
    if is_CoordIRAC_sortredshift:
        coord = 'IRAC' 
        sort  = 'redshift'     
        fn_lst          = [FN_UGNE, FN_LIM, FN_UGNE_NAGN, FN_LIM_NAGN]
        col_sort_lst    = ['z_m_Ugne', 'redshift']*2
        out_label_lst   = [ 'Ugne_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Ugne_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            ]
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   &
                        (df['ra_irac'].notnull())   ]   # only 3 GHz & IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ 'IRAC', 'ra_irac', 'dec_irac', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/'
                        )
    
    # sorted with redshift
    if is_CoordIRAC_sortLIR:
        coord = 'IRAC' 
        sort  = 'LIR'     
        fn_lst          = [FN_UGNE, FN_LIM, FN_UGNE_NAGN, FN_LIM_NAGN]
        col_sort_lst    = ['Ldust_m_Ugne', 'logLIR']*2
        out_label_lst   = [ 'Ugne_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s'%(subset, coord, sort), 
                            'Ugne_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            'Lim_%s_coord%s_sort%s_nAGN'%(subset, coord, sort), 
                            ]
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   &
                        (df['ra_irac'].notnull())   ]   # only 3 GHz & IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ 'IRAC', 'ra_irac', 'dec_irac', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet_IRACDet/'
                        )


def stack_det3GHz(is_sort_3GHz=False, is_sort_redshift=False, is_sort_LIR=False):
    # Coordinate: VLA 3 GHz
    
    # sorted with 3 GHz
    if is_sort_3GHz:        
        fn_lst          = [FN_UGNE, FN_LIM, FN_GAO5SNR, FN_UGNE_NAGN, FN_LIM_NAGN, FN_GAO5SNR_NAGN]
        col_sort_lst    = ['flux_3GHz']*6
        out_label_lst   = [ 'Ugne_coord3GHz_sort3GHz', 'Lim_coord3GHz_sort3GHz', 'Gao5SNR_coord3GHz_sort3GHz',
                            'Ugne_coord3GHz_sort3GHz_nAGN', 'Lim_coord3GHz_sort3GHz_nAGN', 'Gao5SNR_coord3GHz_sort3GHz_nAGN']
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   ]   # only 3 GHz detected sources

            do_stacking(df,
                        stack_lst       = [ '3GHz', 'ra_3GHz', 'dec_3GHz', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet/'
                        )

    # sorted with reshift
    if is_sort_redshift:
        fn_lst          = [FN_UGNE, FN_LIM, FN_UGNE_NAGN, FN_LIM_NAGN]
        col_sort_lst    = ['z_m_Ugne', 'redshift']*2
        out_label_lst   = [ 'Ugne_coord3GHz_sortReshift', 'Lim_coord3GHz_sortReshift',
                            'Ugne_coord3GHz_sortReshift_nAGN', 'Lim_coord3GHz_sortReshift']
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   ]   # only 3 GHz detected sources
            do_stacking(df,
                        stack_lst       = [ '3GHz', 'ra_3GHz', 'dec_3GHz', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet/'
                        )

    # sorted with LIR
    if is_sort_LIR:
        fn_lst          = [FN_UGNE, FN_LIM, FN_UGNE_NAGN, FN_LIM_NAGN]
        col_sort_lst    = ['Ldust_m_Ugne', 'logLIR']*2
        out_label_lst   = [ 'Ugne_coord3GHz_sortLIR', 'Lim_coord3GHz_sortLIR',
                            'Ugne_coord3GHz_sortLIR_nAGN', 'Lim_coord3GHz_sortLIR_nAGN']
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_3GHz'].notnull())   ]   # only IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ 'IRAC', 'ra_3GHz', 'dec_3GHz', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_3GHzDet/'
                                        )


def stack_detIRAC(is_sort_3GHz=False, is_sort_1d4GHz=False, is_sort_IRAC=False, is_sort_redshift=False, is_sort_LIR=False):    
    # Coordinate: Spitzer IRAC
    
    # sorted with 3GHz
    if is_sort_3GHz:
        fn_lst          = [FN_UGNE, FN_LIM, FN_GAO4SNR, FN_UGNE_NAGN, FN_LIM_NAGN, FN_GAO4SNR_NAGN]
        col_sort_lst    = ['flux_3GHz']*6
        out_label_lst   = [ 'Ugne_coordIRAC_sort3GHz', 'Lim_coordIRAC_sort3GHz', 'Gao4SNR_coordIRAC_sort3GHz',
                            'Ugne_coordIRAC_sort3GHz_nAGN', 'Lim_coordIRAC_sort3GHz_nAGN', 'Gao4SNR_coordIRAC_sort3GHz_nAGN']
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_irac'].notnull())   ]   # only IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ 'IRAC', 'ra_irac', 'dec_irac', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg =  '../Data/COSMOS/Image/Stacking/XS_IRACDet/'
                                        )

    # sorted with 1d4GHz
    if is_sort_1d4GHz:
        fn_lst          = [FN_UGNE, FN_LIM, FN_GAO4SNR, FN_UGNE_NAGN, FN_LIM_NAGN, FN_GAO4SNR_NAGN]
        col_sort_lst    = ['Flux_corr_1d4GHz']*6
        out_label_lst   = [ 'Ugne_coordIRAC_sort1d4GHz', 'Lim_coordIRAC_sort1d4GHz', 'Gao4SNR_coordIRAC_sort1d4GHz',
                            'Ugne_coordIRAC_sort1d4GHz_nAGN', 'Lim_coordIRAC_sort1d4GHz_nAGN', 'Gao4SNR_coordIRAC_sort1d4GHz_nAGN']
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_irac'].notnull())   ]   # only IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ 'IRAC', 'ra_irac', 'dec_irac', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg =  '../Data/COSMOS/Image/Stacking/XS_IRACDet/'
                                        )

    # sorted with IRAC
    if is_sort_IRAC:
        fn_lst          = [FN_UGNE, FN_LIM, FN_GAO4SNR, FN_UGNE_NAGN, FN_LIM_NAGN, FN_GAO4SNR_NAGN]
        col_sort_lst    = ['flux_c1_1']*6
        out_label_lst   = [ 'Ugne_coordIRAC_sortIRAC', 'Lim_coordIRAC_sortIRAC', 'Gao4SNR_coordIRAC_sortIRAC',
                            'Ugne_coordIRAC_sortIRAC_nAGN', 'Lim_coordIRAC_sortIRAC_nAGN', 'Gao4SNR_coordIRAC_sortIRAC_nAGN']
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_irac'].notnull())   ]   # only IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ 'IRAC', 'ra_irac', 'dec_irac', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg =  '../Data/COSMOS/Image/Stacking/XS_IRACDet/'
                                        )

    # sorted with reshift
    if is_sort_redshift:
        fn_lst          = [FN_UGNE, FN_LIM, FN_UGNE_NAGN, FN_LIM_NAGN]
        col_sort_lst    = ['z_m_Ugne', 'redshift']*2
        out_label_lst   = [ 'Ugne_coordIRAC_sortRedshift', 'Lim_coordIRAC_sortRedshift',
                            'Ugne_coordIRAC_sortRedshift_nAGN', 'Lim_coordIRAC_sortRedshift_nAGN']
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_irac'].notnull())   ]   # only IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ 'IRAC', 'ra_irac', 'dec_irac', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_IRACDet/'
                                        )
 
    # sorted with LIR
    if is_sort_LIR:
        fn_lst          = [FN_UGNE, FN_LIM, FN_UGNE_NAGN, FN_LIM_NAGN]
        col_sort_lst    = ['Ldust_m_Ugne', 'logLIR']*2
        out_label_lst   = [ 'Ugne_coordIRAC_sortLIR', 'Lim_coordIRAC_sortLIR',
                            'Ugne_coordIRAC_sortLIR_nAGN', 'Lim_coordIRAC_sortLIR_nAGN']
        for i, fn in enumerate(fn_lst):
            df = pd.read_csv(fn+'.csv')
            df = df[    (df['ra_irac'].notnull())   ]   # only IRAC detected sources
            do_stacking(df,
                        stack_lst       = [ 'IRAC', 'ra_irac', 'dec_irac', col_sort_lst[i], 
                                            'mediam', out_label_lst[i]  ],
                        path_saveimg    = '../Data/COSMOS/Image/Stacking/XS_IRACDet/'
                                        )

def do_stacking(df, stack_lst, path_saveimg):
    # 
    '''
    Purpose
        stacking the image based on the input parameters
    ---------------------
    Input Parameter
        df     [DataFrame]: subset catalog for stacking coordinate
        stack_lst   [list]: the list contain the parameters below
                            slt_radec, the coordinate to be used for stacking
                            col_ra, the column name of ra
                            col_dec, the column name of dec
                            col_sort, the column name of sorted column
                            stack_stat, statistic for the stacking, such as medium, mean.
                            out_label, 
    ---------------------
    No Return
    '''  

    # set the stacking parameter
    slt_radec, col_ra, col_dec, col_sort, stack_stat, out_label = stack_lst

    # set the ra, dec, and sort array
    ra_arr      = df[col_ra].to_numpy()
    dec_arr     = df[col_dec].to_numpy()
    sort_arr    = df[col_sort].to_numpy()

    
    # start stacking
    for i, img_filename in enumerate(IMG_FILENAME_LIST):
        stacking_half_width = 88 if i == 0 else 50 # length = 35.35"

        stack = Stacking(img_filename)
        f_img = pyfits.getdata(img_filename)
        f_hd  = pyfits.getheader(img_filename)


        # transfer the wcs to pixel coordinate
        coord_pix_nlst  = ['%s_pix'%(slt_radec), '%s_pix'%(slt_radec)]
        coord_pix_fn    = '%sCoord_pix_%s_img%s.txt'%(PATH_TABLE, slt_radec, IMG_LABEL_LIST[i])
        print(coord_pix_fn)
        coord_pix_arr   = stack.wcs2pix(ra_arr, dec_arr, coord_pix_nlst, is_savetxt=True, outname=coord_pix_fn)
        ra_pix_arr      = coord_pix_arr.T[0]
        dec_pix_arr     = coord_pix_arr.T[1]
        print('Start to stack %d images for %s'%(len(coord_pix_arr), IMG_LABEL_LIST[i]))

        if IS_BINNING:

            # sort the ra and dec based on the sorted array
            ra_pix_sort_2darr   = sort_for_binning(ra_pix_arr , sort_arr, BINNING_ELE_NUM)
            dec_pix_sort_2darr  = sort_for_binning(dec_pix_arr, sort_arr, BINNING_ELE_NUM)
            sort_sort_2darr     = sort_for_binning(sort_arr, sort_arr, BINNING_ELE_NUM)

            bin_gp_num = len(ra_pix_sort_2darr)           
            # create new directory if not exist
            path_data_out = '%s%s_%s_bin%sgp%s/'%(path_saveimg, out_label, stack_stat, BINNING_ELE_NUM, bin_gp_num) 
            try:
                os.makedirs(path_data_out)
            except FileExistsError:
                pass

            # save the sorted and binned csv
            df_sort = df.sort_values(by=[col_sort]).reset_index(drop=True)
            df_sort = drop_unname(df_sort)                
            df_sort.to_csv('%s%s.csv'%(path_data_out, out_label), index=True) 
            
            df_bin = df_sort.groupby(np.arange(len(df_sort))//BINNING_ELE_NUM).median()
            df_bin = drop_unname(df_bin)
            df_bin.to_csv('%s%s_bin.csv'%(path_data_out, out_label), index=True) 

            print('Start to bin the images based on the order of %s.'%(col_sort))
            print('The data is devided to into %d bins, each bin has %d elements.'%(bin_gp_num, BINNING_ELE_NUM))      

            for bin_gp_id in range(bin_gp_num):
                # remove the nan of last binning subsets
                if bin_gp_id != (bin_gp_num-1):
                    coord_pix_each_bin_arr = np.column_stack((  ra_pix_sort_2darr[bin_gp_id], dec_pix_sort_2darr[bin_gp_id]))
                else:
                    nan_num = np.isnan(ra_pix_sort_2darr).sum()
                    coord_pix_each_bin_arr = np.column_stack((  ra_pix_sort_2darr[bin_gp_id][:-nan_num], 
                                                                dec_pix_sort_2darr[bin_gp_id][:-nan_num] ))
                # stacking
                stackging_2d_fn = '%sStacking_img%s_%s_%s_bin%sgp%s_%s.fits'%(path_data_out, 
                                    IMG_LABEL_LIST[i], out_label, stack_stat, BINNING_ELE_NUM, bin_gp_num, bin_gp_id)                    
                stacking_1d, stacking_2d = stack.make_stacking(IMG_RMS_FILENAME_LIST[i], IMG_RMS_LIST[i], \
                    coord_pix_each_bin_arr, stack_stat, stacking_half_width, stackging_2d_fn)
        else:
            try:
                os.makedirs(path_saveimg)
            except FileExistsError:
                pass
            # stack all the elements into one image
            stackging_2d_fn = '%sStacking_%s_%s_%s.fits'%(path_saveimg, 
                                IMG_LABEL_LIST[i], out_label, stack_stat)
            stacking_1d, stacking_2d = stack.make_stacking(IMG_RMS_FILENAME_LIST[i], IMG_RMS_LIST[i], \
                coord_pix_arr, stack_stat, stacking_half_width, stackging_2d_fn, IMG_LABEL_LIST[i], slt_radec)                    

#####################################################################################
# function
#####################################################################################
def drop_unname(df):
    if 'Unnamed: 0' in df.columns:
        df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
    return df

def sort_for_binning(tobe_sort_arr, sort_std_arr, binning_ele_num=20):
    """
    Purpose
        Transfor the world coordinates to the pixel positions of the input FITS file
    ---------------------
    Input Parameter
        tobe_sort_arr   [1D-ndarray]: to
        sort_std_arr    [1D-ndarray]: right ascension (deg in wcs)
        binning_ele_num [1D-ndarray]: declination (deg in wcs)
    ---------------------
    Return
        sorted_arr      [1D-ndarray]: right ascension and declination   (pixel)
    """   
    # sort by the given arr
    sort_ind_arr = np.argsort(sort_std_arr)
    
    # data preprocessing
    fill_nan_num = int(binning_ele_num - len(sort_ind_arr)%binning_ele_num)
    fill_nan_arr = np.zeros(fill_nan_num)
    fill_nan_arr[:] = np.nan
    sort_ind_for_binning_arr = np.concatenate((sort_ind_arr, fill_nan_arr), axis=None)
    
    binning_num = int(len(sort_ind_for_binning_arr)/binning_ele_num) 
    sort_ind_sub_arr = np.array_split(sort_ind_for_binning_arr, binning_num)

    # 
    row_num, column_num = np.shape(sort_ind_sub_arr)
    sorted_arr = np.zeros((row_num, column_num))

    for row in range(row_num):
        if row != row_num-1:
            sorted_arr[row] = tobe_sort_arr[sort_ind_sub_arr[row].astype(int)]
        elif row == row_num-1:
            sorted_last_arr = tobe_sort_arr[sort_ind_sub_arr[-1][:-fill_nan_num].astype(int)]
            sorted_arr[row] = np.concatenate((sorted_last_arr, fill_nan_arr), axis=None)
    
    return sorted_arr

#####################################################################################
class Stacking:

    def __init__(self, filename):
        self.fn     = filename
        self.f      = pyfits.getdata(self.fn)
        self.f_hd   = pyfits.getheader(self.fn)
        self.f_hd_d = {}
        self.set_beam_info()
       
        # convert the image to 2 dimension
        self.f = self.make_image_2d(self.f)

    def make_image_2d(self, fitsfiles):
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

    def make_stacking(self, img_rms_fn, img_rms, coord_pix_ds9_arr, stacking_stat='mean', 
                          stacking_half_width='20', stackging_2d_fn='stacking.fits', 
                          img_label = None, slt_name = None
                          ):
        """
        Purpose
            stack the image based on the input coordinate info
        ---------------------
        Input Parameter
            img_rms_fn          [str]           : filename of the rms image
            img_rms             [float]         : rms value of the original image 
            coord_pix_ds9_arr   [2d-ndarray]    : 2d array contains the ra and dec arr
            stacking_stat       [str]           : statistic method for stacking 
                                                  (e.g., mean, weighted mean, median)
            stacking_half_width [int]           : half side width of the stacking square (pixel)
            stackging_2d_fn     [str]           : filename for the output stacking 2d image
        ---------------------
        Return
            stacking_1d         [float]         : value of the stacking at the postion
            stacking_2d         [2d-ndarray]    : 2d image of the stacking

        """        

        # set the range for stacking (2d array)  
        stacking_width = 2*stacking_half_width + 1      # the width (side length) of stacking area [pixel]

        # set the rms map
        img_rms_f           = self.make_rms_map(img_rms_fn, img_rms)
        img_rms_f           = self.make_image_2d(img_rms_f)
       
        # stacking in 2d
        stacking_num, column_num = np.shape(coord_pix_ds9_arr)
        stacking_1d_arr     = np.zeros(stacking_num)
        stacking_3d_arr     = np.zeros((stacking_num, stacking_width, stacking_width))
        rms_1d_arr          = np.zeros(stacking_num)
        rms_3d_arr          = np.zeros((stacking_num, stacking_width, stacking_width))
        # 1-based (as in the FITS convention, for example coordinates coming from ds9)
        ra_pix_ds9_arr      = coord_pix_ds9_arr.T[0]    # ra    (pixel, 1-based)
        dec_pix_ds9_arr     = coord_pix_ds9_arr.T[1]    # dec   (pixel, 1-based)            

        nan_2d_arr_num = 0
        for i, coord_pix in enumerate(coord_pix_ds9_arr):   

            # 0-based (as in Numpy arrays)
            ra_pix_np           = int(int(ra_pix_ds9_arr[i]) -1)  # ra    (pixel, 0-based)
            dec_pix_np          = int(int(dec_pix_ds9_arr[i])-1)  # dec   (pixel, 0-based)
            ra_pix_np_start     = ra_pix_np  - stacking_half_width
            ra_pix_np_end       = ra_pix_np  + stacking_half_width
            dec_pix_np_start    = dec_pix_np - stacking_half_width
            dec_pix_np_end      = dec_pix_np + stacking_half_width

            # stacking in 1d 
            try:
                stacking_1d_arr[i]  = self.f[dec_pix_np][ra_pix_np]
                rms_1d_arr[i]       = img_rms_f[dec_pix_np][ra_pix_np]
            except IndexError: # some sources are outside of (1.4 GHz) map
                stacking_1d_arr[i]  = np.nan
                rms_1d_arr[i]       = np.nan
                #print('1D', i, ra_pix_np, dec_pix_np)

            # stacking in 2d 
            # split the 2d area
            stacking_2d_arr     = self.f[dec_pix_np_start:dec_pix_np_end+1].T[ra_pix_np_start:ra_pix_np_end+1].T
            rms_2d_arr          = img_rms_f[dec_pix_np_start:dec_pix_np_end+1].T[ra_pix_np_start:ra_pix_np_end+1].T 
           
            dec_len, ra_len = np.shape(stacking_2d_arr)
            if (ra_len != stacking_width) or (dec_len != stacking_width):
                stacking_1d_arr[i]  = np.nan
                rms_1d_arr[i]       = np.nan
                # stacking in 2d, set to np.nan array
                stacking_2d_arr = rms_2d_arr = np.zeros((stacking_width, stacking_width))
                stacking_2d_arr[:] = np.nan
                rms_2d_arr[:] = np.nan
                #print('2D', i, ra_pix_np, dec_pix_np)

            if np.isnan(stacking_2d_arr.sum()):
                nan_2d_arr_num +=1

            # save the stacking image to a 3d array
            stacking_3d_arr[i] = stacking_2d_arr
            rms_3d_arr[i]      = rms_2d_arr


        # stacking for different statistic method
        if stacking_stat == 'mean':
            stacking_1d     = np.nanmean(stacking_1d_arr)
            stacking_2d     = np.nanmean(stacking_3d_arr, axis=0)
        elif stacking_stat == 'weighted_mean':
            stacking_1d     = MyStatistic.nanaverage(stacking_1d_arr, 1/rms_1d_arr, axis=0)
            stacking_2d     = MyStatistic.nanaverage(stacking_3d_arr, 1/rms_3d_arr, axis=0)
        elif stacking_stat == 'mediam':
            stacking_1d     = np.nanmedian(stacking_1d_arr)
            stacking_2d     = np.nanmedian(stacking_3d_arr, axis=0)
        else:
            print('Please select the stacking statistic from: (1) mean (2) weighted_mean (3) mediam.')

        # save the 2d stacking to FITS file
        pyfits.writeto(stackging_2d_fn, stacking_2d, self.f_hd, overwrite=True)

        # calculate the theoretical error and S/R
        stacking_err_theor  = np.divide(img_rms, np.sqrt(stacking_num))
        print('Non-nan num: ', np.sum((~np.isnan(stacking_1d_arr))))
        print('theoretical error = %.2f'%(stacking_err_theor))
        print('theoretical S/N = %.2f'%(stacking_1d/stacking_err_theor))

        return stacking_1d, stacking_2d

    def make_rms_map(self, img_rms_fn, img_rms):
        """
        Set the rms map

        Input Parameter
            img_rms_fn  [str]       : filename of rms image. None if no rmg image
            img_rms     [float]     : rms value of the image
        Return
            img_rms_f   [2d-ndarray]: rms image
        """
        if img_rms_fn is not None:
            img_rms_f   = pyfits.getdata(img_rms_fn)
        else:
            img_rms_f   = self.f*0 + img_rms

        return img_rms_f

    def set_beam_info(self):
        """set the beam information from the header to the f_hd_d dictionaty"""

        self.f_hd_d['pix_deg']          = self.f_hd['CDELT2']                   # 1 pixel length           [deg]
        self.f_hd_d['pix_arcs']         = self.deg2arcs(self.f_hd_d['pix_deg']) # 1 pixel length           [arcs]
        self.f_hd_d['bmaj_deg']         = self.f_hd['BMAJ']                     # major axis of beam size  [deg]
        self.f_hd_d['bmin_deg']         = self.f_hd['BMIN']                     # minor axis of beam size  [deg]
        self.f_hd_d['bpa_deg']          = self.f_hd['BPA']                      # angle of the beam        [deg]
        # major axis of beam size [pixel]
        self.f_hd_d['bmaj_pix']         = self.deg2pixel(self.f_hd_d['bmaj_deg'], self.f_hd_d['pix_deg'] )
        # minor axis of beam size [pixel]  
        self.f_hd_d['bmin_pix']         = self.deg2pixel(self.f_hd_d['bmin_deg'], self.f_hd_d['pix_deg'] )
        # beam area [pixel]
        self.f_hd_d['beam_area_pix']    = np.divide(np.pi * self.f_hd_d['bmaj_pix'] * self.f_hd_d['bmin_pix'], 
                                                    4 * np.log(2))
    @staticmethod
    def deg2arcs(input_deg):
        return input_deg * 36e2

    @staticmethod
    def deg2pixel(input_deg, pix_deg):
        return np.divide(input_deg, pix_deg)

    def wcs2pix(self, ra_wcs_arr, dec_wcs_arr, coord_name_lst, is_savetxt=False, outname='coord_pix.txt'):
        """
        Purpose
            Transfor the world coordinates to the pixel positions of the input FITS file
        ---------------------
        Input Parameter
            fits_fn                [str]: the filename of a FITs file
            ra_wcs_arr      [1D-ndarray]: right ascension (deg in wcs)
            dec_wcs_arr     [1D-ndarray]: declination (deg in wcs)
            coord_name_lst        [list]: 
            is_savetxt            [bool]: 
            outname                [str]: the filename of the coordinate in pixel
        ---------------------
        Return
            coord_pix_arr   [2D-ndarray]: right ascension and declination   (pixel)
        """   

        # loading WCS information from the FITS file
        wcs             = WCS(self.f_hd)
        coord_pix_arr   = np.zeros((len(ra_wcs_arr), 2))
        print(outname)
              
        for i in range(len(ra_wcs_arr)):
            
            # transforms world coordinates to pixel coordinates
            # 1-based (as in the FITS convention, for example coordinates coming from DS9)
            pixcrd                  = wcs.wcs_world2pix([[ra_wcs_arr[i], dec_wcs_arr[i], 0, 0]], 1)[0] 
            coord_pix_arr.T[0][i]   = int(np.around(pixcrd[0]))  # ra_pix,  get the closest pixel by round the number
            coord_pix_arr.T[1][i]   = int(np.around(pixcrd[1]))  # dec_pix, get the closest pixel by round the number

        if is_savetxt:
            # save the coordinate in pixel to the txt file
            columns     = " ".join(coord_name_lst)
            fmt_lst     = ['%d' for i in coord_name_lst]
            fmt_str     = ' '.join(fmt_lst)
            np.savetxt(outname, coord_pix_arr, fmt=fmt_str, header=columns)

        return coord_pix_arr


if __name__ == '__main__':
    main()