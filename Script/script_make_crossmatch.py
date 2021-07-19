#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: script_make_crossmatch.py
Name: Chia-Lin Ko
Date: Jun 09, 2021
Last Modified Date: Jul 20, 2021
------------------------
This program aims to crossmatch tables.
"""
import numpy as np
import pandas as pd
import os

# my own packing
import path
import FuncTableCrossmatch as tc

# Constant
# path
PATH_CATALOG_JCMT       = path.catalog_jcmt
PATH_CATALOG_VLA        = path.catalog_vla
PATH_CATALOG_IRAC       = path.catalog_irac
PATH_CATALOG_MIPS       = path.catalog_mips
PATH_CATALOG_CROSSMATCH = path.catalog_crossmatch


def main():

    fn_dict = set_filename()
    crossmatch_450um(fn_dict)
    crossmatch_1d4GHz(fn_dict)
    # crossmatch_1d4GHz_diffR(fn_dict)
    crossmatch_3GHz(fn_dict)
    crossmatch_MIPS(fn_dict)
    crossmatch_IRAC(fn_dict)

def set_filename():
    fn_dict = {}
    # original catalog
    fn_dict['fn_JCMT_450um_Ugne']  = '%scosmos_jcmt_450um_SEDUgne_2020cat.csv'%(PATH_CATALOG_JCMT)
    fn_dict['fn_JCMT_450um_Lim']   = '%scosmos_jcmt_450um_2020cat.csv'%(PATH_CATALOG_JCMT)
    fn_dict['fn_JCMT_450um_Gao']   = '%scosmos_jcmt_450um_2021cat.csv'%(PATH_CATALOG_JCMT)
    fn_dict['fn_VLA_1d4GHz_dp']    = '%scosmos_vla_1d4GHz_dp_2011cat.csv'%(PATH_CATALOG_VLA)
    fn_dict['fn_VLA_1d4GHz_XS']    = '%scosmos_vla_1d4GHz_XS_2021cat.csv'%(PATH_CATALOG_VLA)
    fn_dict['fn_VLA_3GHz_lp']      = '%scosmos_vla_3GHz_2017cat.csv'%(PATH_CATALOG_VLA)
    fn_dict['fn_VLA_3GHz_AGN']     = '%scosmos_vla_3GHz_multiAGN_2017cat.csv'%(PATH_CATALOG_VLA)
    fn_dict['fn_IRAC']             = '%scosmos_irac_2007cat.csv'%(PATH_CATALOG_IRAC)
    fn_dict['fn_MIPS_LeFloch']     = '%scosmos_mips_24um_LeFloch_2008cat.csv'%(PATH_CATALOG_MIPS)
    fn_dict['fn_MIPS_whwang']      = '%scosmos_mips_24um_whwang_2020cat.csv'%(PATH_CATALOG_MIPS)

    # cross-matched catalog
    fn_dict['fn_match_VLA_3GHzlp_AGN']             = '%scosmos_match_vla_3GHzlp_AGN.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_VLA_3GHzlp_AGN_1d4GHzdp']    = '%scosmos_match_vla_3GHzlp_AGN_1d4GHzdp.csv'%(PATH_CATALOG_CROSSMATCH)

    # Ugne based crossmatch catalogs
    fn_dict['fn_match_Ugne_Lim']                            = '%scosmos_match_450um_Ugne_Lim.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Ugne_Lim_1d4GHzXS']                   = '%scosmos_match_450um_Ugne_Lim_1d4GHzXS.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Ugne_Lim_1d4GHzXS_3GHzlp']            = '%scosmos_match_450um_Ugne_Lim_1d4GHzXS_3GHzlp.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Ugne_Lim_1d4GHzXS_3GHzlp_irac']       = '%scosmos_match_450um_Ugne_Lim_1d4GHzXS_3GHzlp_irac.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Ugne_Lim_1d4GHzXS_3GHzlp_mips']       = '%scosmos_match_450um_Ugne_Lim_1d4GHzXS_3GHzlp_mips.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Ugne_Lim_1d4GHzXS_3GHzlp_mips_irac']  = '%scosmos_match_450um_Ugne_Lim_1d4GHzXS_3GHzlp_mips_irac.csv'%(PATH_CATALOG_CROSSMATCH)

    # Lim based crossmatch catalogs
    fn_dict['fn_match_Lim_1d4GHzXS']                        = '%scosmos_match_450um_Lim_1d4GHzXS.csv'%(PATH_CATALOG_CROSSMATCH)   
    fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp']                 = '%scosmos_match_450um_Lim_1d4GHzXS_3GHzlp.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp_irac']            = '%scosmos_match_450um_Lim_1d4GHzXS_3GHzlp_irac.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp_mips']            = '%scosmos_match_450um_Lim_1d4GHzXS_3GHzlp_mips.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp_mips_irac']       = '%scosmos_match_450um_Lim_1d4GHzXS_3GHzlp_mips_irac.csv'%(PATH_CATALOG_CROSSMATCH)
    
    # Gao based crossmatch catalogs
    fn_dict['fn_match_Gao_Lim']                             = '%scosmos_match_450um_Gao_Lim.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Gao_Lim_1d4GHzXS']                    = '%scosmos_match_450um_Gao_Lim_1d4GHzXS.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Gao_Lim_1d4GHzXS_3GHzlp']             = '%scosmos_match_450um_Gao_Lim_1d4GHzXS_3GHzlp.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Gao_Lim_1d4GHzXS_3GHzlp_irac']        = '%scosmos_match_450um_Gao_Lim_1d4GHzXS_3GHzlp_irac.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Gao_Lim_1d4GHzXS_3GHzlp_mips']        = '%scosmos_match_450um_Gao_Lim_1d4GHzXS_3GHzlp_mips.csv'%(PATH_CATALOG_CROSSMATCH)
    fn_dict['fn_match_Gao_Lim_1d4GHzXS_3GHzlp_mips_irac']   = '%scosmos_match_450um_Gao_Lim_1d4GHzXS_3GHzlp_mips_irac.csv'%(PATH_CATALOG_CROSSMATCH)

    return fn_dict

def crossmatch_450um(fn_dict):

    fn_Ugne     = fn_dict['fn_JCMT_450um_Ugne']
    fn_Lim      = fn_dict['fn_JCMT_450um_Lim']
    fn_Gao      = fn_dict['fn_JCMT_450um_Gao']
    fn_Ugne_Lim = fn_dict['fn_match_Ugne_Lim']
    fn_Gao_Lim  = fn_dict['fn_match_Gao_Lim']

    input_dict = {} # filename, column_ra, column_dec
    input_dict['JCMT_450um_Ugne']   = [fn_Ugne]
    input_dict['JCMT_450um_Gao']    = [fn_Gao]
    input_dict['JCMT_450um_Lim']    = [fn_Lim,   'RA_450Lim',    'DEC_450Lim']

    # cross match with JCMT 450 um, Lim
    do_merge_value(input_dict, f1_key= 'JCMT_450um_Ugne', f2_key = 'JCMT_450um_Lim', 
        f1_value='ID_STUDIES_Ugne', f2_value='ID_STUDIES_Lim', join_type = 'inner', fn_match = fn_Ugne_Lim)
    do_merge_value(input_dict, f1_key= 'JCMT_450um_Gao', f2_key = 'JCMT_450um_Lim', 
        f1_value='ID_STUDIES_Gao', f2_value='ID_STUDIES_Lim', join_type = 'outer', fn_match = fn_Gao_Lim)

def crossmatch_1d4GHz(fn_dict):  
    
    fn_1d4GHzXS         = fn_dict['fn_VLA_1d4GHz_XS']
    fn_Ugne             = fn_dict['fn_match_Ugne_Lim']
    fn_Lim              = fn_dict['fn_JCMT_450um_Lim']
    fn_Gao              = fn_dict['fn_match_Gao_Lim']
    fn_Ugne_1d4GHzXS    = fn_dict['fn_match_Ugne_Lim_1d4GHzXS']
    fn_Lim_1d4GHzXS     = fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp_irac']
    fn_Gao_1d4GHzXS     = fn_dict['fn_match_Gao_Lim_1d4GHzXS']

    input_dict = {} # filename, column_ra, column_dec
    input_dict['VLA_1d4GHz_XS']      = [fn_1d4GHzXS,    'RA_1d4GHz',    'DEC_1d4GHz']
    input_dict['match_Ugne_Lim']     = [fn_Ugne,        'RA_450Lim',    'DEC_450Lim']
    input_dict['JCMT_450um_Lim']     = [fn_Lim,         'RA_450Lim',    'DEC_450Lim']
    input_dict['match_Gao_Lim']      = [fn_Gao,         'RA_450Gao',    'Dec_450Gao']

    do_merge_radec(input_dict, f1_key= 'match_Ugne_Lim', f2_key = 'VLA_1d4GHz_XS', radius=4, fn_bestmatch = fn_Ugne_1d4GHzXS)
    do_merge_radec(input_dict, f1_key= 'JCMT_450um_Lim', f2_key = 'VLA_1d4GHz_XS', radius=4, fn_bestmatch = fn_Lim_1d4GHzXS)
    do_merge_radec(input_dict, f1_key= 'match_Gao_Lim',  f2_key = 'VLA_1d4GHz_XS', radius=4, fn_bestmatch = fn_Gao_1d4GHzXS)

def crossmatch_3GHz(fn_dict):

    fn_3GHzlp               = fn_dict['fn_VLA_3GHz_lp']
    fn_3GHzAGN              = fn_dict['fn_VLA_3GHz_AGN']
    fn_1d4GHzdp             = fn_dict['fn_VLA_1d4GHz_dp']
    fn_3GHzlp_AGN           = fn_dict['fn_match_VLA_3GHzlp_AGN']
    fn_3GHzlp_AGN_1d4GHzdp  = fn_dict['fn_match_VLA_3GHzlp_AGN_1d4GHzdp']
    fn_Ugne         = fn_dict['fn_match_Ugne_Lim_1d4GHzXS']
    fn_Lim          = fn_dict['fn_match_Lim_1d4GHzXS']
    fn_Gao          = fn_dict['fn_match_Gao_Lim_1d4GHzXS']
    fn_Ugne_3GHz    = fn_dict['fn_match_Ugne_Lim_1d4GHzXS_3GHzlp']
    fn_Lim_3GHz     = fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp']
    fn_Gao_3GHz     = fn_dict['fn_match_Gao_Lim_1d4GHzXS_3GHzlp']


    input_dict = {} # filename, column_ra, column_dec
    input_dict['VLA_3GHz_lp']               = [fn_3GHzlp]
    input_dict['VLA_3GHz_AGN']              = [fn_3GHzAGN]
    input_dict['VLA_3GHz_lpAGN']            = [fn_3GHzlp_AGN,           'ra_3GHz',      'dec_3GHz']
    input_dict['VLA_1d4GHz_dp']             = [fn_1d4GHzdp,             'ra_1d4GHz',    'dec_1d4GHz']
    input_dict['VLA_3GHzlpAGN_1d4GHz']      = [fn_3GHzlp_AGN_1d4GHzdp,  'ra_3GHz',      'dec_3GHz']
    input_dict['match_Ugne_Lim_1d4GHzXS']   = [fn_Ugne,                 'RA_450Lim',    'DEC_450Lim']
    input_dict['match_Lim_1d4GHzXS']        = [fn_Lim,                  'RA_450Lim',    'DEC_450Lim']
    input_dict['match_Gao_Lim_1d4GHzXS']    = [fn_Gao,                  'RA_450Gao',    'Dec_450Gao']
    
    do_merge_value(input_dict, f1_key= 'VLA_3GHz_lp', f2_key = 'VLA_3GHz_AGN', 
        f1_value='id_3GHz', f2_value='ID_VLA_3GHz_multiAGN', join_type = 'outer', fn_match = fn_3GHzlp_AGN)
    do_merge_radec(input_dict, f1_key= 'VLA_3GHz_lpAGN', f2_key = 'VLA_1d4GHz_dp', radius=1, fn_bestmatch = fn_3GHzlp_AGN_1d4GHzdp)
    do_merge_radec(input_dict, f1_key= 'match_Ugne_Lim_1d4GHzXS', f2_key = 'VLA_3GHzlpAGN_1d4GHz', radius=4, fn_bestmatch = fn_Ugne_3GHz)
    do_merge_radec(input_dict, f1_key= 'match_Lim_1d4GHzXS', f2_key = 'VLA_3GHzlpAGN_1d4GHz', radius=4, fn_bestmatch = fn_Lim_3GHz)
    do_merge_radec(input_dict, f1_key= 'match_Gao_Lim_1d4GHzXS',  f2_key = 'VLA_3GHzlpAGN_1d4GHz', radius=4, fn_bestmatch = fn_Gao_3GHz)
    set_RadioDet(fn_Ugne_3GHz, fn_Ugne_3GHz)
    set_RadioDet(fn_Lim_3GHz,  fn_Lim_3GHz)
    set_RadioDet(fn_Gao_3GHz,  fn_Gao_3GHz)


def set_RadioDet(fn_in, fn_out):

    df = pd.read_csv(fn_in)
    df['ra_radio']      = np.where(df['ra_3GHz'].notna(), df['ra_3GHz'], df['RA_1d4GHz'])
    df['ra_err_radio']  = np.where(df['ra_3GHz'].notna(), df['ra_err_3GHz'], df['E_RA_1d4GHz'])
    df['dec_radio']     = np.where(df['ra_3GHz'].notna(), df['dec_3GHz'], df['DEC_1d4GHz'])
    df['dec_err_radio'] = np.where(df['ra_3GHz'].notna(), df['dec_err_3GHz'], df['E_DEC_1d4GHz'])
    df['radio_det']     = np.where(df['ra_3GHz'].notna(), '3',
                                   np.where(df['RA_1d4GHz'].notna(), '1.4', np.nan))

    if 'Unnamed: 0' in df.columns:
        df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)

    # save csv
    df.to_csv(fn_out, index=True, header=True)
    print('Save matched catalog %s'%(fn_out))


def crossmatch_MIPS(fn_dict):
    fn_MIPSLeFloch  = fn_dict['fn_MIPS_LeFloch']
    fn_MIPSwhwang   = fn_dict['fn_MIPS_whwang']
    fn_Ugne         = fn_dict['fn_match_Ugne_Lim_1d4GHzXS_3GHzlp']
    fn_Lim          = fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp']
    fn_Gao          = fn_dict['fn_match_Gao_Lim_1d4GHzXS_3GHzlp']
    fn_Ugne_MIPS    = fn_dict['fn_match_Ugne_Lim_1d4GHzXS_3GHzlp_mips']
    fn_Lim_MIPS     = fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp_mips']
    fn_Gao_MIPS     = fn_dict['fn_match_Gao_Lim_1d4GHzXS_3GHzlp_mips']

    input_dict = {} # filename, column_ra, column_dec
    input_dict['MIPS_LeFloch']                              = [fn_MIPSLeFloch,  'RA_mips24_LeFloch',    'DEC_mips24_LeFloch']
    input_dict['MIPS_whwang']                               = [fn_MIPSwhwang,   'RA_mips24_whw',        'DEC_mips24_whw']
    input_dict['match_Ugne_Lim_1d4GHzXS_3GHzlp_for_MIPS']   = [fn_Ugne,         'RA_450Lim',            'DEC_450Lim']
    input_dict['match_Lim_1d4GHzXS_3GHzlp_for_MIPS']        = [fn_Lim,          'RA_450Lim',            'DEC_450Lim']
    input_dict['match_Gao_Lim_1d4GHzXS_3GHzlp_for_MIPS']    = [fn_Gao,          'RA_450Gao',            'Dec_450Gao']
    do_merge_radec(input_dict, f1_key= 'match_Ugne_Lim_1d4GHzXS_3GHzlp_for_MIPS',f2_key = 'MIPS_LeFloch', radius=4, fn_bestmatch = fn_Ugne_MIPS)
    do_merge_radec(input_dict, f1_key= 'match_Lim_1d4GHzXS_3GHzlp_for_MIPS',     f2_key = 'MIPS_LeFloch', radius=4, fn_bestmatch = fn_Lim_MIPS)
    do_merge_radec(input_dict, f1_key= 'match_Gao_Lim_1d4GHzXS_3GHzlp_for_MIPS', f2_key = 'MIPS_LeFloch', radius=4, fn_bestmatch = fn_Gao_MIPS)
    set_RadioNotDet_MipsDet(fn_Ugne_MIPS, fn_Ugne_MIPS)
    set_RadioNotDet_MipsDet(fn_Lim_MIPS,  fn_Lim_MIPS)
    set_RadioNotDet_MipsDet(fn_Gao_MIPS,  fn_Gao_MIPS)

def set_RadioNotDet_MipsDet(fn_in, fn_out):

    df = pd.read_csv(fn_in)
    df['ra_nradio_mips']       = np.where(df['ra_radio'].notna(), np.nan, df['RA_mips24_LeFloch'])
    df['dec_nradio_mips']      = np.where(df['ra_radio'].notna(), np.nan, df['DEC_mips24_LeFloch'])
    if 'Unnamed: 0' in df.columns:
        df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)

    # save csv
    df.to_csv(fn_out, index=True, header=True)
    print('Save matched catalog %s'%(fn_out))

def crossmatch_IRAC(fn_dict):
    
    fn_IRAC             = fn_dict['fn_IRAC']
    # radio detected
    fn_Ugne             = fn_dict['fn_match_Ugne_Lim_1d4GHzXS_3GHzlp']
    fn_Lim              = fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp']
    fn_Gao              = fn_dict['fn_match_Gao_Lim_1d4GHzXS_3GHzlp']
    fn_Ugne_IRAC        = fn_dict['fn_match_Ugne_Lim_1d4GHzXS_3GHzlp_irac']
    fn_Lim_IRAC         = fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp_irac']
    fn_Gao_IRAC         = fn_dict['fn_match_Gao_Lim_1d4GHzXS_3GHzlp_irac']    
    # radio non-detected
    fn_Ugne_MIPS        = fn_dict['fn_match_Ugne_Lim_1d4GHzXS_3GHzlp_mips']
    fn_Lim_MIPS         = fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp_mips']
    fn_Gao_MIPS         = fn_dict['fn_match_Gao_Lim_1d4GHzXS_3GHzlp_mips']
    fn_Ugne_MIPS_IRAC   = fn_dict['fn_match_Ugne_Lim_1d4GHzXS_3GHzlp_mips_irac']
    fn_Lim_MIPS_IRAC    = fn_dict['fn_match_Lim_1d4GHzXS_3GHzlp_mips_irac']
    fn_Gao_MIPS_IRAC    = fn_dict['fn_match_Gao_Lim_1d4GHzXS_3GHzlp_mips_irac']


    input_dict = {} # filename, column_ra, column_dec
    input_dict['IRAC']                                          = [fn_IRAC,         'ra_irac',          'dec_irac']
    input_dict['match_Ugne_Lim_1d4GHzXS_3GHzlp_for_IRAC']       = [fn_Ugne,         'ra_radio',         'dec_radio']
    input_dict['match_Lim_1d4GHzXS_3GHzlp_for_IRAC']            = [fn_Lim,          'ra_radio',         'dec_radio']
    input_dict['match_Gao_Lim_1d4GHzXS_3GHzlp_for_IRAC']        = [fn_Gao,          'ra_radio',         'dec_radio']
    input_dict['match_Ugne_Lim_1d4GHzXS_3GHzlp_MIPS_for_IRAC']  = [fn_Ugne_MIPS,    'ra_nradio_mips',   'dec_nradio_mips']
    input_dict['match_Lim_1d4GHzXS_3GHzlp_MIPS_for_IRAC']       = [fn_Lim_MIPS,     'ra_nradio_mips',   'dec_nradio_mips']
    input_dict['match_Gao_Lim_1d4GHzXS_3GHzlp_MIPS_for_IRAC']   = [fn_Gao_MIPS,     'ra_nradio_mips',   'dec_nradio_mips']
    do_merge_radec(input_dict, f1_key= 'match_Ugne_Lim_1d4GHzXS_3GHzlp_for_IRAC',      f2_key = 'IRAC', radius=1, fn_bestmatch = fn_Ugne_IRAC)
    do_merge_radec(input_dict, f1_key= 'match_Lim_1d4GHzXS_3GHzlp_for_IRAC',           f2_key = 'IRAC', radius=1, fn_bestmatch = fn_Lim_IRAC)
    do_merge_radec(input_dict, f1_key= 'match_Gao_Lim_1d4GHzXS_3GHzlp_for_IRAC',       f2_key = 'IRAC', radius=1, fn_bestmatch = fn_Gao_IRAC)
    do_merge_radec(input_dict, f1_key= 'match_Ugne_Lim_1d4GHzXS_3GHzlp_MIPS_for_IRAC', f2_key = 'IRAC', radius=2, fn_bestmatch = fn_Ugne_MIPS_IRAC)
    do_merge_radec(input_dict, f1_key= 'match_Lim_1d4GHzXS_3GHzlp_MIPS_for_IRAC',      f2_key = 'IRAC', radius=2, fn_bestmatch = fn_Lim_MIPS_IRAC)
    do_merge_radec(input_dict, f1_key= 'match_Gao_Lim_1d4GHzXS_3GHzlp_MIPS_for_IRAC',  f2_key = 'IRAC', radius=2, fn_bestmatch = fn_Gao_MIPS_IRAC)


def crossmatch_1d4GHz_diffR(fn_dict):  
    
    fn_1d4GHzXS = fn_dict['fn_VLA_1d4GHz_XS']
    fn_Ugne     = fn_dict['fn_match_Ugne_Lim']
    fn_Lim      = fn_dict['fn_JCMT_450um_Lim']
    fn_Gao      = fn_dict['fn_match_Gao_Lim']

    input_dict = {} # filename, column_ra, column_dec
    input_dict['VLA_1d4GHz_XS']      = [fn_1d4GHzXS,    'RA_1d4GHz',    'DEC_1d4GHz']
    input_dict['match_Ugne_Lim']     = [fn_Ugne,        'RA_450Lim',    'DEC_450Lim']
    input_dict['JCMT_450um_Lim']     = [fn_Lim,         'RA_450Lim',    'DEC_450Lim']
    input_dict['match_Gao_Lim']      = [fn_Gao,         'RA_450Gao',    'Dec_450Gao']

    # cross match with VLA 1.4 GHz
    r_start = 0.5
    r_end   = 7.5 
    r_step  = 0.5
    radius_arr = np.arange(r_start, r_end+r_step, r_step)

    def set_new_filename(fn_in, r):
        
        # path = '/'.join(fn_in.split('/')[:-1]) + '/'
        path = PATH_CATALOG_CROSSMATCH
        path_new = path + 'diffR_450um_1d4GHzXS/'
        if not os.path.exists(path_new):
            os.makedirs(path_new)
        fn = fn_in.split('/')[-1]
        fn_new = path_new + fn.split('.csv')[0] + '_r%.1f'%(r) + '_1d4GHzXS' + '.csv'

        return fn_new

    for r in radius_arr:

        fn_Ugne_1d4GHzXS = set_new_filename(fn_Ugne,r)
        fn_Lim_1d4GHzXS  = set_new_filename(fn_Lim,r)
        fn_Gao_1d4GHzXS  = set_new_filename(fn_Gao,r)

        do_merge_radec(input_dict, f1_key= 'match_Ugne_Lim', f2_key = 'VLA_1d4GHz_XS', radius=r, fn_bestmatch = fn_Ugne_1d4GHzXS)
        do_merge_radec(input_dict, f1_key= 'JCMT_450um_Lim', f2_key = 'VLA_1d4GHz_XS', radius=r, fn_bestmatch = fn_Lim_1d4GHzXS)
        do_merge_radec(input_dict, f1_key= 'match_Gao_Lim',  f2_key = 'VLA_1d4GHz_XS', radius=r, fn_bestmatch = fn_Gao_1d4GHzXS)

    
def do_merge_value(input_dict, f1_key, f2_key, f1_value, f2_value, join_type='inner', fn_match=None):
    # init parm
    fn_f1 = input_dict[f1_key][0]
    fn_f2 = input_dict[f2_key][0]

    print('-------------------------------------------------')
    print('Start to merge catalog %s and %s'%(f1_key, f2_key))
    df1 = pd.read_csv(fn_f1)
    df2 = pd.read_csv(fn_f2)

    df_m = df1.merge(df2, left_on=f1_value, right_on=f2_value, how=join_type, suffixes=['_1', '_2'])
    print('The match number is: %d'%(df_m.shape[0]))
    if 'Unnamed: 0' in df_m.columns:
        df_m.drop(df_m.columns[df_m.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)

    # save csv
    if fn_match is not None:
        df_m.to_csv(fn_match, index=True, header=True)
        print('Save matched catalog %s'%(fn_match))
 

def do_merge_radec(input_dict, f1_key, f2_key, radius=1, 
        fn_allmatch = None, fn_bestmatch = None, fn_innerjoin = None):

    # init parm
    fn_f1, ra_f1, dec_f1 = input_dict[f1_key]
    fn_f2, ra_f2, dec_f2 = input_dict[f2_key]

    print('-------------------------------------------------')
    print('Start to merge catalog %s and %s'%(f1_key, f2_key))
    df1 = pd.read_csv(fn_f1)
    df2 = pd.read_csv(fn_f2)
    ind, dist = tc.crossmatch_angular(  X1 = df1[[ra_f1, dec_f1]].values, 
                                        X2 = df2[[ra_f2, dec_f2]].values,
                                        max_distance = radius/3600, return_dist = True)
    df_m    = tc.merge_df(df1, df2, ind, dist)      # all match
    df_bm   = tc.best_match(df_m, df_m.columns[0])  # best match
    df_bmi  = tc.inner_join(df_bm)                  # inner join
    print('The all match number is: %d'%(df_m.shape[0]))
    print('The best match number is: %d'%(df_bm.shape[0]))
    print('The inner join best match number is: %d'%(df_bmi.shape[0]))

    # save csv
    if fn_allmatch is not None:
        df_m.to_csv(fn_allmatch, index=True, header=True)
        print('Save all-match catalog %s'%(fn_allmatch))
    if fn_bestmatch is not None:
        df_bm.to_csv(fn_bestmatch, index=True, header=True)
        print('Save best-match catalog %s'%(fn_bestmatch))
    if fn_innerjoin is not None:
        df_bmi.to_csv(fn_innerjoin, index=True, header=True)
        print('Save best-match, inner join catalog %s'%(fn_innerjoin))
    

if __name__ == '__main__':
    main()