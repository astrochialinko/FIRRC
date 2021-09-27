#!/usr/local/anaconda3/envs/astro37/bin/python
"""
File: script_00_remodel_catalog.py
Name: Chia-Lin Ko
Date: Jul 23, 2021
Last Modified Date: Jul 23, 2021
------------------------
This program aims to rename and convert the original catalog 
from the archive websites to csv and fits format for analysis.
"""
import os
from astropy.table import Table

# my own packing
import path

# Constant
# path
PATH_CATALOG_ORG_VLA    = path.catalog_org_vla

PATH_CATALOG_JCMT       = path.catalog_jcmt
PATH_CATALOG_VLA        = path.catalog_vla
PATH_CATALOG_IRAC       = path.catalog_irac
PATH_CATALOG_MIPS       = path.catalog_mips
PATH_CATALOG_CROSSMATCH = path.catalog_crossmatch


def main():


    # VLA 3 GHz COSMOS-XS Survey
    remodel_table(
        fn_in = '%sapjabb77at6_mrt.txt'%(PATH_CATALOG_ORG_VLA+'COSMOS_XS/'),
        fn_out = '%scosmos_vla_3GHz_10GHz_XS_2021cat'%(PATH_CATALOG_VLA),
        added_col_name='_3GHzXS',
        is_save_csv=True,
        is_save_fits=True
        )    


def remodel_table(fn_in, fn_out, added_col_name=None, is_save_csv=True, is_save_fits=True):
    t = Table.read(fn_in, format="ascii.cds")
    if added_col_name is not None:
        rename_column(t, added_col_name)
    if is_save_csv:
        table2csv(t, fn_out+'.csv')
    if is_save_fits:
        table2fits(t, fn_out+'.fits')

def rename_column(astro_table, added_col_name):
    col_names = astro_table.colnames
    for col_name in col_names:
        astro_table.rename_column(col_name, col_name+added_col_name)

def table2csv(astro_table, fn_out):
    astro_table.write(fn_out, format='csv', overwrite=True)

def table2fits(astro_table, fn_out):
    astro_table.write(fn_out, format='fits', overwrite=True)

def table2latex(astro_table, fn_out):
    astro_table.write(fn_out, format='latex', overwrite=True)


if __name__ == '__main__':
    main()