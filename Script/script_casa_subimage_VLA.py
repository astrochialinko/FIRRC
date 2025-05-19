"""
File: script_casa_subimage_vla.py
Name: Chia-Lin Ko
Create Date: Aug 11, 2021
Last Modified Date: Apr 11, 2021
------------------------
This program aims to create a (sub)image from a region of the image
"""
import os
import numpy as np
##########################################################


###### Work Flow #########################################
thesteps = []
step_title = {
               0: 'Create subimages',
               1: 'Export FITS images'
              }

try:
  print 'List of steps to be executed ...', mysteps
  thesteps = mysteps
except:
  print 'global variable mysteps not set.'
if (thesteps==[]):
  thesteps = range(0,len(step_title))
  print 'Executing all steps: ', thesteps
##########################################################



##### Global paramters and keywords ######################
# list
IMG_FILENAME_LIST     = [
  'cosmos_vla_3GHz_2017image_uJy',
  'cosmos_vla_3GHz_2017rms_uJy',
  'cosmos_vla_1d4GHz_XS_2021image_uJy',
  'cosmos_vla_1d4GHz_XS_2021rms_uJy'
                        ]
CATALOG_FILENAME_LIST = [ 
  'Coord_pix_3GHz_slt_nAGN_both_1d4_3GHz',
  #'Coord_pix_3GHz_slt_nAGN_both_1d4_3GHz',
  'Coord_pix_1d4GHz_dp_slt_nAGN_both_1d4_3GHz',
  #'Coord_pix_1d4GHz_lp_slt_nAGN_both_1d4_3GHz'
                        ]

# path
PATH_DATA_IN    = '../Data/COSMOS/Image/VLA/'
PATH_DATA_OUT   = '../Data/COSMOS/Image/VLA/RadioSFG_sub/'
PATH_TABLE      = '../Data/Tables/'

# parm
SUB_WIDTH       = 151 # [pixel]

##########################################################


##### Create subimages ####################################
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i, img_filename in enumerate(IMG_FILENAME_LIST):

    cat_filename  = '%s%s.txt'%(PATH_TABLE, CATALOG_FILENAME_LIST[i])
    radec_list    = np.loadtxt(cat_filename, dtype = str, delimiter=' ')
    lines_num, item_num = np.shape(radec_list)
    print 'Total', lines_num, 'subimages'

    for j in range(lines_num):
      ra      = int(radec_list[j][0])-1    # ra        (pix)
      dec     = int(radec_list[j][1])-1    # dec       (pix)

      Imagename     = '%s%s'%(PATH_DATA_IN, img_filename)
      SubImagename  = '%s%s_sub_%s'%(PATH_DATA_OUT, img_filename, j)
      region        = 'centerbox[[%spix, %spix], [%spix, %spix]]'%(ra, dec, SUB_WIDTH, SUB_WIDTH)
      print 'region = ', region

      # Creat a subimage
      imsubimage(
                  imagename=Imagename+'.image',
                  outfile=SubImagename+'.image',
                  region=region,
                  overwrite=True
                  )

##########################################################


##### Export FITS images #################################
mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i, img_filename in enumerate(IMG_FILENAME_LIST):
    for j in range(lines_num):

      SubImagename  = '%s%s_sub_%s'%(PATH_DATA_OUT, img_filename, j)

      outimages = [
                      SubImagename,
                  ]

      for outcasaimage in outimages:
          exportfits(
                     imagename = outcasaimage+'.image',
                     fitsimage = outcasaimage+'.fits',
                     overwrite = True
                       )

##########################################################

