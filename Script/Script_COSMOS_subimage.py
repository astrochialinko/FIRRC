##########################################################
# 2020.Oct.05
# Using CASA version 5.4.0
#
# Purpose:
#     Create a (sub)image from a region of the image (3 GHz VLA) 
#     
#
# Editor:
#   Chia-Lin Ko
##########################################################
import os
import numpy as np
##########################################################


###### Work Flow #########################################
thesteps = []
step_title = {
               0: 'Import FITS images',
               1: 'Create subimages',
               2: 'Convert unit from Jy/beam to uJy/beam',
               3: 'Change the header',
               4: 'Export FITS images'
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
outnames        = [
                   'vla_3ghz_msmf',
                   'vla_3ghz_msmf.rms'
                  ]

parm_sub        = 'detopt'
parm_unit       = 'uJy'
path_data       = '../Data/COSMOS/Image/'
path_data_orgin = '../Original_data/COSMOS/Image/VLA/'
##########################################################


##### Import FITS images #################################
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i, outname in enumerate(outnames):

    Iimagename   = '%s%s'%(path_data_orgin, outname)

    outimages = [
                    Iimagename,
                ]

    for outcasaimage in outimages:
        importfits(
                   fitsimage = outcasaimage+'.fits',
                   imagename = outcasaimage+'.image',
                   overwrite = True
                     )

##########################################################


##### Creat subimages ####################################
mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i, outname in enumerate(outnames):

    Imagename    = '%s%s.image'%(path_data_orgin, outname)
    SubImagename = '%s%s_sub_%s.image'%(path_data, outname, parm_sub)
    region       = 'box[[12000pix, 15500pix], [17500pix, 21000pix]]'

    # Creat a subimage
    command = 'rm -rf ' + SubImagename
    os.system(command)
    
    imsubimage(
                imagename=Imagename,
                outfile=SubImagename,
                region=region,
                overwrite=True
                )

    # Delect the Iimage (CASA format)
    command = 'rm -rf ' + Imagename
    os.system(command)

##########################################################


##### Convert Unit #######################################
mystep = 2
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i, outname in enumerate(outnames):
    SubImagename     = '%s%s_sub_%s.image'%(path_data, outname, parm_sub)
    SubUnitImagename = '%s%s_sub_%s_%s.image'%(path_data, outname, parm_sub, parm_unit)
    expr             = 'IM0*1e6'            # convert Jy/beam to uJy/beam
    
    # Creat a subimage
    command = 'rm -rf ' + SubUnitImagename
    os.system(command) 

    immath(
            imagename=[SubImagename],
            mode='evalexpr',
            expr=expr,
            outfile=SubUnitImagename
            )

##########################################################


##### Change header ######################################
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i, outname in enumerate(outnames):
    SubUnitImagename = '%s%s_sub_%s_%s.image'%(path_data, outname, parm_sub, parm_unit)
    hdkey   = 'bunit'
    hdvalue = 'uJy/beam'

    imhead(
            imagename=SubUnitImagename,
            mode='put',
            hdkey=hdkey,
            hdvalue=hdvalue
            )


##########################################################


##########################################################


##### Export FITS images #################################
mystep = 4
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i, outname in enumerate(outnames):

    SubUnitImagename   = '%s%s_sub_%s_%s'%(path_data, outname, parm_sub, parm_unit)

    outimages = [
                    SubUnitImagename,
                ]


    for outcasaimage in outimages:
        exportfits(
                   imagename = outcasaimage+'.image',
                   fitsimage = outcasaimage+'.fits',
                   overwrite = True
                     )

    # Delect the Iimage (CASA format)
    # command = 'rm -rf ' + SubIimagename + '.image'
    # os.system(command)

##########################################################

