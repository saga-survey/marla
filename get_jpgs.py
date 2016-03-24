#! /usr/bin/python
##################################################
#  COMPILE SPECTRA FROM ALL SOURCES FOR ALL FLAG 0 HOSTS
#  GAMA:  
#  MMT: 
#  AAT:
#  IMACS:
#  WIYN:  
#
##################################################

import numpy as np

import os
import glob
import pyspherematch as sm
from astropy.table import Table


import urllib
from textwrap import dedent

SAGA_DIR   = os.environ['SAGA_DIR']



############################################################################
def get_jpg(ra,dec,outname):

   # ULR REQUEST FOR SDSS IMAGE CUTOUT
        sdssurl = 'http://skyservice.pha.jhu.edu/DR10/ImgCutout/getjpeg.aspx?ra='+\
                                str(ra)+'&dec='+str(dec)+\
                                '&width=256&height=256&scale=0.4'#&opt=G'

        urllib.urlretrieve(sdssurl, outname)
        return sdssurl

############################################################################
def get_satellites():

	# READ SPECTRA
	file = SAGA_DIR +'/data/saga_spectra_dirty.fits.gz'
	allspec = Table.read(file)


	msk1   = allspec['HOST_SAGA_NAME'] != '' 
	msk2    = allspec['SATS'] == 1
	sats = allspec[msk1&msk2]

	for s in sats:

		outjpg = '{}_{:.1f}.jpg'.format(s['HOST_SAGA_NAME'], s['r'])
		outjpg = str(s['HOST_SAGA_NAME'])+ '_' +str(s['r']).format()  +'.jpg'
		print outjpg,s['RA'],s['DEC']
		get_jpg(s['RA'],s['DEC'],outjpg)



