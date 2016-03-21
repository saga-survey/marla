#! /usr/bin/python
##################################################
#  COMPILE SPECTRA FROM ALL SOURCES FOR ALL FLAG 0 HOSTS
#       1. SDSS 
# 		2. GAMA
#       3. MMT, AAT, IMACS, WIYN 
#
#
#
##################################################

import numpy as np
import os

import matplotlib.pyplot as plt


from astropy import table
from astropy.table import Table

from matplotlib import rc
font = {'family':'serif','size':12}
rc('font',**font)

SAGA_DIR   = os.environ['SAGA_DIR']

############################################################################
def plot_Pclass_specz():


	# READ SPECTRA
	file = SAGA_DIR +'/data/saga_spectra_dirty.fits.gz'
	allspec = Table.read(file)


	# FIND UNIQUE SAGANAMES AND ORDER BY NSATS
	msk1          = allspec['ZQUALITY'] >= 3
	msk2          = allspec['r']  > 17.7


	msk  = msk1 & msk2
	sagaspec     = allspec[msk]



#	plt.subplot(2,1,1)
	plt.plot(sagaspec['PROBABILITY_CLASS1'],sagaspec['SPEC_Z'],'.',markersize=0.8)
	plt.xlabel('Pclass')
	plt.ylabel('Spec z')
	plt.ylim(1e-3,5)
	plt.yscale('log')
	plt.xscale('log')

	plt.show()
	plt.savefig('myfig.pdf')

if __name__ == '__main__':
    create_saga_spectra()

 

 
