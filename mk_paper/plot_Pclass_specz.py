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

from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=9, usetex=True)

SAGA_DIR   = os.environ['SAGA_DIR']

############################################################################
def plot_Pclass_specz():


	# READ SPECTRA
	file = SAGA_DIR +'/data/saga_spectra_dirty.fits.gz'
	allspec = Table.read(file)

	sats = allspec['SATS'] == 1

	# FIND UNIQUE SAGANAMES AND ORDER BY NSATS
	msk1          = allspec['ZQUALITY'] >= 3
	msk2          = allspec['r']  > 17.7


	msk  = msk1 & msk2
	sagaspec     = allspec[msk]
	ssats = sagaspec['SATS'] == 1



#	plt.subplot(2,1,1)
	fig = plt.figure(figsize=(7, 3))
	fig.subplots_adjust(bottom=0.15, top=0.9,left=0.12, right=0.95, wspace=0.3)

	ax = fig.add_subplot(121)
	ax.plot(allspec['PROBABILITY_CLASS1'],allspec['SPEC_Z'],'.',markersize=0.8)
	ax.plot(allspec['PROBABILITY_CLASS1'][sats],allspec['SPEC_Z'][sats],'o',markersize=1)
	ax.set_xlabel('Pclass')
	ax.set_ylabel('Spec z')
	ax.set_ylim(1e-3,5)
	ax.set_yscale('log')
	ax.set_xscale('log')


	ax = fig.add_subplot(122)
	ax.plot(sagaspec['PROBABILITY_CLASS1'],sagaspec['SPEC_Z'],'.',markersize=0.8)
	ax.plot(sagaspec['PROBABILITY_CLASS1'][ssats],sagaspec['SPEC_Z'][ssats],'o',markersize=1)
	ax.set_xlabel('Pclass')
	ax.set_ylabel('Spec z')
	ax.set_ylim(1e-3,5)
	ax.set_yscale('log')
	ax.set_xscale('log')

	plt.show()
	plt.savefig('Pclass.v.specz.pdf')

if __name__ == '__main__':
    create_saga_spectra()

 

 
