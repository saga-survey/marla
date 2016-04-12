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

from astroML.plotting import scatter_contour
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=9, usetex=True)

SAGA_DIR   = os.environ['SAGA_DIR']

############################################################################
def plot_Pclass_specz():


	# READ SPECTRA
	file = SAGA_DIR +'/data/saga_spectra_clean.fits.gz'
	allspec = Table.read(file)



	lowz1 = allspec['SATS'] == 2
	lowz2 = allspec['SPEC_Z'] < 0.02
	lowz=lowz1&lowz2

	sats = allspec['SATS'] == 1

	b = allspec['r'] < 17.7
	f = allspec['r'] >= 17.7

	# FIND UNIQUE SAGANAMES AND ORDER BY NSATS
	

	fig = plt.figure(figsize=(7, 5))
	fig.subplots_adjust(left=0.1, right=0.95, wspace=0.05,
                    bottom=0.1, top=0.95, hspace=0.05)


	thr=500
	##############################
	ax = fig.add_subplot(2, 1, 1)
#	scatter_contour(allspec['PROBABILITY_CLASS1'],allspec['SPEC_Z'], log_counts=True, ax=ax,
#                histogram2d_args=dict(bins=20),
#                plot_args=dict(marker=',', linestyle='none', color='black'),
#                contour_args=dict(cmap=plt.cm.bone))

	ax.plot(allspec['PROBABILITY_CLASS1'][b],allspec['SPEC_Z'][b],'.k',markersize=2,label='All Galaxies')
	ax.plot(allspec['PROBABILITY_CLASS1'][b&lowz],allspec['SPEC_Z'][b&lowz],'.b',markersize=5,label='$0.001 < z < 0.02$')
	ax.plot(allspec['PROBABILITY_CLASS1'][b&sats],allspec['SPEC_Z'][b&sats],'.r',markersize=7,label = 'Satellites')


	ax.set_xlabel('Pclass')
	ax.set_ylabel('Spec z')
	ax.set_ylim(1e-3,1.2)
	ax.set_xlim(1e-10,1.1)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.get_xaxis().set_visible(False)
	ax.text(2e-10,2e-3,'Bright Galaxies, r $<$ 17.7')


	plt.legend(loc=1,fontsize=9)


	##############################
	ax = fig.add_subplot(2, 1, 2)
	ax.plot(allspec['PROBABILITY_CLASS1'][f],allspec['SPEC_Z'][f],'.k',markersize=2,label='All Galaxies')

	ax.plot(allspec['PROBABILITY_CLASS1'][f&lowz],allspec['SPEC_Z'][f&lowz],'.b',markersize=5,label='$0.001 < z < 0.02$')
	ax.plot(allspec['PROBABILITY_CLASS1'][f&sats],allspec['SPEC_Z'][f&sats],'.r',markersize=7,label = 'Satellites')

	ax.set_xlabel('Pclass')
	ax.set_ylabel('Spec z')
	ax.set_ylim(1e-3,1.2)
	ax.set_xlim(1e-10,1.1)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.text(2e-10,2e-3,'Faint Galaxies, r $>$ 17.7')


	plt.show()
	plt.savefig('fig_Pclass.pdf')

if __name__ == '__main__':
    create_saga_spectra()

 

 
