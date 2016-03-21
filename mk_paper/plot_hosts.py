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
from astropy.io import ascii
from FileLoader import GoogleSheets, FitsTable


from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=9, usetex=True)

SAGANAMES  = GoogleSheets('1GJYuhqfKeuJr-IyyGF_NDLb_ezL6zBiX2aeZFHHPr_s', 0, header_start=0)


SAGA_DIR   = os.environ['SAGA_DIR']
SAGA_DROPBOX= os.environ['SAGA_DROPBOX']


############################################################################
def plot_hosts():


	file = SAGA_DROPBOX + 'hosts/host_catalog_all.csv'
	all = ascii.read(file)
	ri_all = all['r'] - all['i']

	file = SAGA_DROPBOX + 'hosts/host_catalog_flag0.csv'
	hosts = ascii.read(file)

	names = SAGANAMES.load()
	m = hosts['NSAID']  & 0
	for n in names:
		msaga = hosts['NSAID'] == n['NSA']
		print np.sum(m)
		m = np.logical_or(m, msaga)
	print hosts['NSAID'][m]


	ri = hosts['r'] - hosts['i']



	fig = plt.figure(figsize=(7, 3))
	fig.subplots_adjust(bottom=0.15, top=0.9,
                    left=0.12, right=0.95, wspace=0.3)


	# first plot the flux distribution
	ax = fig.add_subplot(121)
	ax.plot(all['K_abs'],ri_all, 'o')
	ax.plot(hosts['K_abs'],ri, 'o')
	ax.plot(hosts['K_abs'][m],ri[m], 'ro')

	ax.set_ylim(-0.1, 0.8)
	ax.set_xlim(-24.9,-22.5)
	ax.set_ylabel('(r-i)')
	ax.set_xlabel('$M_K$')
	ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))


# second plot the flux distribution
	ax = fig.add_subplot(122)
	ax.plot(all['r_abs'],ri_all, 'o',label='All Galaxies')
	ax.plot(hosts['r_abs'],ri, 'o',label='MW Hosts')
	ax.plot(hosts['r_abs'][m],ri[m], 'ro',label='This Paper')

	#ax.set_xlim(-0.1, 2.1)
	ax.set_ylim(-0.1,0.8)
	ax.set_ylabel('(r-i)')
	ax.set_xlabel('$M_r$')

	#ax.text(0.04, 0.98, r'${\rm 20\%\ flux\ error}$',
    #    ha='left', va='top', transform=ax.transAxes,
    #    bbox=dict(ec='none', fc='w'))
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels)


#	plt.show()
	plt.savefig('fig_hosts.pdf')

if __name__ == '__main__':
    create_saga_spectra()

 

 
