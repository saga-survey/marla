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


	file = SAGA_DROPBOX + 'hosts/masterlist.csv'
	all = ascii.read(file)
	ri_all = all['r'] - all['i']
	dmod  = 5*np.log10(1e6*all['distance']) - 5.
	K_abs = all['K'] - dmod
	r_abs = all['r'] - dmod

	file = SAGA_DROPBOX + 'hosts/host_catalog_flag0.csv'
	hosts = ascii.read(file)

	names = SAGANAMES.load()
	m = hosts['NSAID']  & 0
	for n in names:
		msaga = hosts['NSAID'] == n['NSA']
		m = np.logical_or(m, msaga)
	print hosts['NSAID'][m]
	print "number of blue hosts = ",np.sum(m)

	ri = hosts['r'] - hosts['i']



	fig = plt.figure(figsize=(8, 3))
	fig.subplots_adjust(bottom=0.15, top=0.9,
                    left=0.12, right=0.95, wspace=0.3)


	# first plot the flux distribution
	ax = fig.add_subplot(121)
	ax.plot(K_abs,ri_all, '.',markersize=4,color='grey')
	ax.plot(hosts['K_abs'],ri, 'bo')
	ax.plot(hosts['K_abs'][m],ri[m], 'ro')

	#MW
	ax.plot([-24.02],[0.33],'y*',markersize=14,color='#ffcc11')


	ax.set_ylim(0.1, 0.7)
	ax.set_xlim(-24.9,-22.5)
	ax.set_ylabel('(r-i)')
	ax.set_xlabel('$M_K$')
	ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))


# second plot the flux distribution
	ax = fig.add_subplot(122)
	ax.plot(r_abs,ri_all, '.',color='grey',markersize=4,label='All Galaxies')
	ax.plot(hosts['r_abs'],ri, 'bo',label='MW Hosts')
	ax.plot(hosts['r_abs'][m],ri[m], 'ro',label='This Paper')

	ax.plot([-20.2],[0.33],'y*',label='Milky Way',color='#ffcc11')
	ax.plot([-20.2],[0.33],'y*',markersize=14,color='#ffcc11')


	ax.set_xlim(-21.5,-19)
	ax.set_ylim(0.1,0.7)
	ax.set_ylabel('(r-i)')
	ax.set_xlabel('$M_r$')

	#ax.text(0.04, 0.98, r'${\rm 20\%\ flux\ error}$',
    #    ha='left', va='top', transform=ax.transAxes,
    #    bbox=dict(ec='none', fc='w'))
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels,fontsize=6)


#	plt.show()
	plt.savefig('fig_hosts.pdf')

if __name__ == '__main__':
    create_saga_spectra()

 

 
