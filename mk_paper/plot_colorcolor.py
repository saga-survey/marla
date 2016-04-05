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
from scipy import integrate


from astropy import table
from astropy.table import Table
from astropy import units as u
from astropy.io import ascii


from astroML.plotting import scatter_contour
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=9, usetex=True)


from utils import get_mcconn_table

SAGA_DIR   = os.environ['SAGA_DIR']

############################################################################
def plot_colorcolor():


# ADD MW SATELLITES?

	# READ SPECTRA
	file = SAGA_DIR +'/data/saga_spectra_dirty.fits.gz'
	allspec = Table.read(file)

	m1 = allspec['ZQUALITY'] > 2
	m2 = allspec['REMOVE'] == -1
	allspec = allspec[m1&m2]

#	msat = allspec['SATS'] == 1
	msat1 = allspec['SPEC_Z'] > 0.005
	msat2 = allspec['SPEC_Z'] < 0.02
	msat = msat1 & msat2
	mname1 = allspec['SATS'] == 1
	mname2 = allspec['HOST_SAGA_NAME'] != ''
	mname = mname1&mname2

	gr = allspec['g'] - allspec['r']
	ri = allspec['r'] - allspec['i']



	fmag = allspec['r'] > 17.7
	bmag = allspec['r'] <= 17.7

	gal = allspec['PHOTPTYPE']  == 3
	star = allspec['PHOTPTYPE'] == 6

	xl = [-0.3,2.0]
	yl=[-0.8,1.3]

	xl1 = gr > xl[0]
	xl2 = gr < xl[1]
	yl1 = ri > yl[0]
	yl2 = ri < yl[1]
	lims = xl1&xl2&yl1&yl2

	fig = plt.figure(figsize=(7, 7))
	fig.subplots_adjust(left=0.1, right=0.95, wspace=0.05,
                    bottom=0.1, top=0.95, hspace=0.05)


	thr = 60
	##############################
	ax = fig.add_subplot(2, 2, 1)
	scatter_contour(gr[gal&bmag&lims], ri[gal&bmag&lims], threshold=thr, log_counts=True, ax=ax,
                histogram2d_args=dict(bins=40),
                plot_args=dict(marker=',', linestyle='none', color='black'),
                contour_args=dict(cmap=plt.cm.bone))


	ax.plot(gr[gal&bmag&msat], ri[gal&bmag&msat],'b.',label = 'low-z')
	ax.plot(gr[gal&bmag&msat&mname], ri[gal&bmag&msat&mname],'r.',label='satellites')
	ax.set_xlim(xl)
	ax.set_ylim(yl)
	ax.text(xl[0]+0.02,yl[1]-0.2,'Bright Galaxies, r $<$ 17.7')
	ax.get_xaxis().set_visible(False)
	plt.axvline(0.8, c='k', ls=':')
	plt.axhline(0.5, c='k', ls=':')
	plt.legend(loc=4,fontsize=9)
	ax.set_ylabel('(r-i)')

	##############################
	ax = fig.add_subplot(2, 2, 2)
	scatter_contour(gr[gal&fmag&lims], ri[gal&fmag&lims], threshold=thr, log_counts=True, ax=ax,
                histogram2d_args=dict(bins=40),
                plot_args=dict(marker=',', linestyle='none', color='black'),
                contour_args=dict(cmap=plt.cm.bone))

	ax.plot(gr[gal&fmag&msat], ri[gal&fmag&msat],'b.')
	ax.plot(gr[gal&fmag&msat&mname], ri[gal&fmag&msat&mname],'r.')
	ax.set_xlim(xl)
	ax.set_ylim(yl)
	ax.text(xl[0]+0.02,yl[1]-0.2,'Faint Galaxies, r $>$ 17.7')
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_visible(False)
	plt.axvline(0.8, c='k', ls=':')
	plt.axhline(0.5, c='k', ls=':')


	##############################
	ax = fig.add_subplot(2, 2, 3)
	scatter_contour(gr[star&bmag&lims], ri[star&bmag&lims], threshold=thr, log_counts=True, ax=ax,
                histogram2d_args=dict(bins=40),
                plot_args=dict(marker=',', linestyle='none', color='black'),
                contour_args=dict(cmap=plt.cm.bone))

	ax.plot(gr[star&bmag&msat], ri[star&bmag&msat],'b.')
	ax.plot(gr[star&bmag&msat&mname], ri[star&bmag&msat&mname],'r.')
	ax.set_xlim(xl)
	ax.set_ylim(yl)
	ax.text(xl[0]+0.03,yl[1]-0.2,'Bright Stars, r $<$ 17.7')
	plt.axvline(0.8, c='k', ls=':')
	plt.axhline(0.5, c='k', ls=':')
	ax.set_xlabel('(g-r)')
	ax.set_ylabel('(r-i)')

	##############################
	ax = fig.add_subplot(2, 2, 4)
	scatter_contour(gr[star&fmag&lims], ri[star&fmag&lims], threshold=thr, log_counts=True, ax=ax,
                histogram2d_args=dict(bins=40),
                plot_args=dict(marker=',', linestyle='none', color='black'),
                contour_args=dict(cmap=plt.cm.bone))

	ax.plot(gr[star&fmag&msat], ri[star&fmag&msat],'b.')
	ax.plot(gr[star&fmag&msat&mname], ri[star&fmag&msat&mname],'r.')
	ax.set_xlim(xl)
	ax.set_ylim(yl)
	ax.text(xl[0]+0.03,yl[1]-0.2,'Faint Stars, r $>$ 17.7')
	ax.get_yaxis().set_visible(False)
	plt.axvline(0.8, c='k', ls=':')
	plt.axhline(0.5, c='k', ls=':')
	ax.set_xlabel('(g-r)')

	plt.show()
	plt.savefig('fig_colorcolor.pdf')


	# DO OBJECTS OUTSIDE OF STRAIGHT COLOR REQUIRES,
	# PASS CRITERIA INCLUDING ERRORS?

	lowz = allspec[msat]

	gr = lowz['g'] - lowz['r']
	ri = lowz['r'] - lowz['i']

	m1 = gr > 0.8
	m2 = ri > 0.5
	for obj in lowz[m1|m2]:
		gmr  = obj['g'] - obj['r'] 
		rmi  = obj['r'] - obj['i'] 
		gmre  = obj['g'] - obj['r'] - 2*obj['g_err'] - 2*obj['r_err']
		rmie  = obj['r'] - obj['i'] - 2*obj['r_err'] - 2*obj['i_err']


#		print gmr,gmre,obj['SATS'],obj['PHOTPTYPE'],obj['SPEC_Z'],obj['REMOVE']
#		print rmi,rmie
		if gmre > 0.8 and rmie > 0.5:
			print obj['SPEC_Z'],obj['RA'],obj['DEC'],obj['r']


	print 'STARS'
	star = lowz['PHOTPTYPE'] ==6
	for obj in lowz[star]:
		print obj['SPEC_Z'],obj['RA'],obj['DEC'],obj['TELNAME'],obj['HOST_NSAID'],obj['SATS'],obj['REMOVE']


#	ax.hexbin(gr[gal&fmag], ri[gal&fmag],extent=[-3,2.5,-3.0,2.5],cmap='Greys',bins='log')


