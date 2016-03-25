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

import  PIL
from PIL import Image

import matplotlib.pyplot as plt


from astropy import table
from astropy.table import Table
from astropy.io import fits


from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=9, usetex=True)



SAGA_DIR   = os.environ['SAGA_DIR']
SAGA_DROPBOX= os.environ['SAGA_DROPBOX']


dir_jpg = '/Users/marlageha/Projects/SAGA/jpegs/'


############################################################################
# PLOT JUST IMAGES SORTED BY HOST
# AND SATELLITE BRIGHTNESS
def plot_images():


    # READ SPECTRA
	file = SAGA_DIR +'/data/saga_spectra_dirty.fits.gz'
	allspec = Table.read(file)

	# FIND UNIQUE SAGANAMES AND ORDER BY NSATS
	msk          = allspec['HOST_SAGA_NAME'] != '' 
	sagaspec     = allspec[msk]
	sorted_hosts = sort_saga_hosts(sagaspec)
	nhost        = np.size(sorted_hosts)/2.
	

	fig = plt.figure(figsize=(5,7))

	# ONE LINE OF IMAGES PER HOST
	nh=0
	for host in sorted_hosts:
		nh=nh+1
		i=0
		
		# SORT SATELLITES ON BRIGHTNESS
		m1 = allspec['HOST_SAGA_NAME'] == host[1]
		m2 = allspec['SATS']  == 1
		sats = allspec[m1&m2]
		sats.sort('r')

		for s in sats:
			i=i+1

			img, t = get_img(s)
			ax = fig.add_subplot(nhost, 9, 9*(nh-1) + i)
			ax.text(0,250,t,color='white',fontsize=2)
			ax.get_xaxis().set_visible(False)
			ax.get_yaxis().set_visible(False)
			ax.imshow(img)

	plt.savefig('fig_images.pdf')


############################################################################
# PLOT IMAGES + SPECTRA
def plot_spectra():


 # READ SPECTRA
	file = SAGA_DIR +'/data/saga_spectra_dirty.fits.gz'
	allspec = Table.read(file)

	# FIND UNIQUE SAGANAMES AND ORDER BY NSATS
	m1       = allspec['HOST_SAGA_NAME'] != '' 
	m2       = allspec['SATS']  == 1
	sats     = allspec[m1&m2]
	nsats    = np.sum(m1&m2)
	sats.sort('r')


	fig = plt.figure(figsize=(5, 7))
#	fig.subplots_adjust(left=0.1, right=0.95, wspace=0.05,
#                    bottom=0.1, top=0.95, hspace=0.05)


	i=0
	for s in sats:
		if s['TELNAME'] == 'MMT':

			img, t    = get_img(s)
			wv, flux  = get_spec(s)


			# plot image
			ax = fig.add_subplot(nsats, 2, i+1)
			ax.get_xaxis().set_visible(False)
			ax.get_yaxis().set_visible(False)

			#plot spectrum
			ax = fig.add_subplot(nsats, 2, i+2)
			ax.plot(wv, flux,'.')
			ax.title=t

			i=i+2

	plt.savefig('fig_spectra.pdf')



def get_img(sat):

	jpgfile = dir_jpg + str(sat['HOST_SAGA_NAME'])+ '_' +str(sat['r']).format()  +'.jpg'
	imgjpg = Image.open(jpgfile)
	dmod = 5.*np.log10(sat['HOST_DIST']*1e6) -5.
	mr = sat['r'] - dmod
	mrs = '{0:.1f}'.format(mr)
	title = str(sat['HOST_SAGA_NAME']) + ' Mr='+mrs

	return imgjpg, title




def get_spec(sat):

	if sat['TELNAME'] == 'MMT':
		specfile = SAGA_DROPBOX + sat['MASKNAME'].replace('zlog','fits.gz')
		hdulist = fits.open(specfile)
		awv = hdulist[0].data
		afl = hdulist[1].data

		slit = int(float(sat['SPECOBJID']))-1
		print slit
		wv   = awv[slit,:]
		fl   = afl[slit,:]
		return wv,fl


def plot_spec_lines(z):


	return redline,bluelines





def sort_saga_hosts(sagaspec):
	"""
	Find unique named SAGA hosts in allspec. 
	Sort names by nsats
	"""
	# FIND UNIQUE SAGA NAMES and CALCULATE NSATS
	unique_hosts = []
	nsats        = []
	for s in sagaspec:
		if s['HOST_SAGA_NAME'] not in unique_hosts: 
			unique_hosts.append(s['HOST_SAGA_NAME'])

			# CALCULATE NSATS FOR GIVEN HOST
			msk1 = sagaspec['HOST_SAGA_NAME'] == s['HOST_SAGA_NAME']
			msk2 = sagaspec['SATS'] == 1
			msk3 = sagaspec['REMOVE'] != 2
			msk = msk1 & msk2 & msk3
 			n = np.sum(msk)

			nsats.append([n,s['HOST_SAGA_NAME']])

#			print  s['HOST_SAGA_NAME'],sagaspec['RA'][msk],sagaspec['DEC'][msk]


	sorted_hosts = sorted(nsats,reverse=True)
	return sorted_hosts
