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
import re

from astropy.io import ascii

from astropy import table
from astropy.table import Table
from tabulate import tabulate
import pyspherematch as sm

from FileLoader import GoogleSheets


SAGA_DIR   = os.environ['SAGA_DIR']
REMOVELIST = GoogleSheets('1Y3nO7VyU4jDiBPawCs8wJQt2s_PIAKRj-HSrmcWeQZo', 1379081675, header_start=1)

############################################################################
def calc_general():


	# CALCULATE NUMBERS FOR MW SIM SECTION
	r =21
	d1 = 0.01*3e5/70.
	d2 = 0.0045*3e5/70.
	dmin  = 5*np.log10(20. *1e6) - 5. #Mpc
	dmax  = 5*np.log10(45. *1e6) - 5.#Mpc

	Mr_min = 21 - dmin
	Mr_max = 21 - dmax
	print "Min/max survey distance ", d1,d2
	print "Min/max satellites detectable",Mr_min,Mr_max


	# CALCULATE NUMBERS FOR SPECTRA SECTION
	file = SAGA_DIR +'/data/saga_spectra_dirty.fits.gz'
	allspec = Table.read(file)

	good1 = allspec['ZQUALITY'] > 2
	good2 = allspec['REMOVE'] == -1
	good = good1 & good2



	# MMT
	mmt = allspec['TELNAME'] == 'MMT'
	print '\\newcommand{\\mmtspec}{'+str(np.sum(mmt&good))+" }"

	# AAT
	aat = allspec['TELNAME'] == 'AAT'
	print '\\newcommand{\\aatspec}{'+str(np.sum(aat&good))+" }"

	# IMACS
	imacs = allspec['TELNAME'] == 'IMACS'
	print '\\newcommand{\\imacsspec}{'+str(np.sum(imacs&good))+" }"

	# WIYN
	wiyn = allspec['TELNAME'] == 'WIYN'
	print '\\newcommand{\\wiynspec}{'+str(np.sum(wiyn&good))+" }"

	print
	all = mmt&aat&imacs&wiyn
	print '\\newcommand{\\allspec}{'+str(np.sum(all&good))+" }"


 
