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


	# READ SPECTRA
	file = SAGA_DIR +'/data/saga_spectra_all.fits.gz'
	allspec = Table.read(file)

	good = allspec['ZQUALITY'] > 2


	# ARE THERE SPECTRA WITH BAD PHOTOMERTY?   
	# IF SO, REMOVE FOR CALCULATIONS BELOW
	msk = allspec['REMOVE'] != -1
	sat = allspec['SATS'] == 1
	print np.sum(msk & good)
	for a in allspec[msk&good&sat]:
		print a['RA'],a['DEC'],a['BINNED1'],a['BAD_COUNTS_ERROR'],a['SATURATED']

	# FIND UNIQUE SAGANAMES AND ORDER BY NSATS
	msk          = allspec['HOST_SAGA_NAME'] != '' 
	sagaspec     = allspec[msk]
	sorted_hosts = sort_saga_hosts(sagaspec)



	# NUMBER OF NON-SDSS SPECTRA



	# NUMBER OF SPECTRA BY TELESCOPE
	mmt = allspec['TELNAME'] == 'MMT'
 	aat = allspec['TELNAME'] == 'AAT'
 	imacs = allspec['TELNAME'] == 'IMACS'
 	wiyn = allspec['TELNAME'] == 'WIYN'

 	print 'Number of good/taken MMT spectra = ',np.sum(mmt & good), np.sum(mmt)
 	print 'Number of good/taken AAT spectra = ',np.sum(aat & good), np.sum(aat)
 	print 'Number of good/taken IMACS spectra = ',np.sum(imacs & good), np.sum(imacs)

	# NUMBER OF STELLAR SPECTRA AND ANY SATELLITES?




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
			msk = msk1 & msk2
			n = np.sum(msk)

			nsats.append([n,s['HOST_SAGA_NAME']])


	sorted_hosts = sorted(nsats,reverse=True)
	return sorted_hosts


###############################################
def tex_table_header():
	table_header = """
		\documentclass[8pt]{article} 
		\usepackage{natbib}
		\usepackage{fancyhdr}
		\usepackage{graphics}
		\pagestyle{fancy}
		\setlength{\ textwidth}{8.5in}
		\setlength{\ textheight}{8.30in}
		\setlength{\ topmargin}{20pt}
		\setlength{\oddsidemargin}{-1in}
		\setlength{\hoffset}{0pt}
		\setlength{\ voffset}{10pt}
		\ begin{document}
		\ begin{center}
		\ begin{table}


	"""
	table_header = re.sub(' ', '', table_header).strip()
#	table_header = re.sub('\s+', ' ', table_header).strip()

	return table_header



def tex_table_footer():
	table_foot = """

		\end{table}
		\end{center}
		\end{document}
	"""
	table_foot = re.sub(' ', '', table_foot).strip()
	return table_foot



if __name__ == '__main__':
    create_saga_spectra()

 

 
