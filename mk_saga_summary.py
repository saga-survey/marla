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

from astropy import table
from astropy.table import Table
from tabulate import tabulate
import pyspherematch as sm

from FileLoader import GoogleSheets




SAGA_DIR   = os.environ['SAGA_DIR']
REMOVELIST = GoogleSheets('1Y3nO7VyU4jDiBPawCs8wJQt2s_PIAKRj-HSrmcWeQZo', 1379081675, header_start=1)

############################################################################
def mk_saga_summary():

	# READ SPECTRA
	file = SAGA_DIR +'/data/saga_spectra_clean.fits.gz'
	allspec = Table.read(file)

	# FIND UNIQUE SAGANAMES
	msk1 = allspec['HOST_SAGA_NAME'] != '' 
	msk2 = allspec['REMOVE'] == -1
	msk = msk1 & msk2

	sagaspec = allspec[msk]


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
	print sorted_hosts


	# CALCULATE STUFF
	data = []
	headers = ['SAGA Name', 'NSAID','RA','Dec','Nsats','MK','Dist']
	for host in sorted_hosts:

		row = []
		msk = sagaspec['HOST_SAGA_NAME'] == host[1]

		row.append(host[1])
		row.append(str(sagaspec['HOST_NSAID'][msk][0]))
		row.append(sagaspec['HOST_RA'][msk][0])
		row.append(sagaspec['HOST_DEC'][msk][0])
		row.append(str(host[0]))
		row.append(str('{:.1f}'.format(sagaspec['HOST_MK'][msk][0])))
		row.append(str(sagaspec['HOST_DIST'][msk][0]))

		data.append(row)

	print tabulate(data, headers, tablefmt="rst")



if __name__ == '__main__':
    create_saga_spectra()

 

 
