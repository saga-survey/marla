#! /usr/bin/python
##################################################
#   TOOLS FOR CREATING SAGA CATALOGS
##################################################
import numpy as np
from astropy.io import ascii
from astropy.io import fits
from astropy import table
import pyfits
import os
import glob
import pyspherematch as sm
from FileLoader import GoogleSheets, FitsTable


# SET-UP DIRECTORIES AND FILES TO BE LOADED
SAGA_DIR = os.getenv('SAGA_DIR', os.curdir)

NSACAT     = FitsTable(os.path.join(SAGA_DIR, 'cats', 'nsa_v0_1_2.fits'))


#####################################################################
def fix_basecats(sagatable):
	""" 
	These are objects that have weird personal problems
	and we need to fix by hand
	"""


	# WRONG REDSHIFT IN THE NSA, but good in SDSS
	m = sagatable['OBJID'] == 1237668367995568266
	sagatable['SPEC_Z'][m] = 0.21068
	sagatable['TELNAME'][m] = 'SDSS'
	sagatable['MASKNAME'][m] = 'SDSS'

	# DON"T BELIEVE THIS NED REDSHIFT, RESET TO -1
	m = sagatable['OBJID'] == 1237667966962434538
	sagatable['SPEC_Z'][m] = -1
	sagatable['ZQUALITY'][m] = -1

	# A SATELLITE WITH A BAD PETRORAD_R
	m=sagatable['OBJID'] == 1237651735757259099
	sagatable['PETRORAD_R'] == 2.97

    # WRONG REDSHIFT IN NSA, BUT GOOD IN SDSS
	m = sagatable['OBJID'] == 1237678881574551723
	sagatable['SPEC_Z'][m] = 1.093277
	sagatable['TELNAME'][m] = 'SDSS'
	sagatable['MASKNAME'][m] = 'SDSS'



    # NSA BELIEVES THE SDSS REDSHIFT, WHICH IS TOO LOW-SN
	m = sagatable['OBJID'] == 1237661356465979704
	if sagatable['TELNAME'][m] == 'NSA':
		sagatable['ZQUALITY'][m] = -1
 	
   # ODYSSET SATELLITE SHRED, BUT IS GOOD
	m = sagatable['OBJID'] == 1237662662147638034
	sagatable['REMOVE'][m] = -1

  # BRIGHT TARGETS FROM PALOMAR  --- NEED TO UPDATE!!
	m = sagatable['OBJID'] == 1237662698115367389
	sagatable['SPEC_Z'][m] = 0.0907
	sagatable['ZQUALITY'][m] = 4
	sagatable['TELNAME'][m] = 'MMT'
	sagatable['MASKNAME'][m] = 'PAL'


	m = sagatable['OBJID'] == 1237679996084486446
	sagatable['SPEC_Z'][m] = 0.0524
	sagatable['ZQUALITY'][m] = 4
	sagatable['TELNAME'][m] = 'MMT'
	sagatable['MASKNAME'][m] = 'PAL'


	return sagatable

######################################################
def fix_nsa():
	""" 
	Fix problems in the v1.0 NSA catalog, write to V2.0
	"""

	nsa = NSACAT.load()

	# NSA HAS WRONG OBJECT, EDIT IT OUT AND RE_WRITE NSA
	# (RA, Dec) = (127.324917502, 25.75292055)
	#
	m = nsa['NSAID'] !=  64408
	nsa = table.Table(nsa[m])

	# WRITE VERSION 1.3
	outfits = os.path.join(SAGA_DIR, 'cats', 'nsa_v0_1_3.fits')
	if os.path.isfile(outfits):
	    os.remove(outfits)
	nsa.write(outfits, format='fits')


def sat_hack(nsaid,sqltable,sagaspec):
	""" 
	Fix remove/sat combo
	"""

	m1 = sagaspec['HOST_NSAID'] == nsaid
	m2 = sagaspec['SATS'] ==1
	m3 = sagaspec['REMOVE'] ==-1

	m = m1&m2&m3
	print 'fix sats by hand',np.sum(m)

	for s in sagaspec[m]:
		msk = sqltable['OBJID'] == s['OBJID']
		sqltable['SATS'][msk] = s['SATS']
		sqltable['REMOVE'][msk] = s['REMOVE']

