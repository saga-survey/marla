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

from astropy.io import ascii
from astropy.io import fits
from astropy import table
from astropy.table import Table
import os
import glob
import pyspherematch as sm

from FileLoader import GoogleSheets
import saga_tools
import read_saga_spectra

SAGA_DIR = os.environ['SAGA_DIR']


# GOOGLE DOCUMENTS URLs
REMOVELIST = GoogleSheets('1Y3nO7VyU4jDiBPawCs8wJQt2s_PIAKRj-HSrmcWeQZo', 1379081675, header_start=1)


############################################################################
def create_allspectra(flagged_obs_hosts=False):


    # RUN EITHER FULL HOST LIST OR JUST FLAG ZERO HOSTS, default to flag zero
    if flagged_obs_hosts:
        sheet = GoogleSheets('1GJYuhqfKeuJr-IyyGF_NDLb_ezL6zBiX2aeZFHHPr_s', 0)
        nsa_col = 'NSA'
    else:
        sheet = GoogleSheets('1b3k2eyFjHFDtmHce1xi6JKuj3ATOWYduTBFftx5oPp8', 448084634)
        nsa_col = 'NSAID'

	hostdata = sheet.load()# READ ALL NON-SDSS SPECTRA AND COMBINE INTO SINGLE TABLE
	gama_table   = read_saga_spectra.read_gama()
	mmt_table    = read_saga_spectra.read_mmt()
	aat_table    = read_saga_spectra.read_aat()
	imacs_table  = read_saga_spectra.read_imacs()
	wiyn_table   = read_saga_spectra.read_wiyn()
	deimos_table = read_saga_spectra.read_deimos()

	#pdb.set_trace()
	sagaspec = table.vstack([gama_table,mmt_table,aat_table,imacs_table,wiyn_table],\
						     metadata_conflicts='silent')



 # FOR EACH FLAG ZERO HOST OR NON-FLAG ZERO OBSERVED HOST
 # READ BASE SDSS FILE
 #  KEEP ONLY SDSS SPECTRA AND MATCHES FROM OTHER SOURCES
	nhost = 0
	for host in hostdata:
	    flag  = host['flag']   # SELECT ONLY NON_FLAG0
	    if flag == 0:

			nra  = host['RA']
			ndec = host['Dec']
			nsaid  = host[nsa_col]# NAME OF BASE SQL FILES
			basefile  = os.path.join(SAGA_DIR, 'hosts', 'base_sql_nsa{0}.fits'.format(nsaid))
			basetable = Table.read(basefile)	
			print basefile


			# KEEP SDSS SPECTRA
			sdss_msk =  basetable['ZQUALITY'] == 4 
			sdss_spec = basetable[sdss_msk]



		   # MATCH GAMA+SAGA IN SDSSS TO GET PHOTOMETRIC PROPERTIES
			id1,id2,d = sm.spherematch(basetable['RA'], basetable['DEC'],sagaspec['RA'], sagaspec['DEC'],1./3600,nnearest=1)
			nmatch = np.size((d > 0.0).nonzero())

			if (nmatch != 0):
				basetable['TELNAME'][id1]    = sagaspec['TELNAME'][id2]
				basetable['MASKNAME'][id1]   = sagaspec['MASKNAME'][id2]
				basetable['ZQUALITY'][id1]   = sagaspec['ZQUALITY'][id2]
				basetable['SPEC_Z'][id1]     = sagaspec['SPEC_Z'][id2]
#				basetable['SPECOBJID'][id1]  = sagaspec['specobjid'][id2]
	#				sdss_spec['spec_repeat'][id1]= sagaspec[id2]['spec_z']


		  # COMBINE SDSS SPECTRA AND SAGA SPECTRA	
			spectable = table.vstack([sdss_spec,basetable[id1]])
			            		  	

		  # KEEP ONLY GALAXIES 	
		  	spectable.columns
			galaxy   = spectable['PHOTPTYPE'] == 3
			spectable = spectable[galaxy]


	      # COMBINE INTO SINGLE ALLSPEC FILE
			if (nhost == 0):
				allspec = spectable
			if (nhost > 0):
				allspec = table.vstack([allspec,spectable])  # APPEND			
			nhost=nhost+1	


	# DOUBLE CHECK TO REMOVE OBJECTS ON REMOVE LIST
	rmv = REMOVELIST.load()
	allspec = saga_tools.rm_removelist_obj(rmv,allspec)


	# INITIALIZE SATS ARRAY (3=primary, 2=lowz, 1=satellite)	
	allspec = saga_tools.fill_sats_array(allspec)

	# CLEAN UP REPEAT ENTRY OF SATELLITES 
#	allspec = saga_tools.repeat_sat_cleanup(allspec)



	# WRITE ALL SPECTRA TAKEN
	file='allspectaken.fits'
	if os.path.isfile(file):
		 os.remove(file)
	print 'writing file: ',file
	allspec.write(file,format='fits')

	# WRITE ALL GOOD SPECTRA
	msk = allspec['ZQUALITY'] >= 3
	allgood=allspec[msk]
	file='allgoodspec.fits'
	if os.path.isfile(file):
		 os.remove(file)
	print 'writing file: ',file
	allgood.write(file,format='fits')



	return allspec 



	 

 

 
