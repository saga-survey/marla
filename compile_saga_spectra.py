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
from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import glob
import pyspherematch as sm


from FileLoader import GoogleSheets, FitsTable
import saga_tools
import read_saga_spectra

SAGA_DIR = os.environ['SAGA_DIR']
SAGA_DROPBOX= os.environ['SAGA_DROPBOX']


# GOOGLE DOCUMENTS URLs
REMOVELIST = GoogleSheets('1Y3nO7VyU4jDiBPawCs8wJQt2s_PIAKRj-HSrmcWeQZo', 1379081675, header_start=1)
#NSACAT     = FitsTable(os.path.join(SAGA_DIR, 'cats', 'nsa_v0_1_3.fits'))


############################################################################
def compile_saga_spectra(flagged_obs_hosts=False):


    # RUN EITHER FULL HOST LIST OR JUST FLAG ZERO HOSTS, default to flag zero
    if flagged_obs_hosts:
        sheet = GoogleSheets('1GJYuhqfKeuJr-IyyGF_NDLb_ezL6zBiX2aeZFHHPr_s', 0)
    else:
        sheet = GoogleSheets('1b3k2eyFjHFDtmHce1xi6JKuj3ATOWYduTBFftx5oPp8', 448084634)
	hostdata     = sheet.load()

	# READ ALL NON-SDSS SPECTRA AND COMBINE INTO SINGLE TABLE
	gama_table   = read_saga_spectra.read_gama()
	mmt_table    = read_saga_spectra.read_mmt()
	aat_table    = read_saga_spectra.read_aat()
	aat_table_mz = read_saga_spectra.read_aat_mz()
	imacs_table  = read_saga_spectra.read_imacs()
	wiyn_table   = read_saga_spectra.read_wiyn()
	deimos_table = read_saga_spectra.read_deimos()


	# LEAST TO MOST IMPORTANT
	sagaspec = table.vstack([wiyn_table,imacs_table,gama_table,deimos_table,aat_table,mmt_table,aat_table_mz],\
						     metadata_conflicts='silent')

	# WRITE FILE WITH ALL SPECTRA
	file  = os.path.join(SAGA_DIR, 'data', 'saga_spectra_raw.fits.gz')	
	write_fits(sagaspec,file)

	# COMBINE MULTIPLE SPECTRA
#	read_saga_spectra.spec_combine(sagaspec)
	spec = sagaspec
	sagaspec = read_saga_spectra.find_uniques(spec)



 # FOR EACH FLAG ZERO HOST OR NON-FLAG ZERO OBSERVED HOST
 # READ BASE SDSS FILE
 #  KEEP ONLY SDSS SPECTRA AND MATCHES FROM OTHER SOURCES
	nhost = 0
	for host in hostdata:
	    flag  = host['flag']   # SELECT ONLY NON_FLAG0#	    
	    if flag == 0:
			nra  = host['RA']
			ndec = host['Dec']
			nsaid  = host['NSAID']# NAME OF BASE SQL FILES
			basefile  = os.path.join(SAGA_DIR, 'base_catalogs', 'base_sql_nsa{0}.fits.gz'.format(nsaid))
			basetable = Table.read(basefile)	
			print nsaid


			# CALCULATE OBJECT DISTANCE FROM HOST
			catsc = SkyCoord(u.Quantity(sagaspec['RA'], u.deg), u.Quantity(sagaspec['DEC'], u.deg))
			hostcoords = SkyCoord(nra*u.deg, ndec*u.deg)
			seps  = catsc.separation(hostcoords)
			rhost = seps.to(u.deg).value
			host_spec_objs = rhost < 1.0


		  # COMBINE SDSS SPECTRA AND SAGA SPECTRA	
			sdss_table   = read_saga_spectra.read_sdss(basetable)
			hostspec     = table.vstack([sdss_table,sagaspec[host_spec_objs]])



			tmp = hostspec
			hostspec = read_saga_spectra.find_uniques(tmp)

			id2,id1,d = sm.spherematch(hostspec['RA'],hostspec['DEC'],basetable['RA'], basetable['DEC'],3./3600,nnearest=1)
			basetable['TELNAME'][id1]    = hostspec['TELNAME'][id2]
			basetable['MASKNAME'][id1]   = hostspec['MASKNAME'][id2]
			basetable['ZQUALITY'][id1]   = hostspec['ZQUALITY'][id2]
			basetable['SPEC_Z'][id1]     = hostspec['SPEC_Z'][id2]
			basetable['SPEC_Z_ERR'][id1]     = hostspec['SPEC_Z_ERR'][id2]
			basetable['SPECOBJID'][id1]  = hostspec['specobjid'][id2]
			basetable['SPEC_REPEAT'][id1]= hostspec['SPEC_REPEAT'][id2]

			m = basetable['SPEC_Z'] != -1

			m2=basetable['TELNAME'][m] == ''
			print 'number with telname = ',np.sum(m2)

	      # COMBINE INTO SINGLE ALLSPEC FILE
			if (nhost == 0):
				allspec = basetable[m]
			if (nhost > 0):
				allspec = table.vstack([allspec,basetable[m]])  # APPEND			
			nhost=nhost+1	



	# INITIALIZE SATS ARRAY (3=primary, 2=lowz, 1=satellite)	
	allspec = saga_tools.fill_sats_array(allspec)


	# CLEAN UP REPEAT ENTRY OF SATELLITES 
	allspec = saga_tools.repeat_sat_cleanup(allspec)


	# WRITE ALL SPECTRA TAKEN
	file  = os.path.join(SAGA_DIR, 'data', 'saga_spectra_dirty.fits.gz')
	write_fits(allspec, file)


	# WRITE ALL GOOD SPECTRA
	# KEEP ALL GOOD SPECTRA WHICH ARE GALAXIES
	zql     = allspec['ZQUALITY'] >= 3 
	rml     = allspec['REMOVE'] == -1 
	galonly = allspec['PHOTPTYPE'] == 3
	clean   = zql & rml & galonly
	allclean = allspec[clean]
	file  = os.path.join(SAGA_DIR, 'data', 'saga_spectra_clean.fits.gz')	
	write_fits(allclean,file)


	# WRITE TEXT FILE OF SATELLITES WHICH CAN BE EDITED
#	file  = os.path.join(SAGA_DIR, 'data', 'saga_satellites.dat')	
#	allsats = allclean['SATS' == 1]
#	write_satellite_file(allsats,file)


	return allspec 


def write_fits(file, filename):
	if os.path.isfile(filename):
		 os.remove(filename)
	print 'writing file: ',filename
	file.write(filename,format='fits')

	 

if __name__ == '__main__':
    create_saga_spectra()

 

 
