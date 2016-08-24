#! /usr/bin/python
##################################################
#  COMPILE SPECTRA FROM ALL SOURCES FOR ALL FLAG 0 HOSTS
#  GAMA:  
#  MMT: 
#  AAT:
#  IMACS:
#  WIYN:  
#
##################################################

import numpy as np

from astropy.io import ascii
from astropy.io import fits
from astropy import table
from astropy.table import Table
import pyfits

import os
import glob
import pyspherematch as sm
from astropy import units as u
from astropy.coordinates import SkyCoord

import pdb

SAGA_DROPBOX= os.environ['SAGA_DROPBOX']
SAGA_DIR = os.environ['SAGA_DIR']

###################################################
# COMBINE MULTIPLE SPECTRA OF INDIVIDUAL OBJECTS
def spec_combine(spec):
	
	# FOR EACH SAGA MEASUREMENT
	for obj in spec:

		# WEIRD WAY TO FIND MATCHING OBJECTS... spherematch was given weird answers
		diff_ra = np.abs(obj['RA']- spec['RA'])
		diff_dec=np.abs(obj['DEC'] - spec['DEC'])
		m1 = diff_ra < (2./3600.)
		m2 = diff_dec < (2./3600.)
		m = m1&m2


		# IF THERE ARE MORE THAN ONE SPECTRA
		if (np.sum(m)>1):

			mspec = spec[m]
			mz = spec['ZQUALITY'] > 2
			mr = spec['RA'] != 0 

			# ARE THERE ONE OR MORE GOOD SPECTRA?
			mgood = mz & mr & m1 & m2
			print np.sum(mgood)

            # IF SO REPLACE NUMBERS INTO FILE
			if np.sum(mgood) == 1:
				spec['TELNAME'][m]  = spec['TELNAME'][mgood]
				spec['MASKNAME'][m] = spec['MASKNAME'][mgood]
				spec['ZQUALITY'][m] = spec['ZQUALITY'][mgood]
				spec['SPEC_Z'][m]   = spec['SPEC_Z'][mgood]
				spec['specobjid'][m]= spec['specobjid'][mgood]
				spec['SPEC_REPEAT'][m]= spec['TELNAME'][mgood]
                # AND SET ALL OTHER MEASUREMENTS TO ZERO
				ra = spec['RA'][mgood]
				spec['RA'][m] = 0
				spec['RA'][mgood] = ra


				

	return spec

def find_uniques(spec_data, remove_imacs=False, nearenough_sep=5*u.arcsec):
	import collections

	# MATCH ON COORDINATES    
	scs = SkyCoord(spec_data['RA'], spec_data['DEC'], unit=u.deg)
	idx1, idx2, sep2d, _ = scs.search_around_sky(scs, nearenough_sep)

	# now contruct the groups from the pairs
	grpdct = {}
	grpi = 0
	for i1, i2 in zip(idx1, idx2):
	    if i1 in grpdct:
	        if i2 in grpdct:
	            # combine the two groups by assigning grp2 items to grp1

	            # this block is by far the slowest part so if the data size grows it should be optimized
	            grp1 = grpdct[i1]
	            grp2 = grpdct[i2]
	            if grp1 != grp2:
	                to_set_to_1 = [i for i, grp in grpdct.iteritems() if grp==grp2]
	                for i in to_set_to_1:
	                    grpdct[i] = grp1
	        else:
	            #add i2 to the group i1 is already in
	            grpdct[i2] = grpdct[i1]
	    else:
	        if i2 in grpdct:
	            #add i1 to the group i2 is already in
	            grpdct[i2] = grpdct[i1]
	        else:
	            # add them both to a new group
	            grpdct[i1] = grpdct[i2] = grpi
	            grpi += 1
	        
	grpnum_to_group_members = collections.defaultdict(list)
	for k, v in grpdct.iteritems():
	    grpnum_to_group_members[v].append(k)
	# convert the members into arrays
	grpnum_to_group_members = {k:np.array(v) for k, v in grpnum_to_group_members.iteritems()}

	# identify which is the "best" spectrum (meaning the first zq=4 spectrum)
	idxs_to_keep = []
	new_repeats = []
	for grpnum, allmembers in grpnum_to_group_members.iteritems():
	    if remove_imacs:
	        members =  allmembers[spec_data['TELNAME'][allmembers]!='IMACS']
	        if len(members)==0:
	            continue
	    else:
	        members = allmembers
	        
	    idxs_to_keep.append(members[np.argsort(spec_data['ZQUALITY'][members])[-1]])
	    new_repeats.append('+'.join(np.unique(spec_data['SPEC_REPEAT'][members])))

	# now build the output table from the input
	unique_objs = spec_data[np.array(idxs_to_keep)]
	del unique_objs['SPEC_REPEAT']
	unique_objs['SPEC_REPEAT'] = new_repeats

	return unique_objs




############################################################################
# READ GAMA CATALOG
def read_gama():

  # READ GAMA FILE, DOWNLOADED FROM ??
	gfile = SAGA_DIR + '/cats/GAMA_SpecObj.fits'
	g    = pyfits.open(gfile)    # BUG: NEED TO READ WITH PYFITS, NOT ASTROPY
	gm   = g[1].data



  # ACCEPT GAMA SPECTRA WITH GOOD QUALITY AND NOT IN SDSS
  	msk = (gm['NQ'] >= 3) & (gm['SURVEY'] != 'SDSS') 
	gama=gm[msk]


  # PLACE HOLDER ARRAYS	
	one = np.ones(gama.size)
	telname = ['GAMA' for i in one]
	spec_repeat = ['GAMA' for i in one]
	gama_specobjid = gama['SPECID']
	gamaname = gama['CATAID'].astype('S80')

	# GAMA GOES UP TO 5, set these to 4 to avoid confusion
	m = gama['nq'] == 5
	gama['nq'][m] = 4



  # CREATE GAMA SPEC TABLE
	gama_table = table.table.Table([gama['RA'], gama['DEC'],gamaname,one,\
								    gama['z'],gama['nq'],telname,spec_repeat], \
	     				            names=('RA', 'DEC', 'MASKNAME','specobjid',\
 		                            'SPEC_Z','ZQUALITY','TELNAME','SPEC_REPEAT'))

	print "GAMA SPECTRA = ",gama.size
	return gama_table



############################################################################
# READ MMT ZLOG CATALOGS
def read_mmt():

	mmt_path  = SAGA_DROPBOX + '/Spectra/Final/MMT/'
	mmt_files = glob.glob(mmt_path+'*zlog')

	n=0
	for mfile in mmt_files:	

		# ACCEPT ALL GOOD SPECTRA
		print mfile
		mdata = ascii.read(mfile,guess=False,format='no_header', delimiter=' ',names=['col1','ra','dec','mag','z','col6','zq','col8','col9','col10','col11'])


		zq  = mdata.field('zq') >= 1 
		obj = mdata.field('mag') != 0 
		msk = zq & obj  
										# ONLY OBJECTS WITH QUALITY GE 3

		mmt = mdata[msk]

		# PLACE HOLDER ARRAYS	
		one = np.ones(len(mmt))
		telname = ['MMT' for i in one]
		spec_repeat = ['MMT' for i in one]
		maskid = [mfile for i in one]
		maskid = [x.split(SAGA_DROPBOX,1)[1] for x in maskid]


	   # CREATE MMT SPEC TABLE
		mmt_table1 = table.table.Table([15.*mmt['ra'], mmt['dec'], maskid, mmt['col8'],\
									    mmt['z'], mmt['zq'], telname,spec_repeat], \
			     				        names=('RA', 'DEC', 'MASKNAME','specobjid',\
			     		                       'SPEC_Z','ZQUALITY','TELNAME','SPEC_REPEAT'))

	   # CREATE OR APPEND TO MMT TABLE	
		if (n==0):  mmt_table = mmt_table1
		if (n > 0): mmt_table = table.vstack([mmt_table,mmt_table1])
		n=n+1
	print "Number of MMT Files = ",n	
	return mmt_table



############################################################################
# READ AAT ZLOG CATALOGS
def read_aat():

	aat_path  = SAGA_DROPBOX + '/Spectra/Final/AAT/'
	aat_files = glob.glob(aat_path+'*zlog')

	n=0
	for afile in aat_files:	

		# ACCEPT ALL GOOD SPECTRA
		adata = ascii.read(afile, data_start=0, delimiter=' ')
		msk = adata.field('col7') >= 1  # ONLY OBJECTS WITH QUALITY GE 1
		aat = adata[msk]

		# PLACE HOLDER ARRAYS	
		one = np.ones(len(aat))
		telname = ['AAT' for i in one]
		spec_repeat = ['AAT' for i in one]
		maskid = [afile for i in one]
		maskid = [x.split(SAGA_DROPBOX,1)[1] for x in maskid]


	   # CREATE MMT SPEC TABLE
		aat_table1 = table.table.Table([aat['col2'], aat['col3'], maskid, aat['col8'],\
									    aat['col5'], aat['col7'], telname,spec_repeat], \
			     				        names=('RA', 'DEC', 'MASKNAME','specobjid',\
			     		                       'SPEC_Z','ZQUALITY','TELNAME','SPEC_REPEAT'))


	   # CREATE OR APPEND TO AAT TABLE	
		if (n==0):  aat_table = aat_table1
		if (n > 0): aat_table = table.vstack([aat_table,aat_table1])
		n=n+1
	print "Number of AAT Files = ",n	
	return aat_table


############################################################################
# READ AAT MZ CATALOGS
def read_aat_mz():

	aat_path  = SAGA_DROPBOX + '/Spectra/Final/AAT/'
	aat_files = glob.glob(aat_path+'*mz')

	n=0
	for afile in aat_files:	

		# ACCEPT ALL GOOD SPECTRA
		adata = ascii.read(afile)
		msk = adata.field('col14') >= 1  # ONLY OBJECTS WITH QUALITY GE 1
		mz = adata[msk]

		# PLACE HOLDER ARRAYS	
		one = np.ones(len(mz))
		telname = ['AAT' for i in one]
		spec_repeat = ['AAT' for i in one]
		maskid = [afile for i in one]
		maskid = [x.split(SAGA_DROPBOX,1)[1] for x in maskid]

		# RA/DEC COMES IN RADIANS!
		ra  = mz['col3']*180/np.pi
		dec = mz['col4']*180/np.pi


	   # CREATE MMT SPEC TABLE
		mz_table1 = table.table.Table([ra, dec, maskid, mz['col1'],\
									    mz['col13'], mz['col14'], telname,spec_repeat],\
			     				        names=('RA', 'DEC', 'MASKNAME','specobjid',\
			     		                       'SPEC_Z','ZQUALITY','TELNAME','SPEC_REPEAT'))


	   # CREATE OR APPEND TO AAT TABLE	
		if (n==0):  aat_table = mz_table1
		if (n > 0): aat_table = table.vstack([aat_table,mz_table1])
		n=n+1
	print "Number of AAT MZ Files = ",n	
	return aat_table



############################################################################
# READ IMACS ZLOG CATALOGS


def read_imacs():

	imacs_path  = SAGA_DROPBOX + '/Spectra/Final/IMACS/'
	imacs_files = glob.glob(imacs_path+'*zlog')

	n=0
	for ifile in imacs_files:	

		# ACCEPT ALL GOOD SPECTRA
		idata = ascii.read(ifile,guess=False,format='no_header', delimiter=' ',names=['specid','ra','dec','col4','z','col6','zq','col8','col9','col10','slitid','col12'])
		msk   = idata.field('zq') >= 1  

		imacs = idata[msk]

		# PLACE HOLDER ARRAYS	
		one = np.ones(len(imacs))
		telname = ['IMACS' for i in one]
		spec_repeat = ['IMACS' for i in one]
		maskid = [ifile for i in one]
		maskid = [x.split(SAGA_DROPBOX,1)[1] for x in maskid]


	   # CREATE MMT SPEC TABLE
		imacs_table1 = table.table.Table([imacs['ra'], imacs['dec'], imacs['slitid'],imacs['col8'],\
									    imacs['z'], imacs['zq'], telname,spec_repeat], \
			     				        names=('RA', 'DEC', 'MASKNAME','specobjid',\
			     		                       'SPEC_Z','ZQUALITY','TELNAME','SPEC_REPEAT'))


	   # CREATE OR APPEND TO imacs TABLE	
		if (n==0):  imacs_table = imacs_table1
		if (n > 0): imacs_table = table.vstack([imacs_table,imacs_table1])
		n=n+1
	print "Number of IMACS Files = ",n	
	return imacs_table




############################################################################
# READ WIYN ZSPEC CATALOGS
def read_wiyn():

	from astropy import units as u
	from astropy.coordinates import SkyCoord

	wiyn_path  = SAGA_DROPBOX + '/Spectra/Final/WIYN/'
	wiyn_files = glob.glob(wiyn_path+'*fits.gz')

	n=0
	for wfile in wiyn_files:	

		# ACCEPT ALL GOOD SPECTRA
		w = fits.getdata(wfile)
		wiyn1 = table.Table(w)	
		msk =  wiyn1['ZQUALITY'] >= 1  # ONLY OBJECTS WITH QUALITY GE 3
		wiyn = wiyn1[msk]

		# PLACE HOLDER ARRAYS	
		one = np.ones(len(wiyn))
		telname = ['WIYN' for i in one]
		spec_repeat = ['WIYN' for i in one]

      # KEEP GERENIC PATH TO DATA FILES (remove specific path)
		maskid = [wfile for i in one]
		maskid = [x.split(SAGA_DROPBOX,1)[1] for x in maskid]

	  # CONVERT WIYN STRING COORDINATES TO DEGREES	
		c=SkyCoord(wiyn['RA'],wiyn['DEC'],frame='icrs',unit=(u.hourangle, u.deg))


	   # CREATE WIYN SPEC TABLE		
		wiyn_table1 = table.table.Table([c.ra.deg, c.dec.deg, maskid, wiyn['FID'],\
									    wiyn['Z'], wiyn['ZQUALITY'], telname,spec_repeat], \
			     				        names=('RA', 'DEC', 'MASKNAME','specobjid',\
			     		                       'SPEC_Z','ZQUALITY','TELNAME','SPEC_REPEAT'))


	   # CREATE OR APPEND TO AAT TABLE	
		if (n==0):  wiyn_table = wiyn_table1
		if (n > 0): wiyn_table = table.vstack([wiyn_table,wiyn_table1])
		n=n+1
	print "Number of WIYN Files = ",n	
	return wiyn_table


def read_sdss(base):

  	msk  = base['ZQUALITY'] == 4
	sdss = base[msk]


  # PLACE HOLDER ARRAYS	

  # CREATE GAMA SPEC TABLE
	sdss_table = table.table.Table([sdss['RA'], sdss['DEC'],sdss['MASKNAME'],sdss['OBJID'],\
								    sdss['SPEC_Z'],sdss['ZQUALITY'],sdss['TELNAME'],sdss['SPEC_REPEAT']], \
	     				            names=('RA', 'DEC', 'MASKNAME','specobjid',\
 		                            'SPEC_Z','ZQUALITY','TELNAME','SPEC_REPEAT'))

	return sdss_table


############################################################################
# CREATE DEIMOS ZSPEC CATALOGS
def read_deimos():

# CREATE DEIMOS SPEC TABL
	#  [ ODYSSEY, ODSSEY-Pen]
	ra   = [247.825839103498, 221.86742, 150.12470]
	dec  = [20.210825313885, -0.28144459, 32.561687]
	mask = ['deimos2014', 'deimos2016-DN1','deimos2016-MD1']
	sid  = [0,0,0]
	v    = [2375/3e5,0.056,1.08]
	zq   = [4,4,4]
	tel  = ['DEIMOS','DEIMOS','DEIMOS']


	deimos_table = table.table.Table([ra, dec, mask, sid,\
									    v, zq, tel, tel], \
			     				        names=('RA', 'DEC', 'MASKNAME','specobjid',\
			     		                       'SPEC_Z','ZQUALITY','TELNAME','SPEC_REPEAT'))
	return deimos_table




