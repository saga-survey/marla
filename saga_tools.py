#! /usr/bin/python
##################################################
#   TOOLS FOR CREATING SAGA CATALOGS
##################################################
import numpy as np
from astropy.io import ascii
from astropy.io import fits
from astropy import table
from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import glob
import pyspherematch as sm


#####################################################################
def photoflags(addlist,sagatable):
	""" 
	SET REMOVE FLAG = 3 for bad SDSS Photo Flags
	SET BACK TO REMOVE = -1 if in by-hand list  
	"""
	binned1   = sagatable['BINNED1'] == 0
	saturated = sagatable['SATURATED'] != 0
	baderr    = sagatable['BAD_COUNTS_ERROR'] != 0 

	flgs = binned1 | saturated | baderr
	sagatable['REMOVE'][flgs] = 3


   # MATCH sql OBJECTS TO THOSE IN GOOGLE DOC ADD LIST
	id1,id2,d = sm.spherematch(sagatable['RA'], sagatable['DEC'],\
		           addlist.field('Targ_RA'), addlist.field('Targ_Dec'),\
		           1./3600,nnearest=1)

	nmatch = np.size((d > 0.0).nonzero())
	print "add list objects = ",nmatch

  # SET REMOVED FLAG BACK TO -1
  	if (nmatch > 0):
		sagatable['REMOVE'][id1] = -1


	return sagatable


#####################################################################
#  USE GOOGLE DOC REMOVE LIST TO SET REMOVE FLAG
#   -1 =  GOOD OBJECT
#    1 =  ON REMOVE LIST, DO NOT USE  (rm_removelist_obj)
#    2 = SHREDDED OBJECT BASED ON NSA (nsa_cleanup)
#    3 = bad photometry flag
#
def rm_removelist_obj(removelist,sagatable):

	
   # MATCH sql OBJECTS TO THOSE IN GOOGLE DOC REMOVE LIST
	id1,id2,d = sm.spherematch(sagatable['RA'], sagatable['DEC'],\
		           removelist.field('Targ_RA'), removelist.field('Targ_Dec'),\
		           1./3600,nnearest=1)

	nmatch = np.size((d > 0.0).nonzero())
	print "remove list objects = ",nmatch

  # SET REMOVED FLAG TO 1
  	if (nmatch > 0):
		sagatable['REMOVE'][id1] = 1

	return sagatable




#####################################################################
#  USE NSA TO CREATE UP REPEATED OR SHREDDED OBJECTS
#  FIND NSA MATCH, THEN USE PETROTH90 RADIUS TO SET REMOVE FLAG
#   -1 = GOOD OBJECT
#    1 =  ON Google REMOVE LIST, DO NOT USE  (rm_removelist_obj)
#    2 = SHREDDED OBJECT BASED ON NSA (nsa_cleanup)
# 	 3 = photoflags
#    4 = SHREDDED OBJECT BASED ON SDSS (z > 0.05)
#
def nsa_cleanup(nsa,sagatable):

   # MATCH NSA to SAGA, BEGIN WITH SMALLER RADIUS
	id2,id1,d = sm.spherematch(nsa['RA'], nsa['DEC'],\
							sagatable['RA'], sagatable['DEC'],\
							5./3600,nnearest = 1)

	nmatch = np.size(id1)

    # FOR EACH UNIQUE NSA MATCH
	for i in range(0,nmatch-1):

		nid = id2[i]
		sid = id1[i]
		nra = nsa['RA'][nid]
		ndec= nsa['DEC'][nid]

	# FIRST PASS USING SPHEREMATCH
		reff_nsa = 3*nsa['PETROTH90'][nid]		
		m1,m2,dd = sm.spherematch(sagatable['RA'], sagatable['DEC'],\
			[nra],[ndec],reff_nsa/3600)


	# FIND POINTS INSIDE ELLIPSE
		a = 2.*nsa['PETROTH90'][nid]
		b = nsa['SERSIC_BA'][nid] * a
		th = (np.pi/180) * (nsa['SERSIC_PHI'][nid]+270)
		
		sth = np.sin(th)
		cth = np.cos(th)

		x = 3600*(sagatable['RA'][m1] - nra)
		y = 3600*(sagatable['DEC'][m1] - ndec)
		tmp1 = (x*cth) - (y*sth)  
		tmp2 = (x*sth) + (y*cth)
		
		rem = (tmp1/a)**2 + (tmp2/b)**2 <= 1.  # true for points inside ellipse
		

      # REMOVE SHREDDED OBJECTS NEAR NSA 
		sagatable['REMOVE'][m1[rem]] = 2     

      # EXPLICITLY ADD NSA PROPERTIES BACK INTO DATABASE
		sagatable['REMOVE'][sid]   = -1   
		sagatable['ZQUALITY'][sid] = 4
		sagatable['TELNAME'][sid]  = 'NSA'
		sagatable['PHOTPTYPE'][sid]= 3
		sagatable['PHOT_SG'][sid]  = 'GALAXY'
		sagatable['RA'][sid]       = nra
		sagatable['DEC'][sid]      = ndec
		sagatable['SPEC_Z'][sid]   = nsa['Z'][nid]
		sagatable['SPEC_Z_WARN'][sid] = 0
		sagatable['MASKNAME'][sid] = nsa['ZSRC'][nid]
		sagatable['OBJ_NSAID'][sid]= nsa['NSAID'][nid]

		# REPLACE PHOTOMETRY
		mag = 22.5 - 2.5*np.log10(nsa['SERSICFLUX'][nid])
		mag[np.isnan(mag)]= -99

		sagatable['u'][sid]= mag[2]
		sagatable['g'][sid]= mag[3]
		sagatable['r'][sid]= mag[4]
		sagatable['i'][sid]= mag[5]
		sagatable['z'][sid]= mag[6]



#		sagatable['expRad_r'][sid] = nsa['SERSIC_TH50'][nid]
#		sagatable['sb_exp_r'][sid] = -99#nsa['PETROTH90']
#		sagatable['petroR90_r'][sid] = nsa['PETROTH90'][nid]#
#		sagatable['petroR50_r'][sid] = nsa['PETROTH50'][nid]
#		sagatable['petroMag_r'][sid] = -99#nsa['PETROTH90']


	return sagatable




#####################################################################
# FOR GALAXIES BEYOND NSA REDSHIFT CUTOFF
# FIND NEARBY SHREDS W/IN REFF AND REMOVE
#
def sdss_cleanup(sagatable):

	# FOR GALAXIES WITH SDSS SPEC BEYOND NSA REDSHIFT CUTOFF
	# AND GOOD MEASURE OF RADIUS
	s1 = sagatable['SPEC_Z'] > 0.05 
	s2 = sagatable['PETRORADERR_R'] > 0.
	s3 = sagatable['PETRORAD_R']/sagatable['PETRORADERR_R'] > 2.  # better than 50% error
	s4 = sagatable['REMOVE'] == -1

	smsk = s1&s2&s3#&s4
	sdss_spec = sagatable[smsk]


	# FOR EACH SDSS WITH REDSHIFT > NSA CUTOFF
	for obj in sdss_spec:

       # USE TWICE REFF
		reff = 1.25*obj['PETRORAD_R']


		catsc = SkyCoord(u.Quantity(sagatable['RA'], u.deg), u.Quantity(sagatable['DEC'], u.deg))
		objcoords = SkyCoord(obj['RA']*u.deg, obj['DEC']*u.deg)
		seps = catsc.separation(objcoords)

		robj = seps.to(u.arcsec).value

		m1 = robj < reff   
		if np.sum(m1) > 3:
			print np.sum(m1),obj['RA'], obj['DEC'],obj['REMOVE']

		# REMOVE SHREDDED OBJECTS NEAR SDSS SPECTRUM 
		sagatable['REMOVE'][m1] = 4  

	# BUT KEEP ORIGINAL OBJECT
	sagatable['REMOVE'][smsk] = -1


	return sagatable


#####################################################################
#
def fill_sats_array(sqltable):
	#  FILL IN THE ARRAY SATS AS:
	#      SATS = 3   if object is SAGA PRIMARY
	#      SATS = 2   if object is low-z, z < 0.05
	#      SATS = 1   if object is SAGA SATELLITE (+/- 200 km/s w/in 300kpc)
	#      SATS = 0   if object is high-z, z > 0.05
	#      SATS = -1  no redshift

	c = 2.997e5	

	# Velocity difference between object and SAGA host
	vdiff = np.abs((c * sqltable['SPEC_Z']) -  sqltable['HOST_VHOST']) 


	# GALAXY ONLY, OBJECTS WITH SPECTRA ONLY
	galcut = (sqltable['SPEC_Z'] != -1) & \
	  	     (sqltable['PHOTPTYPE'] == 3) & \
	  	     (sqltable['ZQUALITY'] > 2)


	# SATS = 0, ALL HIGH-Z (z > 0.05) GALAXIES
	highz  = galcut & \
		    (sqltable['SPEC_Z'] >= 0.05)
	sqltable['SATS'][highz] = 0


	# SATS = 2, ALL LOW-Z (z < 0.05) GALAXIES
	ncut = sqltable['MASKNAME'] != 'ned'  # NED LOWZ VELOCITIES HAVE PROBLEMS
	lowz  = galcut & ncut &\
		    (sqltable['SPEC_Z'] < 0.05)
	sqltable['SATS'][lowz] = 2


	# SATS = 1, SATELLITES!!
	vhost= sqltable['HOST_VHOST']
	robj = sqltable['RHOST_KPC']

	vcut = np.abs(vdiff) < 250
	rcut = (robj < 300) &  (robj > 20)
	ncut = sqltable['MASKNAME'] != 'ned'  # NED LOWZ VELOCITIES HAVE PROBLEMS

	sats = galcut & vcut & rcut & ncut
		      
	sqltable['SATS'][sats] = 1.
	return sqltable


	# SATS = 3,FIND THE PRIMARY GALAXY
#	prim  = galcut & \
#		   vdiff < vhost + 20 & vdiff > vhost - 20 \
#		   sqltable['rhost_kpc'] > 20   




#####################################################################
#  For each satellite, search for other nearby satellites
#  if there are more than one matches, add these extras to repeat spec
#
def repeat_sat_cleanup(sagatable):


	sats = sagatable['SATS'] == 1 #& s['ZQUALITY'] > 2

	for s in sagatable[sats]:
		if s['REMOVE'] == -1:

			rsat    = s['RA']
			dsat    = s['DEC']

			m1,m2,d = sm.spherematch(sagatable['RA'], sagatable['DEC'],\
			                     [rsat],[dsat], 10./3600)

			m = d != 0
			if np.sum(m):
				r = m1[m]
				k = m1[d ==0]
				sagatable['REMOVE'][r] = 3
				sagatable['REMOVE'][k] = -1

				sagatable['SATS'][r] = -99
				sagatable['SATS'][k] = 1


	return sagatable

##################################
def saga_name(names,nsaid):
	"""
	Read Google doc Observed Host List and
	parse SAGA name from file
	"""

	saga_names  = names['SAGA Name']
	saga_nsaids = names['NSA']

	sname = ''
	if nsaid in saga_nsaids: 
		msk = saga_nsaids == nsaid
		sname = saga_names[msk]

	return sname	

#####################################################################
#  ADD EXISTING SAGA SPECTRA TO BASE CATALOGS
#
def add_saga_spec(sagaspec,sqltable):

	
   # MATCH sql OBJECTS TO THOSE IN SAGA SPEC
	id1,id2,d = sm.spherematch(sqltable['RA'], sqltable['DEC'],\
				   sagaspec['RA'], sagaspec['DEC'],\
		           1./3600,nnearest=1)

	nmatch = np.size((d >= 0.0).nonzero())
	print "adding existing SAGA SPECTRA = ",nmatch

  # SET REMOVED FLAG TO 1
  	if (nmatch > 0):
		sqltable['SPEC_Z'][id1]   = sagaspec['SPEC_Z'][id2]
		sqltable['ZQUALITY'][id1] = sagaspec['ZQUALITY'][id2]
		sqltable['TELNAME'][id1]  = sagaspec['TELNAME'][id2]
		sqltable['MASKNAME'][id1] = sagaspec['MASKNAME'][id2]
		sqltable['SPECOBJID'][id1]    = sagaspec['SPECOBJID'][id2]
		sqltable['SPEC_REPEAT'][id1]  = sagaspec['SPEC_REPEAT'][id2]
		sqltable['SATS'][id1]  = sagaspec['SATS'][id2]
		sqltable['REMOVE'][id1]  = sagaspec['REMOVE'][id2]


	return sqltable





