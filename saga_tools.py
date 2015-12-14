#! /usr/bin/python
##################################################
#   TOOLS FOR CREATING SAGA CATALOGS
##################################################
import numpy as np
from astropy.io import ascii
from astropy.io import fits
from astropy import table
import os
import glob
import pyspherematch as sm




#####################################################################
#  USE GOOGLE DOC REMOVE LIST TO SET REMOVE FLAG
#   -1 =  GOOD OBJECT
#    1 =  ON REMOVE LIST, DO NOT USE  (rm_removelist_obj)
#    2 = SHREDDED OBJECT BASED ON NSA (nsa_cleanup)
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
#    1 =  ON REMOVE LIST, DO NOT USE  (rm_removelist_obj)
#    2 = SHREDDED OBJECT BASED ON NSA (nsa_cleanup)
#
def nsa_cleanup(nsa,sagatable):

   # MATCH NSA to SAGA, BEGIN WITH SMALLER RADIUS
	id1,id2,d = sm.spherematch(sagatable['RA'], sagatable['DEC'],\
							   nsa['RA'], nsa['DEC'],\
							   2./3600,nnearest = 1)
	nmatch = np.size((d > 0.0).nonzero())

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
#		sagatable['phot_sg'][sid]  = 3
		sagatable['RA'][sid]       = nra
		sagatable['DEC'][sid]      = ndec
		sagatable['SPEC_Z'][sid]   = nsa['Z'][nid]
		sagatable['SPEC_Z_WARN'][sid] = 0
		sagatable['MASKNAME'][sid] = nsa['ZSRC'][nid]
		sagatable['OBJ_NSAID'][sid]= nsa['NSAID'][nid]

		# REPLACE PHOTOMETRY
		mag = 22.5 - 2.5*np.log10(nsa['NMGY'][nid])
		An  = nsa['EXTINCTION'][nid]

		sagatable['u'][sid]= mag[2] + An[2]     # add back in unextinction corrected
		sagatable['g'][sid]= mag[3] + An[3]  
		sagatable['r'][sid]= mag[4] + An[4]  
		sagatable['i'][sid]= mag[5] + An[5]  
		sagatable['z'][sid]= mag[6] + An[6]  


#		sagatable['expRad_r'][sid] = -99#nsa['PETROTH90']
#		sagatable['sb_exp_r'][sid] = -99#nsa['PETROTH90']
#		sagatable['petroR90_r'][sid] = nsa['PETROTH90'][nid]
#		sagatable['petroR50_r'][sid] = nsa['PETROTH50'][nid]
#		sagatable['petroMag_r'][sid] = -99#nsa['PETROTH90']


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
	vdiff = np.abs((c * sqltable['spec_z']) -  sqltable['HOST_VHOST']) 


	# GALAXY ONLY, OBJECTS WITH SPECTRA ONLY
	galcut = (sqltable['spec_z'] != -1) & \
	  	     (sqltable['phot_sg'] == 3) 


	# SATS = 0, ALL HIGH-Z (z > 0.05) GALAXIES
	highz  = galcut & \
		    (sqltable['spec_z'] >= 0.05)
	sqltable['SATS'][highz] = 0


	# SATS = 2, ALL LOW-Z (z < 0.05) GALAXIES
	ncut = sqltable['MASKNAME'] != 'ned'  # NED LOWZ VELOCITIES HAVE PROBLEMS
	lowz  = galcut & ncut &\
		    (sqltable['spec_z'] < 0.05)
	sqltable['SATS'][lowz] = 2


	# SATS = 1, SATELLITES!!
	vhost= sqltable['HOST_VHOST']
	robj = sqltable['RHOST_KPC']

	vcut = np.abs(vdiff) < 200
	rcut = (robj < 300) &  (robj > 20)
	ncut = sqltable['MASKNAME'] != 'ned'  # NED LOWZ VELOCITIES HAVE PROBLEMS
	sats = galcut & vcut & rcut & ncut
		      
	sqltable['SATS'][sats] = 1
	return sqltable


	# SATS = 3,FIND THE PRIMARY GALAXY
#	prim  = galcut & \
#		   vdiff < vhost + 20 & vdiff > vhost - 20 \
#		   sqltable['rhost_kpc'] > 20   




#####################################################################
# Take input list and remove objects lised in Google Doc Remove list
#  For each satellite, search for other nearby satellites
#  if there are more than one matches, add these extras to repeat spec
#
def repeat_sat_cleanup(sagatable):



	isat = np.where(sagatable['SATS'] == 1)[0]
	print isat.size

	for i in isat:
		if sagatable['REMOVE'][i] == -1:

			rsat    = sagatable['ra'][i]
			dsat    = sagatable['dec'][i]

			m1,m2,d = sm.spherematch(sagatable['ra'], sagatable['dec'],\
			                     [rsat],[dsat], 10./3600)
			
#			sagatable['REMOVE'][m1] = 3
#			sagatable['REMOVE'][i]  = -1


	return sagatable

##################################
def saga_name(names,nsaid):
	"""Read Google doc Observed Host List and
	parse SAGA name from file"""

	saga_names  = names['SAGA Name']
	saga_nsaids = names['NSA']

	sname = ''
	if nsaid in saga_nsaids: 
		msk = saga_nsaids == nsaid
		sname = saga_names[msk]

	return sname	




