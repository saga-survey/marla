#! /usr/bin/python
##################################################
#  CREATE BASE SDSS CATALOGS
#  
#  DOWNLOAD SQL FILE FROM SDSS (using download_host_sqlfile.py)
#  READ SQL FILE, CLEAN USING NSA, REMOVE TARGETS FROM REMOVE LIST
# 
# FOR EACH SQL FILE:
#      1. READ SQL FILES 
#      2. ADD COLUMNS FOR HOST PROPERTIES
#      3. ADD COLUMNS FOR SPECTRIOSCOPY AND REMOVE LIST
#	   4. CLEAN PHOTOMETRY USING NSA and REMOVELIST
#      5. WRITE FILES BASE_SQL_NSA*  
#
##################################################
__all__ = ['run_hostlist','create_base_catalog']


import numpy as np

from astropy.io import ascii
from astropy.io import fits
from astropy import table
from astropy.table import Table
import os

from FileLoader import GoogleSheets, FitsTable
import pyspherematch as sm

SAGA_DIR = os.getenv('SAGADIR', os.curdir)
remove_list = GoogleSheets('1Y3nO7VyU4jDiBPawCs8wJQt2s_PIAKRj-HSrmcWeQZo', 1379081675, header_start=1)
nsa_catalog = FitsTable(os.path.join(SAGA_DIR, '/cats/nsa_v0_1_2.fits'))


def run_hostlist():	

  # READ HOST LIST FROM GOOGLE DOCS
  hostdata = GoogleSheets('1b3k2eyFjHFDtmHce1xi6JKuj3ATOWYduTBFftx5oPp8', 448084634).load()
  nsa_col = 'NSAID'
 

  # FOR EACH HOST, READ SQL AND CREATE BASE CATALOGS
  for host in hostdata:
    nid = host[nsa_col]
    create_base_catalog(nid, host)





def create_base_catalog(nsaid,host):

  # READ SQL FILE 
  sqlfile = 'sql_nsa' + str(nsaid) + '.fits'
#    sqlfile = path + '/hosts/sql_nsa' + str(nsaid) + '.fits'
  sqltable = Table.read(sqlfile)	

	
  # GET BASIC HOST PARAMETERS FROM GOOGLE HOST LIST  
  hostra   = host['RA']  #RA
  hostdec  = host['Dec']  #DEC
  hostdist = host['distance']  #DISTANCE
  hostv    = host['vhelio']  #VHELIO
  hostMK   = host['K'] #M_K
  hostflag = host['flag'] #HOST FLAG


  # SET VIRIAL RADIUS = 300 kpc, EXCLUDE 20 kpc AROUND HOST
  rkpc = 300.
  rvir = (180./np.pi)* np.arcsin(0.3/hostdist)
  rgal = (180./np.pi)* np.arcsin(0.02/hostdist)

  cdec  = np.cos(np.radians(hostdec))
  dra   = sqltable['RA'] - hostra
  ddec  = sqltable['DEC'] - hostdec
  rhost_arcm = 60*((dra * cdec) ** 2 + ddec ** 2) ** 0.5
  rhost_kpc = 1000.*hostdist*np.sin((np.pi/180.)*rhost_arcm/60.)

  # PLACE HOLDER ARRAYS	
  one = np.ones(rhost_arcm.size)
  intone = one.astype(int)
  strone = ['    ' for i in one]
  strtwo = ['                                                ' for i in one]

  # ADD EXTRA COLUMNS - 
  sqltable.add_column(table.Column(name='HOST_RA',    data=hostra*one))
  sqltable.add_column(table.Column(name='HOST_DEC',   data=hostdec*one))
  sqltable.add_column(table.Column(name='HOST_DIST',  data=hostdist*one))
  sqltable.add_column(table.Column(name='HOST_VHOST', data=hostv*one))
  sqltable.add_column(table.Column(name='HOST_MK',    data=hostMK*one))
  sqltable.add_column(table.Column(name='HOST_NSAID', data=nsaid*intone))
  sqltable.add_column(table.Column(name='HOST_FLAG',  data=hostflag*intone))
  sqltable.add_column(table.Column(name='HOST_SAGA_NAME', data=strtwo))

  # ADD EXTRA COLUMNS - OBJECT
  sqltable.add_column(table.Column(name='RHOST_ARCM', data=rhost_arcm*one))
  sqltable.add_column(table.Column(name='RHOST_KPC',  data=rhost_kpc*one))
  sqltable.add_column(table.Column(name='OBJ_NSAID', data=-1*intone))
  sqltable.add_column(table.Column(name='SATS', data=-1*intone))

  sqltable.add_column(table.Column(name='PCLASS_1', data=-1*one))
  sqltable.add_column(table.Column(name='PCLASS_1_WISE', data=-1*one))

  # SHOULD THIS OBJECT BE REMOVED FROM SDSS CATALOG?
  # -1  =  Good source
  #  1  =  On remove list
  #  2  =  Overlaps with NSA GALAXY 
  sqltable.add_column(table.Column(name='REMOVE', data=-1*intone))

  # ADD EXTRA COLUMNS - SPECTRA
  sqltable.add_column(table.Column(name='TELNAME', data=strone))
  sqltable.add_column(table.Column(name='MASKNAME', data=strtwo))
  sqltable.add_column(table.Column(name='ZQUALITY', data=-1*one))
  sqltable.add_column(table.Column(name='SPEC_REPEAT', data=strtwo))


  # ADD IN BEN'S PREDICTIONS
#  id1,id2,d = sm.spherematch(sqltable['RA'], sqltable['DEC'],ML['RA'], ML['DEC'],1./3600)
#  sqltable['PCLASS_1'][id1] = ML['PROBABILITY_CLASS_1'][id2]	
#  sqltable['PCLASS_1_WISE'][id1] = ML['PROBABILITY_CLASS_1_WISE'][id2]	


  # REPLACE WISE NUMBERS
  wise=0
  if wise:
    wbasefile = path + '/hosts/wise/base_sql_nsa' +  str(nsaid) + '_nw1.fits'
    w = fits.getdata(wbasefile)
    wbasetable = table.Table(w)	
    id1,id2,d = sm.spherematch(sqltable['ra'], sqltable['dec'],wbasetable['RA'], wbasetable['DEC'],\
    							1./3600)
    print 'Read WISE catalog: ',wbasefile
    print id1.size,sqltable['ra'].size
  	
    sqltable['w1'][id1]    = wbasetable['W1'][id2]#
    sqltable['w1err'][id1] = wbasetable['W1ERR'][id2]
    sqltable['w2'][id1]    = wbasetable['W2'][id2]
    sqltable['w2err'][id1] = wbasetable['W2ERR'][id2]
    sqltable['w1'][np.isnan(sqltable['w1'])]=9999
    sqltable['w1err'][np.isnan(sqltable['w1err'])]=9999



  # INITALIZE SDSS SPECTRAL ENTRIE
  s0  = sqltable['SPEC_Z'] != -1
  s1  = sqltable['SPEC_Z_WARN'] ==  0
  msk = s0 & s1 
  sqltable['TELNAME'][msk] = 'SDSS'
  sqltable['MASKNAME'][msk] = 'SDSS'
  sqltable['SPEC_REPEAT'][msk] = 'SDSS'
  sqltable['ZQUALITY'][msk] = 4


  # IF THIS IS A SAGA HOST, SET SAGA NAME
#  sqltable['HOST_SAGA_NAME'] = saga.saga_name(nsaid)

  # SET REMOVE FLAGS
#  sqltable = saga.rm_removelist_obj(removelist, sqltable)	


  # CLEAN USING NSAID
#  sqltable = saga.nsa_cleanup(nsa,sqltable)

  # WRITE FITS FILE
  write_base_fits(nsaid,sqltable)



def write_base_fits(nsaid,sqltable):
  outfits = path + '/hosts/base_sql_nsa' +  str(nsaid) + '.fits'

  if os.path.isfile(outfits):
    os.remove(outfits)

  sqltable.write(outfits, format='fits')




if __name__ == '__main__':
    run_query(False)
    run_query(True)


	      




