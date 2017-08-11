#! /usr/bin/python
##################################################
#  CREATE BASE SDSS CATALOGS
#
#  DOWNLOAD SQL FILE FROM SDSS (using download_host_sqlfile.py)
#  READ SQL FILE, CLEAN USING NSA, REMOVE TARGETS FROM REMOVE LIST
#
#  FOR EACH SQL FILE:
#    1. READ SQL FILES
#    2. ADD COLUMNS FOR HOST PROPERTIES
#    3. ADD COLUMNS FOR SPECTRIOSCOPY AND REMOVE LIST
#    4. CLEAN PHOTOMETRY USING NSA and REMOVELIST
#    5. WRITE FILES BASE_SQL_NSA*
#
##################################################
__all__ = ['run_hostlist','create_base_catalog']

import os
import numpy as np
from astropy.table import Table, Column
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u

import saga_tools
import fix_byhand

import pyspherematch as sm
from FileLoader import GoogleSheets, FitsTable


# SET-UP DIRECTORIES AND FILES TO BE LOADED
#SAGA_DIR    = os.getenv('SAGA_DIR', os.curdir)
SAGA_DIR     = os.environ['SAGA_DIR']
SAGA_DROPBOX = os.environ['SAGA_DROPBOX']

REMOVELIST = GoogleSheets('1Y3nO7VyU4jDiBPawCs8wJQt2s_PIAKRj-HSrmcWeQZo', 1379081675, header_start=1)
ADDLIST    = GoogleSheets('1Y3nO7VyU4jDiBPawCs8wJQt2s_PIAKRj-HSrmcWeQZo', 286645731, header_start=1)
SAGANAMES  = GoogleSheets('1GJYuhqfKeuJr-IyyGF_NDLb_ezL6zBiX2aeZFHHPr_s', 0, header_start=0)
NSACAT     = FitsTable(os.path.join(SAGA_DIR, 'cats', 'nsa_v0_1_3.fits'))
SAGACAT    = FitsTable(os.path.join(SAGA_DIR, 'data', 'saga_spectra_dirty.fits.gz'))
CLEAN      =  FitsTable(os.path.join(SAGA_DIR, 'data', 'saga_spectra_clean.fits.gz'))

MLFILE = '/Users/marlageha/Dropbox/SAGA/data/SAGA.objid.noclean.PROBS_WISE_PROBS.sept28.fits.gz'
ML = Table.read(MLFILE)

##################################  
def run_hostlist(nowise=False,noML=False,nosaga=False,flagged_obs_hosts=False):
    """
    For each host in hostlist, create base catalog
    default is to read offline Lang WISE catalogs

    Parameters
    ----------
    nowise : bool, optional.  Turn off WISE catalog read
    """

    # RUN EITHER FULL HOST LIST OR JSUT FLAG ZERO HOSTS, default to flag zero
    sheet = GoogleSheets('1b3k2eyFjHFDtmHce1xi6JKuj3ATOWYduTBFftx5oPp8', 448084634)
    flag0data = sheet.load()

    file      = SAGA_DROPBOX + '/hosts/submaster_all.ecsv.gz'
    submaster = ascii.read(file, format='ecsv')

##    print googleflag0.columns
    flag0 = np.in1d(submaster['NSAID'],flag0data['NSAID'])
    hostdata = submaster[flag0]

    # FOR EACH HOST, READ SQL AND CREATE BASE CATALOGS
    for host in hostdata:
        nid = host['NSAID']
        catalog = create_base_catalog(nid, host,nowise,noML,nosaga)
        write_base_fits(nid, catalog)



##################################  
def run_single_host(nid,nowise=False,noML=False,nosaga=False):
    """
    Run base catalog for single host

    ACCEPTS ONLY FLAG_ZERO HOSTS FOR MOMENT!

    Parameters
    ----------
    nowise : bool, optional.  Turn off WISE catalog read
    """

#    sheet = GoogleSheets('1b3k2eyFjHFDtmHce1xi6JKuj3ATOWYduTBFftx5oPp8', 448084634)
#    hostdata = sheet.load()

    file      = SAGA_DROPBOX + '/hosts/submaster_all.ecsv.gz'
    hostdata = ascii.read(file, format='ecsv')


    msk = hostdata['NSAID'] == nid
    if np.sum(msk) == 0:

        print 'NO HOST FOUND'
        catalog = create_base_catalog(nid, hostdata[0],nowise,noML,nosaga)
        write_base_fits(nid, catalog)
    else:
        catalog = create_base_catalog(nid, hostdata[msk],nowise,noML,nosaga)
        write_base_fits(nid, catalog)



##################################  
def run_named_hosts(nowise=False,noML=False,nosaga=False):
    """
    Run base catalog for SAGA named_hosts host

    ACCEPTS ONLY FLAG_ZERO HOSTS FOR MOMENT!

    Parameters
    ----------
    nowise : bool, optional.  Turn off WISE catalog read
    """

#    sheet = GoogleSheets('1b3k2eyFjHFDtmHce1xi6JKuj3ATOWYduTBFftx5oPp8', 448084634)
#    hostdata = sheet.load()

    file      = SAGA_DROPBOX + '/hosts/submaster_all.ecsv.gz'
    hostdata = ascii.read(file, format='ecsv')

    sheet = GoogleSheets('1GJYuhqfKeuJr-IyyGF_NDLb_ezL6zBiX2aeZFHHPr_s', 0)
    named_hosts = sheet.load()


    for sagahost in named_hosts:
        nid = sagahost['NSA']
        msk = hostdata['NSAID'] == nid
        if np.sum(msk) == 0:
            print nid#,sagahost['SAGA Name']
            print 'NO HOST FOUND'
            print
        else:
            catalog = create_base_catalog(nid, hostdata[msk],nowise,noML,nosaga)
            write_base_fits(nid, catalog)


##################################
def create_base_catalog(nsaid, host,nowise,noML,nosaga):
    """
    Create single base catalog from SQL request
    with value-added quantities    
    """

    # READ SQL FILE
    sqlfile  = os.path.join(SAGA_DIR, 'base_catalogs', 'sql_nsa{0}.fits.gz'.format(nsaid))
    sqltable = Table.read(sqlfile)

    # GET BASIC HOST PARAMETERS FROM GOOGLE HOST LIST
    hostra   = host['RA']        #RA
    hostdec  = host['Dec']       #DEC
    hostdist = host['distance']  #DISTANCE
    hostv    = host['vhelio']    #VHELIO
    hostMK   = host['M_K'] 
    hostMr   = host['M_r'] 
    hostMg   = host['M_g'] 

#    hostflag = host['flag']      #HOST FLAG



    # CALCULATE OBJECT DISTANCE FROM HOST
    catsc = SkyCoord(u.Quantity(sqltable['RA'], u.deg), u.Quantity(sqltable['DEC'], u.deg))
    hostcoords = SkyCoord(hostra, hostdec,unit='deg')
    seps = catsc.separation(hostcoords)

    rhost_arcm = seps.to(u.arcmin).value
    rhost_kpc  = 1000.*hostdist*np.sin(np.deg2rad(rhost_arcm/60.))



    # ADD EXTRA COLUMNS -
    size = len(sqltable)
    cols = [_filled_column('W1',-1.,size),
            _filled_column('W1ERR',-1.,size),
            _filled_column('W2',-1.,size),
            _filled_column('W2ERR',-1.,size),
            _filled_column('HOST_RA', hostra, size), #HOST
            _filled_column('HOST_DEC', hostdec, size),
            _filled_column('HOST_DIST', hostdist, size),
            _filled_column('HOST_VHOST', hostv, size),
            _filled_column('HOST_MK', hostMK, size),
            _filled_column('HOST_MR', hostMr, size),
            _filled_column('HOST_MG', hostMg, size),
            _filled_column('HOST_NSAID', nsaid, size),
            _filled_column('HOST_SAGA_NAME', ' '*48, size),
            _filled_column('HOST_NGC_NAME', ' '*48, size),
            Column(rhost_arcm, 'RHOST_ARCM'), #OBJECT
            Column(rhost_kpc, 'RHOST_KPC'),
            _filled_column('OBJ_NSAID', -1, size),
            _filled_column('SATS', -1, size),
            _filled_column('PROBABILITY_CLASS1', -1., size),
            _filled_column('RESCALED_PROBABILITY_CLASS1', -1., size),
            _filled_column('REMOVE', -1, size), # -1 = Good source; 1 = On remove list; 2 = Overlaps with NSA GALAXY
            _filled_column('TELNAME', ' '*6, size), #SPECTRA
            _filled_column('MASKNAME', ' '*48, size),
            _filled_column('ZQUALITY', -1, size),
            _filled_column('SPEC_CLASS', ' '*2, size),
            _filled_column('SPECOBJID', ' '*48L, size),
            _filled_column('SPEC_REPEAT', ' '*48, size),
            _filled_column('SPEC_SN', -99., size),
            _filled_column('SPEC_HA_EW', -99., size),
            _filled_column('SPEC_HA_EWERR', -99., size)]


    sqltable.add_columns(cols)
    del cols


    # ADD WISE NUMBERS
    if not nowise:
        wbasefile = os.path.join(SAGA_DIR, 'unwise', 'unwise_{0}.fits.gz'.format(nsaid))
        print 'Read WISE catalog: ', wbasefile

        wbasetable = FitsTable(wbasefile).load()
        wbasetable['W1_MAG'][np.isnan(wbasetable['W1_MAG'])]=-1
        wbasetable['W1_MAG_ERR'][np.isnan(wbasetable['W1_MAG_ERR'])]=-1
        wbasetable['W2_MAG'][np.isnan(wbasetable['W2_MAG'])]=-1
        wbasetable['W2_MAG_ERR'][np.isnan(wbasetable['W2_MAG_ERR'])]=-1
        id1, id2, d = sm.spherematch(sqltable['RA'], sqltable['DEC'], wbasetable['RA'], wbasetable['DEC'], 1./3600)
        sqltable['W1'][id1]    = wbasetable['W1_MAG'][id2]
        sqltable['W1ERR'][id1] = wbasetable['W1_MAG_ERR'][id2]
        sqltable['W2'][id1]    = wbasetable['W2_MAG'][id2]
        sqltable['W2ERR'][id1] = wbasetable['W2_MAG_ERR'][id2]



    # ADD ML PREDICTIONS USING OLD FILE
    if not noML:
        id1, id2, d = sm.spherematch(sqltable['RA'], sqltable['DEC'], ML['RA'], ML['DEC'], 1./3600)
        sqltable['PROBABILITY_CLASS1'][id1]    = ML['PROBABILITY_CLASS_1'][id2]





    # INITALIZE SDSS SPECTRAL ENTRIES
    msk = sqltable['SPEC_Z'] != -1
    sqltable['TELNAME'][msk] = 'SDSS'
    sqltable['MASKNAME'][msk] = 'SDSS'
    sqltable['SPEC_REPEAT'][msk] = 'SDSS'
    sqltable['ZQUALITY'][msk] = 4

    # SET BAD SDSSS SPECTROSCOPY TO ZQ = -1
    # THERE ARE A FEW RECOVERABLE SPECTRA WITH SPEC_Z_WARN = 4, 
    # BUT WILL NOT RECOVER THEM HERE
    msk_badsdss = sqltable['SPEC_Z_WARN'] != 0
    sqltable['ZQUALITY'][msk & msk_badsdss] = 1


    # IF THIS IS A SAGA HOST, SET SAGA NAME
    names = SAGANAMES.load()
    sname = ''
    ngc   = '' 
    if nsaid in names['NSA']: 
        msk = names['NSA'] == nsaid
        sname = names['SAGA'][msk]
        ngc = names['NGC'][msk]
    sqltable['HOST_SAGA_NAME'] = sname
    sqltable['HOST_NGC_NAME'] = ngc


    # ADD EXISTING SAGA SPECTRA TO FILE
    sagaspec = SAGACAT.load()
    if not nosaga:
        sqltable = saga_tools.add_saga_spec(sagaspec,sqltable)

    # SET REMOVE FLAGS
    rmv = REMOVELIST.load()
    sqltable = saga_tools.rm_removelist_obj(rmv, sqltable)


    # CLEAN AND DE-SHRED USING NSAID
    nsa = NSACAT.load()
    print 'Cleaning on NSA catalog'
    sqltable = saga_tools.nsa_cleanup(nsa, sqltable)


    # CLEAN AND DE_SHRED HIGHER REDSHIFT OBJECTS
    print 'Cleaning on SDSS catalog'
    sqltable = saga_tools.sdss_cleanup(sqltable)


    # REMOVE OBJECTS WITH BAD PHOTO-FLAGS, ADD OBJECTS BACK IN BY HAND
    addlst = ADDLIST.load()
    print 'Cleaning on PhotoFlags'
    sqltable = saga_tools.photoflags(addlst,sqltable)




    # FIX MAJOR PROBLEM SPECTRA BY HAND
    sqltable = fix_byhand.fix_basecats(sqltable)

    # FIX SATS ID BY HAND
    clean_sagaspec = CLEAN.load()
    sqltablle = fix_byhand.sat_hack(nsaid,sqltable,clean_sagaspec)

    print
    return sqltable


##################################
def _filled_column(name, fill_value, size):
    """
    Tool to allow for large strings
    """
    return Column([fill_value]*int(size), name)


##################################
def write_base_fits(nsaid, sqltable):
    outfits = os.path.join(SAGA_DIR, 'base_catalogs', 'base_sql_nsa{0}.fits.gz'.format(nsaid))
    if os.path.isfile(outfits):
        os.remove(outfits)
    sqltable.write(outfits, format='fits')


if __name__ == '__main__':
    run_hostlist()
