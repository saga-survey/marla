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
import pyspherematch as sm
from FileLoader import GoogleSheets, FitsTable

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
        catalog = create_base_catalog(nid, host)
        write_base_fits(nid, catalog)


def _filled_column(name, fill_value, size):
    return Column([fill_value]*int(size), name)


def create_base_catalog(nsaid,host):
    # READ SQL FILE
    sqlfile = 'sql_nsa{0}.fits'.format(nsaid)
    #sqlfile = os.path.join(SAGA_DIR, 'hosts', sqlfile)
    sqltable = Table.read(sqlfile)

    # GET BASIC HOST PARAMETERS FROM GOOGLE HOST LIST
    hostra   = host['RA']  #RA
    hostdec  = host['Dec']  #DEC
    hostdist = host['distance']  #DISTANCE
    hostv    = host['vhelio']  #VHELIO
    hostMK   = host['K'] #M_K
    hostflag = host['flag'] #HOST FLAG

    # SET VIRIAL RADIUS = 300 kpc, EXCLUDE 20 kpc AROUND HOST
    #rkpc = 300.
    #rvir = np.rad2deg(np.arcsin(0.3/hostdist))
    #rgal = np.rad2deg(np.arcsin(0.02/hostdist))

    cdec  = np.cos(np.deg2rad(hostdec))
    dra   = sqltable['RA'] - hostra
    ddec  = sqltable['DEC'] - hostdec
    rhost_arcm = 60*((dra * cdec) ** 2 + ddec ** 2) ** 0.5
    rhost_kpc = 1000.*hostdist*np.sin(np.deg2rad(rhost_arcm/60.))

    # ADD EXTRA COLUMNS -
    size = len(sqltable)
    cols = [_filled_column('HOST_RA', hostra, size), #HOST
            _filled_column('HOST_DEC', hostdec, size),
            _filled_column('HOST_DIST', hostdist, size),
            _filled_column('HOST_VHOST', hostv, size),
            _filled_column('HOST_MK', hostMK, size),
            _filled_column('HOST_NSAID', nsaid, size),
            _filled_column('HOST_FLAG', hostflag, size),
            _filled_column('HOST_SAGA_NAME', ' '*48, size),
            Column(rhost_arcm, 'RHOST_ARCM'), #OBJECT
            Column(rhost_kpc, 'RHOST_KPC'),
            _filled_column('OBJ_NSAID', -1, size),
            _filled_column('SATS', -1, size),
            _filled_column('PCLASS_1', -1, size),
            _filled_column('PCLASS_1_WISE', -1, size),
            _filled_column('REMOVE', -1, size), # -1 = Good source; 1 = On remove list; 2 = Overlaps with NSA GALAXY
            _filled_column('TELNAME', ' '*4, size), #SPECTRA
            _filled_column('MASKNAME', ' '*48, size),
            _filled_column('ZQUALITY', -1, size),
            _filled_column('SPEC_REPEAT', ' '*48, size)]
    sqltable.add_columns(cols)

    # ADD IN BEN'S PREDICTIONS
    #id1, id2, d = sm.spherematch(sqltable['RA'], sqltable['DEC'], ML['RA'], ML['DEC'], 1./3600)
    #sqltable['PCLASS_1'][id1] = ML['PROBABILITY_CLASS_1'][id2]
    #sqltable['PCLASS_1_WISE'][id1] = ML['PROBABILITY_CLASS_1_WISE'][id2]

    # REPLACE WISE NUMBERS
    wise = 0
    if wise:
        wbasefile = os.path.join(SAGA_DIR, 'hosts', 'wise', 'base_sql_nsa{0}_nw1.fits',format(nsaid))
        wbasetable = FitsTable(wbasefile).load()
        id1, id2, d = sm.spherematch(sqltable['RA'], sqltable['DEC'], wbasetable['RA'], wbasetable['DEC'], 1./3600)
        print 'Read WISE catalog: ', wbasefile
        print id1.size, sqltable['ra'].size

        #TODO: probably need to add columns to sqltable
        sqltable['w1'][id1]    = wbasetable['W1'][id2]
        sqltable['w1err'][id1] = wbasetable['W1ERR'][id2]
        sqltable['w2'][id1]    = wbasetable['W2'][id2]
        sqltable['w2err'][id1] = wbasetable['W2ERR'][id2]
        sqltable['w1'][np.isnan(sqltable['w1'])] = 9999
        sqltable['w1err'][np.isnan(sqltable['w1err'])] = 9999

    # INITALIZE SDSS SPECTRAL ENTRIE
    s0  = sqltable['SPEC_Z'] != -1
    s1  = sqltable['SPEC_Z_WARN'] ==  0
    msk = s0 & s1
    sqltable['TELNAME'][msk] = 'SDSS'
    sqltable['MASKNAME'][msk] = 'SDSS'
    sqltable['SPEC_REPEAT'][msk] = 'SDSS'
    sqltable['ZQUALITY'][msk] = 4

    # IF THIS IS A SAGA HOST, SET SAGA NAME
    #sqltable['HOST_SAGA_NAME'] = saga.saga_name(nsaid)

    # SET REMOVE FLAGS
    #sqltable = saga.rm_removelist_obj(remove_list, sqltable)

    # CLEAN USING NSAID
    #sqltable = saga.nsa_cleanup(nsa_catalog, sqltable)

    return sqltable


def write_base_fits(nsaid,sqltable):
    outfits = os.path.join(SAGA_DIR, 'hosts', 'base_sql_nsa{0}.fits'.format(nsaid))
    if os.path.isfile(outfits):
        os.remove(outfits)
    sqltable.write(outfits, format='fits')


if __name__ == '__main__':
    run_hostlist()







