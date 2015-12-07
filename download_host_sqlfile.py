#! /usr/bin/python
##################################################
# READ host list, print sql query
# 
# USES Dan FM's
#    https://github.com/dfm/casjobs
#
# TO RUN CASJOBS, NEED TO GET AN ACCOUNT FROM:
#    http://skyserver.sdss3.org/CasJobs/CreateAccount.aspx
#
# THEN EDIT YOUR .BASH_PROFILE
#    export CASJOBS_WSID='2090870927'   # get your WSID from site above
#    export CASJOBS_PW='my password'
#
# MG 7/2014
##################################################
__all__ = ['run_query', 'run_casjob', 'construct_sdss_query']

import time
import os
import re
from astropy.table import Table
from astropy import units as u
from casjobs import CasJobs


def get_google_csv_url(key, gid):
    return 'https://docs.google.com/spreadsheets/d/{key}/export?format=csv&gid={gid}'.format(**locals())


def run_query(flagged_obs_hosts=False):
    """
    Send casjob queries for each host

    Parameters
    ----------
    flagged_obs_hosts : bool, optional

    Returns
    -------
    Downloads host files into sql files located in SAGADIR/hosts

    """

    # RUN EITHER FULL HOST LIST OR JSUT FLAG ZERO HOSTS
    if flagged_obs_hosts:
        url = get_google_csv_url('1GJYuhqfKeuJr-IyyGF_NDLb_ezL6zBiX2aeZFHHPr_s', 0)
        nsa_col = 'NSA'
    else:
        url = get_google_csv_url('1b3k2eyFjHFDtmHce1xi6JKuj3ATOWYduTBFftx5oPp8', 448084634)
        nsa_col = 'NSAID'
    
    # READ HOST LIST FROM GOOGLE DOCS
    hostdata = Table.read(url, format='ascii.csv')


    # FOR EACH HOST, DOWNLOAD SQL QUERY
    for host in hostdata:
        if flagged_obs_hosts and not host['Host Flag']:
            continue
        nid = host[nsa_col]
        if nid == -1:
            continue
        print nid

        # CREATE SQL QUERY
        qry = construct_sdss_query(nid, host['RA'], host['Dec'])
        print qry
        run_casjob(qry, 'sql_nsa{}'.format(nid))



def run_casjob(query, outname):
    """
    Run single casjob

    Parameters
    ----------
    query : output from construct_sdss_query

    Returns
    -------
    Downloads casjob output

    """

    if not all(k in os.environ for k in ('CASJOBS_WSID', 'CASJOBS_PW')):
        raise ValueError('You are not setup to run casjobs')
    
    SAGA_DIR = os.getenv('SAGADIR', os.curdir)
    
    # USES POST 
    cjob = CasJobs(base_url='http://skyserver.sdss.org/casjobs/services/jobs.asmx', request_type='POST')
    outfits = os.path.join(SAGA_DIR, outname + '.fits')


    # IF FILE DOESN"T ALREADY EXIST, SUBMIT JOB TO CAS
    if not os.path.isfile(outfits):
        job_id = cjob.submit(query, context='DR10')
        code = None
        while code != 5:
            time.sleep(60)
            code, status = cjob.status(job_id)
            print code, status

        print 'downloading', outfits
        cjob.request_and_get_output(outname, 'FITS', outfits)
        cjob.drop_table(outname)


def construct_sdss_query(outname, ra, dec, radius=1.0):
    """
    Generates the query to send to the SDSS to get the full SDSS catalog around
    a target.

    Parameters
    ----------
    outname : string
    ra : `Quantity` or float
        The center/host RA (iconn degrees if float)
    dec : `Quantity` or float
        The center/host Dec (in degrees if float)
    radius : `Quantity` or float
        The radius to search out to (in degrees if float)
    
    Returns
    -------
    query : str
        The SQL query to send to the SDSS skyserver
    """

    query_template = """
    SELECT  p.objId  as objID,
    p.ra, p.dec, p.type, dbo.fPhotoTypeN(p.type) as phot_sg, p.flags, p.specObjID, 
    p.modelMag_u as u, p.modelMag_g as g, p.modelMag_r as r,p.modelMag_i as i,p.modelMag_z as z,
    p.modelMagErr_u as u_err, p.modelMagErr_g as g_err, p.modelMagErr_r as r_err,p.modelMagErr_i as i_err,p.modelMagErr_z as z_err,
    p.psfMag_u as psf_u, p.psfMag_g as psf_g, p.psfMag_r as psf_r, p.psfMag_i as psf_i, p.psfMag_z as psf_z,
    p.extinction_u as Au, p.extinction_g as Ag, p.extinction_r as Ar, p.extinction_i as Ai, p.extinction_z as Az,
    p.fibermag_r, p.fiber2mag_r,
    p.expRad_r, p.expMag_r + 2.5*log10(2*PI()*p.expRad_r*p.expRad_r + 1e-20) as sb_exp_r,
    p.petroR50_r, p.petroR90_r,p.petroMag_r,
    p.petroMag_r + 2.5*log10(2*PI()*p.petroR50_r*p.petroR50_r) as sb_petro_r,
    ISNULL(w.j_m_2mass,9999) as J, ISNULL(w.j_msig_2mass,9999) as Jerr, 
    ISNULL(w.H_m_2mass,9999) as H, ISNULL(w.h_msig_2mass,9999) as Herr, 
    ISNULL(w.k_m_2mass,9999) as K, ISNULL(w.k_msig_2mass,9999) as Kerr,
    ISNULL(w.w1mpro,9999) as w1, ISNULL(w.w1sigmpro,9999) as w1err, 
    ISNULL(w.w2mpro,9999) as w2, ISNULL(w.w2sigmpro,9999) as w2err,
    ISNULL(s.z, -1) as spec_z, ISNULL(s.zErr, -1) as spec_z_err, ISNULL(s.zWarning, -1) as spec_z_warn, 
    ISNULL(pz.z,-1) as photoz,ISNULL(pz.zerr,-1) as photoz_err

    FROM dbo.fGetNearbyObjEq({ra}, {dec}, {radarcmin}) n, PhotoPrimary p
    INTO mydb.sql_nsa{outname}
    LEFT JOIN SpecObj s ON p.specObjID = s.specObjID
    LEFT JOIN PHOTOZ  pz ON p.ObjID = pz.ObjID
    LEFT join WISE_XMATCH as wx on p.objid = wx.sdss_objid
    LEFT join wise_ALLSKY as w on  wx.wise_cntr = w.cntr
    WHERE n.objID = p.objID AND (flags & dbo.fPhotoFlags('BINNED1')) != 0
        AND (flags & dbo.fPhotoFlags('SATURATED')) = 0
        AND (flags & dbo.fPhotoFlags('BAD_COUNTS_ERROR')) = 0 
    """

    if isinstance(ra, u.Quantity):
        ra = ra.to(u.deg).value

    if isinstance(dec, u.Quantity):
        dec = dec.to(u.deg).value

    if isinstance(radius, u.Quantity):
        radarcmin = radius.to(u.arcmin).value
    else:
        radarcmin = (radius*u.deg).to(u.arcmin).value

    # ``**locals()`` means "use the local variable names to fill the template"
    q = query_template.format(**locals())
    q = re.sub('\s+', ' ', q).strip()
    return q


if __name__ == '__main__':
    run_query(False)
    run_query(True)
