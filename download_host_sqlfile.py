#! /usr/bin/python
##################################################
# READ host list, print sql query
# 
#  > import saga
#  > saga.()
#
# USES Dan FM's
# https://github.com/dfm/casjobs
#
# TO RUN CASJOBS, NEED TO GET AN ACCOUNT FROM:
#    http://skyserver.sdss3.org/CasJobs/CreateAccount.aspx
#
# THEN EDIT YOUR .BASH_PROFILE
#    CASJOBS_WSID='2090870927'   # get your WSID from site above
#    export CASJOBS_WSID
#
#    CASJOBS_PW = 'my password'
#    export CASJOBS_PW
#
# MG 7/2014
##################################################

import numpy as np
from astropy.io import ascii
from astropy import units as u
import requests
import os


path = os.environ['SAGADIR']


# GOOGLE DOCUMENTS URLs
SAGA_URLS = {'Host_noflags':'https://docs.google.com/a/yale.edu/spreadsheets/d/1b3k2eyFjHFDtmHce1xi6JKuj3ATOWYduTBFftx5oPp8/export?format=csv&gid=448084634',
                 'Obs_Hosts': 'https://docs.google.com/spreadsheet/pub?key=0AggNS3_oqq91dHdVZ1J3cVZrRWVjcmkyeUF1WjVNQ1E&output=csv&gid=0',
                 'Remove':'http://docs.google.com/spreadsheets/d/1Y3nO7VyU4jDiBPawCs8wJQt2s_PIAKRj-HSrmcWeQZo/export?format=csv&gid=1379081675'}


#################################################################
def run_query():
    """
    Send casjob queries for each host

    Parameters
    ----------
    None  


    Returns
    -------
    Downloads Flag0 host files into sql files located in SAGADIR/hosts

    """

	# READ HOST LIST FROM GOOGLE DOCS
    csvurl = SAGA_URLS.get('Host_noflags', None)
    res = requests.get(csvurl)
    csv = res.content.split('\n')[1:]
    res.close() 
    hostdata = ascii.read(csv)


	# FOR EACH HOST, DOWNLOAD SQL QUERY
    for host in hostdata:
		nid  = host['col1']
		print nid
		if nid != -1:

			nra  = host['col2']
			ndec = host['col3']

	      # CREATE SQL QUERY
			qry = construct_sdss_query(nid,nra,ndec, 1.0)
			print
			print qry


		  # RUN CASJOB.  NEED TO DEFINE CASJOBS_WSID and CASJOBS_PW in PATH
			cjob = casjobs.CasJobs(base_url='http://skyserver.sdss3.org/casjobs/services/jobs.asmx')

			out = 'sql_nsa'+str(nid)
			outfits = out+'.fits'

		  # IF FILE DOESN"T ALREADY EXIST, SUBMIT JOB TO CAS
			if not os.path.isfile(outfits):
		   		job_id = cjob.submit(qry, context='DR10')

			  # WAIT FOR JOB TO FINISH
				time.sleep(90)
			   	code = 0 
			   	code, status  = cjob.status(job_id)
			   	print code, status
				while (code != 5):
					time.sleep(60)
					code, status  = cjob.status(job_id)
					print code, status


			   # DOWNLOAD FILE AND REMOVE FROM MYDB LISTING (SAVE SPACE)
				print 'downloading ',outfits
				cjob.request_and_get_output(out, 'FITS', outfits)
				cjob.drop_table(out)


#################################################################

def run_query_nonFlag0():
    """
    RUN SQL QUERY FOR NON-FLAG0 hosts
    Based on non-flag0 hosts in OBS_HOST google docs

    """

    # READ HOST LIST FROM GOOGLE DOCS
    csvurl = SAGA_URLS.get('Obs_Hosts', None)
    res = requests.get(csvurl)
    csv = res.content.split('\n')[1:]
    res.close() 
    hostdata = ascii.read(csv)



    # FOR EACH HOST, DOWNLOAD SQL QUERY
    for host in hostdata:
        flag  = host['col8']   # SELECT ONLY NON_FLAG0
        print flag
        if flag != 0:

            nra  = host['col4']
            ndec = host['col5']
            nid  = host['col3']

          # CREATE SQL QUERY
            qry = construct_sdss_query(nid,nra,ndec, 1.0)
            print
            print qry


          # RUN CASJOB.  NEED TO DEFINE CASJOBS_WSID and CASJOBS_PW in PATH
            cjob = casjobs.CasJobs(base_url='http://skyserver.sdss3.org/casjobs/services/jobs.asmx')

            out = 'sql_nsa'+str(nid)
            outfits = out+'.fits'

          # IF FILE DOESN"T ALREADY EXIST, SUBMIT JOB TO CAS
            if not os.path.isfile(outfits):
                job_id = cjob.submit(qry, context='DR10')

              # WAIT FOR JOB TO FINISH
                time.sleep(90)
                code = 0 
                code, status  = cjob.status(job_id)
                print code, status
                while (code != 5):
                    time.sleep(60)
                    code, status  = cjob.status(job_id)
                    print code, status


               # DOWNLOAD FILE AND REMOVE FROM MYDB LISTING (SAVE SPACE)
                print 'downloading ',outfits
                cjob.request_and_get_output(out, 'FITS', outfits)
                cjob.drop_table(out)



#################################################################

def construct_sdss_query(outname, ra, dec, radius=1*u.deg):
    """
    Generates the query to send to the SDSS to get the full SDSS catalog around
    a target.

    Parameters
    ----------
    ra : `Quantity` or float
        The center/host RA (iconn degrees if float)
    dec : `Quantity` or float
        The center/host Dec (in degrees if float)
    radius : `Quantity` or float
        The radius to search out to (in degrees if float)
    into : str or None
        The name of the table to construct in your `mydb` if you want to use
        this with CasJobs, or None to have no "into" in the SQL. This also
        adjust other parts of the query a little to work with CasJobs instead
        of the direct query.
    magcut : 2-tuple or None
        if not None, adds a magnitude cutoff.  Should be a 2-tuple
        ('magname', faintlimit). Ignored if None.


    Returns
    -------
    query : str
        The SQL query to send to the SDSS skyserver


    """
    from textwrap import dedent



    query_template = dedent("""
    SELECT  p.objId  as objID,
    p.ra, p.dec, p.type as phot_sg, p.flags, p.specObjID, 
    p.modelMag_u as u, p.modelMag_g as g, p.modelMag_r as r,p.modelMag_i as i,p.modelMag_z as z,
    p.modelMagErr_u as u_err, p.modelMagErr_g as g_err, p.modelMagErr_r as r_err,p.modelMagErr_i as i_err,p.modelMagErr_z as z_err,
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
    """)



    if isinstance(ra, float):
        ra = ra*u.deg
    ra = ra.to(u.deg).value

    if isinstance(dec, float):
        dec = dec*u.deg
    dec = dec.to(u.deg).value

    if isinstance(radius, float):
        radius = radius*u.deg
    radarcmin = radius.to(u.arcmin).value

    # ``**locals()`` means "use the local variable names to fill the template"
    return query_template.format(**locals())




