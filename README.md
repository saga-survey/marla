# marla
The directory contains scripts to create SAGA base catalogs for a list
of SAGA host galaxies.


--  There are three Google docs used in the scripts below and must be
maintained to run the base catalogs:

	SAGA_HOST_LISTS:  contains the full SAGA host list and the Flag
	Zero list.
	Remove_Lists/targets:  contains a list of objects which have been flagged by hand as junk.
	SAGA_Hosts+Satellites:  Official file of SAGA common names



-- Directory and file setup:
	1.   Set environment variable SAGADIR to top level directory where
         SAGA data will live.

	2. Within SAGADIR, create two directories /hosts to house the sql
    and base catalogs.   /cats to hold catalogs.

	3.  Download NSA catalogs and place in SAGADIR/cats
	
	4.  Casjobs:    To run casjobs, need to setup an account at:
         http://skyserver.sdss3.org/CasJobs/CreateAccount.aspx
	     Then edit your .bash_profile
              >  export CASJOBS_WSID='2090870927'   # get your WSID from site above
	          >export CASJOBS_PW='my password'


-- Download sql files.   Downloads *all* objects within 1 degree of all hosts in host list
	> import download_host_sqlfile
	> download_host_sqlfile.run_query()


--  Create value-added SAGA base catalogs.
	> import create_base_catalogs
	> create_base_catalogs.run_query(nowise=1)

	If you want to include Lang WISE data, need to download _nw files
    into SAGA/nw_hosts/
