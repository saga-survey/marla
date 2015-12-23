# marla
The directory contains scripts to create SAGA base catalogs for the a
host list.

### Google docs
There are three google docs used in the scripts and these must be maintained to run the base catalogs:

	SAGA_HOST_LISTS:  Contains the full SAGA host list and the Flag Zero list.
	REMOVE_LIST:          Contains a list of objects which have been flagged by hand as junk.
	SAGA_HOSTS+SATELLITES:  Official file of SAGA common names

### Directories and File Setup:
These local directories and files are require to run scripts.

--  Environment variables:
    SAGA_DIR:           Top level directory where SAGA data will live on local disk
	SAGA_DROPBOX:  Local directory for DropBox/SAGA

-- Create Local Directories:
	SAGA_DIR/hosts:    Houses the sql downloads and base files
	SAGA_DIR/cats:      Local copy of NSA catalog


-- Download NSA catalog (nsa_v0_1_2.fits) and place in SAGADIR/cats
	
	4.  Casjobs:    To run casjobs, need to setup an account at:
         http://skyserver.sdss3.org/CasJobs/CreateAccount.aspx
	     Then edit your .bash_profile
              >  export CASJOBS_WSID='2090870927'   # get your WSID from site above
	          >export CASJOBS_PW='my password'


## Download sql files.
SDSS downloads of *all* objects within 1 degree of all hosts in host list
	> import download_host_sqlfile
	> download_host_sqlfile.run_query()


--  Create value-added SAGA base catalogs.
	> import create_base_catalogs
	> create_base_catalogs.run_query(nowise=1)

	If you want to include Lang WISE data, need to download _nw files
    into SAGA/nw_hosts/


-- Compile all SAGA spectra
     - download GAMA spectro catalog DR2: http://www.gama-survey.org/dr2/data/cat/SpecCat/v08/
