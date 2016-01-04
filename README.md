# marla
This directory contains scripts to create SAGA base catalogs and
compile SAGA spectra into a single file.

#### Google docs:
There are three google docs used in the scripts and these must be maintained to run the base catalogs:

	SAGA_HOST_LISTS:  Contains the full SAGA host list and the Flag Zero list.
	REMOVE_LIST:  Contains a list of objects which have been flagged by hand as junk.
	SAGA_HOSTS+SATELLITES:  Official file of SAGA common names

#### Directories and File Setup:
These local directories and files are require to run scripts.

--  Environment variables:

	SAGA_DIR:           Top level directory where SAGA data will live on local disk
	SAGA_DROPBOX:  Local directory for DropBox/SAGA

-- Create local directories and download data files, [GAMA DR2](http://www.gama-survey.org/dr2/data/cat/SpecCat/v08/), [NSA catalog](http://nsatlas.org/data):

	SAGA_DIR/hosts:    Houses the sql downloads and base files
	SAGA_DIR/cats:      Local copy of outside data catalogs
	Download NSA catalog (nsa_v0_1_2.fits) and place in SAGADIR/cats
	Download GAMA catalog (GAMA_SpecObj.fits) and place in SAGADIR/cats



-- Casjobs:  [http://skyserver.sdss3.org/CasJobs/CreateAccount.aspx](http://skyserver.sdss3.org/CasJobs/CreateAccount.aspx)

	To run casjobs, need to setup an account, then edit your .bash_profile as:
              > export CASJOBS_WSID='2090870927'   # get your WSID from site above
	          > export CASJOBS_PW='my password'


#### Download SQL files:
SDSS downloads of *all* objects within 1 degree of all hosts in host list

	> import download_host_sqlfile
	> download_host_sqlfile.run_query()

This typically takes a few minutes per host to download from the SDSS
casjobs site.   

#### Create SAGA Base Catalogs
Create value-added SAGA base catalogs.

	> import create_base_catalogs
	> create_base_catalogs.run_query(nowise=1)

If you want to include Lang WISE data, need to download _nw files into SAGA/nw_hosts/, otherwise, set nowise = 1

To create a single host catalog (not run full host list):

	> create_base_catalogs.run_single_host(nsaid)

####Compile all SAGA spectra
Compile all SDSS, GAMA and SAGA-aquired spectra, write to file

	> import compile_saga_spectra
   > compile_saga_spectra
