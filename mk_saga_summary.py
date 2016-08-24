#! /usr/bin/python
##################################################
#  COMPILE SPECTRA FROM ALL SOURCES FOR ALL FLAG 0 HOSTS
#       1. SDSS 
# 		2. GAMA
#       3. MMT, AAT, IMACS, WIYN 
#
#
#
##################################################

import numpy as np
import os
import re

from astropy.io import ascii

from astropy import table
from astropy.table import Table
from tabulate import tabulate
import pyspherematch as sm

from FileLoader import GoogleSheets


SAGA_DIR   = os.environ['SAGA_DIR']
REMOVELIST = GoogleSheets('1Y3nO7VyU4jDiBPawCs8wJQt2s_PIAKRj-HSrmcWeQZo', 1379081675, header_start=1)

############################################################################
def mk_saga_summary():


	# READ SPECTRA
	file = SAGA_DIR +'/data/saga_spectra_dirty.fits.gz'
	allspec = Table.read(file)


	# FIND UNIQUE SAGANAMES AND ORDER BY NSATS
	msk          = allspec['HOST_SAGA_NAME'] != '' 
	sagaspec     = allspec[msk]
	sorted_hosts = sort_saga_hosts(sagaspec)



	# FOR EACH HOST CALCULATE STUFF
	data = []
#	headers = ['SAGA Name', 'NSAID','RA','Dec','Nsats','MK','N_obj','gri (r<21)','gri (r< 20.5)','ML_{0.01,r< 21}','ML_{0.1,r< 21}']
	headers = ['Name', 'NSAID','RA','Dec','Nsats','N_corr','MK','N_obj','gri (r< 20.75)','ML_0.1']
	for host in sorted_hosts:

		print host[1], host[0]
		msk1   = sagaspec['HOST_SAGA_NAME'] == host[1] 
		msk2   = sagaspec['ZQUALITY'] >= 3
		msk3   = sagaspec['SATS'] != 3  # remove host

		msk   = msk1 & msk2 & msk3
		spec  = sagaspec[msk]
		nsaid = sagaspec['HOST_NSAID'][msk][0]
		dist = sagaspec['HOST_DIST'][msk][0]


		# READ BASE CATALOG
		basefile  = os.path.join(SAGA_DIR, 'base_catalogs','base_sql_nsa{0}.fits.gz'.format(nsaid))
		basetable = Table.read(basefile)	





		b_err1 = basetable['r_err'] < 10
#		b_err2 = basetable['r'] > 18
		b_rerr = b_err1 #& b_err2     #  this cut only below 18mag


		b_rmv  = basetable['REMOVE'] == -1
		b_rvir = basetable['RHOST_KPC'] <= 300
		b_fib  = basetable['FIBERMAG_R'] <= 23
		b_gal = basetable['PHOT_SG'] == 'GALAXY'

		b_r21  = basetable['r'] - basetable['EXTINCTION_R'] <= 20.75
		b_r205 = basetable['r'] - basetable['EXTINCTION_R'] <= 20.75



		nobj_all  =  np.sum(b_rmv)  # all objects
		nobj_rvir =  np.sum(b_rmv & b_rvir & b_gal & b_fib)
		nobj_gal  =  np.sum(b_rmv & b_rvir & b_gal &b_r21 & b_fib)
		print nobj_all, nobj_rvir, nobj_gal

		gri21   = gri_completeness(basetable[b_rmv &  b_rerr &b_rvir & b_gal & b_r21 & b_fib],sagaspec[msk])

		ML0001 = ML_completeness(basetable[b_rmv & b_rerr & b_rvir & b_gal & b_r21 & b_fib], sagaspec[msk], 0.0001)
		ML001 = ML_completeness(basetable[b_rmv & b_rerr & b_rvir & b_gal & b_r21 & b_fib], sagaspec[msk], 0.001)
		ML01 = ML_completeness(basetable[b_rmv & b_rerr & b_rvir & b_gal & b_r21 & b_fib], sagaspec[msk], 0.01)
		ML1 = ML_completeness(basetable[b_rmv & b_rerr & b_rvir & b_gal & b_r21 & b_fib], sagaspec[msk], 0.1)


		# ALSO DETERMINE NUMBER OF SATELLITE PASSING CRITERIA
		b_sats = basetable['SATS'] ==1
		nsats_corrected = b_sats & b_rmv & b_rerr & b_rvir & b_gal & b_r21
		ncorr = np.sum(nsats_corrected)



		row = []
		row.append(host[1])
		row.append(str(nsaid))
		row.append(sagaspec['HOST_RA'][msk][0])
		row.append(sagaspec['HOST_DEC'][msk][0])
		row.append(str(host[0]))
		row.append(ncorr)
		dmod =5*np.log10(sagaspec['HOST_DIST'][0]*1e6) - 5.
		mk = sagaspec['HOST_MK'][msk][0] 
		row.append(str('{:.1f}'.format(mk)))
		row.append(str('{:.1f}'.format(mk)))
		row.append(str(nobj_gal))
		row.append(gri21)


#		row.append(ML001)
		row.append(ML01)




		data.append(row)


	# WRITE OUT TABLE
	print tabulate(data, headers, tablefmt="rst")

	f = open('data_ML.tex', 'wb')
	f.write(tex_table_header())
	f.write('\n')
	f.write(tabulate(data, headers, tablefmt="latex"))
	f.write('\n')
	f.write(tex_table_footer())



def gri_completeness(base, spec):



	gmag = base['g'] - base['EXTINCTION_G']
	rmag = base['r'] - base['EXTINCTION_R']
	imag = base['i'] - base['EXTINCTION_I']

	gr = gmag - rmag
	ri = rmag - imag

	grerr = np.sqrt(base['g_err']**2 + base['r_err']**2)
	rierr = np.sqrt(base['r_err']**2 + base['i_err']**2)


	cgmr = gr - 2.*grerr
	crmi = ri - 2.*rierr


	msk1 = cgmr < 0.85
	msk2 = crmi < 0.55
	grimsk = msk1&msk2
	ngri = np.sum(msk1&msk2)


	id1,id2,d = sm.spherematch(base['RA'][grimsk], base['DEC'][grimsk],spec['RA'], spec['DEC'],1./3600,nnearest=1)
	ngri_spec = np.size(d)

	p = 100.*ngri_spec/ngri
	
	compl = '{:02.0f}% ({}/{})'.format(p,ngri_spec,ngri)
	return compl




def ML_completeness(base, spec, prob):

	ML_msk = base['PROBABILITY_CLASS1'] > prob
	nML = np.sum(ML_msk)


	id1,id2,d = sm.spherematch(base['RA'][ML_msk], base['DEC'][ML_msk],spec['RA'], spec['DEC'],1./3600,nnearest=1)
	nML_spec = np.size(d)

	p=0
	if (nML != 0): p = 100.*nML_spec/nML
	
	compl = '{:02.0f}% ({}/{})'.format(p,nML_spec,nML)
	return compl



def sort_saga_hosts(sagaspec):
	"""
	Find unique named SAGA hosts in allspec. 
	Sort names by nsats
	"""

	# FIND UNIQUE SAGA NAMES and CALCULATE NSATS
	unique_hosts = []
	nsats        = []
	for s in sagaspec:
		if s['HOST_SAGA_NAME'] not in unique_hosts: 
			unique_hosts.append(s['HOST_SAGA_NAME'])

			# CALCULATE NSATS FOR GIVEN HOST
			msk1 = sagaspec['HOST_SAGA_NAME'] == s['HOST_SAGA_NAME']
			msk2 = sagaspec['SATS'] == 1
			msk3 = sagaspec['REMOVE'] == -1
			msk = msk1 & msk2 & msk3
 			n = np.sum(msk)

			nsats.append([n,s['HOST_SAGA_NAME']])





#			print  s['HOST_SAGA_NAME'],sagaspec['RA'][msk],sagaspec['DEC'][msk]


	sorted_hosts = sorted(nsats,reverse=True)
	return sorted_hosts


###############################################
def tex_table_header():
	table_header = """
		\documentclass[8pt]{article} 
		\usepackage{natbib}
		\usepackage{fancyhdr}
		\usepackage{graphics}
		\pagestyle{fancy}
		\setlength{\ textwidth}{8.5in}
		\setlength{\ textheight}{8.30in}
		\setlength{\ topmargin}{20pt}
		\setlength{\oddsidemargin}{-1in}
		\setlength{\hoffset}{0pt}
		\setlength{\ voffset}{10pt}
		\ begin{document}
		\ begin{center}
		\ begin{table}


	"""
	table_header = re.sub(' ', '', table_header).strip()
#	table_header = re.sub('\s+', ' ', table_header).strip()

	return table_header



def tex_table_footer():
	table_foot = """

		\end{table}
		\end{center}
		\end{document}
	"""
	table_foot = re.sub(' ', '', table_foot).strip()
	return table_foot



if __name__ == '__main__':
    create_saga_spectra()

 

 
