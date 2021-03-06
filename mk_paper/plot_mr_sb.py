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
	file = SAGA_DIR +'/data/saga_spectra_clean.fits.gz'
	allspec = Table.read(file)


	# FIND UNIQUE SAGANAMES AND ORDER BY NSATS
	msk          = allspec['HOST_SAGA_NAME'] != '' 
	sagaspec     = allspec[msk]
	sorted_hosts = sort_saga_hosts(sagaspec)



	# FOR EACH HOST CALCULATE STUFF
	data = []
	headers = ['SAGA Name', 'NSAID','RA','Dec','Nsats','MK','N_obj','gri (r<21)','gri (r< 20.5)']
	for host in sorted_hosts:

		print host[1]
		msk1   = sagaspec['HOST_SAGA_NAME'] == host[1] 
		msk2   = sagaspec['ZQUALITY'] >= 3
		msk   = msk1 & msk2
		spec  = sagaspec[msk]
		nsaid = sagaspec['HOST_NSAID'][msk][0]

		# READ BASE CATALOG
		basefile  = os.path.join(SAGA_DIR, 'hosts','base_sql_nsa{0}.fits.gz'.format(nsaid))
		basetable = Table.read(basefile)	

#		b_crap  = basetable['REMOVE'] != -1
		b_rmv  = basetable['REMOVE'] == -1
		b_rvir = basetable['RHOST_KPC'] < 300
		b_gal = basetable['PHOT_SG'] == 'GALAXY'
		b_r21  = basetable['r'] < 21.0
		b_r205 = basetable['r'] < 20.5



#		n_crap    =  np.sum(b_crap)  # all objects
		nobj_all  =  np.sum(b_rmv)  # all objects
		nobj_rvir =  np.sum(b_rmv & b_rvir)
		nobj_gal  =  np.sum(b_rmv & b_rvir & b_gal &b_r21)
		print nobj_all, nobj_rvir, nobj_gal

		gri205  = gri_completeness(basetable[b_rmv & b_rvir & b_gal & b_r205],sagaspec[msk])
		gri21   = gri_completeness(basetable[b_rmv & b_rvir & b_gal & b_r21],sagaspec[msk])



		row = []
		row.append(host[1])
		row.append(str(nsaid))
		row.append(sagaspec['HOST_RA'][msk][0])
		row.append(sagaspec['HOST_DEC'][msk][0])
		row.append(str(host[0]))
		row.append(str('{:.1f}'.format(sagaspec['HOST_MK'][msk][0])))
		row.append(str(nobj_gal))
		row.append(gri21)
		row.append(gri205)


#		ML1 = ML_completeness(0.1)
#		ML2 = ML_completeness(0.5)


		data.append(row)


	# WRITE OUT TABLE
	print tabulate(data, headers, tablefmt="rst")

	f = open('data_new.tex', 'wb')
	f.write(tex_table_header())
	f.write('\n')
	f.write(tabulate(data, headers, tablefmt="latex"))
	f.write('\n')
	f.write(tex_table_footer())



def gri_completeness(base, spec):

	gmr  = base['g'] - base['r'] - 2*base['g_err'] - 2*base['r_err']
	rmi  = base['r'] - base['i'] - 2*base['r_err'] - 2*base['i_err']

	msk1 = gmr < 1.0
	msk2 = rmi < 0.5
	grimsk = msk1&msk2
	ngri = np.sum(msk1&msk2)


	id1,id2,d = sm.spherematch(base['RA'][grimsk], base['DEC'][grimsk],spec['RA'], spec['DEC'],1./3600,nnearest=1)
	ngri_spec = np.size(d)

	p = 100.*ngri_spec/ngri
	
	compl = '{:02.0f}% ({}/{})'.format(p,ngri_spec,ngri)
	return compl




def ML_completeness(base, spec, prob):

	ML_msk = base['PCLASS_1'] > prob
	nML = np.sum(ML_msk)


	id1,id2,d = sm.spherematch(base['RA'][ML_msk], base['DEC'][ML_msk],spec['RA'], spec['DEC'],1./3600,nnearest=1)
	nML_spec = np.size(d)

	p = 100.*nML_spec/nML
	
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
			msk = msk1 & msk2
			n = np.sum(msk)

			nsats.append([n,s['HOST_SAGA_NAME']])


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

 

 
