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

import matplotlib.pyplot as plt
from scipy import integrate


from astropy import table
from astropy.table import Table
from astropy import units as u
from astropy.io import ascii



from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=9, usetex=True)


from utils import get_mcconn_table

SAGA_DIR   = os.environ['SAGA_DIR']

############################################################################
def plot_mwsats():


	# READ SPECTRA
	file = SAGA_DIR +'/data/saga_spectra_dirty.fits.gz'
	allspec = Table.read(file)

	m = allspec['SATS'] == 1
	sats = allspec[m]

   # MCCONNACHIE DATA
	mcconndat = get_mcconn_table()
	sample = mcconndat['distance']<0
	sample[mcconndat['name']=='LMC'] = True
	sample[mcconndat['name']=='SMC'] = True
	sample[mcconndat['name']=='Fornax'] = True
	sample[mcconndat['name']=='Leo I'] = True
	ldat = mcconndat[sample]


	sagadists = [5, 20, 30, 40] * u.Mpc
	distmods = 5*np.log10(sagadists/(10*u.pc))

	data = ldat[ldat['Vabs']<-10]
	data = data[np.argsort(data['Vabs'])]

	fig = plt.figure(figsize=(7, 3))
	fig.subplots_adjust(bottom=0.15, top=0.9,left=0.12, right=0.95, wspace=0.3)


    # PLOT SIZES
	plt.subplot(1,2,1)
	for gal in data:
		plt.plot(sagadists.value, ((gal['rh_phys']/sagadists.to(u.kpc).value)*u.rad).to(u.arcsec).value,'-o',label=gal['name'])
	plt.plot(sats['HOST_DIST'],sats['PETRORAD_R'],'r.',label='SAGA sats')

	plt.ylim(0,20)
	plt.xlim(15,42)
	plt.xlabel('D [Mpc]')
	plt.ylabel('Radius [arcsec]')

	plt.axhline(1.5, c='k', ls='--')

    # PLOT SAGA DATA

	######################	
	plt.subplot(1,2,2)


	Vmr=0.3
	rad = 3*u.arcsec/2. #3" diam = fibermag

	cc = []
	for _ in range(5):
	    cc.extend(plt.rcParams['axes.color_cycle'])

	#first do plummer
	for i, gal in enumerate(data):
	    fracs = np.array([frac_inside(gal, rad, dist, 'plummer') for dist in sagadists])
	    fibmags = gal_mag_frac(gal, fracs, sagadists) - Vmr
	    totmags = gal_mag_frac(gal, np.ones_like(fracs), sagadists) - Vmr
	    
	    plt.plot(fibmags, totmags,'-o',color=cc[i],label=gal['name'])
	    
#	        plt.text(fibmags[fi], totmags[ti], gal['name'],color=cc[i])
	#plt.legend(loc=0)


	  # PLOT SAGA DATA
	plt.plot(sats['FIBERMAG_R'],sats['r'],'r.',label='SAGA sats')

	plt.legend(loc=2,fontsize=6)

	for i, gal in enumerate(data):
	    fracs = np.array([frac_inside(gal, rad, dist, 'exp') for dist in sagadists])
	    fibmags = gal_mag_frac(gal, fracs, sagadists) - Vmr
	    totmags = gal_mag_frac(gal, np.ones_like(fracs), sagadists) - Vmr
	    
	    plt.plot(fibmags, totmags,'--o',color=cc[i],label=gal['name'])

	    
#	x1t1 = [min(fiblim[0], totlim[0]), max(fiblim[1], totlim[1])]
#	plt.fill_between(x1t1,x1t1, totlim[1],facecolor=[0.6]*3)
	plt.axhline(21, c='k', ls=':')
	plt.axvline(23, c='k', ls=':')


	    
	plt.xlabel('Approximate fibermag\_r')
	plt.ylabel('r')
	plt.title('Solid=Plummer, Dashed=Exponential')
	plt.xlim(16,24) #tweak for text
	plt.ylim(12,22)






	plt.show()
	plt.savefig('fig_mwsats.pdf')

def plummer_prof(r, totflux=1, Reff=1):
	a = Reff*(2.**(2./3)-1.)**(-0.5)
	return 2*a*(3 * totflux / (4*np.pi*a**3))*(1+(r/a)**2)**-2.5


def exp_prof(r, totflux=1, Reff=1):
    from scipy.special import lambertw
    #the parenthetical below is non-trivial, but it's ~1.67835
    rs = Reff / -(lambertw(-1/2/np.exp(1),-1).real+1.)
    rs = Reff/1.67835
    print rs
    return totflux*np.exp(-r/rs)/np.pi/2./rs/rs


def frac_inside(gal, radius, dist, prof):
    """
    Compute the magnitude inside the given angular `radius` for the
    table entry `gal` at a distance of `dist`, assuming a profile
    `prof` ('plummer' or 'exp')
    """    
    if prof == 'plummer':
        prof_func = plummer_prof
    elif prof == 'exp':
        prof_func = exp_prof
    else:
        raise ValueError('invalid profile type '+ str(prof))
        
    rh_at_distance = (u.arcmin*gal['rh'] * (gal['distance']*u.kpc)/dist).to(u.arcmin)
    return integrate.quad(lambda x:prof_func(x, 1, rh_at_distance.value)*x*np.pi*2, 0, radius.to(u.arcmin).value)[0]

def gal_mag_frac(gal, frac, dist):
    """
    determine the V magnitude for `frac` of the flux from `gal` at distance `dist`
    """
    flux = 10**((gal['Vabs'])/-2.5)
    fluxd = flux * ((10*u.pc/dist)**2).decompose().value
    return -2.5*np.log10(frac*fluxd)



if __name__ == '__main__':
    create_saga_spectra()

 

 
