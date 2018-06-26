import numpy as np
import argparse
import ellc

import gatspy
from gatspy import datasets, periodic
from gatspy.periodic import LombScargleMultiband, LombScargle, LombScargleFast, LombScargleMultibandFast

import astropy
from astropy import units, constants

import csv
import multiprocessing, logging

import pandas as pd

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

#OpSim (in OpSimRun1) - http://ops2.lsst.org/docs/current/architecture.html
import sqlite3
from astropy.coordinates import SkyCoord, ICRS

#for Katie's code
import scipy.special as ss
import astropy.stats as astroStats
from scipy.interpolate import interp1d
import scipy.stats
import time
import os

import GxSampleThinDisklogitAMG as GxSample

global do_plot, verbose, do_parallel, filters, sigma_sys, totaltime, cadence
do_plot = False 
verbose = False
do_parallel = False 
sigma_sys = 0.005  #systematic photometric error
years = 10.
totaltime = 365.* years
cadence = 3.

#filters = ['u_', 'g_', 'r_', 'i_', 'z_']
filters = ['u_', 'g_', 'r_', 'i_', 'z_', 'y_']

#This is copied from Katie's GxRealizationThinDisk.py, but here we've created a method
def LSSTsim(Ncores, Nbin):

	def paramTransform(dat):
		datMin = min(dat)-0.0001
		datMax = max(dat)+0.0001
		datZeroed = dat-datMin

		datTransformed = datZeroed/(datMax-datMin)
		return datTransformed

	# SET TIME TO TRACK COMPUTATION TIME
	##############################################################################
	start_time = time.time()
	#np.random.seed()

	# CONSTANTS
	##############################################################################
	G = 6.67384* 10.**-11.0
	c = 2.99792458* 10.**8.0
	parsec = 3.08567758* 10.**16.
	Rsun = 6.955* 10.**8.
	Msun = 1.9891* 10.**30.
	day = 86400.0
	rsun_in_au = 215.0954
	day_in_year = 365.242
	sec_in_day = 86400.0
	sec_in_hour = 3600.0
	hrs_in_day = 24.0
	sec_in_year = 3.15569*10**7.0
	Tobs = 3.15569*10**7.0
	geo_mass = G/c**2
	m_in_AU = 1.496*10**11.0


	
	mTotDisk = 2.15*10**10
	binID = '0012'

	##############################################################################
	# STELLAR TYPES - KW
	#
	#   0 - deeply or fully convective low mass MS star
	#   1 - Main Sequence star
	#   2 - Hertzsprung Gap
	#   3 - First Giant Branch
	#   4 - Core Helium Burning
	#   5 - First Asymptotic Giant Branch
	#   6 - Second Asymptotic Giant Branch
	#   7 - Main Sequence Naked Helium star
	#   8 - Hertzsprung Gap Naked Helium star
	#   9 - Giant Branch Naked Helium star
	#  10 - Helium White Dwarf
	#  11 - Carbon/Oxygen White Dwarf
	#  12 - Oxygen/Neon White Dwarf
	#  13 - Neutron Star
	#  14 - Black Hole
	#  15 - Massless Supernova
	##############################################################################


	# LOAD THE FIXED POPULATION DATA
	##############################################################################
	############## UNITS ################
	#                                   #
	# mass: Msun, porb: year, sep: rsun #
	#                                   #
	#####################################
			
	dts = {'names':('binNum','tBornBin','Time','tBorn1','tBorn2','commonEnv','id1','id2',
			'm1','m2','m1Init','m2Init','Lum1','Lum2','rad1','rad2',
			'T1','T2','massc1','massc2','radc1','radc2','menv1','menv2',
			'renv1','renv2','spin1','spin2','rrol1','rrol2','porb','sep','ecc'),
		'formats':('i','f','f','f','f','i','i','i',
			'f','f','f','f','f','f','f','f',
			'f','f','f','f','f','f','f','f',
			'f','f','f','f','f','f','f','f','f')}
			
	#FixedPop = np.loadtxt('fixedData_'+binID+'.dat', delimiter = ',', dtype=dts)
	FixedPop = pd.read_hdf('../input/dat_ThinDisk_12_0_12_0.h5', key='bcm')
	FixedPopLog = np.loadtxt('../input/fixedPopLogCm_'+binID+'.dat', delimiter = ',')
			
	# COMPUTE THE NUMBER AT PRESENT DAY NORMALIZED BY TOTAL MASS OF THE GX COMPONENT
	##############################################################################
	mTotFixed = sum(FixedPopLog[:,2])
	nPop = int(len(FixedPop)*mTotDisk/mTotFixed)
	print('The number of binaries in the Gx for: '+str(binID)+' is: '+str(nPop))
		
	# TRANSFORM THE FIXED POP DATA TO HAVE LIMITS [0,1] &
	# COMPUTE THE BINWIDTH TO USE FOR THE KDE SAMPLE; SEE KNUTH_BIN_WIDTH IN ASTROPY
	##############################################################################

	# UNITS: 
	# MASS [MSUN], ORBITAL PERIOD [LOG10(YEARS)], LUMINOSITIES [LSUN], RADII [RSUN]
	#FixedPop['m1'] = FixedPop['mass_1'] #or maybe some more efficient way
	###
	#print (FixedPop['m1'])

	FixedPop['m1'] = FixedPop['mass_1']
	#print (FixedPop['m1'])
	FixedPop['m2'] = FixedPop['mass_2']
	FixedPop['Lum1'] = FixedPop['lumin_1']
	FixedPop['Lum2'] = FixedPop['lumin_2']
	FixedPop['rad1'] = FixedPop['rad_1']
	FixedPop['rad2'] = FixedPop['rad_2']

	m1Trans = ss.logit(paramTransform(FixedPop['m1']))
	bwM1 = astroStats.scott_bin_width(m1Trans)
			
	m2Trans = ss.logit(paramTransform(FixedPop['m2']))
	bwM2 = astroStats.scott_bin_width(m2Trans)
			
	porbTrans = ss.logit(paramTransform(np.log10(FixedPop['porb'])))
	bwPorb = astroStats.scott_bin_width(porbTrans)
			
	Lum1Trans = ss.logit(paramTransform(FixedPop['Lum1']))
	bwLum1 = astroStats.scott_bin_width(FixedPop['Lum1'])

	Lum2Trans = ss.logit(paramTransform(FixedPop['Lum2']))
	bwLum2 = astroStats.scott_bin_width(FixedPop['Lum2'])
				
	# The eccentricity is already transformed, but only fit KDE to ecc if ecc!=0.0
	eIndex, = np.where(FixedPop['ecc']>1e-2)
	if len(eIndex) > 50:

		eccTrans = FixedPop['ecc']
		for jj in eccTrans.keys():
			if eccTrans[jj] > 0.999:
				eccTrans[jj] = 0.999
			elif eccTrans[jj] < 1e-4:
				eccTrans[jj] = 1e-4
		eccTrans = ss.logit(eccTrans)
		bwEcc = astroStats.scott_bin_width(eccTrans)
	else:
		bwEcc = 100.0

	rad1Trans = ss.logit(paramTransform(FixedPop['rad1']))
	bwRad1 = astroStats.scott_bin_width(rad1Trans)

	rad2Trans = ss.logit(paramTransform(FixedPop['rad2']))
	bwRad2 = astroStats.scott_bin_width(rad2Trans)
	#print(bwEcc,bwPorb,bwM1,bwM2,bwLum1,bwLum2,bwRad1,bwRad2)
	popBw = min(bwEcc,bwPorb,bwM1,bwM2,bwLum1,bwLum2,bwRad1,bwRad2)
			
	# GENERATE THE DATA LIST DEPENDING ON THE TYPE OF COMPACT BINARY TYPE
	##############################################################################
	type1Save = 1
	if type1Save < 14 and len(eIndex)>50:
		print('both bright stars and eccentric')
		datList = np.array((m1Trans, m2Trans, porbTrans, eccTrans, rad1Trans, rad2Trans, Lum1Trans, Lum2Trans))
	elif type1Save < 14 and len(eIndex)<50:
		print('both bright stars and circular')
		datList = np.array((m1Trans, m2Trans, porbTrans, rad1Trans, rad2Trans, Lum1Trans, Lum2Trans))
					
	# GENERATE THE KDE FOR THE DATA LIST
	##############################################################################\
	print(popBw)
	print(datList)
	for i,x in enumerate(datList):
		print(len(x))
		f = plt.figure()
		plt.plot(x)
		f.savefig(str(i)+".png")

	sampleKernel = scipy.stats.gaussian_kde(datList)#, bw_method=popBw)
				
			
	# CALL THE MONTE CARLO GALAXY SAMPLE CODE
	##############################################################################
	print('nSample: '+str(Nbin))
	output = multiprocessing.Queue()
	processes = [multiprocessing.Process(target = GxSample.GxSample, \
					args   = (x, Ncores, FixedPop, sampleKernel,\
					binID, Nbin, popBw, len(eIndex), Tobs, output)) \
					for x in range(Ncores)]
	for p in processes:
		p.start()

	for p in processes:
		p.join()

	nSRC = [output.get() for p in processes]
			
	print('The number of sources in pop '+binID+' is ', nSRC)

	gxDatTot = []
	for kk in range(Ncores):
		print(kk)
		if os.path.getsize('gxRealization_'+str(kk)+'_'+str(binID)+'.dat') > 0:
			gxReal = np.loadtxt('gxRealization_'+str(kk)+'_'+str(binID)+'.dat', delimiter = ',')
		else:
			gxReal = []
		if len(gxReal)>0:
			gxDatTot.append(gxReal)
		os.remove('gxRealization_'+str(kk)+'_'+str(binID)+'.dat')
	gxDatSave = np.vstack(gxDatTot)
	print(np.shape(gxDatSave))

	return gxDatSave
		
def getFieldID(db, myRA, myDEC, deglim = 3.5/2.):
	#uses RA/Dec (from galactic coordinates) to return locatiom's fieldID according to OpSim
	#field-of-view == 3.5-degree diameter (also returned with fieldFov key)

	cursor = db.cursor()
	cursor.execute("SELECT fieldid, fieldra, fielddec FROM field")
	c = np.array(cursor.fetchall())
	RA = c[:,1]
	Dec = c[:,2]
	dbCoord = SkyCoord(ra = RA*units.degree, dec = Dec*units.degree, frame='icrs')
	inCoord = SkyCoord(ra = myRA*units.degree, dec = myDEC*units.degree, frame='icrs')

	imin, sep2d, dist3d = inCoord.match_to_catalog_sky(dbCoord)

	dbID = (c[imin,0]).astype('int') 

	mask = np.where(sep2d.to(units.degree).value > deglim)

	#this check apparently isn't necessary because it looks like the entire sky is covered with fieldIDs, but I suppose some of these fieldIDs don't have any observation dates (in the northern hemisphere)
	if (len(mask[0]) > 0):
		print (mask[0])
		print("WARNING: coordinate outside LSST FOV", myRA[mask], myDec[mask])
		dbID[mask] = -999

	return dbID

def getSummaryCursor(db, dbID):
	#gets SQlite cursor to pull information from OpSim
	print("getting Summary cursor ...")
	cursor = db.cursor()
	cursor.execute("SELECT expDate, fieldid, filter FROM summary") #WHERE FieldID IN dbID")
	c = np.array(cursor.fetchall()) #NOTE: this takes a LONG time
	print("have cursor.")
	## have some reservations about how I pulled as there is fieldID in summary AND in field apparently ##
	return c

def getexpDate(c, dbID, filtin):
	#matches FieldID to existing OpSim database ID and matches observation filters to get dates (in seconds since the start of the 
	# survey)
	date = c[:,0]
	FieldID = c[:,1]
	filt = c[:,2]
	# print("filtin = ", filtin)
	# print("db stuff", dbID, FieldID, date, filt)
	# print("dbID = ", dbID)
	# print("FieldID = ", FieldID)
	# print("date = ", date)
	# print("filt = ", filt)

	# posID = np.where(FieldID == str(dbID))
	# print("posID = ", posID)
	# posFilt = np.where(filt == str(filtin[:-1]))
	# print("posFilt = ", posFilt)
	posIDFilt = np.where(np.logical_and(FieldID == str(dbID), filt == filtin[:-1]))
	print("posIDFilt = ", posIDFilt)

	#sort_idx = dbID.argsort()
	#OpSimdates1 = sort_idx[np.searchsorted(dbID,FieldID,sorter = sort_idx)]

	OpSimdates = posIDFilt[0]

	# OpSimdates = np.nonzero(np.in1d(np.nonzero(np.in1d(FieldID, str(dbID))), np.nonzero(np.in1d(filt, filtin[:-1]))))
	# if OpSimdates == np.nan:
	if (len(OpSimdates) < 1):
		return [None]
	else:


		#OpSimdates = np.where(np.logical_and(np.array(FieldID) == np.array(dbID), np.array(filt) == filtin))
		# print ('FieldID == dbID = ', np.nonzero(np.in1d(FieldID, dbID)))
		# print ('filt == filtin = ', np.nonzero(np.in1d(filt, filtin)))
		print('OpSimdates =', OpSimdates)
		#print('sort_idx OpSim = ', OpSimdates1)
		#print ('filt == filtin = ', np.where (np.array(filt)==np.array(filtin)))
		dates = np.array([float(d) for d in date[OpSimdates] ])/86400. #converting seconds to days
		return dates

def getRad(m):
	#use stellar mass to get stellar radius (not necessary, read out by Katie's work)
	if (m > 1):  #*units.solMass
		eta = 0.57
	else:
		eta = 0.8
	return (m)**eta #* units.solRad (units.solMass)

def calcTeff(L, R):
	#use stellar radius and stellar luminosity to get the star's effective temperature
	logTeff = 3.762 + 0.25*np.log10(L) - 0.5*np.log10(R) 
	return 10.**logTeff
	
def calclogg(m, L, T):
	#use stellar mass, luminosity, and effective temperature to get log(gravity)
	return np.log10(m) + 4.*np.log10(T) - np.log10(L) - 10.6071
	
def getLum(m):
	#use stellar mass to return stellar luminosit (not necessary, read out by Katie's work)
	if (m<0.43):
		cons = 0.23
		coeff = 2.3
	if m>0.43 and m<2.0:
		cons = 1
		coeff = 4
	if m>2.0 and m<20.0:
		cons = 1.5
		coeff = 3.5
	else:
		cons= 3200
		coeff = 1
	return cons*(m**coeff) 

def getSig2Rand(filt, magnitude):
	#returns 2 sigma random error baseed on the pass band (y-values may be wonky - need to check for seeing and 
	# against others)
	X = 1. #function of distance??
	t_vis = 30. #seconds
	if filt == 'u_':  
		gamma = 0.037
		seeing = 0.77
		m_sky = 22.9
		C_m = 22.92
		k_m = 0.451
		#m_5 = 23.68
	if filt == 'g_':
		gamma = 0.038
		seeing = 0.73
		m_sky = 22.3
		C_m = 24.29
		k_m = 0.163
		#m_5 = 24.89
	if filt == 'r_':
		gamma = 0.039
		seeing = 0.70
		m_sky = 21.2
		C_m = 24.33
		k_m = 0.087
		#m_5 = 24.43
	if filt == 'i_':
		gamma = 0.039
		seeing = 0.67
		m_sky = 20.5
		C_m = 24.20
		k_m = 0.065
		#m_5 = 24.00
	if filt == 'z_':
		gamma = 0.040
		seeing = 0.65
		m_sky = 19.6
		C_m = 24.07
		k_m = 0.043
		#m_5 = 24.45
	if filt == 'y_':
		gamma = 0.0039
		seeing = 0.65 #not sure where this is from - not in Ivezic; still the z value
		m_sky = 18.61
		C_m = 23.73
		k_m = 0.170
		#m_5 = 24.45 #not sure where this is from - not in Ivezic; still the z value
		#from Ivezic et al 2008 - https://arxiv.org/pdf/0805.2366.pdf - Table 2 (p26)

	m_5 = C_m + (0.50*(m_sky - 21.)) + (2.50*np.log10(0.7/seeing)) + (1.25*np.log10(t_vis/30.)) - (k_m*(X-1.))
	return (0.04 - gamma)*(10**(0.4*(magnitude - m_5))) + gamma*((10**(0.4*(magnitude - m_5)))**2)*(magnitude**2)


def ellc_gatspy_sim(j, t_zero, period, a, q, f_c, f_s, ld_1, ld_2, Teff1, Teff2, logg1, logg2, M_H, R_1, R_2, incl, sbratio, shape_1, shape_2, appmag, bnum, n_band, n_base, return_dict, plot_dict, db, dbID_OpSim):
	#this is the general simulation - ellc light curves and gatspy periodograms
	#global filters, sigma_sys, totaltime, cadence

	print("here in ellc_gatspy_sim")
	totalt = []
	totalmag = []
	totaldmag = []
	totalfilts = []
	

	if (verbose):
		print (j, 't_zero = ', t_zero)
		print (j, 'period = ', period)
		print (j, 'semimajor = ', a)
		print (j, 'massrat = ', q)
		print (j, 'f_c = ', f_c)
		print (j, 'f_s = ', f_s)
		print (j, 'ld_1 = ', ld_1)
		print (j, 'ld_2 = ', ld_2)

		print (j, 'radius_1 = ', R_1)
		print (j, 'radius_2 = ', R_2)
		print (j, 'incl = ', incl)
		print (j, 'sbratio = ', sbratio)
		print (j, 'Teff1 = ', Teff1)
		print (j, 'Teff2 = ', Teff2)

	if (do_plot):
		# for plotting the periodograms
		for_plotting = dict()
		pds = np.linspace(0.2, 2.*period, 10000)
		for_plotting['periods'] = pds
		for_plotting['seed'] = bnum 

	#set the random seed
	#np.random.seed(seed = seed)


	delta_mag = 0. 

	if (opsim):
		cursor = getSummaryCursor(db, dbID_OpSim)

	for i, filt in enumerate(filters):

		if (opsim):
			print("I'm running code with OpSim!")
			time = getexpDate(cursor, dbID_OpSim, filt)
			nobs = len(time)
			print("time_OpSim = ", time)
		else:
			print("I'm NOT running code with OpSim!")
			nobs = int(round(totaltime / (cadence * len(filters)))) 
			time = np.sort(totaltime * np.random.random(size=nobs))
			#time.sort() #maybe not needed and taking time
		drng = max(time) - min(time)

		# if (time_OpSim[0] != None):
		if (time[0] != None):

			filtellc = filt
			if (filt == 'y_'):
				filtellc = 'z_' #because we don't have limb darkening for y_
			##########################
			#Can we find limb darkening coefficients for y band??  (Then we wouldn't need this trick)
			##########################

			ldy_filt = ellc.ldy.LimbGravityDarkeningCoeffs(filtellc)
			a1_1, a2_1, a3_1, a4_1, y = ldy_filt(Teff1, logg1, M_H)
			a1_2, a2_2, a3_2, a4_2, y = ldy_filt(Teff2, logg2, M_H)
			ldc_1 = [a1_1, a2_1, a3_1, a4_1] 
			ldc_2 = [a1_2, a2_2, a3_2, a4_2]

			lc = ellc.lc(time,t_zero=t_zero, period=period, a=a, q=q,
				f_c=f_c, f_s=f_s,
				ld_1=ld_1,ldc_1=ldc_1,ld_2=ld_2,ldc_2=ldc_2,
				radius_1=R_1, radius_2=R_2,incl=incl,sbratio=sbratio, 
				shape_1=shape_1, shape_2=shape_2, grid_1='default',grid_2='default') 

			if (verbose):
				print (j, 'ldc_1 = ', ldc_1)
				print (j, 'ldc_2 = ', ldc_2)    
				print (j, 'time = ', time[0:10])
				print (j, 'lc = ', lc[0:10])

			if (min(lc) >= 0):
				magn = -2.5*np.log10(lc)
				magnitude = appmag + magn   
				#Ivezic 2008, https://arxiv.org/pdf/0805.2366.pdf , Table 2
				sigma2_rand = getSig2Rand(filt, magnitude)   #random photometric error
				sigma = ((sigma_sys**2.) + (sigma2_rand))**(1./2.)

				#now add the uncertaintly onto the magnitude
				#magnitude_obs = magnitude
				magnitude_obs = [np.random.normal(loc=x, scale=sig) for (x,sig) in zip(magnitude, sigma)]

				delta_mag = max([delta_mag, abs(min(magnitude_obs) - max(magnitude_obs))])

				if (verbose):   
					print(j, 'magn = ', magn[0:10])
					print(j, 'magnitude = ', magnitude[0:10])  
					print(j, 'magnitude_obs = ', magnitude_obs[0:10])  
					print(j, "sigma2_rand = ", sigma2_rand[0:10])
					print(j, "sigma = ", sigma[0:10])
					print(j, "delta_mag = ", delta_mag)

				t = np.array(time)
				totalt.extend(t)
				mag = np.array(magnitude_obs)
				totalmag.extend(mag)
				dmag = np.array(sigma)
				totaldmag.extend(dmag)
				filts = np.array(filt)
				specfilt = nobs*[filt]
				totalfilts.extend(specfilt)

				
				#model = LombScargle(fit_period = True)
				model = LombScargleFast(fit_period = True)
				model.optimizer.period_range = (0.2, drng)
				model.fit(t, mag, dmag)
				LSS = model.best_period
				print ('LSS running')
				return_dict[j] = return_dict[j] + [LSS]
				print ('LSS in return_dict')
				if (verbose): 
					print(j, "LSS = ",LSS)


				if (do_plot):
					if (i == 0):
						for_plotting['phase_obs'] = dict()
						for_plotting['mag_obs'] = dict()
						for_plotting['mag'] = dict()
						for_plotting['scores'] = dict()
						for_plotting['LSS'] = dict()


					phase_obs = np.array([(tt % period)/period for tt in time])
					scores = model.score(pds)
					for_plotting['phase_obs'][filt] = phase_obs
					for_plotting['mag_obs'][filt] = magnitude_obs
					for_plotting['mag'][filt] = magnitude
					for_plotting['scores'][filt] = scores
					for_plotting['LSS'][filt] = LSS

			else:
				return_dict[j] = return_dict[j] + [-999.]
				if (verbose):
					print(j, "bad fluxes", filt)
					#raise Exception('stopping')




def afromP(m1, m2, P):
	#returns the semimajor axis from the period and stellar masses
	return (((P**2.) * constants.G * (m1 + m2) / (4*np.pi**2.))**(1./3.)).decompose().to(units.AU)

def Eggleton_RL(q,a):
#Eggleton (1983) Roche Lobe
#assuming synchronous rotation
#but taking the separation at pericenter
	return a*0.49*q**(2./3.)/(0.6*q**(2./3.) + np.log(1. + q**(1./3.)))

def make_plots(plotdict):
	#colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
	#print([ matplotlib.colors.to_hex(c) for c in colors])
	colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

	f, ax = plt.subplots(len(filters)+1, 2)

	pds = plotdict['periods']
	LSM = plotdict['LSM']
	period = plotdict['pvalue']
	P_multi = plotdict['P_multi']
	seed = plotdict['seed']
	print("plotting line ",seed)
	for ii,filt in enumerate(filters):
		phase_obs = plotdict['phase_obs'][filt]
		mag_obs = plotdict['mag_obs'][filt]
		mag = plotdict['mag'][filt]
		scores = plotdict['scores'][filt]
		LSS = plotdict['LSS'][filt]

		sx = np.argsort(phase_obs)
		ax[ii][0].plot(phase_obs[sx], np.array(mag_obs)[sx], 'o', mfc='none', mec = colors[ii])
		ax[ii][0].plot(phase_obs[sx], np.array(mag)[sx], color = "black")
		ax[ii][0].set_ylim(ax[ii][0].get_ylim()[::-1])
		ax[ii][0].set_xlim(0,1)
		ax[ii][0].set_ylabel(filt)
		ax[ii][0].set_xticklabels([])

		ax[ii][1].plot(pds, scores, color = colors[ii])
		ax[ii][1].plot([LSS,LSS],[0,1], color = "dimgray", lw = 2)
		ax[ii][1].plot([period,period],[0,1],'--', color = "black")
		ax[ii][1].set_xlim(0, 2.*period)
		ax[ii][1].set_ylim(0, max(scores))
		ax[ii][1].set_xticklabels([])
		ax[ii][1].set_yticklabels([])

	plt.locator_params(axis='y', nticks=2)
	ii = len(filters)
	ax[ii][1].plot(pds, P_multi, color = colors[ii])
	ax[ii][1].set_xlim(0, 2.*period)
	ax[ii][1].set_ylim(0, max(P_multi))
	ax[ii][1].plot([period,period],[0,1],'--', color = "black")
	ax[ii][1].plot([LSM,LSM],[0,1], ':', color = "dimgray")

	f.subplots_adjust(hspace=0.1, wspace=0.1)
	f.delaxes(ax[ii][0])
	ax[ii-1][0].set_xlabel("phase")
	ax[ii][1].set_xlabel("period (days)")
	ax[ii][1].set_yticklabels([])

	f.savefig("lc_gatspy_fig_"+str(seed).rjust(10,'0')+".png", bbox_inches='tight')



def define_args(parser):
	parser.add_argument("-n", "--n_procs", type=int, help="Number of proceses [0]")
	parser.add_argument("-c", "--n_bin", type=int, help="Number of binaries [100000]")
	parser.add_argument("-p", "--plots", action='store_true', help="Set to create plots")
	parser.add_argument("-v", "--verbose", action='store_true', help="Set to show verbose output")
	parser.add_argument("-i", "--input_file", type=str, help="input file name")
	parser.add_argument("-o", "--output_file", type=str, help="output file name")
	parser.add_argument("-a", "--n_band", type=int, help="Nterms_band input for gatspy [2]")
	parser.add_argument("-b", "--n_base", type=int, help="Nterms_base input for gatspy [2]")
	parser.add_argument("-s", "--seed", type=int, help="random seed []")
	parser.add_argument("-l", "--opsim", action='store_true', help="set to run LSST OpSim, else run nobs =")

if __name__ == "__main__":
	#global do_plot, verbose, do_parallel, filters, sigma_sys, totaltime, cadence

	db = sqlite3.connect('../db/minion_1016_sqlite.db')


	#https://docs.python.org/2/howto/argparse.html
	n_appmag_error = 0
	n_incl_error = 0
	n_period_error = 0
	n_R_error = 0
	n_totalrun = 0

	parser = argparse.ArgumentParser()
	define_args(parser)
	args = parser.parse_args()
	#to print out the options that were selected (probably some way to use this to quickly assign args)
	opts = vars(args)
	options = { k : opts[k] for k in opts if opts[k] != None }
	print(options)

	ofile = 'output_file.csv'
	infile ='/projects/p30137/LSSTmsPop/gxRealization_0012_0.dat'##

	n_bin = 100000
	n_band = 2
	n_base = 2  
	n_procs = 0
	if (args.n_procs is not None):
		n_procs = args.n_procs
		if (n_procs > 1):
			do_parallel = True
	if (n_procs < 1):
		n_procs = 1
		do_parallel = False 

	if (args.n_bin is not None):
		n_bin = args.n_bin

	if (args.input_file is not None):
		infile = args.input_file
		
	if (args.output_file is not None):
		ofile = args.output_file

	if (args.n_band is not None):
		n_band = args.n_band
	if (args.n_base is not None):
		n_base = args.n_base

	do_plot = args.plots
	verbose = args.verbose
	opsim = args.opsim

	#set the random seed
	if (args.seed is not None):
		np.random.seed(seed = args.seed)
	else:
		np.random.seed()

#####################
#run Katie's code
	gxDat = LSSTsim(n_procs, n_bin)
	print(gxDat)
##############################

	manager = multiprocessing.Manager()
	return_dict = manager.dict()
		
	plot_dict = manager.dict()


	#logger = multiprocessing.log_to_stderr()
	##logger.setLevel(logging.DEBUG)
	#logger.setLevel(multiprocessing.SUBDEBUG)



	csvfile = open(ofile, 'wt')
	csvwriter = csv.writer(csvfile, delimiter=',',
								quotechar='|', quoting=csv.QUOTE_MINIMAL)
	csvwriter.writerow(['PERIOD', 'MASS_1', 'MASS_2', 'RADIUS_1', 'RADIUS_2', 'a', 'INCLINATION', 'MIN_INCLINATION', 'xGx', 'yGx', 'zGx', 'dist_kpc', 'eccentricity', 'max(app_magnitude)', 'appmag_error', 'inclination_error', 'period_error', 'radius_error', 'u_LSS_PERIOD', 'g_LSS_PERIOD', 'r_LSS_PERIOD', 'i_LSS_PERIOD', 'z_LSS_PERIOD', 'LSM_PERIOD', 'delta_mag', 'chi2', 'mean(dmag)'])
	

	jobs = []
	j = 0

	for i, line in enumerate(gxDat):
		m_1 = line[0]
		m_2 = line[1]
		porb = line[2]
		ecc = line[3]
		r_1 = line[4]
		r_2 = line[5]
		L_1 = line[6]
		L_2 = line[7]
		inc = line[12]
		om = line[13]
		dist_kpc = line[11]
		xGx = line[8] #parsec
		yGx = line[9] #pc
		zGx = line[10] #pc
		coord = SkyCoord(x=xGx, y=yGx, z=zGx, unit='pc', representation='cartesian', frame='galactocentric')
		ra_dec_dist = coord.icrs
		ra_dec_dist_array = ra_dec_dist.to_string().split(" ")

		myRA = float(ra_dec_dist_array[0])
		myDEC = float(ra_dec_dist_array[1])

		period = 10.** (porb)
		t_zero = np.random.random() * period

		dbID_OpSim = getFieldID(db, myRA, myDEC, deglim = 3.5/2.)
				
		dist_solrad = dist_kpc *44334448006.9 
				
		appmag = -26.74 - 2.5*np.log10(((L_1 + L_2)/1.)*((4.84814e-9/dist_kpc)**2))
			
		sbratio = L_2/L_1
		a = afromP(m_1*units.solMass, m_2*units.solMass, period*units.day).to(units.solRad).value

		R_1 = (r_1/a)
		R_2 = (r_2/a)

		R_1e = r_1/Eggleton_RL(m_1/m_2, a * (1. - ecc))
		R_2e = r_2/Eggleton_RL(m_2/m_1, a * (1. - ecc))

		M_H = 0
		Teff1 = min(50000., max(3500., calcTeff(L_1, r_1)))
		Teff2 = min(50000., max(3500., calcTeff(L_2, r_2)))
		logg1 = min(5., max(0., calclogg(m_1, L_1, Teff1)))
		logg2 = min(5., max(0., calclogg(m_2, L_2, Teff2)))
				
		ld_1 = 'claret'
		ld_2 = 'claret'
			
		incl = (inc * 180)/ np.pi
		min_incl = (90. - 2.*(np.arctan((((r_1 + r_2)/(a)))) * (180/np.pi)))
				
		runthis = True
		n_totalrun += 1
		if appmag <= 11. or appmag>= 24.:
			appmag_error = 1
			n_appmag_error += 1 
			runthis = False
		else:
			appmag_error = 0
				
		if incl <= min_incl:
			incl_error = 1
			n_incl_error += 1
			runthis = False
		else:
			incl_error = 0
					
		if period >= 2 * totaltime:
			period_error = 1
			n_period_error +=1
			runthis = False
		else:
			period_error = 0
				
		if R_1 <= 0 or R_1 >=1 or R_1e <= 0 or R_1e >=1 or R_2 <=0 or R_2 >= 1 or R_2e <= 0 or R_2e >=1:
			R_error = 1
			n_R_error += 1
			runthis = False
		else:
			R_error = 0
		

		
		q = m_2/m_1
		shape_1 = 'sphere'
		shape_2 = 'sphere'

		f_c = np.sqrt(ecc)*np.cos(om*np.pi/180.)
		f_s = np.sqrt(ecc)*np.sin(om*np.pi/180.)
				#removed delta_mag from output, now removing max(appmag), min(appmag) in favor of appmag because appmag is float rather than array???#
		output = [period, m_1, m_2, R_1, R_2, a, incl, min_incl, xGx, yGx, zGx, dist_kpc, ecc, appmag, appmag_error, incl_error, period_error, R_error]
			
				#print(i, len(jobs), appmag_error, incl_error, period_error, R_1_error, R_2_error)#, r_1, a, ecc, R_1, R_1e, R_2, R_2e)

		if (runthis):
			return_dict[j] = output
			plot_dict[j] = dict()
			p = multiprocessing.Process(target=ellc_gatspy_sim, args=(j, t_zero, period, a, q, f_c, f_s, ld_1, ld_2, Teff1, Teff2, logg1, logg2, M_H, R_1, R_2, incl, sbratio, shape_1, shape_2, appmag, i, n_band, n_base, return_dict, plot_dict, db, dbID_OpSim))
			jobs.append(p)
			j += 1
		else:
			for filt in filters:
				output.append(-999) #LSS
			output.append(-999) #LSM
			csvwriter.writerow(output)

			############################## 
		if (len(jobs) == n_procs or (i >= n_bin and len(jobs) > 0)):
			#print ('n_appmag_errorGx = ', n_appmag_errorGx)
			print ('n_appmag_error = ',  n_appmag_error)
			#print ('n_incl_errorGx = ', n_incl_errorGx)
			print ('n_incl_error = ', n_incl_error)
			#print ('n_period_errorGx = ', n_period_errorGx)
			print ('n_period_error = ', n_period_error)
			#print ('n_R_errorGx = ', n_R_errorGx)
			print ('n_R_error = ', n_R_error)
			#print ('n_totalrunGx = ', n_totalrunGx)
			print ('n_totalrun = ', n_totalrun
				)
				### makes printed 1/0 redundant - possible we want to remove ###
			if (do_parallel):
				for k,job in enumerate(jobs):
					print("starting job ",k)
					x = job.start()
				for k,job in enumerate(jobs):
					print("joining job ",k)
					job.join()
			else:
				for k,job in enumerate(jobs):
					print("running job",k)
					job.run()

			for j in range(n_procs):
				csvwriter.writerow(return_dict[j])
				csvfile.flush()
				if (do_plot):
					if ('LSM' in plot_dict[j]):
						 make_plots(plot_dict[j])


			 #raise Exception('stopping')
			jobs = []
			j = 0

		if (i > n_bin):
		   break

	csvfile.close()


