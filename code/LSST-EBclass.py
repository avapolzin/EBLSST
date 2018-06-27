### I need to pull out gatspy into it's own function for easier plug-and-play
### WE NEED TO DOUBLE CHECK getsig2rand

#python libraries
import numpy as np
import argparse
import time
import os
from astropy import units, constants
import astropy.stats as astroStats
from astropy.coordinates import SkyCoord, ICRS
import csv
import multiprocessing, logging
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
#OpSim (in OpSimRun1) - http://ops2.lsst.org/docs/current/architecture.html
import sqlite3
import scipy.special as ss
from scipy.interpolate import interp1d
import scipy.stats

#3rd party codes
import ellc
import gatspy
from gatspy import datasets, periodic
from gatspy.periodic import LombScargleMultiband, LombScargle, LombScargleFast, LombScargleMultibandFast
import GxSampleThinDisklogitAMG as GxSample

class LSST-EBclass(object):

    def __init__(self, *args,**kwargs):

    	self.do_plot = False 
		self.verbose = False
		self.do_parallel = False 
		self.sigma_sys = 0.005  #systematic photometric error
		self.years = 10.
		self.totaltime = 365.* years
		self.cadence = 3.

		#self.filters = ['u_', 'g_', 'r_', 'i_', 'z_']
		self.filters = ['u_', 'g_', 'r_', 'i_', 'z_', 'y_']

		self.n_bin = 100000
		self.n_band = 2
		self.n_base = 2  
		self.n_cores = 0

		self.ofile = 'output_file.csv' #output file name
		self.GalaxyFile ='../input/dat_ThinDisk_12_0_12_0.h5' #for Katie's model


		self.plot_dict = None
		self.return_dict = None

		self db = sqlite3.connect('../db/minion_1016_sqlite.db')

		self.parser = argparse.ArgumentParser()


	#Some approximate function for deriving stellar parameters
	def getRad(self, m):
		#use stellar mass to get stellar radius (not necessary, read out by Katie's work)
		if (m > 1):  #*units.solMass
			eta = 0.57
		else:
			eta = 0.8
		return (m)**eta #* units.solRad (units.solMass)

	def calcTeff(self, L, R):
		#use stellar radius and stellar luminosity to get the star's effective temperature
		logTeff = 3.762 + 0.25*np.log10(L) - 0.5*np.log10(R) 
		return 10.**logTeff
	
	def calclogg(self,m, L, T):
		#use stellar mass, luminosity, and effective temperature to get log(gravity)
		return np.log10(m) + 4.*np.log10(T) - np.log10(L) - 10.6071
	
	def getLum(self, m):
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

	def afromP(self, m1, m2, P):
		#returns the semimajor axis from the period and stellar masses
		return (((P**2.) * constants.G * (m1 + m2) / (4*np.pi**2.))**(1./3.)).decompose().to(units.AU)

	def Eggleton_RL(self, q,a):
	#Eggleton (1983) Roche Lobe
	#assuming synchronous rotation
	#but taking the separation at pericenter
		return a*0.49*q**(2./3.)/(0.6*q**(2./3.) + np.log(1. + q**(1./3.)))

	#This is copied from Katie's GxRealizationThinDisk.py, but here we've created a method
	#This will draw the binaries from the Galaxy realization
	def LSSTsim(self):

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
		FixedPop = pd.read_hdf(self.GalaxyFile, key='bcm')
		FixedPopLog = np.loadtxt('fixedPopLogCm_'+binID+'.dat', delimiter = ',')
				
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
		##############################################################################
		if (self.verbose):
			print(popBw)
			print(datList)

		sampleKernel = scipy.stats.gaussian_kde(datList)#, bw_method=popBw)
					
				
		# CALL THE MONTE CARLO GALAXY SAMPLE CODE
		##############################################################################
		print('nSample: '+str(self.n_bin))
		output = multiprocessing.Queue()
		processes = [multiprocessing.Process(target = GxSample.GxSample, \
						args = (x, self.n_cores, FixedPop, sampleKernel,\
						binID, self.n_bin, popBw, len(eIndex), Tobs, output)) \
						for x in range(self.n_cores)]
		for p in processes:
			p.start()

		for p in processes:
			p.join()

		nSRC = [output.get() for p in processes]
				
		print('The number of sources in pop '+binID+' is ', nSRC)

		gxDatTot = []
		for kk in range(self.n_cores):
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

	#database manipulation
	def getSummaryCursor(self):
		#gets SQlite cursor to pull information from OpSim
		print("getting Summary cursor ...")
		cursor = self.db.cursor()
		cursor.execute("SELECT fieldid, fieldra, fielddec, expDate, filter FROM summary") 
		self.cursor = np.array(cursor.fetchall()) #NOTE: this takes a LONG time
		print("have cursor.")


	def getFieldID(self, myRA, myDEC, deglim = 3.5/2.):
		#uses RA/Dec (from galactic coordinates) to return locatiom's fieldID according to OpSim
		#field-of-view == 3.5-degree diameter (also returned with fieldFov key)

		cursor = self.db.cursor()
		RA = self.cursor[:,1]
		Dec = self.cursor[:,2]
		dbCoord = SkyCoord(ra = RA*units.degree, dec = Dec*units.degree, frame='icrs')
		inCoord = SkyCoord(ra = myRA*units.degree, dec = myDEC*units.degree, frame='icrs')

		imin, sep2d, dist3d = inCoord.match_to_catalog_sky(dbCoord)

		dbID = (self.cursor[imin,0]).astype('int') 

		mask = np.where(sep2d.to(units.degree).value > deglim)

		#this check apparently isn't necessary because it looks like the entire sky is covered with fieldIDs, but I suppose some of these fieldIDs don't have any observation dates (in the northern hemisphere)
		if (len(mask[0]) > 0):
			print (mask[0])
			print("WARNING: coordinate outside LSST FOV", myRA[mask], myDec[mask])
			dbID[mask] = -999

		return dbID

	def getexpDate(self, dbID, filtin):
		#matches FieldID to existing OpSim database ID and matches observation filters to get dates (in seconds since the start of the 
		# survey)
		date = c[:,3]
		FieldID = c[:,0]
		filt = c[:,4]

		posIDFilt = np.where(np.logical_and(FieldID == str(dbID), filt == filtin[:-1]))
		if (self.verbose):
			print("posIDFilt = ", posIDFilt)

		OpSimdates = posIDFilt[0]

		if (len(OpSimdates) < 1):
			return [None]
		else:
			if (self.verbose):
				print('OpSimdates =', OpSimdates)
			dates = np.array([float(d) for d in date[OpSimdates] ])/86400. #converting seconds to days
			return dates

	def getSig2Rand(self, filt, magnitude):
		#returns 2 sigma random error based on the pass band (y-values may be wonky - need to check for seeing and 
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


	def make_plots(self):
		#colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
		#print([ matplotlib.colors.to_hex(c) for c in colors])
		colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

		f, ax = plt.subplots(len(self.filters)+1, 2)

		pds = self.plot_dict['periods']
		LSM = self.plot_dict['LSM']
		period = self.plot_dict['pvalue']
		P_multi = self.plot_dict['P_multi']
		seed = self.plot_dict['seed']
		print("plotting line ",seed)
		for ii,filt in enumerate(self.filters):
			phase_obs = self.plot_dict['phase_obs'][filt]
			mag_obs = self.plot_dict['mag_obs'][filt]
			mag = self.plot_dict['mag'][filt]
			scores = self.plot_dict['scores'][filt]
			LSS = self.plot_dict['LSS'][filt]

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
		ii = len(self.filters)
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


	def define_args(self):
		self.parser.add_argument("-n", "--n_cores", 	type=int, help="Number of cores [0]")
		self.parser.add_argument("-c", "--n_bin", 		type=int, help="Number of binaries [100000]")
		self.parser.add_argument("-i", "--gal_file", 	type=str, help="Galaxy input file name")
		self.parser.add_argument("-o", "--output_file", type=str, help="output file name")
		self.parser.add_argument("-a", "--n_band", 		type=int, help="Nterms_band input for gatspy [2]")
		self.parser.add_argument("-b", "--n_base", 		type=int, help="Nterms_base input for gatspy [2]")
		self.parser.add_argument("-s", "--seed", 		type=int, help="random seed []")
		self.parser.add_argument("-p", "--plots", 		action='store_true', help="Set to create plots")
		self.parser.add_argument("-v", "--verbose", 	action='store_true', help="Set to show verbose output")
		self.parser.add_argument("-l", "--opsim", 		action='store_true', help="set to run LSST OpSim, else run nobs =")


	def apply_args(self):

		args = self.parser.parse_args()
		#to print out the options that were selected (probably some way to use this to quickly assign args)
		opts = vars(args)
		options = { k : opts[k] for k in opts if opts[k] != None }
		print(options)

		if (args.n_cores is not None):
			self.n_cores = args.n_cores
			if (n_cores > 1):
				self.do_parallel = True
		if (n_cores < 1):
			self.n_cores = 1
			self.do_parallel = False 

		if (args.n_bin is not None):
			self.n_bin = args.n_bin

		if (args.gal_file is not None):
			self.GalaxyFile = args.gal_file
			
		if (args.output_file is not None):
			self.ofile = args.output_file

		if (args.n_band is not None):
			self.n_band = args.n_band
		if (args.n_base is not None):
			self.n_base = args.n_base

		self.do_plot = args.plots
		self.verbose = args.verbose
		self.do_opsim = args.do_opsim

		#set the random seed
		if (args.seed is not None):
			np.random.seed(seed = args.seed)
		else:
			np.random.seed()

#
	#run everything
	def runAll(self):

		self.define_args()
		self.apply_args()

		#get the summary cursor
		getSummaryCursor()