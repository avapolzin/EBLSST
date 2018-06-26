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

		self.Ncores = 1
		self.Nbin = 1000

		self.ofile = 'output_file.csv' #output file name
		self.GalaxyFile ='../input/dat_ThinDisk_12_0_12_0.h5' #for Katie's model


		self db = sqlite3.connect('../db/minion_1016_sqlite.db')


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
		print('nSample: '+str(self.Nbin))
		output = multiprocessing.Queue()
		processes = [multiprocessing.Process(target = GxSample.GxSample, \
						args   = (x, self.Ncores, FixedPop, sampleKernel,\
						binID, self.Nbin, popBw, len(eIndex), Tobs, output)) \
						for x in range(self.Ncores)]
		for p in processes:
			p.start()

		for p in processes:
			p.join()

		nSRC = [output.get() for p in processes]
				
		print('The number of sources in pop '+binID+' is ', nSRC)

		gxDatTot = []
		for kk in range(self.Ncores):
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



	#run everything
	def runAll(self):

		#get the summary cursor
		getSummaryCursor()