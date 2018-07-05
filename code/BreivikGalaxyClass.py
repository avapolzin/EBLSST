import math
import scipy.special as ss
import scipy.stats
import multiprocessing 
import numpy as np
import os
import time
import pandas as pd
import astropy.stats as astroStats

#This is copied from Katie's GxRealizationThinDisk.py, but here we've created a Class
#This will draw the binaries from the Galaxy realization
class BreivikGalaxyClass(object):

	def __init__(self, *args,**kwargs):
		self.verbose = False
		self.GalaxyFile ='../input/dat_ThinDisk_12_0_12_0.h5' #for Katie's model
		self.GalaxyFileLogPrefix ='../input/fixedPopLogCm_'

		self.n_bin = 100000
		self.n_cores = 1
		self.popID = '0012'
		self.seed = None

	def GxSample(self, x, pop, sampleKernel, bw, nEcc, Tobs, output):
		def untransform(dat, datSet):
			datMin = min(datSet)
			datMax = max(datSet)
			
			datUntransformed = dat*(datMax-datMin)
			datUnzeroed = datUntransformed+datMin
			
			return datUnzeroed  

	  
		# CONSTANTS
		##############################################################################
		G = 6.67384*math.pow(10, -11.0)
		c = 2.99792458*math.pow(10, 8.0)
		parsec = 3.08567758*math.pow(10, 16)
		Rsun = 6.955*math.pow(10, 8)
		Msun = 1.9891*math.pow(10,30)
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
		
	  
		##### UNITS EXPECTED FOR POPULATION ###################
		# mass: Msun, orbital period: days, Tobs: seconds     #
		#######################################################
		
		# seed the random generator
		########################################
		if (self.seed == None):
			np.random.seed()
		else:
			np.random.seed(seed = self.seed)

		# solar coordinates in the galaxy: in parsecs from 
		# (Chaper 8 of Galactic Structure and stellar Pops book) Yoshii (2013)
		############################################################################
		x_sun = 8000.0
		y_sun = 0.0
		z_sun = 25
		
		# sample the number of binaries given by nBin
		#################################################
		power = []
		freq = []
		n = 0
		nWrite = 0
		eccBin = 0
		
		# open the file to save the gx data
		#################################################
		gxFile = 'gxRealization_'+str(x)+'_'+str(self.popID)+'.dat'
		binDat = [] 
		
		nSample = int(self.n_bin/float(self.n_cores))
		dataSample = sampleKernel.resample(nSample)
			
		# check to see if the population has BHs or is ecc
		###########################################################
		   
		m1T = ss.expit(dataSample[0,:])
		m2T = ss.expit(dataSample[1,:])
		porbT = ss.expit(dataSample[2,:])
		eccT = ss.expit(dataSample[3,:])
		rad1T = ss.expit(dataSample[4,:])
		rad2T = ss.expit(dataSample[5,:])
		Lum1T = ss.expit(dataSample[6,:])
		Lum2T = ss.expit(dataSample[7,:])
			
		m1 = untransform(m1T, pop['m1'])
		m2 = untransform(m2T, pop['m2'])
		porb = untransform(porbT, np.log10(pop['porb']))
		ecc = eccT             
		ii=0
		for e in ecc:
			if e < 0:
				ecc[ii] = abs(e)
			elif e > 1:
				ecc[ii] = 1-(ecc-1)
			ii+=1
		rad1 = 10**(untransform(rad1T, pop['rad1']))
		rad2 = 10**(untransform(rad2T, pop['rad2']))
		Lum1 = 10**(untransform(Lum1T, pop['Lum1']))
		Lum2 = 10**(untransform(Lum2T, pop['Lum2']))
		
		# COMPUTE THE POSITION AND ORIENTATION OF EACH BINARY
		##############################################################################
		# First assign the position relative to the galactic center
		# Assign the radial position of the binary
		norm_r = 1/2.5
		a_0_r = np.random.uniform(0, 1.0, len(m1))
		r = -2.5*np.log(1-a_0_r)
		# Assign the azimuthal position of the star
		phi = np.random.uniform(0, 2*np.pi, len(m1))
		# Assign the z position of the star with r as a parameter
		norm_zr = 0.71023
		a_0_zr = np.random.uniform(0, 1.0, len(m1))
		z = 1/1.42*np.arctanh(a_0_zr/(0.704)-1)
		# convert to cartesian and parsecs
		xGX = r*np.cos(phi)*1000.0
		yGX = r*np.sin(phi)*1000.0
		zGX = z*1000.0
		# compute the distance to Earth/LISA/us in kiloparsecs
		dist = ((xGX-x_sun)**2+(yGX-y_sun)**2+(zGX-z_sun)**2)**(1/2.0)
		dist_kpc = dist/1000.0
			
		inc = np.arccos(2.*np.random.uniform(0,1.0,len(m1)) - 1.)
		OMEGA = np.random.uniform(0,2*math.pi,len(m1))
		omega = np.random.uniform(0,2*math.pi,len(m1))

		binDat = np.vstack((m1, m2, porb, ecc, rad1, rad2, Lum1, Lum2, xGX, yGX, zGX, dist_kpc, inc, OMEGA, omega)).T
		radTotAU = (rad1+rad2)/rsun_in_au
		radAng = radTotAU/dist
		binEclipseIndex, = np.where(radAng>inc*4.8e-6)
	 
		np.savetxt(gxFile, binDat, delimiter = ',')     
				
		#gxFile.close() 
		output.put(np.shape(binDat)) 


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
				
		FixedPop = pd.read_hdf(self.GalaxyFile, key='bcm')
		FixedPopLog = np.loadtxt(self.GalaxyFileLogPrefix+self.popID+'.dat', delimiter = ',')
				
		# COMPUTE THE NUMBER AT PRESENT DAY NORMALIZED BY TOTAL MASS OF THE GX COMPONENT
		##############################################################################
		mTotFixed = sum(FixedPopLog[:,2])
		nPop = int(len(FixedPop)*mTotDisk/mTotFixed)
		print('The number of binaries in the Gx for: '+str(self.popID)+' is: '+str(nPop))
			
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
		processes = [multiprocessing.Process(target = self.GxSample, \
						args = (x, FixedPop, sampleKernel, popBw, len(eIndex), Tobs, output)) \
						for x in range(self.n_cores)]
		for p in processes:
			p.start()

		for p in processes:
			p.join()

		nSRC = [output.get() for p in processes]
				
		print('The number of sources in pop '+self.popID+' is ', nSRC)

		gxDatTot = []
		for kk in range(self.n_cores):
			if os.path.getsize('gxRealization_'+str(kk)+'_'+str(self.popID)+'.dat') > 0:
				gxReal = np.loadtxt('gxRealization_'+str(kk)+'_'+str(self.popID)+'.dat', delimiter = ',')
			else:
				gxReal = []
			if len(gxReal)>0:
				gxDatTot.append(gxReal)
			os.remove('gxRealization_'+str(kk)+'_'+str(self.popID)+'.dat')
		gxDatSave = np.vstack(gxDatTot)

		return gxDatSave
