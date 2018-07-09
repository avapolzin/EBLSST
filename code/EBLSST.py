##Have Adam double check the conversion from bolometric to apparent magnitude

import math
import scipy.special as ss
import scipy.stats
from scipy.interpolate import interp1d
import multiprocessing 
import logging
import numpy as np
import os
import time
import pandas as pd
import csv
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
#OpSim (in OpSimRun1) - http://ops2.lsst.org/docs/current/architecture.html
import sqlite3

import astropy.stats as astroStats
from astropy import units, constants
from astropy.coordinates import SkyCoord, ICRS

#3rd party codes
import ellc
import gatspy
from gatspy import datasets, periodic
from gatspy.periodic import LombScargleMultiband, LombScargle, LombScargleFast, LombScargleMultibandFast

#only using a small portion of vespa, to get the A_V value, but NOTE vespa also can access TRILEGAL galaxy model...
import vespa
#extinction will allow me to convert A_V to any wavelength.  Not sure which reference is best.  I will use ccm89, for now. 
#import extinction
#could use this instead, seems to be linked more closely to astropy : https://dust-extinction.readthedocs.io/en/latest/index.html
#pip install git+https://github.com/karllark/dust_extinction.git
from dust_extinction.parameter_averages import F99

filters = ['u_', 'g_', 'r_', 'i_', 'z_', 'y_']


class EclipsingBinary(object):
	def __init__(self, *args,**kwargs):

		self.MbolSun = 4.73 #as "recommended" for Flower's bolometric correction in http://iopscience.iop.org/article/10.1088/0004-6256/140/5/1158/pdf

		#these is defined by the user
		self.m1 = None #*units.solMass
		self.m2 = None #*units.solMass
		self.r1 = None #*units.solRad
		self.r2 = None #*units.solRad
		self.L1 = None #*units.solLum
		self.L2 = None #*units.solLum
		self.period = None #*units.day 
		self.eccentricity = None
		self.omega = None
		self.inclination = None
		self.t_zero = None
		self.dist = None #*units.kpc
		self.xGx = None #*units.parsec
		self.yGx = None #*units.parsec
		self.zGx = None #*units.parsec
		self.lineNum = 0
		self.verbose = False
		self.RV = 3.1

		#from https://www.lsst.org/scientists/keynumbers
		#in nm
		self.wavelength = {
			'u_': (324. + 395.)/2.,
			'g_': (405. + 552.)/2.,
			'r_': (552. + 691.)/2.,
			'i_': (691. + 818.)/2.,
			'z_': (818. + 921.)/2.,
			'y_': (922. + 997. )/2.
		}

		#these will be calculated after calling self.initialize()
		self.RL1 = None
		self.RL2 = None
		self.T1 = None
		self.T2 = None
		self.T12 = None
		self.L1 = None
		self.L2 = None
		self.g1 = None
		self.g2 = None
		self.a = None 
		self.q = None
		self.f_c = None
		self.f_s = None
		self.R_1 = None
		self.R_2 = None
		self.sbratio = None
		self.RA = None
		self.Dec = None
		self.Mbol = None
		self.AV = None
		self.appMagMean = dict()
		self.appMagMeanAll = None
		self.absMagMean = dict()
		self.Ared = dict()
		self.BC = dict()

		#for light curves
		self.SED = None
		self.filters = filters
		self.M_H = 0
		self.ld_1 = 'claret'
		self.ld_2 = 'claret'
		self.grid_1 = 'default'
		self.grid_2 = 'default'
		self.shape_1 = 'sphere'
		self.shape_2 = 'sphere'
		self.sigma_sys = 0.005  #systematic photometric error
		self.obsDates = dict()
		self.appMag = dict()
		self.appMagObs = dict()
		self.appMagObsErr = dict()
		self.deltaMag = dict()
		self.maxDeltaMag = 0.
		self.doOpSim = True
		self.fieldCursor = None
		self.summaryCursor = None
		self.observable = True
		self.appmag_failed = 0
		self.incl_failed = 0
		self.period_failed = 0
		self.radius_failed = 0
		self.OpSimID = None
		self.years = 10.
		self.totaltime = 365.* self.years
		self.cadence = 3.
		self.Nfilters = 6.
		self.nobs = 0
		#this is for the magnitude uncertainties
		self.sigmaDict = {
			'u_': {
				'gamma'	: 0.037,
				'seeing': 0.77,
				'm_sky'	: 22.9,
				'C_m' 	: 22.92,
				'k_m' 	: 0.451},
			'g_': {
				'gamma'	: 0.038,
				'seeing': 0.73,
				'm_sky'	: 22.3,
				'C_m'	: 24.29,
				'k_m'	:0.163},
			'r_': {
				'gamma'	: 0.039,
				'seeing': 0.70,
				'm_sky'	: 21.2,
				'C_m'	: 24.33,
				'k_m'	: 0.087},
			'i_': {
				'gamma'	: 0.039,
				'seeing': 0.67,
				'm_sky'	: 20.5,
				'C_m'	: 24.20,
				'k_m'	: 0.065},
			'z_': {
				'gamma'	: 0.040,
				'seeing': 0.65,
				'm_sky'	: 19.6,
				'C_m'	: 24.07,
				'k_m'	: 0.043},
			'y_': {
		#from Ivezic et al 2008 - https://arxiv.org/pdf/0805.2366.pdf - Table 2 (p26)
				'gamma'	: 0.0039,
				'seeing': 0.65, #not sure where this is from - not in Ivezic; still the z value
				'm_sky'	: 18.61,
				'C_m'	: 23.73,
				'k_m'	: 0.170}
		}


		#set within the "driver" code, for gatspy
		self.LSS = dict()
		self.LSSmodel = dict()
		self.LSM = -999.
		self.LSMmodel = None

		self.seed = None




	def getFlowerBCV(self, Teff):
		#from http://iopscience.iop.org/article/10.1088/0004-6256/140/5/1158/pdf
		#which updates/corrects from Flower, P. J. 1996, ApJ, 469, 355
		lt = np.log10(Teff)
		a = 0.
		b = 0.
		c = 0.
		d = 0.
		e = 0.
		f = 0.
		if (lt < 3.7):
			a = -0.190537291496456*10.**5.
			b =  0.155144866764412*10.**5.
			c = -0.421278819301717*10.**4.
			d =  0.381476328422343*10.**3.
		if (lt >= 3.7 and lt < 3.9):
			a = -0.370510203809015*10.**5.
			b =  0.385672629965804*10.**5.
			c = -0.150651486316025*10.**5.
			d =  0.261724637119416*10.**4.
			e = -0.170623810323864*10.**3.
		if (lt >= 3.9):
			a = -0.118115450538963*10.**6.
			b =  0.137145973583929*10.**6.
			c = -0.636233812100225*10.**5.
			d =  0.147412923562646*10.**5.
			e = -0.170587278406872*10.**4.
			f =  0.788731721804990*10.**2.


		BCV = a + b*lt + c*lt**2. + d*lt**3. + e*lt**4 + f*lt**5.

		return BCV

	#Some approximate function for deriving stellar parameters
	def getRad(self, m):
		#(not needed with Katie's model, but included here in case needed later)
		#use stellar mass to get stellar radius (not necessary, read out by Katie's work)
		if (m > 1):  #*units.solMass
			eta = 0.57
		else:
			eta = 0.8
		return (m)**eta #* units.solRad (units.solMass)

	def getTeff(self, L, R):
		#use stellar radius and stellar luminosity to get the star's effective temperature
		logTeff = 3.762 + 0.25*np.log10(L) - 0.5*np.log10(R) 
		return 10.**logTeff
	
	def getlogg(self,m, L, T):
		#use stellar mass, luminosity, and effective temperature to get log(gravity)
		return np.log10(m) + 4.*np.log10(T) - np.log10(L) - 10.6071
	
	def getLum(self, m):
		#(not needed with Katie's model, but included here in case needed later)
		#use stellar mass to return stellar luminosity (not necessary, read out by Katie's work)
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

	def getafromP(self, m1, m2, P):
		#returns the semimajor axis from the period and stellar masses
		return (((P**2.) * constants.G * (m1 + m2) / (4*np.pi**2.))**(1./3.)).decompose().to(units.AU)

	def Eggleton_RL(self, q,a):
		#Eggleton (1983) Roche Lobe
		#assuming synchronous rotation
		#but taking the separation at pericenter
		return a*0.49*q**(2./3.)/(0.6*q**(2./3.) + np.log(1. + q**(1./3.)))


	def setLightCurve(self, filt, t_vis=30., X=1.):
		def getSig2Rand(filt, magnitude):
			#returns 2 sigma random error based on the pass band (y-values may be wonky - need to check for seeing and 
			# against others)
			#X = 1. #function of distance??
			#t_vis = 30. #seconds

			m_5 = self.sigmaDict[filt]['C_m'] + (0.50*(self.sigmaDict[filt]['m_sky'] - 21.)) + (2.50*np.log10(0.7/self.sigmaDict[filt]['seeing'])) + (1.25*np.log10(t_vis/30.)) - (self.sigmaDict[filt]['k_m']*(X-1.))
			return (0.04 - self.sigmaDict[filt]['gamma'])*(10**(0.4*(magnitude - m_5))) + self.sigmaDict[filt]['gamma']*((10**(0.4*(magnitude - m_5)))**2)*(magnitude**2)

		self.appMagObs[filt] = [None]
		self.deltaMag[filt] = [None]

		#in case the user did not initialize
		if (self.T1 == None):
			self.initialize()

		#limb darkenning
		##########################
		#Can we find limb darkening coefficients for y band??  (Then we wouldn't need this trick)
		##########################
		filtellc = filt
		if (filt == 'y_'):
			filtellc = 'z_' #because we don't have limb darkening for y_
		ldy_filt = ellc.ldy.LimbGravityDarkeningCoeffs(filtellc)
		a1_1, a2_1, a3_1, a4_1, y = ldy_filt(self.T1, self.g1, self.M_H)
		a1_2, a2_2, a3_2, a4_2, y = ldy_filt(self.T2, self.g2, self.M_H)
		ldc_1 = [a1_1, a2_1, a3_1, a4_1] 
		ldc_2 = [a1_2, a2_2, a3_2, a4_2]

		#light curve
		lc = ellc.lc(self.obsDates[filt], ldc_1=ldc_1, ldc_2=ldc_2, 
			t_zero=self.t_zero, period=self.period, a=self.a, q=self.q,
			f_c=self.f_c, f_s=self.f_s, ld_1=self.ld_1,  ld_2=self.ld_2,
			radius_1=self.R_1, radius_2=self.R_2, incl=self.inclination, sbratio=self.sbratio, 
			shape_1=self.shape_1, shape_2=self.shape_2, grid_1=self.grid_1,grid_2=self.grid_2) 
		if (min(lc) > 0):

			magn = -2.5*np.log10(lc)
			self.appMag[filt] = self.appMagMean[filt] + magn   
			#Ivezic 2008, https://arxiv.org/pdf/0805.2366.pdf , Table 2
			sigma2_rand = getSig2Rand(filt, self.appMag[filt])   #random photometric error
			self.appMagObsErr[filt] = ((self.sigma_sys**2.) + (sigma2_rand))**(1./2.)

			#now add the uncertainty onto the magnitude
			self.appMagObs[filt] = np.array([np.random.normal(loc=x, scale=sig) for (x,sig) in zip(self.appMag[filt], self.appMagObsErr[filt])])

			self.deltaMag[filt] = abs(min(self.appMagObs[filt]) - max(self.appMagObs[filt]))

	#For OpSim database
	def getFieldID(self, myRA, myDEC, deglim = 3.5/2.):
		#uses RA/Dec (from galactic coordinates) to return locatiom's fieldID according to OpSim
		#field-of-view == 3.5-degree diameter (also returned with fieldFov key)

		RA = self.fieldCursor[:,1].astype(float)
		Dec = self.fieldCursor[:,2].astype(float)
		dbCoord = SkyCoord(ra = RA*units.degree, dec = Dec*units.degree, frame='icrs')
		inCoord = SkyCoord(ra = myRA*units.degree, dec = myDEC*units.degree, frame='icrs')

		imin, sep2d, dist3d = inCoord.match_to_catalog_sky(dbCoord)

		dbID = (self.fieldCursor[imin,0]).astype('int') 

		mask = np.where(sep2d.to(units.degree).value > deglim)
		#this check apparently isn't necessary because it looks like the entire sky is covered with fieldIDs, but I suppose some of these fieldIDs don't have any observation dates (in the northern hemisphere)
		if (len(mask[0]) > 0):
			print(mask[0])
			print("WARNING: coordinate outside LSST FOV", myRA[mask], myDec[mask])
			dbID[mask] = -999

		if (self.verbose):
			print("have Field ID", dbID)

		return dbID

	def getOpSimDates(self, filtin):
		#matches FieldID to existing OpSim database ID and matches observation filters to get dates (in seconds since the start of the 
		# survey)
		FieldID = self.summaryCursor[:,0].astype('int')
		date = self.summaryCursor[:,1].astype('float')
		filt = self.summaryCursor[:,2]

		posIDFilt = np.where(np.logical_and(FieldID == self.OpSimID, filt == filtin[:-1]))

		if (self.verbose):
			print("posIDFilt = ", posIDFilt, filtin)

		OpSimdates = posIDFilt[0]

		if (len(OpSimdates) < 1):
			return [None]
		else:
			if (self.verbose):
				print('OpSimdates =', OpSimdates)
			dates = np.array([float(d) for d in date[OpSimdates] ])/86400. #converting seconds to days
			return dates

	def checkIfObservable(self):

		if (self.appMagMeanAll <= 11. or self.appMagMeanAll >= 24.):
			self.appmag_failed = 1
		
		incl = self.inclination * 180/np.pi
		min_incl = 90. - np.arctan2(self.r1 + self.r2, 2.*self.a)*180./np.pi
		if (incl <= min_incl):
			self.incl_failed = 1
		
		if (self.period >= 2. * self.totaltime):
			self.period_failed = 1
				
		if (self.R_1 <= 0 or self.R_1 >=1 or self.R_2 <=0 or self.R_2 >= 1 or self.R_1e >=1 or self.R_2e >=1):
			self.radius_failed = 1
			
		if (self.radius_failed or self.period_failed or self.incl_failed or self.appmag_failed):
			self.observable = False

	def initializeSeed(self):
		if (self.seed == None):
			np.random.seed()
		else:
			np.random.seed(seed = self.seed)

	def initialize(self):

		#should I initialize the seed here?  
		#No I am initializing it in LSSTEBworker.py
		#self.initializeSeed()

		self.q = self.m2/self.m1
		self.T1 = min(50000., max(3500., self.getTeff(self.L1, self.r1)))
		self.T2 = min(50000., max(3500., self.getTeff(self.L2, self.r2)))
		self.g1 = min(5., max(0., self.getlogg(self.m1, self.L1, self.T1)))
		self.g2 = min(5., max(0., self.getlogg(self.m2, self.L2, self.T2)))
		self.a = self.getafromP(self.m1*units.solMass, self.m2*units.solMass, self.period*units.day).to(units.solRad).value
		self.f_c = np.sqrt(self.eccentricity)*np.cos(self.omega*np.pi/180.)
		self.f_s = np.sqrt(self.eccentricity)*np.sin(self.omega*np.pi/180.)
		self.R_1 = (self.r1/self.a)
		self.R_2 = (self.r2/self.a)
		self.sbratio = self.L2/self.L1
		self.R_1e = self.r1/self.Eggleton_RL(self.m1/self.m2, self.a * (1. - self.eccentricity))
		self.R_2e = self.r2/self.Eggleton_RL(self.m2/self.m1, self.a * (1. - self.eccentricity))

		#estimate a combined Teff value, as I do in the N-body codes (but where does this comes from?)
		logLb = np.log10(self.L1 + self.L2)
		logRb = 0.5*np.log10(self.r1**2. + self.r2**2.)
		self.T12 = 10.**(3.762 + 0.25*logLb - 0.5*logRb)
		#print(self.L1, self.L2, self.T1, self.T2, self.T12)

		coord = SkyCoord(x=self.xGx, y=self.yGx, z=self.zGx, unit='pc', representation='cartesian', frame='galactocentric')
		self.RA = coord.icrs.ra.to(units.deg).value
		self.Dec = coord.icrs.dec.to(units.deg).value

		#one option for getting the exinction
		self.AV = vespa.stars.extinction.get_AV_infinity(self.RA, self.Dec, frame='icrs')
		ext = F99(Rv=self.RV)
		for f in self.filters:
			#print(extinction.fitzpatrick99(np.array([self.wavelength[f]*10.]), self.AV, self.RV, unit='aa')[0] , ext(self.wavelength[f]*units.nm)*self.AV)
			#self.Ared[f] = extinction.fitzpatrick99(np.array([self.wavelength[f]*10.]), self.AV, self.RV, unit='aa')[0] #or ccm89
			self.Ared[f] = ext(self.wavelength[f]*units.nm)*self.AV

		self.Mbol = self.MbolSun - 2.5*np.log10(self.L1 + self.L2)

		######################
		#Now I'm assuming that the magnitude is bolometric -- but we should account for redenning in each filter
		#######################
		# self.absMagMean = self.Mbol
		# self.appMagMean = self.absMagMean + 5.*np.log10(self.dist*100.)  #multiplying by 1000 to get to parsec units

		self.appMagMeanAll = 0.
		for f in self.filters:
			#BCV = self.getFlowerBCV(self.T12)
			self.SED.getBCf(self.T12*units.K) #now, how do I use this??
			self.absMagMean[f] = self.Mbol# - self.BC[f]
			self.appMagMean[f] = self.absMagMean[f] + 5.*np.log10(self.dist*100.) + self.Ared[f]  #multiplying by 1000 to get to parsec units
			self.LSS[f] = -999.
			self.appMagMeanAll += self.appMagMean[f]
		self.appMagMeanAll /= len(self.filters)

		#check if we can observe this (not accounting for the location in galaxy)
		self.checkIfObservable()

		#if we're using OpSim, then get the field ID
		#get the field ID number from OpSim where this binary would be observed
		if (self.doOpSim and self.observable):
			self.OpSimID = self.getFieldID(self.RA, self.Dec)

	def observe(self, filt):

		#get the observation dates
		if (self.doOpSim):
			self.obsDates[filt] = self.getOpSimDates(filt)
		else:
			nobs = int(round(self.totaltime / (self.cadence * self.Nfilters)))
			self.obsDates[filt] = np.sort(self.totaltime * np.random.random(size=nobs))

		self.nobs += len(self.obsDates[filt])
		#get the light curve, and related information
		self.setLightCurve(filt)

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################


#This is copied from Katie's GxRealizationThinDisk.py, but here we've created a Class
#This will draw the binaries from the Galaxy realization
class BreivikGalaxy(object):

	def __init__(self, *args,**kwargs):
		self.verbose = False
		self.GalaxyFile ='../input/dat_ThinDisk_12_0_12_0.h5' #for Katie's model
		self.GalaxyFileLogPrefix ='../input/fixedPopLogCm_'

		self.n_bin = 100000
		self.n_cores = 4
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
		np.random.seed()

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


###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################


class SED(object):
	def __init__(self, *args,**kwargs):
		self.filterFilesRoot = '../input/filters/'
		#['u_band_Response.dat','g_band_Response.dat','r_band_Response.dat','i_band_Response.dat','z_band_Response.dat','y_band_Response.dat']
		self.filters = filters
		self.filterThroughput = dict()

		self.BCf = dict()

	def readFilters(self):
		for f in self.filters:
			# #https://github.com/lsst-pst/syseng_throughputs/tree/master/components/camera/filters
			# fname = self.filterFilesRoot + f + 'band_Response.dat'
			# df = pd.read_csv(fname, delim_whitespace=True, header=None, names=['w','t'])

			#https://github.com/lsst/throughputs/tree/master/baseline
			fname = self.filterFilesRoot + 'filter_'+f[0]+'.dat' 
			df = pd.read_csv(fname, delim_whitespace=True, header=None, names=['w','t'], skiprows = 6)
			self.filterThroughput[f] = {'w':df['w'].values*units.nm, 't':df['t'].values}

	def getBCf(self, T, dw = 0.1, wmin =10, wmax = 10000):

		#could improve this with a Kurucz model
		def blackbody(w, T):
			#erg/s/cm^2/AA/steradian
			Bl = 2.*constants.h*constants.c**2./w**5. / (np.exp( constants.h*constants.c / (w*constants.k_B*T)) -1.)
			return Bl.decompose()

		w = np.arange(wmin, wmax, dw)*units.nm
		f = blackbody(w,T).value
		#norm = np.sum(f)*dw
		norm = max(f)
		fn = f/norm
		
		ftot = np.sum(fn)

		for f in self.filters:
			ft = np.interp(w, self.filterThroughput[f]['w'], self.filterThroughput[f]['t'])
			# plt.plot(w,ft,'--')
			tot = np.sum(fn*ft)
			self.BCf[f] = tot/ftot

		#print(self.BCf)

		# plt.plot(w,fn)
		# for f in self.filters:
		# 	plt.plot(self.filterThroughput[f]['w'], self.filterThroughput[f]['t'])
		# plt.xlim(100, 2000)
		# plt.show()

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

class LSSTEBworker(object):

	def __init__(self, *args,**kwargs):

		#NOTE: these need to be defined on the command line.  The default, if not defined, will be False
		self.do_plot = False 
		self.verbose = False
		self.doOpSim = False

		self.useFast = True
		self.doLSM = True
		self.do_parallel = False 

		self.years = 10.
		self.totaltime = 365.* self.years
		self.cadence = 3.

		self.filters = filters

		self.n_bin = 100000
		self.n_band = 2
		self.n_base = 2  
		self.n_cores = 1

		self.ofile = 'output_file.csv' #output file name
		self.dbFile = '../db/minion_1016_sqlite.db' #for the OpSim database
		self.db = None
		self.cursor = None

		#dictionaries -- could be handled by the multiprocessing manager, redefined in driver
		self.return_dict = dict()

		self.csvwriter = None #will hold the csvwriter object

		#some counters
		self.n_totalrun = 0
		self.n_appmag_failed = 0
		self.n_incl_failed = 0
		self.n_period_failed = 0
		self.n_radius_failed = 0


		self.seed = None

		self.SED = None

	#database manipulation
	def getCursors(self):
		#gets SQlite cursor to pull information from OpSim
		self.db = sqlite3.connect(self.dbFile)
		cursor = self.db.cursor()

		cursor.execute("SELECT fieldid, expDate, filter FROM summary") 
		self.summaryCursor = np.array(cursor.fetchall()) #NOTE: this takes a LONG time
		print("have summary cursor.")

		cursor.execute("SELECT fieldid, fieldra, fielddec FROM field")
		self.fieldCursor = np.array(cursor.fetchall()) #NOTE: this takes a LONG time
		print("have field cursor.")



	def make_gatspy_plots(self, j):
		EB = self.return_dict[j]

		#colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
		#print([ matplotlib.colors.to_hex(c) for c in colors])
		colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

		f, ax = plt.subplots(len(self.filters)+1, 2)

		LSM = EB.LSM
		period = EB.period
		pds = np.linspace(0.2, 2.*period, 10000)
		for ii,filt in enumerate(self.filters):

			drng = max(EB.obsDates[filt]) - min(EB.obsDates[filt])

			phase_obs = np.array([(tt % period)/period for tt in EB.obsDates[filt]])
			scores = EB.LSSmodel[filt].score(pds)
			mag_obs = EB.appMagObs[filt]
			mag = EB.appMag[filt]
			LSS = EB.LSS[filt]

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

		if (self.doLSM):
			plt.locator_params(axis='y', nticks=2)
			P_multi = EB.LSMmodel.periodogram(pds)
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

	def run_ellc_gatspy(self, j):
		#this is the general simulation - ellc light curves and gatspy periodograms

		EB = self.return_dict[j]

		#for the multiband gatspy fit
		allObsDates = np.array([])
		allAppMagObs = np.array([])
		allAppMagObsErr = np.array([])
		allObsFilters = np.array([])

		if (self.verbose):
			print("in run_ellc_gatspy")


		for i, filt in enumerate(self.filters):

			#observe the EB (get dates, create the light curve for this filter)
			EB.observe(filt)
			EB.LSS[filt] = -999.

			if (EB.obsDates[filt][0] != None and min(EB.appMagObs[filt]) > 0):

				#run gatspy for this filter
				drng = max(EB.obsDates[filt]) - min(EB.obsDates[filt])
				#print("filter, nobs", filt, len(EB.obsDates[filt]))
				if (self.useFast and len(EB.obsDates[filt]) > 50):
					model = LombScargleFast(fit_period = True)
				else:
					model = LombScargle(fit_period = True)
				model.optimizer.period_range = (0.2, drng)
				model.fit(EB.obsDates[filt], EB.appMagObs[filt], EB.appMagObsErr[filt])
				EB.LSS[filt] = model.best_period
				EB.LSSmodel[filt] = model
				EB.maxDeltaMag = max(EB.deltaMag[filt], EB.maxDeltaMag)

				#to use for the multiband fit
				allObsDates = np.append(allObsDates, EB.obsDates[filt])
				allAppMagObs = np.append(allAppMagObs, EB.appMagObs[filt])
				allAppMagObsErr = np.append(allAppMagObsErr, EB.appMagObsErr[filt])
				allObsFilters = np.append(allObsFilters, np.full(len(EB.obsDates[filt]), filt))

				if (self.verbose): 
					print(j, 'filter = ', filt)  
					print(j, 'obsDates = ', EB.obsDates[filt][0:10])
					print(j, 'appMagObs = ', EB.appMagObs[filt][0:10])
					print(j, 'delta_mag = ', EB.deltaMag[filt])
					print(j, 'LSS = ',EB.LSS[filt])

		if (len(allObsDates) > 0 and self.doLSM): 
			drng = max(allObsDates) - min(allObsDates)
			if (self.useFast and len(allObsDates) > 50*len(self.filters)):
				model = LombScargleMultibandFast(fit_period = True)
			else:
				model = LombScargleMultiband(Nterms_band=self.n_band, Nterms_base=self.n_base, fit_period = True)
			model.optimizer.period_range = (0.2, drng)
			model.fit(allObsDates, allAppMagObs, allAppMagObsErr, allObsFilters)
			EB.LSM = model.best_period
			EB.LSMmodel = model
			if (self.verbose): 
				print(j, 'LSM =', EB.LSM)


		#not sure if I need to do this...
		self.return_dict[j] = EB




	def getEB(self, line, i):
		EB = EclipsingBinary()

		# EB.seed = self.seed + i
		EB.initializeSeed()
		EB.SED = self.SED

		#solar units
		EB.m1 = line[0]
		EB.m2 = line[1]
		EB.r1 = line[4]
		EB.r2 = line[5]
		EB.L1 = line[6]
		EB.L2 = line[7]
		EB.period = 10.**line[2] #days
		EB.eccentricity = line[3]
		EB.inclination = line[12] #radians
		EB.omega = line[13] #radians
		EB.dist = line[11] #kpc
		#pc
		EB.xGx = line[8] 
		EB.yGx = line[9] 
		EB.zGx = line[10] 

		EB.t_zero = np.random.random() * EB.period

		#for observations
		EB.doOpSim = self.doOpSim
		EB.years = self.years
		EB.totaltime = self.totaltime 
		EB.cadence= self.cadence 
		EB.Nfilters = len(self.filters)
		EB.verbose = self.verbose
		if (self.doOpSim):
			EB.summaryCursor = self.summaryCursor
			EB.fieldCursor = self.fieldCursor
		EB.initialize()

		#some counters for how many EBs we could potentially observe with LSST
		self.n_totalrun += 1
		self.n_appmag_failed += EB.appmag_failed
		self.n_incl_failed += EB.incl_failed
		self.n_period_failed += EB.period_failed
		self.n_radius_failed += EB.radius_failed
			
		return EB


	def writeOutputLine(self, EB, header = False):
		if (header):
			self.csvwriter.writerow(['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'RA', 'Dec', 'd', 'nobs','appMagMean', 'maxDeltaMag', 'mag_failure', 'incl_failure', 'period_failure', 'radius_failure', 'u_LSS_PERIOD', 'g_LSS_PERIOD', 'r_LSS_PERIOD', 'i_LSS_PERIOD', 'z_LSS_PERIOD', 'y_LSS_PERIOD','LSM_PERIOD'])

		else:
			output = [EB.period, EB.m1, EB.m2, EB.r1, EB.r2, EB.eccentricity, EB.inclination, EB.RA, EB.Dec, EB.dist, EB.nobs, EB.appMagMeanAll, EB.maxDeltaMag, EB.appmag_failed, EB.incl_failed, EB.period_failed, EB.radius_failed]

			#this is for gatspt
			for filt in self.filters:
				output.append(EB.LSS[filt]) 
			output.append(EB.LSM) 
			self.csvwriter.writerow(output)	

	def initialize(self):
		if (self.seed == None):
			np.random.seed()
		else:
			np.random.seed(seed = self.seed)

		self.SED = SED()
		self.SED.readFilters()
