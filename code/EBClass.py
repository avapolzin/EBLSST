import numpy as np
from astropy import units, constants
from astropy.coordinates import SkyCoord, ICRS
import ellc

class EBClass(object):
	def __init__(self, *args,**kwargs):

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

		#these will be calculated after calling self.initialize()
		self.RL1 = None
		self.RL2 = None
		self.T1 = None
		self.T2 = None
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
		self.appMagMean = None
		self.absMagMean = None
		self.Mbol = None


		#for light curves
		self.filters = None
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
		self.LSSTcursor = None
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
			self.appMag[filt] = self.appMagMean + magn   
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

		RA = self.LSSTcursor[:,1].astype(float)
		Dec = self.LSSTcursor[:,2].astype(float)
		dbCoord = SkyCoord(ra = RA*units.degree, dec = Dec*units.degree, frame='icrs')
		inCoord = SkyCoord(ra = myRA*units.degree, dec = myDEC*units.degree, frame='icrs')

		imin, sep2d, dist3d = inCoord.match_to_catalog_sky(dbCoord)

		dbID = (self.LSSTcursor[imin,0]).astype('int') 

		mask = np.where(sep2d.to(units.degree).value > deglim)

		#this check apparently isn't necessary because it looks like the entire sky is covered with fieldIDs, but I suppose some of these fieldIDs don't have any observation dates (in the northern hemisphere)
		if (len(mask[0]) > 0):
			print (mask[0])
			print("WARNING: coordinate outside LSST FOV", myRA[mask], myDec[mask])
			dbID[mask] = -999

		return dbID

	def getOpSimDates(self, filtin):
		#matches FieldID to existing OpSim database ID and matches observation filters to get dates (in seconds since the start of the 
		# survey)
		date = c[:,3]
		FieldID = c[:,0]
		filt = c[:,4]

		posIDFilt = np.where(np.logical_and(FieldID == str(self.OpSimID), filt == filtin[:-1]))
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

	def checkIfObservable(self):

		if (self.appMagMean <= 11. or self.appMagMean>= 24.):
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

	def initialize(self):
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

		coord = SkyCoord(x=self.xGx, y=self.yGx, z=self.zGx, unit='pc', representation='cartesian', frame='galactocentric')
		self.RA = coord.icrs.ra.to(units.deg).value
		self.Dec = coord.icrs.dec.to(units.deg).value

		MbolSun = 4.74
		self.Mbol = MbolSun - 2.5*np.log10(self.L1 + self.L2)

		######################
		#Now I'm assuming that the magnitude is bolometric -- but we should account for redenning in each filter
		#######################
		self.absMagMean = self.Mbol
		self.appMagMean = self.absMagMean + 5.*np.log10(self.dist*100.)  #multiplying by 1000 to get to parsec units

		#if we're using OpSim, then get the field ID
		#get the field ID number from OpSim where this binary would be observed
		if (self.doOpSim):
			self.OpSimID = self.getFieldID(self.RA, self.Dec)

		for f in self.filters:
			self.LSS[f] = -999.

		#check if we can observe this (not accounting for the location in galaxy)
		self.checkIfObservable()


	def observe(self, filt):

		#get the observation dates
		if (self.doOpSim):
			self.obsDates[filt] = getOpSimDates(filt)
		else:
			nobs = int(round(self.totaltime / (self.cadence * self.Nfilters)))
			self.obsDates[filt] = np.sort(self.totaltime * np.random.random(size=nobs))

		#get the light curve, and related information
		self.setLightCurve(filt)

