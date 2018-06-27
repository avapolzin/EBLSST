import numpy as np
from astropy import units, constants
from astropy.coordinates import SkyCoord, ICRS
import ellc


class EBclass(object):
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

		self.M_H = 0
		self.ld_1 = 'claret'
		self.ld_2 = 'claret'
		self.grid_1 = 'default'
		self.grid_2 = 'default'
		self.shape_1 = 'sphere'
		self.shape_2 = 'sphere'

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
		self.appMag = None
		self.absMag = None
		self.Mbol = None

		#will be set within the LSST-EBclass
		self.observable = True
		self.appmag_failed = 0
		self.incl_failed = 0
		self.period_failed = 0
		self.radius_failed = 0
		self.OpSimID = None


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
		self.Dec = oord.icrs.dec.to(units.deg).value

		MbolSun = 4.74
		self.Mbol = MbolSun - 2.5*np.log10(self.L1 + self.L2)

		######################
		#Now I'm assuming that the magnitude is bolometric -- but we should account for redenning in each filter
		#######################
		self.absMag = self.Mbol
		self.appMag = self.absMag + 5.*np.log10(self.dist*100.)  #multiplying by 1000 to get to parsec units


	def LightCurve(self, obsDates, filter):
		#in case the user did not initialize
		if (self.T1 == None):
			self.initialize()

		#limb darkenning
		ldy_filt = ellc.ldy.LimbGravityDarkeningCoeffs(filter)
		a1_1, a2_1, a3_1, a4_1, y = ldy_filt(self.T1, self.g1, self.M_H)
		a1_2, a2_2, a3_2, a4_2, y = ldy_filt(self.T2, self.g2, self.M_H)
		ldc_1 = [a1_1, a2_1, a3_1, a4_1] 
		ldc_2 = [a1_2, a2_2, a3_2, a4_2]

		#light curve
		lc = ellc.lc(obsDates, ldc_1=ldc_1, ldc_2=ldc_2, 
			t_zero=self.t_zero, period=self.period, a=self.a, q=self.q,
			f_c=self.f_c, f_s=self.f_s, ld_1=self.ld_1,  ld_2=self.ld_2,
			radius_1=self.R_1, radius_2=self.R_2, incl=self.inclination, sbratio=self.sbratio, 
			shape_1=self.shape_1, shape_2=self.shape_2, grid_1=self.grid_1,grid_2=self.grid_2) 

		return lc
