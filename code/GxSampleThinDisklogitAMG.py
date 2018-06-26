
# Code: GxSampleThinDisk.py
# Version: 1
# Version changes: SAMPLE GALACTIC POPULATION ACCORDING TO SPECIFIED FLAGS
#
# Edited on:  27 MAR 2017	

##############################################################################
#  IMPORT ALL NECESSARY PYTHON PACKAGES
##############################################################################
import math
import scipy.special as ss
import scipy.stats
import multiprocessing as mp
import numpy as np
import os


##############################################################################
#  DEFINE THE FUNCTIONS YOU WILL USE!
##############################################################################
def GxSample(x, nCores, pop, sampleKernel, popFlag, nBin, bw, nEcc, Tobs, output):
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
	gxFile = 'gxRealization_'+str(x)+'_'+str(popFlag)+'.dat'
	binDat = [] 
	
	nSample = int(nBin/float(nCores))
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
	print('you should see me!'       ) 
	
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
		
	inc = np.random.uniform(0,math.pi/2,len(m1))
	OMEGA = np.random.uniform(0,2*math.pi,len(m1))
	omega = np.random.uniform(0,2*math.pi,len(m1))

	binDat = np.vstack((m1, m2, porb, ecc, rad1, rad2, Lum1, Lum2, xGX, yGX, zGX, dist_kpc, inc, OMEGA, omega)).T
	radTotAU = (rad1+rad2)/rsun_in_au
	radAng = radTotAU/dist
	binEclipseIndex, = np.where(radAng>inc*4.8e-6)
 
	print(len(binEclipseIndex), len(binDat))

	np.savetxt(gxFile, binDat, delimiter = ',')     
			
	#gxFile.close() 
	output.put(np.shape(binDat)) 
		
