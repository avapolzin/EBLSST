### I need to pull out gatspy into it's own function for easier plug-and-play
### WE NEED TO DOUBLE CHECK getsig2rand

#python libraries
import numpy as np
import argparse
import time
import os
from astropy import units, constants
import astropy.stats as astroStats
import csv
import multiprocessing, logging
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
#OpSim (in OpSimRun1) - http://ops2.lsst.org/docs/current/architecture.html
import sqlite3
from scipy.interpolate import interp1d
import scipy.stats

#3rd party codes
import ellc
import gatspy
from gatspy import datasets, periodic
from gatspy.periodic import LombScargleMultiband, LombScargle, LombScargleFast, LombScargleMultibandFast

#our codes
from BreivikGalaxyClass import BreivikGalaxyClass as galaxy

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
		self.n_cores = 1

		self.ofile = 'output_file.csv' #output file name
		self.GalaxyFile ='../input/dat_ThinDisk_12_0_12_0.h5' #for Katie's model
		self db = sqlite3.connect('../db/minion_1016_sqlite.db')


		#dictionaries handled by the multiprocessing manager
		self.manager = multiprocessing.Manager()
		self.return_dict = self.manager.dict()
		self.plot_dict = self.manager.dict()	

		self.parser = argparse.ArgumentParser()


		#some counters
		self.n_totalrun = 0
		self.n_appmag_failed = 0
		self.n_incl_failed = 0
		self.n_period_failed = 0
		self.n_radius_failed = 0



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

	def run_ellc_gatspy(self, j, return_dict, plot_dict):
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


		for i, filt in enumerate(self.filters):

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
				#create the light curve for this filter




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


			if (len(totalmag) > 0): 
				t = np.array(totalt)
				mag = np.array(totalmag)
				dmag = np.array(totaldmag)
				filts = np.array(totalfilts)
				drng = max(t) - min(t)
				print("drng",drng)

				model = LombScargleMultiband(Nterms_band=n_band, Nterms_base=n_base, fit_period = True)
				#model = LombScargleMultibandFast(Nterms=2, fit_period = True) #this is not working in parallel for some reason
				model.optimizer.period_range = (0.2, drng)
				model.fit(t, mag, dmag, filts)
				LSM = model.best_period
				print ('LSM running')
				return_dict[j] = return_dict[j] + [LSM]
				#n_totalrun += 1
				print ('LSM in return_dict')
				return_dict[j] = return_dict[j] + [delta_mag]
			else:
				return_dict[j] = return_dict[j] + [-999.]
				return_dict[j] = return_dict[j] + [-999.]


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


	def checkIfObservable(self, EB):

		self.n_totalrun += 1
	
		if (EB.appmag <= 11. or EB.appmag>= 24.):
			EB.appmag_failed = 1
			self.n_appmag_failed += 1 
		
		incl = EB.inclination * 180/np.pi
		min_incl = 90. - np.arctan2(EB.r1 + EB.r2, 2.*EB.a)*180./np.pi
		if (incl <= min_incl):
			EB.incl_failed = 1
			self.n_incl_error += 1
		
		if (EB.period >= 2. * self.totaltime):
			EB.period_failed = 1
			self.n_period_error +=1
				
		if (EB.R_1 <= 0 or EB.R_1 >=1 or EB.R_2 <=0 or EB.R_2 >= 1 or EB.R_1e >=1 or EB.R_2e >=1):
			EB.radius_failed = 1
			self.n_radius_failed += 1
			
		if (EB.radius_failed or EB.period_failed or EB.incl_failed or EB.appmag_failed):
			EB.observable = False


	def getEB(self, EB, line):
		EB = EBclass()

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

		EB.t_zero = np.random.random() * period

		EB.initialize()

		self.checkIfObservable(EB)

		print(EB)

		return EB



	#run everything
	def runAll(self):

		self.define_args()
		self.apply_args()

		#get the summary cursor
		if (self.do_opsim):
			self.getSummaryCursor()

		#run Katie's code to generate the binaries
		g = Galaxy()
		gxDat = g.LSSTsim()
		print(gxDat)

		#https://docs.python.org/2/howto/argparse.html
		n_appmag_error = 0
		n_incl_error = 0
		n_period_error = 0
		n_R_error = 0
		n_totalrun = 0


		#if we want logging
		#logger = multiprocessing.log_to_stderr()
		##logger.setLevel(logging.DEBUG)
		#logger.setLevel(multiprocessing.SUBDEBUG)



		csvfile = open(ofile, 'wt')
		csvwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
		csvwriter.writerow(['PERIOD', 'MASS_1', 'MASS_2', 'RADIUS_1', 'RADIUS_2', 'a', 'INCLINATION', 'MIN_INCLINATION', 'xGx', 'yGx', 'zGx', 'dist_kpc', 'eccentricity', 'max(app_magnitude)', 'appmag_error', 'inclination_error', 'period_error', 'radius_error', 'u_LSS_PERIOD', 'g_LSS_PERIOD', 'r_LSS_PERIOD', 'i_LSS_PERIOD', 'z_LSS_PERIOD', 'LSM_PERIOD', 'delta_mag', 'chi2', 'mean(dmag)'])
		

		jobs = []
		j = 0

		for i, line in enumerate(gxDat):

			#define the binary parameters
			EB = getEB(line)

			#get the field ID number from OpSim where this binary would be observed
			if (self.do_opsim):
				EB.OpSimID = self.getFieldID(EB.RA, EB.Dec)

			output = [period, m_1, m_2, R_1, R_2, a, incl, min_incl, xGx, yGx, zGx, dist_kpc, ecc, appmag, appmag_error, incl_error, period_error, R_error]
				
			if (EB.observable):
				return_dict[j] = EB
				plot_dict[j] = dict()
				p = multiprocessing.Process(target=run_ellc_gatspy, args=(j, return_dict, plot_dict))
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


