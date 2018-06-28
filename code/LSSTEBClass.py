### I need to pull out gatspy into it's own function for easier plug-and-play
### WE NEED TO DOUBLE CHECK getsig2rand

#python libraries
import numpy as np
import argparse
import os
from astropy import units, constants
import csv
import multiprocessing, logging
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

#class that defines the eclipsing binaries
from EBClass import EBClass

class LSSTEBClass(object):

	def __init__(self, *args,**kwargs):

		#NOTE: these need to be defined on the command line.  The default, if not defined, will be False
		self.do_plot = False 
		self.verbose = False
		self.doOpSim = False

		self.do_parallel = False 
		self.years = 10.
		self.totaltime = 365.* self.years
		self.cadence = 3.

		#self.filters = ['u_', 'g_', 'r_', 'i_', 'z_']
		self.filters = ['u_', 'g_', 'r_', 'i_', 'z_', 'y_']

		self.n_bin = 100000
		self.n_band = 2
		self.n_base = 2  
		self.n_cores = 1

		self.ofile = 'output_file.csv' #output file name
		self.GalaxyFile ='../input/dat_ThinDisk_12_0_12_0.h5' #for Katie's model
		self.db = sqlite3.connect('../db/minion_1016_sqlite.db')
		self.cursor = None

		#dictionaries handled by the multiprocessing manager
		self.manager = multiprocessing.Manager()
		self.return_dict = self.manager.dict()

		self.parser = argparse.ArgumentParser()
		self.csvwriter = None #will hold the csvwriter object

		#some counters
		self.n_totalrun = 0
		self.n_appmag_failed = 0
		self.n_incl_failed = 0
		self.n_period_failed = 0
		self.n_radius_failed = 0



	#database manipulation
	def getSummaryCursor(self):
		#gets SQlite cursor to pull information from OpSim
		cursor = self.db.cursor()
		cursor.execute("SELECT fieldid, fieldra, fielddec, expDate, filter FROM summary") 
		self.cursor = np.array(cursor.fetchall()) #NOTE: this takes a LONG time
		print("have cursor.")



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
				#model = LombScargle(fit_period = True)
				model = LombScargleFast(fit_period = True)
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
					print(j, 'delta_mag = ', EB.delta_mag[filt])
					print(j, 'LSS = ',EB.LSS[filt])

		if (len(allObsDates) > 0): 
			drng = max(allObsDates) - min(allObsDates)
			model = LombScargleMultiband(Nterms_band=self.n_band, Nterms_base=self.n_base, fit_period = True)
			model.optimizer.period_range = (0.2, drng)
			model.fit(allObsDates, allAppMagObs, allAppMagObsErr, allObsFilters)
			EB.LSM = model.best_period
			EB.LSMmodel = model
			if (self.verbose): 
				print(j, 'LSM =', EB.LSM)


		#not sure if I need to do this...
		self.return_dict[j] = EB

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
		#https://docs.python.org/2/howto/argparse.html
		args = self.parser.parse_args()
		#to print out the options that were selected (probably some way to use this to quickly assign args)
		opts = vars(args)
		options = { k : opts[k] for k in opts if opts[k] != None }
		print(options)

		if (args.n_cores is not None):
			self.n_cores = args.n_cores
			if (n_cores > 1):
				self.do_parallel = True
		if (self.n_cores < 1):
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
		self.doOpSim = args.opsim

		#set the random seed
		if (args.seed is not None):
			np.random.seed(seed = args.seed)
		else:
			np.random.seed()





	def getEB(self, line, i):
		EB = EBClass()


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
		EB.filters = self.filters
		EB.doOpSim = self.doOpSim
		EB.LSSTcursor = self.cursor
		EB.years = self.years
		EB.totaltime = self.totaltime 
		EB.cadence= self.cadence 
		EB.Nfilters = len(self.filters)
		EB.initialize()

		#some counters for how many EBs we could potentially observe with LSST
		self.n_totalrun += 1
		self.n_appmag_failed += EB.appmag_failed
		self.n_incl_failed += EB.incl_failed
		self.n_period_failed += EB.period_failed
		self.n_radius_failed += EB.radius_failed
			
		return EB


	def writeOutputLine(self, EB):
		output = [EB.period, EB.m1, EB.m2, EB.r1, EB.r2, EB.eccentricity, EB.inclination, EB.RA, EB.Dec, EB.dist, EB.appMagMean, EB.maxDeltaMag, EB.appmag_failed, EB.incl_failed, EB.period_failed, EB.radius_failed]

		#this is for gatspt
		for filt in self.filters:
			output.append(EB.LSS[filt]) 
		output.append(EB.LSM) 
		self.csvwriter.writerow(output)	
