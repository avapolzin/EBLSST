from LSSTEBClass import LSSTEBClass
from BreivikGalaxyClass import BreivikGalaxyClass
import multiprocessing, logging
import csv

if __name__ == "__main__":

	worker = LSSTEBClass()

	#check for command-line arguments
	worker.define_args()
	worker.apply_args()

	#get the summary cursor
	if (worker.doOpSim):
		print('Getting OpSim cursor...')
		worker.getSummaryCursor()

	#run Katie's code to generate the binaries
	g = BreivikGalaxyClass()
	gxDat = g.LSSTsim()

	#if we want logging
	#logger = multiprocessing.log_to_stderr()
	##logger.setLevel(logging.DEBUG)
	#logger.setLevel(multiprocessing.SUBDEBUG)

	#set up the output file
	csvfile = open(worker.ofile, 'wt')	
	worker.csvwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

	#write header
	worker.csvwriter.writerow(['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'RA', 'Dec', 'd', 'appMagMean', 'maxDeltaMag', 'mag_failure', 'incl_failure', 'period_failure', 'radius_failure', 'u_LSS_PERIOD', 'g_LSS_PERIOD', 'r_LSS_PERIOD', 'i_LSS_PERIOD', 'z_LSS_PERIOD', 'y_LSS_PERIOD','LSM_PERIOD', 'delta_mag', 'chi2', 'mean(dmag)'])
	

	jobs = []
	j = 0
	for i, line in enumerate(gxDat):

		#define the binary parameters
		EB = worker.getEB(line, i)
		EB.lineNum = i

		if (EB.observable):
			worker.return_dict[j] = EB
			p = multiprocessing.Process(target=worker.run_ellc_gatspy, args=(j,))
			jobs.append(p)
			j += 1
		else:
			worker.writeOutputLine(EB)


		if (len(jobs) == worker.n_cores or (i >= worker.n_bin and len(jobs) > 0)):
			#could print this to look at progress
			# print ('n_appmag_failure = ',  worker.n_appmag_failure)
			# print ('n_incl_failure = ', worker.n_incl_failure)
			# print ('n_period_failure = ', worker.n_period_failure)
			# print ('n_R_failure = ', worker.n_radius_failure)
			# print ('n_totalrun = ', worker.n_totalrun)

			#run the jobs
			if (worker.do_parallel):
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

			#print the output
			for j in range(n_procs):
				worker.writeOutputLine(j)
				csvfile.flush()
				if (do_plot):
					if ('LSM' in plot_dict[j]):
						 worker.make_gatspy_plots(j)


			 #raise Exception('stopping')
			jobs = []
			j = 0


	csvfile.close()


