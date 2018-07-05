#!/software/anaconda3.6/bin/python

from LSSTEBClass import LSSTEBClass
from BreivikGalaxyClass import BreivikGalaxyClass
import multiprocessing, logging
import csv
import argparse
import numpy as np
from mpi4py import MPI


def define_args():
	parser = argparse.ArgumentParser()

	parser.add_argument("-n", "--n_cores", 		type=int, help="Number of cores [0]")
	parser.add_argument("-c", "--n_bin", 		type=int, help="Number of binaries per process [100000]")
	parser.add_argument("-i", "--gal_file", 	type=str, help="Galaxy input file name")
	parser.add_argument("-o", "--output_file", 	type=str, help="output file name")
	parser.add_argument("-a", "--n_band", 		type=int, help="Nterms_band input for gatspy [2]")
	parser.add_argument("-b", "--n_base", 		type=int, help="Nterms_base input for gatspy [2]")
	parser.add_argument("-s", "--seed", 		type=int, help="random seed []")
	parser.add_argument("-p", "--plots", 		action='store_true', help="Set to create plots")
	parser.add_argument("-v", "--verbose", 		action='store_true', help="Set to show verbose output")
	parser.add_argument("-l", "--opsim", 		action='store_true', help="set to run LSST OpSim, else run nobs =")

	#https://docs.python.org/2/howto/argparse.html
	args = parser.parse_args()
	#to print out the options that were selected (probably some way to use this to quickly assign args)
	opts = vars(args)
	options = { k : opts[k] for k in opts if opts[k] != None }
	print(options)

	return args

def apply_args(worker, args):

	if (args.n_cores is not None):
		worker.n_cores = args.n_cores
	if (worker.n_cores > 1):
		worker.do_parallel = True
	else:	
		worker.n_cores = 1
		worker.do_parallel = False 

	if (args.n_bin is not None):
		worker.n_bin = args.n_bin

	if (args.gal_file is not None):
		worker.GalaxyFile = args.gal_file
		
	if (args.output_file is not None):
		worker.ofile = args.output_file

	if (args.n_band is not None):
		worker.n_band = args.n_band
	if (args.n_base is not None):
		worker.n_base = args.n_base

	worker.do_plot = args.plots
	worker.verbose = args.verbose
	worker.doOpSim = args.opsim

	#set the random seed
	if (args.seed is not None):
		worker.seed = args.seed




if __name__ == "__main__":

	comm = MPI.COMM_WORLD
	size = comm.Get_size()
	rank = comm.Get_rank()

	sendbuf = None
	root = 0

	args = define_args()
	if (args.n_bin == None):
		args.n_bin = 2

	n_binr = args.n_bin
	n_bins = n_binr*size
	nfields = 15 #this is the number of fields returned from Katie's KDE
	# sendbuf = np.empty((n_bins, nfields), dtype='float64')
	# recvbuf = np.empty((n_binr, nfields), dtype='float64')
	sendbuf = np.empty((size, n_binr*nfields), dtype='float64')
	recvbuf = np.empty(n_binr*nfields, dtype='float64')

	#only the rank 0 process will generate the galactic model
	if (rank == 0):

		#only do Katie's model once (though it's not a huge time sync...)
		print("in rank 0")

		#Katie's code to generate the binaries
		g = BreivikGalaxyClass()
		g.n_bin = n_bins

		#define the correct paths to the input files and db
		g.GalaxyFile ='/projects/p30137/ageller/EB-LSST/input/dat_ThinDisk_12_0_12_0.h5' #for Katie's model
		g.GalaxyFileLogPrefix ='/projects/p30137/ageller/EB-LSST/input/fixedPopLogCm_'

		#now get the binaries
		gxDat = g.LSSTsim()

		print("reshaping to send to other processes")
		sendbuf = np.reshape(gxDat, (size, n_binr*nfields))


	#scatter to the all of the processes
	comm.Scatter(sendbuf,recvbuf, root=root) 
	#now reshape again to get back to the right format
	gxDat = np.reshape(recvbuf, (n_binr, nfields))


	#Our LSST EB class to use gatspy and ellc
	worker = LSSTEBClass()
	#check for command-line arguments
	apply_args(worker, args)	
	if (worker.seed == None):
		worker.seed = 1234
	worker.seed += rank
	worker.initializeSeed() #right now this just sets the random seed
	#worker.doLSM = False

	#get the summary cursor for OpSim, if necessary
	if (worker.doOpSim):
		worker.dbFile = '/projects/p30137/ageller/EB-LSST/db/minion_1016_sqlite.db' #for the OpSim database	
		print('Getting OpSim cursors...')
		worker.getCursors()


	#set up the output file
	worker.ofile = 'output_files/'+str(rank).zfill(3) + worker.ofile
	csvfile = open(worker.ofile, 'wt')	
	worker.csvwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
	#write header
	worker.writeOutputLine(None, header=True)

	j = 0 #this is used with multiprocessing to keep track of the return_dict, but not needed here
	for i, line in enumerate(gxDat):
		print(rank, i)
		line = gxDat[i]

		#define the binary parameters
		EB = worker.getEB(line, i)
		EB.lineNum = i

		if (EB.observable):
			worker.return_dict[j] = EB
			worker.run_ellc_gatspy(j)
			EB = worker.return_dict[j]

		print(rank, EB.period)
		worker.writeOutputLine(EB)
		csvfile.flush()




	csvfile.close()


